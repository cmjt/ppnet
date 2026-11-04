// attach_sim.cpp
//
// Stateful simulator of attachment marks (see R/attachment.R). Old perps join
// an event independently with p_v = 1 - exp(-u_v m w_v / W), where u_v is the
// perp's random effect (1 without one) and W = sum of w over old perps, without
// the random effect. Equivalently, draw H ~ Poisson(m sum_v u_v w_v / W) hits and give
// each hit to perp v with probability u_v w_v / sum_u u_u w_u; a perp joins
// if it gets at least one hit. The decay exp(-beta_e t) is common to all
// weights and cancels from both ratios, so each perp's weight
//   w_v = (d_v + delta) exp(theta x_v + theta_paths y_v) exp(beta_e (last_v - t_ref))
// changes only when its own state does. w and u w are kept in Fenwick trees,
// so each event costs O((m + changes) log n_perps). Random numbers come from
// R, so set.seed() makes runs reproducible.
//
// With 4-cycles (theta_c) joins at an event are drawn from the autologistic
// model whose full conditionals the likelihood uses,
//   logit P(v joins | rest) = log(expm1(u_v q_v)) + theta_c h(z_v),
// with z_v the shared events (within the window) of v with the other old
// joiners: by Gibbs sampling, starting from the independent draw and
// sweeping over every old perp `n_sweeps` times. That costs O(n_sweeps *
// n_perps) per event.

#include <Rcpp.h>
#include <cmath>
#include <functional>
#include <queue>
#include <unordered_map>
#include <vector>
using namespace Rcpp;

namespace {

struct Fenwick {
  std::vector<double> tree, value;
  // Sum of the first k values.
  double prefix(int k) const {
    double s = 0.0;
    for (int j = k; j > 0; j -= j & (-j)) s += tree[j - 1];
    return s;
  }
  // A new node n covers values (n - lowbit(n), n], which can include earlier
  // ones, so it starts from their sum rather than from zero.
  void push(double v) {
    const int n = (int)tree.size() + 1;
    tree.push_back(v + prefix(n - 1) - prefix(n - (n & (-n))));
    value.push_back(v);
  }
  void set(int i, double v) {
    const double d = v - value[i];
    value[i] = v;
    for (int j = i + 1; j <= (int)tree.size(); j += j & (-j)) tree[j - 1] += d;
  }
  void rebuild() {
    std::fill(tree.begin(), tree.end(), 0.0);
    for (int i = 0; i < (int)value.size(); ++i) {
      for (int j = i + 1; j <= (int)tree.size(); j += j & (-j)) tree[j - 1] += value[i];
    }
  }
  double total() const { return prefix((int)tree.size()); }
  // Smallest index whose prefix sum exceeds target.
  int find(double target) const {
    int pos = 0, step = 1;
    while (step * 2 <= (int)tree.size()) step *= 2;
    for (; step > 0; step /= 2) {
      if (pos + step <= (int)tree.size() && tree[pos + step - 1] <= target) {
        pos += step;
        target -= tree[pos - 1];
      }
    }
    return std::min(pos, (int)tree.size() - 1);
  }
};

struct AttachSim {
  // parameters
  double m, beta_e, delta, theta, theta_p, window, lambda_new, lambda_new0, size_new,
    random_effect_shape;
  double theta_c;
  bool use_degree, use_cof, use_paths, use_cyc, log1p_scale, negbin, min_one;
  int n_sweeps;
  // 4-cycles: per perp, partners and the times of their shared events
  std::vector<std::unordered_map<int, std::vector<double> > > shared;
  // perp state: x distinct co-offenders, y 3-paths (co-offenders with
  // multiplicity), both over events within the window
  std::vector<double> u, last;
  std::vector<int> deg, x, y;
  // events within the window that still give 3-paths: (time, perps)
  std::queue<std::pair<double, std::vector<int> > > recent;
  Fenwick s, s0;  // u w for sampling, w for W
  double t_ref;
  std::unordered_map<long long, double> tie;
  typedef std::pair<double, long long> Expiry;
  std::priority_queue<Expiry, std::vector<Expiry>, std::greater<Expiry> > expiries;
  std::vector<int> stamp;
  int stamp_id;

  double weight(int v) const {
    auto h = [&](int n) { return log1p_scale ? std::log1p((double)n) : (double)n; };
    const double eta = (use_cof ? theta * h(x[v]) : 0.0) + (use_paths ? theta_p * h(y[v]) : 0.0);
    return (use_degree ? deg[v] + delta : 1.0) * std::exp(eta) *
      std::exp(beta_e * (last[v] - t_ref));
  }
  void update(int v) {
    const double w = weight(v);
    s0.set(v, w);
    s.set(v, u[v] * w);
  }
  static long long key(int a, int b) {
    return ((long long)std::min(a, b) << 32) | (long long)std::max(a, b);
  }
  double draw_count(double mu) const {
    if (mu <= 0.0) return 0.0;
    return negbin ? R::rnbinom(size_new, size_new / (size_new + mu)) : R::rpois(mu);
  }
};

}  // namespace

// [[Rcpp::export]]
SEXP attach_sim_new(List params, bool use_degree, bool use_cof, bool use_paths, bool use_cyc,
                    bool log1p_scale, double window, bool negbin, bool min_one,
                    double random_effect_shape, double t0, int n_sweeps = 10) {
  AttachSim* sim = new AttachSim();
  auto get = [&](const char* nm, double dflt) {
    return params.containsElementNamed(nm) ? as<double>(params[nm]) : dflt;
  };
  sim->m = get("m", NA_REAL);
  sim->beta_e = get("beta_edges", NA_REAL);
  sim->delta = get("delta", 0.0);
  sim->theta = get("theta", 0.0);
  sim->theta_p = get("theta_paths", 0.0);
  sim->theta_c = get("theta_cycles", 0.0);
  sim->use_cyc = use_cyc;
  sim->n_sweeps = n_sweeps;
  sim->lambda_new = get("lambda_new", NA_REAL);
  sim->lambda_new0 = get("lambda_new0", 0.0);
  sim->size_new = get("size_new", 1.0);
  sim->use_degree = use_degree;
  sim->use_cof = use_cof;
  sim->use_paths = use_paths;
  sim->log1p_scale = log1p_scale;
  sim->window = window;
  sim->negbin = negbin;
  sim->min_one = min_one;
  sim->random_effect_shape = random_effect_shape;
  sim->t_ref = t0;
  sim->stamp_id = 0;
  return XPtr<AttachSim>(sim, true);
}

// One event at time t: returns the 0-based indices of old perps that joined
// and the number of new perps (which get the next indices).
// [[Rcpp::export]]
List attach_sim_step(SEXP ptr, double t) {
  XPtr<AttachSim> xp(ptr);
  AttachSim& S = *xp;
  const int n_old = S.u.size();

  // Keep exp(beta_e (last - t_ref)) bounded by moving the reference time.
  if (S.beta_e * (t - S.t_ref) > 300.0) {
    S.t_ref = t;
    for (int v = 0; v < n_old; ++v) {
      S.s0.value[v] = S.weight(v);
      S.s.value[v] = S.u[v] * S.s0.value[v];
    }
    S.s0.rebuild();
    S.s.rebuild();
  }

  std::vector<int> changed;
  ++S.stamp_id;
  auto touch = [&](int v) {
    if (S.stamp[v] != S.stamp_id) {
      S.stamp[v] = S.stamp_id;
      changed.push_back(v);
    }
  };

  // Co-offender ties whose latest shared event is more than `window` ago.
  while (!S.expiries.empty() && S.expiries.top().first < t) {
    const AttachSim::Expiry top = S.expiries.top();
    S.expiries.pop();
    auto it = S.tie.find(top.second);
    if (it == S.tie.end() || it->second != top.first) continue;
    const int lo = (int)(top.second >> 32), hi = (int)(top.second & 0xffffffffLL);
    --S.x[lo];
    --S.x[hi];
    touch(lo);
    touch(hi);
    S.tie.erase(it);
  }
  // Events more than `window` ago no longer give 3-paths.
  while (!S.recent.empty() && S.recent.front().first + S.window < t) {
    const std::vector<int>& pj = S.recent.front().second;
    for (size_t a = 0; a < pj.size(); ++a) {
      S.y[pj[a]] -= (int)pj.size() - 1;
      touch(pj[a]);
    }
    S.recent.pop();
  }
  for (size_t a = 0; a < changed.size(); ++a) S.update(changed[a]);

  // Old perps: Poisson(m sum u w / W) hits spread in proportion to u w.
  std::vector<int> joined;
  if (n_old > 0) {
    const double W = S.s.total(), W0 = S.s0.total();
    const int hits = W > 0.0 && W0 > 0.0 ? (int)R::rpois(S.m * W / W0) : 0;
    ++S.stamp_id;
    for (int h = 0; h < hits; ++h) {
      const int v = S.s.find(R::unif_rand() * W);
      if (S.stamp[v] != S.stamp_id) {
        S.stamp[v] = S.stamp_id;
        joined.push_back(v);
      }
    }
  }

  // 4-cycles: Gibbs sweeps over the old perps from the independent draw.
  if (S.use_cyc && n_old > 0) {
    const double W0 = S.s0.total(), cut = t - S.window;
    std::vector<char> in(n_old, 0);
    std::vector<double> z(n_old, 0.0);
    // Shared events of u with each partner within the window (pruned).
    auto toggle = [&](int u, double sign) {
      for (auto it = S.shared[u].begin(); it != S.shared[u].end();) {
        std::vector<double>& tt = it->second;
        size_t drop = 0;
        while (drop < tt.size() && tt[drop] < cut) ++drop;
        if (drop > 0) tt.erase(tt.begin(), tt.begin() + drop);
        if (tt.empty()) {
          it = S.shared[u].erase(it);
          continue;
        }
        if (it->first < n_old) z[it->first] += sign * (double)tt.size();
        ++it;
      }
    };
    for (size_t a = 0; a < joined.size(); ++a) {
      in[joined[a]] = 1;
      toggle(joined[a], 1.0);
    }
    auto h = [&](double n) { return S.log1p_scale ? std::log1p(n) : n; };
    for (int sweep = 0; sweep < S.n_sweeps; ++sweep) {
      for (int v = 0; v < n_old; ++v) {
        const double q = S.u[v] * S.m * S.s0.value[v] / W0;
        const double x = std::log(std::expm1(q)) + S.theta_c * h(z[v]);
        const bool y = R::unif_rand() < 1.0 / (1.0 + std::exp(-x));
        if (y != (bool)in[v]) {
          in[v] = y;
          toggle(v, y ? 1.0 : -1.0);
        }
      }
    }
    joined.clear();
    for (int v = 0; v < n_old; ++v) if (in[v]) joined.push_back(v);
  }

  // New perps, given how many old ones joined.
  const bool shifted = S.min_one && joined.empty();
  const int n_new = (int)shifted + (int)S.draw_count(shifted ? S.lambda_new0 : S.lambda_new);
  for (int a = 0; a < n_new; ++a) {
    S.u.push_back(std::isfinite(S.random_effect_shape)
                    ? R::rgamma(S.random_effect_shape, 1.0 / S.random_effect_shape) : 1.0);
    S.last.push_back(t);
    S.deg.push_back(0);
    S.x.push_back(0);
    S.y.push_back(0);
    S.stamp.push_back(-1);
    if (S.use_cyc) S.shared.emplace_back();
    S.s.push(0.0);
    S.s0.push(0.0);
  }

  // Attendance, recency and ties among everyone at the event.
  std::vector<int> perps(joined);
  for (int a = 0; a < n_new; ++a) perps.push_back(n_old + a);
  changed.clear();
  ++S.stamp_id;
  for (size_t a = 0; a < perps.size(); ++a) {
    const int v = perps[a];
    ++S.deg[v];
    S.y[v] += (int)perps.size() - 1;
    S.last[v] = t;
    touch(v);
    for (size_t b = a + 1; b < perps.size(); ++b) {
      const int w = perps[b];
      const long long k = AttachSim::key(v, w);
      if (S.tie.find(k) == S.tie.end()) {
        ++S.x[v];
        ++S.x[w];
        touch(w);
      }
      S.tie[k] = t + S.window;
      if (std::isfinite(S.window)) S.expiries.push(AttachSim::Expiry(t + S.window, k));
    }
  }
  if (std::isfinite(S.window) && perps.size() > 1) S.recent.push(std::make_pair(t, perps));
  if (S.use_cyc) {
    for (size_t a = 0; a < perps.size(); ++a) {
      for (size_t b = a + 1; b < perps.size(); ++b) {
        S.shared[perps[a]][perps[b]].push_back(t);
        S.shared[perps[b]][perps[a]].push_back(t);
      }
    }
  }
  for (size_t a = 0; a < changed.size(); ++a) S.update(changed[a]);

  return List::create(_["joined"] = wrap(joined), _["n_new"] = n_new);
}
