// attach_sim.cpp
//
// Stateful simulator of attachment marks (see R/attachment.R). Old perps join
// an event independently with p_v = 1 - exp(-m w_v / W). Equivalently, draw
// H ~ Poisson(m) hits and give each hit to perp v with probability w_v / W;
// a perp joins if it gets at least one hit. The decay exp(-beta_e t) is common
// to all weights and cancels from w_v / W, so each perp's sampling weight
//   s_v = u_v (d_v + delta) exp(theta x_v) exp(beta_e (last_v - t_ref))
// changes only when its own state does, and is kept in a Fenwick tree. Each
// event then costs O((m + changes) log n_perps). Random numbers come from R,
// so set.seed() makes runs reproducible.

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
  double m, beta_e, delta, theta, window, lambda_new, lambda_new0, size_new, frailty_shape;
  bool use_degree, use_cof, log1p_scale, negbin, min_one;
  // perp state
  std::vector<double> u, last;
  std::vector<int> deg, x;
  Fenwick s;
  double t_ref;
  std::unordered_map<long long, double> tie;
  typedef std::pair<double, long long> Expiry;
  std::priority_queue<Expiry, std::vector<Expiry>, std::greater<Expiry> > expiries;
  std::vector<int> stamp;
  int stamp_id;

  double weight(int v) const {
    const double h = log1p_scale ? std::log1p((double)x[v]) : (double)x[v];
    return u[v] * (use_degree ? deg[v] + delta : 1.0) * (use_cof ? std::exp(theta * h) : 1.0) *
      std::exp(beta_e * (last[v] - t_ref));
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
SEXP attach_sim_new(List params, bool use_degree, bool use_cof, bool log1p_scale,
                    double window, bool negbin, bool min_one, double frailty_shape,
                    double t0) {
  AttachSim* sim = new AttachSim();
  auto get = [&](const char* nm, double dflt) {
    return params.containsElementNamed(nm) ? as<double>(params[nm]) : dflt;
  };
  sim->m = get("m", NA_REAL);
  sim->beta_e = get("beta_edges", NA_REAL);
  sim->delta = get("delta", 0.0);
  sim->theta = get("theta", 0.0);
  sim->lambda_new = get("lambda_new", NA_REAL);
  sim->lambda_new0 = get("lambda_new0", 0.0);
  sim->size_new = get("size_new", 1.0);
  sim->use_degree = use_degree;
  sim->use_cof = use_cof;
  sim->log1p_scale = log1p_scale;
  sim->window = window;
  sim->negbin = negbin;
  sim->min_one = min_one;
  sim->frailty_shape = frailty_shape;
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
    for (int v = 0; v < n_old; ++v) S.s.value[v] = S.weight(v);
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
  for (size_t a = 0; a < changed.size(); ++a) S.s.set(changed[a], S.weight(changed[a]));

  // Old perps: Poisson(m) hits spread in proportion to the weights.
  std::vector<int> joined;
  if (n_old > 0) {
    const double W = S.s.total();
    const int hits = W > 0.0 ? (int)R::rpois(S.m) : 0;
    ++S.stamp_id;
    for (int h = 0; h < hits; ++h) {
      const int v = S.s.find(R::unif_rand() * W);
      if (S.stamp[v] != S.stamp_id) {
        S.stamp[v] = S.stamp_id;
        joined.push_back(v);
      }
    }
  }

  // New perps, given how many old ones joined.
  const bool shifted = S.min_one && joined.empty();
  const int n_new = (int)shifted + (int)S.draw_count(shifted ? S.lambda_new0 : S.lambda_new);
  for (int a = 0; a < n_new; ++a) {
    S.u.push_back(std::isfinite(S.frailty_shape)
                    ? R::rgamma(S.frailty_shape, 1.0 / S.frailty_shape) : 1.0);
    S.last.push_back(t);
    S.deg.push_back(0);
    S.x.push_back(0);
    S.stamp.push_back(-1);
    S.s.push(0.0);
  }

  // Attendance, recency and ties among everyone at the event.
  std::vector<int> perps(joined);
  for (int a = 0; a < n_new; ++a) perps.push_back(n_old + a);
  changed.clear();
  ++S.stamp_id;
  for (size_t a = 0; a < perps.size(); ++a) {
    const int v = perps[a];
    ++S.deg[v];
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
  for (size_t a = 0; a < changed.size(); ++a) S.s.set(changed[a], S.weight(changed[a]));

  return List::create(_["joined"] = wrap(joined), _["n_new"] = n_new);
}
