// attach_walk.cpp
//
// Attachment of old perps to events in a bipartite event-perp history. At
// event k each old perp v has weight
//   w_v = c_v * exp(-beta * (t_k - last_v)),   c_v = (d_v + delta) * exp(theta * x_v)
// (d_v: events attended, last_v: time of the latest one, x_v: co-offender
// feature) and joins independently with p_v = 1 - exp(-q_v), q_v = m w_v / W,
// W = sum of w over all old perps. Because sum_v q_v = m, the event's log-
// likelihood is
//   sum_{v joined} log(expm1(q_v)) - m,
// so only W needs every perp, and W is kept as a running sum that decays
// between events and changes only when a perp's state changes.

#include <Rcpp.h>
#include <cmath>
#include <functional>
#include <queue>
#include <unordered_map>
#include <vector>
using namespace Rcpp;

// Parameter-free pass: the state changes of every perp, and which old perps
// joined each event. State snapshots (node, d, last, x) at index k are
// applied, in order, before event k is scored; x counts distinct co-offenders
// sharing an event no more than `window` before.
// [[Rcpp::export]]
List attach_history_cpp(NumericVector t, IntegerVector node_ptr,
                        IntegerVector edge_ptr, IntegerVector edge_i,
                        IntegerVector edge_j, IntegerVector focal,
                        int n_nodes, double window) {
  const int n_ev = t.size();
  std::vector<int> born(n_nodes, -1), deg(n_nodes, 0), x(n_nodes, 0), stamp(n_nodes, -1);
  std::vector<double> last(n_nodes, 0.0);

  std::vector<int> upd_idx, upd_node, upd_deg, upd_x, obs_idx, obs_node;
  std::vector<double> upd_last;

  // Active co-offender ties, keyed by node pair, holding their expiry time.
  std::unordered_map<long long, double> tie;
  typedef std::pair<double, long long> Expiry;
  std::priority_queue<Expiry, std::vector<Expiry>, std::greater<Expiry> > expiries;
  const bool finite_window = std::isfinite(window);

  int stamp_id = 0;
  std::vector<int> changed;
  auto touch = [&](int v) {
    if (stamp[v] != stamp_id) {
      stamp[v] = stamp_id;
      changed.push_back(v);
    }
  };
  auto record = [&](int idx) {
    for (size_t a = 0; a < changed.size(); ++a) {
      const int v = changed[a];
      upd_idx.push_back(idx);
      upd_node.push_back(v);
      upd_deg.push_back(deg[v]);
      upd_last.push_back(last[v]);
      upd_x.push_back(x[v]);
    }
    changed.clear();
    ++stamp_id;
  };

  for (int k = 0; k < n_ev; ++k) {
    for (int v = node_ptr[k]; v < node_ptr[k + 1]; ++v) born[v] = k;

    // Ties whose latest shared event is more than `window` before t_k.
    while (!expiries.empty() && expiries.top().first < t[k]) {
      const Expiry top = expiries.top();
      expiries.pop();
      auto it = tie.find(top.second);
      if (it == tie.end() || it->second != top.first) continue;  // renewed since
      const int lo = top.second / n_nodes, hi = top.second % n_nodes;
      --x[lo];
      --x[hi];
      touch(lo);
      touch(hi);
      tie.erase(it);
    }
    record(k);

    std::vector<int> perps;
    for (int e = edge_ptr[k]; e < edge_ptr[k + 1]; ++e) {
      if (edge_i[e] == focal[k]) perps.push_back(edge_j[e]);
      else if (edge_j[e] == focal[k]) perps.push_back(edge_i[e]);
    }
    for (size_t a = 0; a < perps.size(); ++a) {
      if (born[perps[a]] < k) {
        obs_idx.push_back(k);
        obs_node.push_back(perps[a]);
      }
    }

    // After the event: attendance, recency and new or renewed ties.
    for (size_t a = 0; a < perps.size(); ++a) {
      const int v = perps[a];
      ++deg[v];
      last[v] = t[k];
      touch(v);
      for (size_t b = a + 1; b < perps.size(); ++b) {
        const int u = perps[b];
        const long long key = (long long)std::min(u, v) * n_nodes + std::max(u, v);
        const double expiry = t[k] + window;
        if (tie.find(key) == tie.end()) {
          ++x[u];
          ++x[v];
          touch(u);
        }
        tie[key] = expiry;
        if (finite_window) expiries.push(Expiry(expiry, key));
      }
    }
    if (k + 1 < n_ev) record(k + 1);
    else changed.clear();
  }

  auto ptr = [&](const std::vector<int>& idx) {
    IntegerVector p(n_ev + 1, 0);
    for (size_t a = 0; a < idx.size(); ++a) ++p[idx[a] + 1];
    for (int k = 0; k < n_ev; ++k) p[k + 1] += p[k];
    return p;
  };

  return List::create(
    _["upd_ptr"] = ptr(upd_idx), _["upd_node"] = wrap(upd_node),
    _["upd_deg"] = wrap(upd_deg), _["upd_last"] = wrap(upd_last),
    _["upd_x"] = wrap(upd_x),
    _["obs_ptr"] = ptr(obs_idx), _["obs_node"] = wrap(obs_node)
  );
}

// Attachment log-likelihood and its gradient with respect to
// (m, beta_edges, delta, theta). upd_feat is the co-offender feature for each
// snapshot (already transformed); use_degree drops the (d + delta) factor.
// [[Rcpp::export]]
NumericVector attach_ll_cpp(NumericVector t, IntegerVector upd_ptr,
                            IntegerVector upd_node, IntegerVector upd_deg,
                            NumericVector upd_last, NumericVector upd_feat,
                            IntegerVector obs_ptr, IntegerVector obs_node,
                            int n_nodes, double m, double beta, double delta,
                            double theta, bool use_degree) {
  const int n_ev = t.size();
  const NumericVector fail = NumericVector::create(R_NegInf, NA_REAL, NA_REAL, NA_REAL, NA_REAL);

  std::vector<char> active(n_nodes, 0);
  std::vector<double> c(n_nodes, 0.0), last(n_nodes, 0.0), feat(n_nodes, 0.0), dd(n_nodes, 0.0);
  std::vector<int> active_nodes;

  // W and the sums its derivatives need, all at time tc:
  //   W  = sum c E,  Wd = sum c / (d + delta) E,  Wt = sum c x E,  A = sum c a E
  // with E = exp(-beta a) and a = tc - last.
  double W = 0.0, Wd = 0.0, Wt = 0.0, A = 0.0, tc = t.size() > 0 ? t[0] : 0.0;
  auto contrib = [&](int v, double sign) {
    const double a = tc - last[v];
    const double cE = c[v] * std::exp(-beta * a);
    W += sign * cE;
    Wd += sign * (use_degree ? cE / dd[v] : 0.0);
    Wt += sign * cE * feat[v];
    A += sign * cE * a;
  };
  auto recompute = [&]() {
    W = Wd = Wt = A = 0.0;
    for (size_t a = 0; a < active_nodes.size(); ++a) contrib(active_nodes[a], 1.0);
  };

  double ll = 0.0, g_m = 0.0, g_beta = 0.0, g_delta = 0.0, g_theta = 0.0;
  for (int k = 0; k < n_ev; ++k) {
    // Decay to t_k; A follows the same recursion as the Hawkes R term.
    const double dt = t[k] - tc;
    const double e = std::exp(-beta * dt);
    A = e * (A + dt * W);
    W *= e;
    Wd *= e;
    Wt *= e;
    tc = t[k];

    for (int u = upd_ptr[k]; u < upd_ptr[k + 1]; ++u) {
      const int v = upd_node[u];
      if (active[v]) {
        contrib(v, -1.0);
      } else {
        active[v] = 1;
        active_nodes.push_back(v);
      }
      dd[v] = upd_deg[u] + delta;
      last[v] = upd_last[u];
      feat[v] = upd_feat[u];
      c[v] = (use_degree ? dd[v] : 1.0) * std::exp(theta * feat[v]);
      contrib(v, 1.0);
    }
    // Bound the drift of the running sums from repeated add/remove.
    if (k % 512 == 511) recompute();

    if (active_nodes.empty()) {
      if (obs_ptr[k + 1] > obs_ptr[k]) return fail;
      continue;
    }
    if (!std::isfinite(W) || W <= 0.0) return fail;

    const double Gd = Wd / W, Gt = Wt / W, Gb = -A / W;
    ll -= m;
    g_m -= 1.0;
    for (int o = obs_ptr[k]; o < obs_ptr[k + 1]; ++o) {
      const int v = obs_node[o];
      const double a = tc - last[v];
      const double q = m * c[v] * std::exp(-beta * a) / W;
      if (!(q > 0.0) || !std::isfinite(q)) return fail;
      ll += std::log(std::expm1(q));
      const double rq = q / (-std::expm1(-q));  // q * dlog(expm1(q)) / dq
      g_m += rq / m;
      g_beta += rq * (-a - Gb);
      if (use_degree) g_delta += rq * (1.0 / dd[v] - Gd);
      g_theta += rq * (feat[v] - Gt);
    }
  }
  return NumericVector::create(ll, g_m, g_beta, g_delta, g_theta);
}
