// attach_walk.cpp
//
// Attachment of old perps to events in a bipartite event-perp history. At
// event k each old perp v has weight
//   w_v = c_v * exp(-beta * (t_k - last_v)),   c_v = (d_v + delta) * exp(theta' s_v)
// (d_v: events attended, last_v: time of the latest one, s_v: change
// statistics from v's own history, e.g. co-offender ties or 3-paths) and joins independently with p_v = 1 - exp(-q_v), q_v = m w_v / W,
// W = sum of w over all old perps. Because sum_v q_v = m, the event's log-
// likelihood is
//   sum_{v joined} log(expm1(q_v)) - m,
// so only W needs every perp, and W is kept as a running sum that decays
// between events and changes only when a perp's state changes.
//
// With random_effect_var > 0 each perp has its own u_v ~ Gamma(1 / random_effect_var,
// 1 / random_effect_var) and joins with p_v = 1 - exp(-u_v q_v), where W leaves
// u out. u_v is integrated out perp by perp (see the end of attach_ll_cpp).

#include <Rcpp.h>
#include <cmath>
#include <functional>
#include <queue>
#include <unordered_map>
#include <vector>
using namespace Rcpp;

// Parameter-free pass: the state changes of every perp, and which old perps
// joined each event. State snapshots (node, d, last, x, y) at index k are
// applied, in order, before event k is scored. Over v's events no more than
// `window` before, x counts distinct co-offenders and y counts 3-paths
// k - v - j - u completed by v joining event k, i.e. the 2-paths v - j - u:
// co-offenders with multiplicity, sum over v's events j of (n_j - 1).
//
// With cycles = true it also lists, for each event k, the old perps v whose
// 4-cycle change statistic is positive: z_vk = number of 4-cycles
// k - v - j - u - k closed by v joining k, one for each other old perp u at
// k and event j within `window` that u and v both attended. These are the
// past co-offenders of the old perps at k (cyc_node, cyc_z, and cyc_y = 1 if
// v joined); every other old perp has z = 0.
// [[Rcpp::export]]
List attach_history_cpp(NumericVector t, IntegerVector node_ptr,
                        IntegerVector edge_ptr, IntegerVector edge_i,
                        IntegerVector edge_j, IntegerVector focal,
                        int n_nodes, double window, bool cycles = false) {
  const int n_ev = t.size();
  std::vector<int> born(n_nodes, -1), deg(n_nodes, 0), x(n_nodes, 0), y(n_nodes, 0),
    stamp(n_nodes, -1);
  std::vector<double> last(n_nodes, 0.0);

  std::vector<int> upd_idx, upd_node, upd_deg, upd_x, upd_y, obs_idx, obs_node;
  // 4-cycles: per perp, its partners and the times of their shared events.
  std::vector<std::unordered_map<int, std::vector<double> > > shared(cycles ? n_nodes : 0);
  std::vector<int> cyc_idx, cyc_node, cyc_z, cyc_y, zc(cycles ? n_nodes : 0, 0),
    joined_now(cycles ? n_nodes : 0, 0);
  std::vector<int> zc_list;
  std::vector<double> upd_last;

  // Active co-offender ties, keyed by node pair, holding their expiry time.
  std::unordered_map<long long, double> tie;
  typedef std::pair<double, long long> Expiry;
  std::priority_queue<Expiry, std::vector<Expiry>, std::greater<Expiry> > expiries;
  const bool finite_window = std::isfinite(window);
  // Events still inside the window, for expiring their 3-paths.
  std::vector<std::vector<int> > ev_perps(finite_window ? n_ev : 0);
  std::queue<int> ev_queue;

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
      upd_y.push_back(y[v]);
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
    // Events more than `window` before t_k no longer give 3-paths.
    while (!ev_queue.empty() && t[ev_queue.front()] + window < t[k]) {
      const std::vector<int>& pj = ev_perps[ev_queue.front()];
      for (size_t a = 0; a < pj.size(); ++a) {
        y[pj[a]] -= (int)pj.size() - 1;
        touch(pj[a]);
      }
      ev_perps[ev_queue.front()].clear();
      ev_queue.pop();
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

    if (cycles) {
      // Shared events of each old joiner u with each partner v within the
      // window, dropping older ones (and empty partners) as they are seen.
      const double cut = t[k] - window;
      for (size_t a = 0; a < perps.size(); ++a) {
        if (born[perps[a]] < k) joined_now[perps[a]] = 1;
      }
      for (size_t a = 0; a < perps.size(); ++a) {
        const int u = perps[a];
        if (born[u] >= k) continue;
        for (auto it = shared[u].begin(); it != shared[u].end();) {
          std::vector<double>& tt = it->second;
          size_t drop = 0;
          while (drop < tt.size() && tt[drop] < cut) ++drop;
          if (drop > 0) tt.erase(tt.begin(), tt.begin() + drop);
          if (tt.empty()) {
            it = shared[u].erase(it);
            continue;
          }
          const int v = it->first;
          if (zc[v] == 0) zc_list.push_back(v);
          zc[v] += (int)tt.size();
          ++it;
        }
      }
      for (size_t a = 0; a < zc_list.size(); ++a) {
        const int v = zc_list[a];
        cyc_idx.push_back(k);
        cyc_node.push_back(v);
        cyc_z.push_back(zc[v]);
        cyc_y.push_back(joined_now[v]);
        zc[v] = 0;
      }
      zc_list.clear();
      for (size_t a = 0; a < perps.size(); ++a) joined_now[perps[a]] = 0;
      for (size_t a = 0; a < perps.size(); ++a) {
        for (size_t b = a + 1; b < perps.size(); ++b) {
          shared[perps[a]][perps[b]].push_back(t[k]);
          shared[perps[b]][perps[a]].push_back(t[k]);
        }
      }
    }

    // After the event: attendance, recency and new or renewed ties.
    for (size_t a = 0; a < perps.size(); ++a) {
      const int v = perps[a];
      ++deg[v];
      y[v] += (int)perps.size() - 1;
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
    if (finite_window && perps.size() > 1) {
      ev_perps[k] = perps;
      ev_queue.push(k);
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
    _["upd_x"] = wrap(upd_x), _["upd_y"] = wrap(upd_y),
    _["obs_ptr"] = ptr(obs_idx), _["obs_node"] = wrap(obs_node),
    _["cyc_ptr"] = ptr(cyc_idx), _["cyc_node"] = wrap(cyc_node),
    _["cyc_z"] = wrap(cyc_z), _["cyc_y"] = wrap(cyc_y)
  );
}

// Attachment log-likelihood and its gradient with respect to
// (m, beta_edges, delta, theta_1..theta_P, random_effect_var). upd_feat holds
// the P change statistics for each snapshot (one column each, already
// transformed); use_degree drops the
// (d + delta) factor. gl_x and gl_w are a Gauss-Legendre rule on [-1, 1],
// used panel by panel for the random-effect integrals.
//
// 4-cycles (theta_c, from the cyc_* lists of attach_history_cpp, cyc_z
// already transformed) make joins at an event depend on each other. Each
// old perp v is then scored by its conditional given the others at event k,
//   logit P(v joins | rest) = log(expm1(q_v)) + theta_c z_vk,
// the full conditionals of a per-event autologistic model, so this is a
// pseudo-likelihood; at theta_c = 0 it is the likelihood above. Perps with
// z_vk = 0 keep their independent terms, so only the cyc_* entries change.
// With per_event, the gradient of each event's contribution is returned as
// a row of `scores`, for sandwich standard errors.
//
// Returns value = (ll, d/dm, d/dbeta, d/ddelta, d/dtheta_1..P,
// d/drandom_effect_var, d/dtheta_c).
// [[Rcpp::export]]
List attach_ll_cpp(NumericVector t, IntegerVector upd_ptr,
                            IntegerVector upd_node, IntegerVector upd_deg,
                            NumericVector upd_last, NumericMatrix upd_feat,
                            IntegerVector obs_ptr, IntegerVector obs_node,
                            int n_nodes, double m, double beta, double delta,
                            NumericVector theta, bool use_degree, double random_effect_var,
                            NumericVector gl_x, NumericVector gl_w,
                            IntegerVector cyc_ptr, IntegerVector cyc_node,
                            NumericVector cyc_z, IntegerVector cyc_y, double theta_c,
                            bool per_event) {
  const int n_ev = t.size();
  const int P = theta.size();
  const int NG = 5 + P;  // gradient entries
  if (upd_feat.ncol() != P) stop("upd_feat needs one column per theta.");
  const bool has_cyc = cyc_ptr.size() == n_ev + 1;
  NumericVector fail_v(NG + 1, NA_REAL);
  fail_v[0] = R_NegInf;
  const List fail = List::create(_["value"] = fail_v, _["scores"] = R_NilValue);
  // Below 1e-8 the random-effect term is numerically the model without it.
  const bool has_re = random_effect_var > 1e-8;
  if (has_re && has_cyc) stop("4-cycles cannot be combined with the random effect.");
  if (has_re && per_event) stop("Per-event scores are not available with the random effect.");

  std::vector<char> active(n_nodes, 0);
  // feat[P * v + j]: perp v's statistic j.
  std::vector<double> c(n_nodes, 0.0), last(n_nodes, 0.0), feat((size_t)P * n_nodes, 0.0),
    dd(n_nodes, 0.0);
  std::vector<int> active_nodes;

  // W and the sums its derivatives need, all at time tc:
  //   W  = sum c E,  Wd = sum c / (d + delta) E,  Wt_j = sum c s_j E,  A = sum c a E
  // with E = exp(-beta a) and a = tc - last.
  double W = 0.0, Wd = 0.0, A = 0.0, tc = t.size() > 0 ? t[0] : 0.0;
  std::vector<double> Wt(P, 0.0);
  auto contrib = [&](int v, double sign) {
    const double a = tc - last[v];
    const double cE = c[v] * std::exp(-beta * a);
    W += sign * cE;
    Wd += sign * (use_degree ? cE / dd[v] : 0.0);
    for (int j = 0; j < P; ++j) Wt[j] += sign * cE * feat[(size_t)P * v + j];
    A += sign * cE * a;
  };
  auto recompute = [&]() {
    W = Wd = A = 0.0;
    std::fill(Wt.begin(), Wt.end(), 0.0);
    for (size_t a = 0; a < active_nodes.size(); ++a) contrib(active_nodes[a], 1.0);
  };

  // Random effect: per-event quantities r_k = exp(-beta (t_k - t_ref)) / W_k and
  // their products with the gradient terms, within epochs over which t_ref
  // is fixed. Over a stretch [s, e) of unchanged state, perp v's
  //   sum_k q_vk = m c_v exp(-beta (t_ref - last_v)) sum_{s <= k < e} r_k,
  // taken from suffix sums of r once the walk is done. r_k tends to fall
  // over an epoch, so suffix sums keep late stretches accurate where prefix
  // sums would lose them to cancellation.
  // Suffix-sum columns: R0, RD, RB, RA, then one per statistic.
  enum { R0, RD, RB, RA, RT };
  const int NR = RT + P;
  std::vector<double> suf, ep_tref;
  std::vector<int> ep_of, ep_begin;
  // Random effect: each join's perp, q and d log q / d (beta, delta, theta);
  // ob_dt holds P values per join.
  std::vector<int> ob_v;
  std::vector<double> ob_q, ob_db, ob_dd, ob_dt;
  std::vector<int> seg_v, seg_s, seg_e, seg_start;
  std::vector<double> seg_c, seg_last, seg_feat, seg_dd;
  if (has_re) {
    suf.assign((size_t)NR * n_ev, 0.0);
    ep_of.assign(n_ev, 0);
    seg_start.assign(n_nodes, 0);
  }
  auto close_segment = [&](int v, int s, int e) {
    if (e <= s) return;
    seg_v.push_back(v);
    seg_s.push_back(s);
    seg_e.push_back(e);
    seg_c.push_back(c[v]);
    seg_last.push_back(last[v]);
    for (int j = 0; j < P; ++j) seg_feat.push_back(feat[(size_t)P * v + j]);
    seg_dd.push_back(dd[v]);
  };

  double ll = 0.0, g_m = 0.0, g_beta = 0.0, g_delta = 0.0, g_fv = 0.0, g_c = 0.0;
  std::vector<double> g_theta(P, 0.0);
  // Per-event scores: the change in each gradient entry over event k.
  NumericMatrix scores(per_event ? n_ev : 0, per_event ? NG : 0);
  std::vector<double> gprev(NG, 0.0);
  auto grad_now = [&](std::vector<double>& g) {
    g[0] = g_m; g[1] = g_beta; g[2] = g_delta;
    for (int j = 0; j < P; ++j) g[3 + j] = g_theta[j];
    g[3 + P] = g_fv; g[4 + P] = g_c;
  };
  std::vector<double> gcur(NG, 0.0);
  auto log1pexp = [](double x) { return x > 0 ? x + std::log1p(std::exp(-x)) : std::log1p(std::exp(x)); };
  for (int k = 0; k < n_ev; ++k) {
    // Decay to t_k; A follows the same recursion as the Hawkes R term.
    const double dt = t[k] - tc;
    const double e = std::exp(-beta * dt);
    A = e * (A + dt * W);
    W *= e;
    Wd *= e;
    for (int j = 0; j < P; ++j) Wt[j] *= e;
    tc = t[k];

    for (int u = upd_ptr[k]; u < upd_ptr[k + 1]; ++u) {
      const int v = upd_node[u];
      if (has_re && active[v]) close_segment(v, seg_start[v], k);
      if (has_re) seg_start[v] = k;
      if (active[v]) {
        contrib(v, -1.0);
      } else {
        active[v] = 1;
        active_nodes.push_back(v);
      }
      dd[v] = upd_deg[u] + delta;
      last[v] = upd_last[u];
      double eta = 0.0;
      for (int j = 0; j < P; ++j) {
        feat[(size_t)P * v + j] = upd_feat(u, j);
        eta += theta[j] * upd_feat(u, j);
      }
      c[v] = (use_degree ? dd[v] : 1.0) * std::exp(eta);
      contrib(v, 1.0);
    }
    // Bound the drift of the running sums from repeated add/remove.
    if (k % 512 == 511) recompute();

    if (active_nodes.empty()) {
      if (obs_ptr[k + 1] > obs_ptr[k]) return fail;
      continue;
    }
    if (!std::isfinite(W) || W <= 0.0) return fail;

    const double Gd = Wd / W, Gb = -A / W;
    std::vector<double> Gt(P);
    for (int j = 0; j < P; ++j) Gt[j] = Wt[j] / W;
    if (has_re) {
      // The -sum_v q_vk = -m term moves into each perp's marginal below.
      if (ep_tref.empty() || beta * (tc - ep_tref.back()) > 300.0) {
        ep_tref.push_back(tc);
        ep_begin.push_back(k);
      }
      const int E = (int)ep_tref.size() - 1;
      ep_of[k] = E;
      const double r = std::exp(-beta * (tc - ep_tref[E]) - std::log(W));
      double* sk = &suf[(size_t)NR * k];
      sk[R0] = r;
      sk[RD] = r * Gd;
      sk[RB] = r * Gb;
      sk[RA] = r * (tc - ep_tref[E]);
      for (int j = 0; j < P; ++j) sk[RT + j] = r * Gt[j];
    } else {
      ll -= m;
      g_m -= 1.0;
    }
    for (int o = obs_ptr[k]; o < obs_ptr[k + 1]; ++o) {
      const int v = obs_node[o];
      const double a = tc - last[v];
      const double q = m * c[v] * std::exp(-beta * a) / W;
      if (!(q > 0.0) || !std::isfinite(q)) return fail;
      if (has_re) {
        ob_v.push_back(v);
        ob_q.push_back(q);
        ob_db.push_back(-a - Gb);
        ob_dd.push_back(use_degree ? 1.0 / dd[v] - Gd : 0.0);
        for (int j = 0; j < P; ++j) ob_dt.push_back(feat[(size_t)P * v + j] - Gt[j]);
        continue;
      }
      ll += std::log(std::expm1(q));
      const double rq = q / (-std::expm1(-q));  // q * dlog(expm1(q)) / dq
      g_m += rq / m;
      g_beta += rq * (-a - Gb);
      if (use_degree) g_delta += rq * (1.0 / dd[v] - Gd);
      for (int j = 0; j < P; ++j) g_theta[j] += rq * (feat[(size_t)P * v + j] - Gt[j]);
    }

    // 4-cycles: swap the independent terms of perps with z > 0 for their
    // conditionals given the rest of the event.
    if (has_cyc) {
      for (int o = cyc_ptr[k]; o < cyc_ptr[k + 1]; ++o) {
        const int v = cyc_node[o];
        const bool y = cyc_y[o] == 1;
        const double a = tc - last[v];
        const double q = m * c[v] * std::exp(-beta * a) / W;
        if (!(q > 0.0) || !std::isfinite(q)) return fail;
        const double eta = std::log(std::expm1(q)), x = eta + theta_c * cyc_z[o];
        const double cond = y ? -log1pexp(-x) : -log1pexp(x);
        const double indep = y ? eta - q : -q;
        ll += cond - indep;
        const double sig = 1.0 / (1.0 + std::exp(-x));
        g_c += ((y ? 1.0 : 0.0) - sig) * cyc_z[o];
        // d(cond - indep) / d log q
        const double dq = q * (((y ? 1.0 : 0.0) - sig) / (-std::expm1(-q)) -
                               (y ? 1.0 / std::expm1(q) : -1.0));
        g_m += dq / m;
        g_beta += dq * (-a - Gb);
        if (use_degree) g_delta += dq * (1.0 / dd[v] - Gd);
        for (int j = 0; j < P; ++j) g_theta[j] += dq * (feat[(size_t)P * v + j] - Gt[j]);
      }
    }
    if (per_event) {
      grad_now(gcur);
      for (int j = 0; j < NG; ++j) {
        scores(k, j) = gcur[j] - gprev[j];
        gprev[j] = gcur[j];
      }
    }
  }

  if (has_re) {
    // Q_v = sum of q_vk over every event where v is a candidate, and its
    // derivatives, from the suffix sums.
    for (size_t i = 0; i < active_nodes.size(); ++i) {
      close_segment(active_nodes[i], seg_start[active_nodes[i]], n_ev);
    }
    const int n_ep = ep_begin.size();
    auto ep_end = [&](int E) { return E + 1 < n_ep ? ep_begin[E + 1] : n_ev; };
    for (int E = 0; E < n_ep; ++E) {
      for (int k = ep_end(E) - 2; k >= ep_begin[E]; --k) {
        for (int j = 0; j < NR; ++j) suf[(size_t)NR * k + j] += suf[(size_t)NR * (k + 1) + j];
      }
    }

    std::vector<double> Qv(n_nodes, 0.0), dQ_beta(n_nodes, 0.0), dQ_delta(n_nodes, 0.0),
      dQ_theta((size_t)P * n_nodes, 0.0);
    std::vector<double> S(NR);
    for (size_t i = 0; i < seg_v.size(); ++i) {
      const int v = seg_v[i], s = seg_s[i], e = seg_e[i];
      for (int E = ep_of[s]; E <= ep_of[e - 1]; ++E) {
        const double gap = ep_tref[E] - seg_last[i];
        const double F = m * seg_c[i] * std::exp(-beta * gap);
        if (F == 0.0) break;  // later epochs are smaller still
        const int lo = std::max(s, ep_begin[E]), hi = std::min(e, ep_end(E));
        for (int j = 0; j < NR; ++j) {
          S[j] = suf[(size_t)NR * lo + j] - (hi < ep_end(E) ? suf[(size_t)NR * hi + j] : 0.0);
        }
        const double Q = F * S[R0];
        Qv[v] += Q;
        for (int j = 0; j < P; ++j) {
          dQ_theta[(size_t)P * v + j] += seg_feat[(size_t)P * i + j] * Q - F * S[RT + j];
        }
        if (use_degree) dQ_delta[v] += Q / seg_dd[i] - F * S[RD];
        dQ_beta[v] += -F * S[RA] - gap * Q - F * S[RB];
      }
    }

    // Joins grouped by perp.
    std::vector<int> ob_ptr(n_nodes + 1, 0), ob_ord(ob_v.size());
    for (size_t o = 0; o < ob_v.size(); ++o) ++ob_ptr[ob_v[o] + 1];
    for (int v = 0; v < n_nodes; ++v) ob_ptr[v + 1] += ob_ptr[v];
    {
      std::vector<int> pos(ob_ptr.begin(), ob_ptr.end() - 1);
      for (size_t o = 0; o < ob_v.size(); ++o) ob_ord[pos[ob_v[o]]++] = o;
    }

    // Perp v's likelihood given u is
    //   prod_{joins} (1 - exp(-u q_k)) * exp(-u Q'),  Q' = Q_v - sum_{joins} q_k,
    // integrated against u ~ Gamma(a, a), a = 1 / random_effect_var. Perps with no
    // returns have the closed form a^a / (a + Q')^a. For the others, in
    // s = log u the log integrand
    //   h(s) = sum_joins log(1 - exp(-u q_k)) - u Q' - a (e^s - 1 - s)
    // (plus the constant log_prior) is concave, so its mode is found by Newton, the range where h is within
    // 40 of its maximum is bracketed, and the integral over that range is
    // taken by composite Gauss-Legendre. Derivatives are posterior means of
    // the derivatives of the log integrand.
    // The prior's log density in s, a log a - lgamma(a) + a s - a e^s, is
    // split as log_prior - a (e^s - 1 - s) so that the terms of order a
    // cancel analytically; log_prior = a log a - a - lgamma(a) and
    // lpsi = log a - digamma(a) use Stirling series for large a.
    const double a = 1.0 / random_effect_var;
    double log_prior, lpsi;
    if (a >= 15.0) {
      const double a2 = a * a;
      log_prior = 0.5 * std::log(a / (2.0 * M_PI)) -
        (1.0 / 12.0 - (1.0 / 360.0 - 1.0 / (1260.0 * a2)) / a2) / a;
      lpsi = 0.5 / a + (1.0 / 12.0 - (1.0 / 120.0 - 1.0 / (252.0 * a2)) / a2) / a2;
    } else {
      log_prior = a * std::log(a) - a - std::lgamma(a);
      lpsi = std::log(a) - R::digamma(a);
    }
    const int n_gl = gl_x.size();
    std::vector<double> sv, uv, lp;
    std::vector<int> keep;  // nodes with non-negligible posterior weight
    // phi(x) = x exp(-x) / (1 - exp(-x)) = d log(1 - exp(-x)) / d log x
    auto phi = [](double x) { return x * std::exp(-x) / (-std::expm1(-x)); };

    for (size_t i = 0; i < active_nodes.size(); ++i) {
      const int v = active_nodes[i];
      const int o0 = ob_ptr[v], o1 = ob_ptr[v + 1], n = o1 - o0;
      double Qj = 0.0;
      for (int o = o0; o < o1; ++o) Qj += ob_q[ob_ord[o]];
      const double Qp = std::max(Qv[v] - Qj, 0.0), B = a + Qp;

      // Eu = E[u] and d_a = d/da of log L_v, both under v's posterior.
      double Eu, d_a;
      keep.clear();
      if (n == 0) {
        ll += -a * std::log1p(Qp / a);
        Eu = a / B;
        d_a = Qp / B - std::log1p(Qp / a);
      } else {
        auto h = [&](double sl) {
          const double u = std::exp(sl);
          double r = -a * (std::expm1(sl) - sl) - u * Qp;
          for (int o = o0; o < o1; ++o) r += std::log(-std::expm1(-u * ob_q[ob_ord[o]]));
          return r;
        };
        auto h12 = [&](double sl, double& d1, double& d2) {
          const double u = std::exp(sl);
          d1 = -a * std::expm1(sl) - u * Qp;
          d2 = -u * B;
          for (int o = o0; o < o1; ++o) {
            const double x = u * ob_q[ob_ord[o]], f = phi(x);
            d1 += f;
            d2 += f * (1.0 - x / (-std::expm1(-x)));
          }
        };
        // Newton from the Poisson-hit mode, halving steps that lower h.
        double sl = std::log((a + n) / (a + Qv[v])), hs = h(sl), d1, d2;
        for (int it = 0; it < 200; ++it) {
          h12(sl, d1, d2);
          double step = -d1 / d2;
          if (!std::isfinite(step)) return fail;
          double cand = sl + step, hc = h(cand);
          while (!(hc >= hs) && std::abs(step) > 1e-14) {
            step /= 2.0;
            cand = sl + step;
            hc = h(cand);
          }
          sl = cand;
          hs = hc;
          if (std::abs(step) < 1e-10) break;
        }
        h12(sl, d1, d2);
        if (!(d2 < 0.0) || !std::isfinite(hs)) return fail;
        const double sd = 1.0 / std::sqrt(-d2);
        // Panels on each side of the mode end at sd * (1, 2, 4, ...), out to
        // where h has fallen by 40 (concavity makes the doubling safe), so
        // nodes are dense near the mode and sparse in the tails.
        sv.clear();
        lp.clear();
        double lmax = R_NegInf;
        for (int side = -1; side <= 1; side += 2) {
          double inner = 0.0, outer = sd;
          for (int it = 0; it < 60; ++it) {
            const double width = outer - inner;
            for (int g = 0; g < n_gl; ++g) {
              const double sg = sl + side * (inner + 0.5 * width * (gl_x[g] + 1.0));
              sv.push_back(sg);
              lp.push_back(std::log(0.5 * width * gl_w[g]) + h(sg));
              lmax = std::max(lmax, lp.back());
            }
            if (h(sl + side * outer) < hs - 40.0) break;
            inner = outer;
            outer *= 2.0;
          }
        }
        const int G = sv.size();
        uv.resize(G);
        double tot = 0.0;
        for (int g = 0; g < G; ++g) tot += (lp[g] = std::exp(lp[g] - lmax));
        ll += lmax + std::log(tot);
        // lp now holds the posterior weights.
        Eu = 0.0;
        d_a = lpsi;
        for (int g = 0; g < G; ++g) {
          lp[g] /= tot;
          if (lp[g] < 1e-17) continue;
          uv[g] = std::exp(sv[g]);
          Eu += lp[g] * uv[g];
          d_a -= lp[g] * (std::expm1(sv[g]) - sv[g]);
          keep.push_back(g);
        }
        ll += log_prior;
      }

      // d/dpsi of the log integrand: sum_joins phi(u q) dlog q - u dQ', with
      // dQ' = dQ_v - sum_joins q dlog q.
      double jq_b = 0.0, jq_d = 0.0;
      std::vector<double> jq_t(P, 0.0);
      for (int o = o0; o < o1; ++o) {
        const int oo = ob_ord[o];
        double Ephi = 0.0;
        for (size_t j = 0; j < keep.size(); ++j) Ephi += lp[keep[j]] * phi(uv[keep[j]] * ob_q[oo]);
        g_m += Ephi / m;
        g_beta += Ephi * ob_db[oo];
        g_delta += Ephi * ob_dd[oo];
        jq_b += ob_q[oo] * ob_db[oo];
        jq_d += ob_q[oo] * ob_dd[oo];
        for (int j = 0; j < P; ++j) {
          g_theta[j] += Ephi * ob_dt[(size_t)P * oo + j];
          jq_t[j] += ob_q[oo] * ob_dt[(size_t)P * oo + j];
        }
      }
      g_m -= Eu * Qp / m;
      g_beta -= Eu * (dQ_beta[v] - jq_b);
      g_delta -= Eu * (dQ_delta[v] - jq_d);
      for (int j = 0; j < P; ++j) g_theta[j] -= Eu * (dQ_theta[(size_t)P * v + j] - jq_t[j]);
      g_fv -= a * a * d_a;
    }
    if (!std::isfinite(ll)) return fail;
  }
  NumericVector out(NG + 1);
  out[0] = ll;
  grad_now(gcur);
  for (int j = 0; j < NG; ++j) out[1 + j] = gcur[j];
  return List::create(_["value"] = out, _["scores"] = per_event ? (SEXP)scores : R_NilValue);
}
