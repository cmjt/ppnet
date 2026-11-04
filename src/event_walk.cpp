// event_walk.cpp
//
// Array-based walks through an event history. Nodes and edges are 0-based
// integer indices, grouped by event in time order: event k owns nodes
// node_ptr[k] .. node_ptr[k + 1] - 1 and edges edge_ptr[k] .. edge_ptr[k + 1] - 1.
// Only per-node state is held, so memory is O(nodes + edges) however long the
// history is.

#include <Rcpp.h>
#include <cmath>
#include <vector>
using namespace Rcpp;

// Attachment mark log-likelihood (BA and BA-bip).
//
// At event k every candidate node born before k gets weight
//   w_v = (deg_v + delta) * exp(-beta_edges * (t_k - birth_v))
// and is attached independently with p_v = w_v / sum(w), clamped to
// [eps, 1 - eps]. deg_v counts each incident edge twice, matching
// sna::degree() on an undirected network, which the network-based
// evaluators use.
//
// A candidate is observed as attached at event k if it shares an edge of
// event k with focal[k], or with any node when focal[k] < 0 (BA).
//
// Returns c(loglik, d loglik / d beta_edges). With a_v the age and
// abar = sum(w * a) / sum(w), d log p_v / d beta_edges = abar - a_v; clamped
// probabilities contribute no gradient.
// [[Rcpp::export]]
NumericVector mark_ll_attach_cpp(NumericVector t, IntegerVector node_ptr,
                          LogicalVector node_cand, NumericVector node_birth,
                          IntegerVector edge_ptr, IntegerVector edge_i,
                          IntegerVector edge_j, IntegerVector focal,
                          double beta_edges, double delta, double eps) {
  const int n_ev = t.size();
  const int n_nodes = node_birth.size();

  std::vector<int> deg(n_nodes, 0), stamp(n_nodes, -1), cand;
  std::vector<double> w;
  cand.reserve(n_nodes);
  w.reserve(n_nodes);

  double ll = 0.0, grad = 0.0;
  const NumericVector neg_inf = NumericVector::create(R_NegInf, NA_REAL);
  for (int k = 0; k < n_ev; ++k) {
    for (int e = edge_ptr[k]; e < edge_ptr[k + 1]; ++e) {
      const int i = edge_i[e], j = edge_j[e];
      if (focal[k] < 0) {
        stamp[i] = k;
        stamp[j] = k;
      } else if (i == focal[k]) {
        stamp[j] = k;
      } else if (j == focal[k]) {
        stamp[i] = k;
      }
    }

    const int n_cand = cand.size();
    if (n_cand > 0) {
      w.resize(n_cand);
      double W = 0.0, WA = 0.0;
      for (int c = 0; c < n_cand; ++c) {
        const int v = cand[c];
        const double age = t[k] - node_birth[v];
        w[c] = (2.0 * deg[v] + delta) * std::exp(-beta_edges * age);
        W += w[c];
        WA += w[c] * age;
      }
      if (!std::isfinite(W) || W <= 0.0) return neg_inf;
      const double abar = WA / W;

      for (int c = 0; c < n_cand; ++c) {
        const int v = cand[c];
        const double p_raw = w[c] / W;
        const double p = std::min(std::max(p_raw, eps), 1.0 - eps);
        const bool clamped = p != p_raw;
        const double dlogp = abar - (t[k] - node_birth[v]);
        if (stamp[v] == k) {
          ll += std::log(p);
          if (!clamped) grad += dlogp;
        } else {
          ll += std::log1p(-p);
          if (!clamped) grad -= p * dlogp / (1.0 - p);
        }
      }
    }

    // Nodes arrive before edges, as in net_add_event().
    for (int v = node_ptr[k]; v < node_ptr[k + 1]; ++v) {
      if (node_cand[v]) cand.push_back(v);
    }
    for (int e = edge_ptr[k]; e < edge_ptr[k + 1]; ++e) {
      ++deg[edge_i[e]];
      ++deg[edge_j[e]];
    }
  }
  return NumericVector::create(ll, grad);
}

// Walk a bipartite event-perp history: the event's perps are the neighbours
// of focal[k] (its event node) among edges of event k. For each event, links
// holds the distinct earlier events, more than `min_gap` and no more than
// `window` before it, that share at least one perp with it (in the order
// first met), and n_perps / n_old count its perps and those born earlier.
static void link_walk(NumericVector t, IntegerVector node_ptr, IntegerVector edge_ptr,
                      IntegerVector edge_i, IntegerVector edge_j, IntegerVector focal,
                      int n_nodes, double window, double min_gap,
                      std::vector<std::vector<int> >& links,
                      std::vector<int>& n_perps, std::vector<int>& n_old) {
  const int n_ev = focal.size();
  std::vector<int> born(n_nodes, -1), seen(n_ev, -1);
  std::vector<std::vector<int> > history(n_nodes);
  links.assign(n_ev, std::vector<int>());
  n_perps.assign(n_ev, 0);
  n_old.assign(n_ev, 0);

  for (int k = 0; k < n_ev; ++k) {
    for (int v = node_ptr[k]; v < node_ptr[k + 1]; ++v) born[v] = k;

    std::vector<int> perps;
    for (int e = edge_ptr[k]; e < edge_ptr[k + 1]; ++e) {
      if (edge_i[e] == focal[k]) perps.push_back(edge_j[e]);
      else if (edge_j[e] == focal[k]) perps.push_back(edge_i[e]);
    }

    for (size_t a = 0; a < perps.size(); ++a) {
      const int v = perps[a];
      if (born[v] < k) ++n_old[k];
      const std::vector<int>& h = history[v];
      for (size_t b = 0; b < h.size(); ++b) {
        const double gap = t[k] - t[h[b]];
        if (seen[h[b]] != k && gap <= window && gap > min_gap) {
          seen[h[b]] = k;
          links[k].push_back(h[b]);
        }
      }
    }
    for (size_t a = 0; a < perps.size(); ++a) history[perps[a]].push_back(k);
    n_perps[k] = perps.size();
  }
}

// Per-event summaries for a bipartite event-perp history; connectedness is
// the number of linked earlier events (see link_walk()).
// [[Rcpp::export]]
DataFrame event_summary_cpp(NumericVector t, IntegerVector node_ptr, IntegerVector edge_ptr,
                            IntegerVector edge_i, IntegerVector edge_j,
                            IntegerVector focal, int n_nodes, double window,
                            double min_gap = 0.0) {
  const int n_ev = focal.size();
  std::vector<std::vector<int> > links;
  std::vector<int> np, no;
  link_walk(t, node_ptr, edge_ptr, edge_i, edge_j, focal, n_nodes, window, min_gap,
            links, np, no);

  IntegerVector n_perps(n_ev), n_new(n_ev), n_old(n_ev), connectedness(n_ev);
  for (int k = 0; k < n_ev; ++k) {
    n_perps[k] = np[k];
    n_old[k] = no[k];
    n_new[k] = np[k] - no[k];
    connectedness[k] = links[k].size();
  }
  return DataFrame::create(_["n_perps"] = n_perps, _["n_new"] = n_new,
                           _["n_old"] = n_old, _["connectedness"] = connectedness);
}

// Connectedness updates for path-dependent excitation. Event k links to the
// earlier events in links[k]; each of those gains one to its connectedness
// at t_k. Returns, 0-based, upd_ptr (event k's updates are
// upd_ptr[k] .. upd_ptr[k + 1] - 1), upd_j (the event updated) and upd_c (its
// connectedness just before the update), with c0 the connectedness at birth.
// [[Rcpp::export]]
List connectedness_updates_cpp(NumericVector t, IntegerVector node_ptr, IntegerVector edge_ptr,
                               IntegerVector edge_i, IntegerVector edge_j,
                               IntegerVector focal, int n_nodes, double window,
                               double min_gap = 0.0) {
  const int n_ev = focal.size();
  std::vector<std::vector<int> > links;
  std::vector<int> np, no;
  link_walk(t, node_ptr, edge_ptr, edge_i, edge_j, focal, n_nodes, window, min_gap,
            links, np, no);

  IntegerVector c0(n_ev), upd_ptr(n_ev + 1);
  std::vector<int> cnt(n_ev, 0), uj, uc;
  for (int k = 0; k < n_ev; ++k) {
    upd_ptr[k] = uj.size();
    for (size_t a = 0; a < links[k].size(); ++a) {
      const int j = links[k][a];
      uj.push_back(j);
      uc.push_back(cnt[j]);
      ++cnt[j];
    }
    c0[k] = links[k].size();
    cnt[k] = c0[k];
  }
  upd_ptr[n_ev] = uj.size();
  return List::create(_["c0"] = c0, _["upd_ptr"] = upd_ptr,
                      _["upd_j"] = IntegerVector(uj.begin(), uj.end()),
                      _["upd_c"] = IntegerVector(uc.begin(), uc.end()));
}
