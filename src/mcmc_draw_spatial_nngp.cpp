#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector mcmc_draw_spatial_nngp_cpp(
  NumericVector spatial_effect,
  double sigma2,
  List neighbors,
  List neighbors_inverse,
  List pos_in_neighbors,
  List B_s,
  NumericVector F_s,
  NumericVector covariate_sums,
  NumericVector resid_sums
) {
  const int N = spatial_effect.size();
  RNGScope scope;

  for (int s0 = 0; s0 < N; ++s0) {
    const int s = s0 + 1;

    double sum_B_F_inv_B = 0.0;
    double sum_B_F_inv_a = 0.0;

    IntegerVector inv_tt = neighbors_inverse[s0];
    for (int j = 0; j < inv_tt.size(); ++j) {
      const int tt = inv_tt[j];
      const int tt0 = tt - 1;

      NumericVector B_tt = B_s[tt0];
      IntegerVector pos_vec = pos_in_neighbors[tt0];

      const int pos_s = pos_vec[s0];
      const double B_tts = B_tt[pos_s - 1];
      const double F_tt  = F_s[tt0];

      double a_tt_s = spatial_effect[tt0];

      IntegerVector nb_tt = neighbors[tt0];
      for (int k = 0; k < nb_tt.size(); ++k) {
        const int l = nb_tt[k];
        if (l == s) continue;
        const int l0 = l - 1;
        const int pos_l = pos_vec[l0];
        a_tt_s -= B_tt[pos_l - 1] * spatial_effect[l0];
      }

      sum_B_F_inv_B += (B_tts * B_tts) / F_tt;
      sum_B_F_inv_a += B_tts * (1.0 / F_tt) * a_tt_s;
    }

    const double V_s =
      1.0 / ((1.0 / sigma2) * covariate_sums[s0] +
             (1.0 / F_s[s0]) +
             sum_B_F_inv_B);

    double mu_s = resid_sums[s0] + sum_B_F_inv_a;

    IntegerVector nb_s = neighbors[s0];
    if (nb_s.size() > 0) {
      NumericVector B_ss = B_s[s0];
      double dot = 0.0;
      for (int k = 0; k < nb_s.size(); ++k) {
        dot += B_ss[k] * spatial_effect[nb_s[k] - 1];
      }
      mu_s += (1.0 / F_s[s0]) * dot;
    }

    spatial_effect[s0] = R::rnorm(V_s * mu_s, std::sqrt(V_s));
  }

  return spatial_effect;
}
