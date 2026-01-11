#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List get_F_and_B_cpp(
  Rcpp::List kernels,
  double tau,
  Rcpp::List neighbors,
  Rcpp::List kernels_partial_inverse
) {
  const int S = kernels.size();

  NumericVector F_s(S);
  List B_s(S);

  for (int s = 0; s < S; ++s) {
    NumericMatrix K = kernels[s];

    SEXP nb_sexp = neighbors[s];
    int nb_len = 0;
    if (nb_sexp != R_NilValue) {
      nb_len = Rf_length(nb_sexp);
    }

    if (nb_len > 0) {
      NumericMatrix Kpi = kernels_partial_inverse[s];

      NumericVector B(nb_len);
      for (int j = 0; j < nb_len; ++j) {
        double acc = 0.0;
        for (int k = 0; k < nb_len; ++k) {
          acc += K(0, k + 1) * Kpi(k, j);
        }
        B[j] = acc;
      }
      B_s[s] = B;

      double dot = 0.0;
      for (int k = 0; k < nb_len; ++k) {
        dot += B[k] * K(k + 1, 0);
      }
      F_s[s] = tau * (K(0, 0) - dot);

    } else {
      B_s[s] = NumericVector(0);

      F_s[s] = tau * K(0, 0);
    }
  }

  return List::create(
    Named("B") = B_s,
    Named("F") = F_s
  );
}
