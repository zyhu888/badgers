// [[Rcpp::depends(RcppArmadillo)]]

#define ARMA_WARN_LEVEL 1
#include <RcppArmadillo.h>
#include <queue>
using namespace Rcpp;

// Bit-exact Rcpp port of utility_v2.R (verified with identical() against the
// R implementation, on both macOS and the Linux deployment target).
//
// Must compile with -ffp-contract=off (e.g. Sys.setenv(PKG_CXXFLAGS =
// "-ffp-contract=off") before sourceCpp()). Without it the compiler fuses
// expressions like a*a - b*b into one rounding step, which rounds
// differently from R's separate multiply/subtract steps.
//
// Manual ports of R's sum()/mean()/rowMeans() use `long double` accumulators
// to match R's own internal LDOUBLE accumulation (src/main/summary.c,
// array.c). On ARM64 (long double == double) this makes no difference, but
// on x86_64 long double is genuine 80-bit precision and this is required
// for bit-exactness.
#ifdef __clang__
#pragma clang fp contract(off)
#endif

// MH acceptance-ratio density for the alpha update.
// [[Rcpp::export]]
double p_density_mmh_cpp(double alpha_j_prop, double alpha_j, double tau_j,
                          NumericVector alpha, NumericMatrix phi,
                          int j, double nu, double mu) {
  double lambda_j_prop = std::exp(alpha_j_prop);
  double lambda_j = std::exp(alpha_j);

  int p = alpha.size();
  int j0 = j - 1; // R's 1-based index -> 0-based

  // comp1 = (lambda_j_prop^2)/(lambda_j^2)
  double comp1 = (lambda_j_prop * lambda_j_prop) / (lambda_j * lambda_j);

  // comp2 = exp(-(lambda_j_prop^2 - lambda_j^2) * tau_j^2 / 2)
  double comp2 = std::exp(-(lambda_j_prop * lambda_j_prop - lambda_j * lambda_j) *
                           (tau_j * tau_j) / 2.0);

  // vec = alpha + alpha_aug - 2*mu, where alpha_aug is alpha with element j0
  // replaced by alpha_j_prop; for k != j0, alpha[k]+alpha_aug[k] == 2*alpha[k]
  // exactly (IEEE754 doubling), so alpha_aug is never materialized.
  double two_mu = 2.0 * mu;
  std::vector<double> vec(p);
  for (int k = 0; k < p; ++k) {
    vec[k] = (k == j0) ? (alpha[j0] + alpha_j_prop - two_mu)
                       : (alpha[k] + alpha[k] - two_mu);
  }

  // comp3 = exp((-1/(2*nu)) * (alpha_j_prop - alpha_j) * (phi_j %*% vec))
  // phi_j = phi[j, ] is read directly from phi's column-major storage with
  // stride p (no row copy needed).
  arma::blas_int bn = p, binc_phi = p, binc_vec = 1;
  double quad = arma::ddot_(&bn, phi.begin() + j0, &binc_phi, vec.data(), &binc_vec);
  double comp3 = std::exp((-1.0 / (2.0 * nu)) * (alpha_j_prop - alpha_j) * quad);

  return comp1 * comp2 * comp3;
}

// BLAS/LAPACK-backed helpers matching R's %*%/crossprod/tcrossprod/chol2inv/
// solve() bit-for-bit. Armadillo's own dot()/operator*/dgetri diverge from R
// on small matrices, so these call the underlying BLAS/LAPACK routines
// directly instead of going through Armadillo's higher-level API.
namespace badgers_blas {

// A %*% x
inline arma::vec matvec(const arma::mat &A, const arma::vec &x) {
  arma::vec y(A.n_rows, arma::fill::zeros);
  char trans = 'N';
  arma::blas_int m = A.n_rows, n = A.n_cols, lda = A.n_rows, incx = 1, incy = 1;
  double alpha = 1.0, beta = 0.0;
  if (m == 0 || n == 0) return y;
  arma::dgemv_(&trans, &m, &n, &alpha, const_cast<double*>(A.memptr()), &lda,
               const_cast<double*>(x.memptr()), &incx, &beta, y.memptr(), &incy,
               (arma::blas_len)1);
  return y;
}

// tcrossprod(x) == x %*% t(x)
inline arma::mat outer(const arma::vec &x) {
  arma::blas_int n = x.n_elem, k = 1;
  arma::mat M(n, n, arma::fill::zeros);
  if (n == 0) return M;
  char transa = 'N', transb = 'T';
  double alpha = 1.0, beta = 0.0;
  arma::dgemm_(&transa, &transb, &n, &n, &k, &alpha, const_cast<double*>(x.memptr()), &n,
               const_cast<double*>(x.memptr()), &n, &beta, M.memptr(), &n,
               (arma::blas_len)1, (arma::blas_len)1);
  return M;
}

// sum(x * y); long double accumulator matches R's sum() (see mean_twopass below)
inline double seqsum_prod(const arma::vec &x, const arma::vec &y) {
  long double s = 0.0L;
  for (arma::uword k = 0; k < x.n_elem; ++k) {
    double term = x[k] * y[k];
    s += term;
  }
  return (double) s;
}

// crossprod(X) == t(X) %*% X
inline arma::mat crossprod(const arma::mat &X) {
  arma::blas_int n = X.n_cols, k = X.n_rows, lda = X.n_rows, ldc = X.n_cols;
  arma::mat C(n, n, arma::fill::zeros);
  if (n == 0 || k == 0) return C;
  char uplo = 'U', trans = 'T';
  double alpha = 1.0, beta = 0.0;
  arma::dsyrk_(&uplo, &trans, &n, &k, &alpha, const_cast<double*>(X.memptr()), &lda,
               &beta, C.memptr(), &ldc, (arma::blas_len)1, (arma::blas_len)1);
  for (int a = 0; a < n; ++a)
    for (int b = 0; b < a; ++b)
      C(a, b) = C(b, a); // dsyrk only fills the upper triangle
  return C;
}

// tcrossprod(A, B) == A %*% t(B)
inline arma::mat tcrossprod(const arma::mat &A, const arma::mat &B) {
  arma::blas_int m = A.n_rows, n = B.n_rows, k = A.n_cols;
  arma::blas_int lda = A.n_rows, ldb = B.n_rows, ldc = A.n_rows;
  arma::mat C(m, n, arma::fill::zeros);
  if (m == 0 || n == 0 || k == 0) return C;
  char transa = 'N', transb = 'T';
  double alpha = 1.0, beta = 0.0;
  arma::dgemm_(&transa, &transb, &m, &n, &k, &alpha, const_cast<double*>(A.memptr()), &lda,
               const_cast<double*>(B.memptr()), &ldb, &beta, C.memptr(), &ldc,
               (arma::blas_len)1, (arma::blas_len)1);
  return C;
}

// chol2inv(R): inverse of A = R'R from upper-triangular chol factor R
inline arma::mat chol2inv(const arma::mat &R_chol) {
  arma::mat A = R_chol; // dpotri overwrites in place
  int n = A.n_rows;
  if (n == 0) return A;
  char uplo = 'U';
  arma::blas_int bn = n, lda = n, info = 0;
  arma::lapack::potri(&uplo, &bn, A.memptr(), &lda, &info);
  if (info != 0) stop("LAPACK failed to invert a Cholesky factor (info=%d).", info);
  for (int a = 0; a < n; ++a)
    for (int b = 0; b < a; ++b)
      A(a, b) = A(b, a); // dpotri only fills the upper triangle
  return A;
}

// solve(A): general square-matrix inverse, via dgesv (LU + solve against
// identity) -- NOT dgetrf+dgetri, which diverges from R's solve()
inline arma::mat solve_general(const arma::mat &A_in) {
  arma::mat A = A_in; // dgesv overwrites A with its LU factors
  int n = A.n_rows;
  arma::mat B(n, n, arma::fill::eye);
  if (n == 0) return B;
  arma::blas_int bn = n, nrhs = n, lda = n, ldb = n, info = 0;
  std::vector<arma::blas_int> ipiv(n);
  arma::lapack::gesv(&bn, &nrhs, A.memptr(), &lda, ipiv.data(), B.memptr(), &ldb, &info);
  if (info != 0) stop("LAPACK failed to solve a linear system (info=%d).", info);
  return B;
}

} // namespace badgers_blas

// R-statistics-function ports needed for setup and post-processing.
namespace badgers_stats {

// mean(x): R's two-pass compensated algorithm (sum, then a correction pass
// over the residuals), not sum/n. Long double accumulators match R's
// real_mean() (src/main/summary.c). mean(numeric(0)) is NaN, which 0.0/0
// naturally produces here.
inline double mean_twopass(const std::vector<double> &vals) {
  int m = (int) vals.size();
  long double s = 0.0L;
  for (double v : vals) s += v;
  s /= m;
  if (R_FINITE((double) s)) {
    long double t = 0.0L;
    for (double v : vals) t += (v - s);
    s += t / m;
  }
  return (double) s;
}

// rowMeans(M): plain sum/n per row, no two-pass correction (unlike mean()),
// but still a long double accumulator (R's do_colsum in array.c).
inline double rowmean_naive(const arma::mat &M, int row) {
  long double s = 0.0L;
  for (arma::uword c = 0; c < M.n_cols; ++c) s += M(row, c);
  return (double) (s / M.n_cols);
}

// quantile(x, probs), R's default type=7, ported from stats:::quantile.default:
// index=(n-1)*p+1 (1-based), lo=floor(index), hi=ceiling(index) (not lo+1),
// qs=x[lo], and only when index>lo and x[hi]!=qs: qs=(1-h)*qs+h*x[hi].
inline double quantile7(std::vector<double> x, double p) {
  std::sort(x.begin(), x.end());
  int n = (int) x.size();
  double index = (n - 1) * p + 1.0;
  double lo = std::floor(index);
  double hi = std::ceil(index);
  int lo_i = (int) lo - 1, hi_i = (int) hi - 1;
  double qs = x[lo_i];
  if (index > lo && x[hi_i] != qs) {
    double h = index - lo;
    qs = (1.0 - h) * qs + h * x[hi_i];
  }
  return qs;
}

// All-pairs shortest-path distances for an unweighted undirected graph
// (replaces igraph::distances(), avoiding that dependency).
inline arma::mat bfs_distances(const arma::Mat<int> &adj) {
  int p = adj.n_rows;
  arma::mat dist(p, p, arma::fill::value(arma::datum::inf));
  std::vector<std::vector<int>> nbrs(p);
  for (int a = 0; a < p; ++a)
    for (int b = 0; b < p; ++b)
      if (a != b && adj(a, b) != 0) nbrs[a].push_back(b);

  for (int src = 0; src < p; ++src) {
    std::vector<int> d(p, -1);
    d[src] = 0;
    std::queue<int> q;
    q.push(src);
    while (!q.empty()) {
      int u = q.front(); q.pop();
      for (int v : nbrs[u]) if (d[v] == -1) { d[v] = d[u] + 1; q.push(v); }
    }
    for (int v = 0; v < p; ++v) dist(src, v) = (d[v] == -1) ? arma::datum::inf : (double) d[v];
  }
  return dist;
}

} // namespace badgers_stats

// V2 block Gibbs update for Gamma with edge-specific and diagonal penalties.
//
// `tau1` is dropped from the signature: in R it is scratch space passed in
// as a fresh zeros(p,p) and fully overwritten before use, so its input value
// never affects the output; allocated as zeros internally instead.
// [[Rcpp::export]]
List AdaptiveGLasso_cpp(NumericMatrix S_, double n, int p,
                         NumericMatrix gamma_, NumericMatrix inv_gamma_,
                         double a_iota, NumericVector b_iota_edge, double b_iota_diag,
                         IntegerVector upperind, IntegerVector lowerind,
                         IntegerMatrix ind_noi_all) {
  arma::mat S = as<arma::mat>(S_);
  arma::mat gamma = as<arma::mat>(gamma_);
  arma::mat inv_gamma = as<arma::mat>(inv_gamma_);

  int m_edge = upperind.size(); // length(upperind) == p*(p-1)/2
  double iota_shape_post = a_iota + 1.0;

  // gamma_adjust = pmax(abs(gamma[upperind]), 1e-12), read from the original
  // (pre-update) gamma -- upperind is a 1-based column-major linear index.
  std::vector<double> gamma_adjust(m_edge);
  for (int i = 0; i < m_edge; ++i) {
    double v = std::fabs(gamma.memptr()[upperind[i] - 1]);
    gamma_adjust[i] = std::max(v, 1e-12);
  }

  // Sample one penalty for each undirected off-diagonal pair.
  NumericVector iota_edge(m_edge);
  for (int i = 0; i < m_edge; ++i) {
    double rate = gamma_adjust[i] + b_iota_edge[i];
    iota_edge[i] = R::rgamma(iota_shape_post, 1.0 / rate);
  }

  // Sample one penalty for each diagonal element.
  NumericVector iota_diag(p);
  for (int i = 0; i < p; ++i) {
    double rate = b_iota_diag + gamma(i, i) / 2.0;
    iota_diag[i] = R::rgamma(iota_shape_post, 1.0 / rate);
  }

  // Sample the off-diagonal scale augmentation.
  std::vector<double> mu_prime(m_edge), shape_ig(m_edge);
  for (int i = 0; i < m_edge; ++i) {
    mu_prime[i] = std::min(iota_edge[i] / gamma_adjust[i], 1e12);
    shape_ig[i] = iota_edge[i] * iota_edge[i];
  }

  // rinvgauss fast path (statmod::rinvgauss): all rnorm() first, then all runif()
  std::vector<double> ig_phi(m_edge); // dispersion = 1/shape
  for (int i = 0; i < m_edge; ++i) ig_phi[i] = 1.0 / shape_ig[i];

  std::vector<double> Y(m_edge), X1(m_edge);
  for (int i = 0; i < m_edge; ++i) {
    double z = R::norm_rand();
    Y[i] = z * z;
  }
  for (int i = 0; i < m_edge; ++i) {
    double Yphi = Y[i] * ig_phi[i] * mu_prime[i];
    if (Yphi == 0.0) {
      X1[i] = 1.0;
    } else if (Yphi > 5e5) {
      X1[i] = 1.0 / Yphi;
    } else {
      X1[i] = 1.0 + Yphi / 2.0 * (1.0 - std::sqrt(1.0 + 4.0 / Yphi));
    }
  }
  std::vector<double> tau_temp(m_edge);
  for (int i = 0; i < m_edge; ++i) {
    double u = R::unif_rand();
    bool firstroot = u < 1.0 / (1.0 + X1[i]);
    double ri = firstroot ? X1[i] : 1.0 / X1[i];
    double ig_val = mu_prime[i] * ri;
    tau_temp[i] = std::max(1e-12, 1.0 / ig_val);
  }

  arma::mat tau1(p, p, arma::fill::zeros);
  for (int i = 0; i < m_edge; ++i) {
    tau1.memptr()[upperind[i] - 1] = tau_temp[i];
    tau1.memptr()[lowerind[i] - 1] = tau_temp[i];
  }

  // Update Gamma column by column.
  for (int jj = 1; jj <= p; ++jj) {
    int j0 = jj - 1;

    arma::uvec ind_noi(p - 1);
    for (int k = 0; k < p - 1; ++k) ind_noi[k] = (arma::uword)(ind_noi_all(k, j0) - 1);

    arma::vec tau_col(p - 1);
    for (int k = 0; k < p - 1; ++k) tau_col[k] = tau1(ind_noi[k], j0);

    arma::mat ig_sub = inv_gamma.submat(ind_noi, ind_noi);
    arma::vec ig_col(p - 1);
    for (int k = 0; k < p - 1; ++k) ig_col[k] = inv_gamma(ind_noi[k], j0);
    double ig_jj = inv_gamma(j0, j0);

    arma::mat outer_ig = badgers_blas::outer(ig_col);
    arma::mat inv_gamma_11(p - 1, p - 1);
    for (int a = 0; a < p - 1; ++a)
      for (int b = 0; b < p - 1; ++b)
        inv_gamma_11(a, b) = ig_sub(a, b) - outer_ig(a, b) / ig_jj;

    double s_plus_iota = S(j0, j0) + iota_diag[j0];
    arma::mat gamma_i(p - 1, p - 1);
    for (int a = 0; a < p - 1; ++a)
      for (int b = 0; b < p - 1; ++b)
        gamma_i(a, b) = s_plus_iota * inv_gamma_11(a, b);
    for (int k = 0; k < p - 1; ++k) gamma_i(k, k) += 1.0 / tau_col[k];

    arma::mat gamma_i_sym(p - 1, p - 1);
    for (int a = 0; a < p - 1; ++a)
      for (int b = 0; b < p - 1; ++b)
        gamma_i_sym(a, b) = (gamma_i(a, b) + gamma_i(b, a)) / 2.0;

    arma::mat R_chol = arma::chol(gamma_i_sym); // upper triangular

    arma::vec Scol(p - 1);
    for (int k = 0; k < p - 1; ++k) Scol[k] = S(ind_noi[k], j0);
    arma::vec fw = arma::solve(arma::trimatl(R_chol.t()), Scol);
    arma::vec mu_gamma_i = -arma::solve(arma::trimatu(R_chol), fw);

    arma::vec z(p - 1);
    for (int k = 0; k < p - 1; ++k) z[k] = R::norm_rand();
    arma::vec beta_i = mu_gamma_i + arma::solve(arma::trimatu(R_chol), z);

    for (int k = 0; k < p - 1; ++k) {
      gamma(ind_noi[k], j0) = beta_i[k];
      gamma(j0, ind_noi[k]) = beta_i[k];
    }

    double gamma_indiv = R::rgamma(n / 2.0 + 1.0, 1.0 / (s_plus_iota / 2.0));

    arma::vec inv_gamma_11_beta = badgers_blas::matvec(inv_gamma_11, beta_i);
    double quad = badgers_blas::seqsum_prod(beta_i, inv_gamma_11_beta);
    gamma(j0, j0) = gamma_indiv + quad;

    arma::mat outer_igb = badgers_blas::outer(inv_gamma_11_beta);
    for (int a = 0; a < p - 1; ++a) {
      for (int b = 0; b < p - 1; ++b) {
        inv_gamma(ind_noi[a], ind_noi[b]) = inv_gamma_11(a, b) + outer_igb(a, b) / gamma_indiv;
      }
    }
    for (int k = 0; k < p - 1; ++k) {
      double v = -inv_gamma_11_beta[k] / gamma_indiv;
      inv_gamma(ind_noi[k], j0) = v;
      inv_gamma(j0, ind_noi[k]) = v;
    }
    inv_gamma(j0, j0) = 1.0 / gamma_indiv;
  }

  return List::create(_["gamma"] = wrap(gamma), _["inv_gamma"] = wrap(inv_gamma),
                       _["iota_edge"] = iota_edge, _["iota_diag"] = iota_diag);
}

// rinvgamma() note: provided by both the invgamma and nimble packages with
// different `scale` conventions; the driver script's library() load order
// resolves it to nimble's (mean = scale/(shape-1)), ported here as
// 1/R::rgamma(shape, 1/scale).

// ---- sigma^2 update (Gibbs) ----
// y_Xbeta = ytrain-(xtrain%*%beta)
// sigma2_scale_post = b_sigma + (sum(y_Xbeta^2) + sum(beta^2/tau^2))/2
// sigma2 = rinvgamma(1, shape = sigma2_shape_post, scale = sigma2_scale_post)
// [[Rcpp::export]]
double update_sigma2_cpp(NumericVector ytrain, NumericMatrix xtrain, NumericVector beta,
                          NumericVector tau, double b_sigma, double sigma2_shape_post) {
  int ntrain = xtrain.nrow();
  int p = xtrain.ncol();

  arma::mat X = as<arma::mat>(xtrain);
  arma::vec b = as<arma::vec>(beta);
  arma::vec Xbeta = badgers_blas::matvec(X, b);

  long double sumsq_resid = 0.0L;
  for (int i = 0; i < ntrain; ++i) {
    double r = ytrain[i] - Xbeta[i];
    sumsq_resid += r * r;
  }

  long double sumsq_beta_tau = 0.0L;
  for (int k = 0; k < p; ++k) {
    double term = (beta[k] * beta[k]) / (tau[k] * tau[k]);
    sumsq_beta_tau += term;
  }

  // Each sum() is cast to double individually before adding, matching R's
  // two separate sum() calls added at the interpreter level.
  double sigma2_scale_post = b_sigma + ((double) sumsq_resid + (double) sumsq_beta_tau) / 2.0;
  return 1.0 / R::rgamma(sigma2_shape_post, 1.0 / sigma2_scale_post);
}

// ---- alpha update (p sequential MH steps) ----
// alpha is mutated in place across j, so later j's use the just-updated
// alpha vector. Draws exactly 2 random numbers per j (rnorm then runif)
// unconditionally, matching R's `if (runif(1) < prob)`.
// [[Rcpp::export]]
List update_alpha_cpp(NumericVector alpha, NumericVector tau, NumericMatrix phi,
                       NumericVector alpha_sd_prop, double nu, double mu) {
  int p = alpha.size();
  NumericVector alpha_new = clone(alpha);
  IntegerVector accepted(p);

  for (int jj = 1; jj <= p; ++jj) {
    int j0 = jj - 1;
    double tau_j = tau[j0];
    double alpha_j = alpha_new[j0];
    double alpha_j_prop = R::rnorm(alpha_j, alpha_sd_prop[j0]);

    double comp = p_density_mmh_cpp(alpha_j_prop, alpha_j, tau_j, alpha_new, phi, jj, nu, mu);
    double u = R::unif_rand();
    bool accept;
    if (std::isfinite(comp) && comp > 0.0) {
      accept = u < std::min(comp, 1.0);
    } else {
      long double delta = (long double) alpha_j_prop - alpha_j;
      long double quad = 0.0L;
      for (int h = 0; h < p; ++h) {
        long double midpoint = (long double) alpha_new[h] - mu;
        if (h == j0) midpoint += delta / 2.0L;
        quad += (long double) phi(j0, h) * midpoint;
      }
      // Combine log(tau) and alpha before exponentiating, including on
      // platforms where long double has the same range as double.
      long double log_tau = std::log((long double) tau_j);
      long double base = 2.0L * ((long double) alpha_j + log_tau);
      long double penalty;
      if (delta == 0.0L) {
        penalty = 0.0L;
      } else if (std::fabs(delta) < 0.5L) {
        penalty = std::exp(base) * std::expm1(2.0L * delta);
      } else {
        penalty = std::exp(2.0L * ((long double) alpha_j_prop + log_tau)) -
          std::exp(base);
      }
      long double log_ratio = 2.0L * delta - penalty / 2.0L - delta * quad / nu;
      accept = std::log((long double) u) < log_ratio;
    }
    if (accept) {
      alpha_new[j0] = alpha_j_prop;
      accepted[j0] = 1;
    }
  }

  return List::create(_["alpha"] = alpha_new, _["accepted"] = accepted);
}

// ---- beta update (Gibbs) ----
// A_inv = crossprod(xtrain); diag(A_inv) += 1/tau^2
// A_inv_chol = chol(A_inv/sigma2); A = chol2inv(A_inv_chol)
// beta_mean_post = tcrossprod(A,xtrain)%*%ytrain/sigma2
// beta = beta_mean_post + rmnorm_chol(1, 0, A_inv_chol, prec_param=TRUE)
//
// tcrossprod(A,xtrain)%*%ytrain/sigma2 is computed in exactly that order
// (%*% binds tighter than / in R), not regrouped into A%*%(t(xtrain)%*%ytrain)
// -- matrix multiplication is not bit-exact-associative.
//
// rmnorm_chol(1, 0, R, prec_param=TRUE) is backsolve(R, rnorm(p)) (per
// nimble's C source), ported as arma::solve(trimatu(...)). Caveat: on macOS,
// nimble.so links Apple's Accelerate framework as well as R's own BLAS, and
// its rmnorm_chol resolves to Accelerate for the triangular solve -- about 1
// ULP off from this port (which matches R's own backsolve()/BLAS) on that
// platform. Confirmed a non-issue on the Linux deployment target, where
// nimble uses the same BLAS as R.
// [[Rcpp::export]]
NumericVector update_beta_cpp(NumericMatrix xtrain, NumericVector ytrain,
                               NumericVector tau, double sigma2) {
  int p = xtrain.ncol();
  arma::mat X = as<arma::mat>(xtrain);

  arma::mat A_inv = badgers_blas::crossprod(X);
  for (int k = 0; k < p; ++k) A_inv(k, k) += 1.0 / (tau[k] * tau[k]);

  arma::mat scaled = A_inv / sigma2;
  arma::mat A_inv_chol = arma::chol(scaled);

  arma::mat A = badgers_blas::chol2inv(A_inv_chol);

  arma::mat T = badgers_blas::tcrossprod(A, X);
  arma::vec yv = as<arma::vec>(ytrain);
  arma::vec beta_mean_post = badgers_blas::matvec(T, yv);
  for (int k = 0; k < p; ++k) beta_mean_post[k] /= sigma2;

  arma::vec z(p);
  for (int k = 0; k < p; ++k) z[k] = R::norm_rand();
  arma::vec draw = arma::solve(arma::trimatu(A_inv_chol), z);

  NumericVector beta(p);
  for (int k = 0; k < p; ++k) beta[k] = beta_mean_post[k] + draw[k];
  return beta;
}

// ---- tau update (p independent Gibbs draws) ----
// 1/tau_j^2 ~ IN(lambda_j*sigma/|beta_j|, lambda_j^2)
// [[Rcpp::export]]
NumericVector update_tau_cpp(NumericVector beta, NumericVector alpha, double sigma2) {
  int p = beta.size();
  NumericVector tau(p);
  double sqrt_sigma2 = std::sqrt(sigma2);

  for (int j = 0; j < p; ++j) {
    double lambda_j = std::exp(alpha[j]);
    double mean_param = lambda_j * sqrt_sigma2 / std::fabs(beta[j]);
    double shape_param = lambda_j * lambda_j;

    double z = R::norm_rand();
    double Y = z * z;
    double phi = 1.0 / shape_param;
    double Yphi = Y * phi * mean_param;
    double X1;
    if (Yphi == 0.0) {
      X1 = 1.0;
    } else if (Yphi > 5e5) {
      X1 = 1.0 / Yphi;
    } else {
      X1 = 1.0 + Yphi / 2.0 * (1.0 - std::sqrt(1.0 + 4.0 / Yphi));
    }
    double u = R::unif_rand();
    bool firstroot = u < 1.0 / (1.0 + X1);
    double ri = firstroot ? X1 : 1.0 / X1;
    if (beta[j] == 0.0) {
      // At beta=0, tau^2 | beta has Gamma(1/2, rate=lambda^2/2).
      // Use the already drawn normal; retain the uniform RNG consumption.
      tau[j] = std::fabs(z) / lambda_j;
      continue;
    }
    double ig_draw = mean_param * ri;

    double tau2 = 1.0 / ig_draw;
    tau[j] = std::sqrt(tau2);
  }
  return tau;
}

// ---- Phi update (Gibbs, per-column with 0/1/>1-neighbor branches) ----
// B is A_Phi = eta*(11^T + zeta*I) + (alpha-mu)(alpha-mu)^T/nu
//
// diag(zeta,p) is a p x p matrix with zeta on the diagonal and 0
// elsewhere (R's diag(x, nrow) form, not the diag(x)-is-identity-of-size-x
// case), so eta*(diag(zeta,p)+1) has eta*(zeta+1) on the diagonal and
// eta off it.
//
// The non_zero_phi draw in the >1-neighbor branch uses R's own
// backsolve()/forwardsolve() directly (not rmnorm_chol), so this block has
// no Accelerate-related caveat.
// [[Rcpp::export]]
List update_phi_cpp(NumericMatrix phi_, NumericMatrix inv_phi_,
                     NumericVector alpha, double mu,
                     double eta, double zeta, double nu) {
  int p = alpha.size();
  arma::mat phi = as<arma::mat>(phi_);
  arma::mat inv_phi = as<arma::mat>(inv_phi_);

  arma::vec alpha_mu(p);
  for (int k = 0; k < p; ++k) alpha_mu[k] = alpha[k] - mu;

  arma::mat outer_am = badgers_blas::outer(alpha_mu);
  arma::mat B(p, p);
  for (int a = 0; a < p; ++a) {
    for (int b = 0; b < p; ++b) {
      double base = (a == b) ? eta * (zeta + 1.0) : eta;
      B(a, b) = base + outer_am(a, b) / nu;
    }
  }

  double gamma_shape = 1.0 + (eta * (1.0 + zeta) + 1.0) / 2.0;
  int m_full = p - 1;

  for (int jj = 1; jj <= p; ++jj) {
    int j0 = jj - 1;
    double B_jj = B(j0, j0);

    arma::uvec not_j(m_full);
    { int c = 0; for (int k = 0; k < p; ++k) if (k != j0) not_j[c++] = k; }

    double inv_phi_22 = inv_phi(j0, j0);
    arma::vec inv_phi_12(m_full);
    for (int k = 0; k < m_full; ++k) inv_phi_12[k] = inv_phi(not_j[k], j0);
    arma::mat inv_phi_11_sub = inv_phi.submat(not_j, not_j);

    arma::mat outer_ip12 = badgers_blas::outer(inv_phi_12);
    arma::mat P_11(m_full, m_full);
    for (int a = 0; a < m_full; ++a)
      for (int b = 0; b < m_full; ++b)
        P_11(a, b) = inv_phi_11_sub(a, b) - outer_ip12(a, b) / inv_phi_22;

    std::vector<arma::uword> nz_ind;
    for (int k = 0; k < p; ++k)
      if (k != j0 && phi(k, j0) != 0.0) nz_ind.push_back((arma::uword)k);
    int m = (int)nz_ind.size();

    double xi = 0.0, diag_phi = 0.0;

    if (m != 0) {
      arma::uvec nzu(nz_ind);
      arma::vec non_zero_phi(m);

      if (m > 1) {
        arma::mat ip_nz_nz = inv_phi.submat(nzu, nzu);
        arma::vec ip_nz_j(m);
        for (int k = 0; k < m; ++k) ip_nz_j[k] = inv_phi(nzu[k], j0);
        arma::mat outer_ipnzj = badgers_blas::outer(ip_nz_j);
        arma::mat prec_nz(m, m);
        for (int a = 0; a < m; ++a)
          for (int b = 0; b < m; ++b)
            prec_nz(a, b) = ip_nz_nz(a, b) - outer_ipnzj(a, b) / inv_phi_22;

        arma::mat prec_nz_sqrt = arma::chol(prec_nz); // upper triangular

        double sd_mid = 1.0 / std::sqrt(B_jj);
        arma::vec mid(m);
        for (int k = 0; k < m; ++k) mid[k] = R::rnorm(0.0, sd_mid);

        double neg_inv_Bjj = -1.0 / B_jj;
        arma::vec pt_mean_o1(m);
        for (int k = 0; k < m; ++k) pt_mean_o1[k] = neg_inv_Bjj * B(nzu[k], j0);

        arma::vec pt_mean_o2 = arma::solve(arma::trimatl(prec_nz_sqrt.t()), pt_mean_o1);
        arma::vec mean_nz = arma::solve(arma::trimatu(prec_nz_sqrt), pt_mean_o2);
        arma::vec bs_mid = arma::solve(arma::trimatu(prec_nz_sqrt), mid);
        for (int k = 0; k < m; ++k) non_zero_phi[k] = bs_mid[k] + mean_nz[k];

        xi = R::rgamma(gamma_shape, 1.0 / (B_jj / 2.0));

        arma::vec partial = badgers_blas::matvec(prec_nz_sqrt, non_zero_phi);
        arma::blas_int bm = m, binc = 1;
        double cp = arma::ddot_(&bm, partial.memptr(), &binc, partial.memptr(), &binc);
        diag_phi = xi + cp;

      } else { // m == 1
        arma::uword nzi = nzu[0];
        double ip_nz_j = inv_phi(nzi, j0);
        double prec_nz = inv_phi(nzi, nzi) - (ip_nz_j * ip_nz_j) / inv_phi_22;

        double mean_nz = (-1.0 / (B_jj * prec_nz)) * B(nzi, j0);
        double sd_nz = 1.0 / std::sqrt(prec_nz * B_jj);
        non_zero_phi[0] = R::rnorm(mean_nz, sd_nz);

        xi = R::rgamma(gamma_shape, 1.0 / (B_jj / 2.0));
        diag_phi = xi + (non_zero_phi[0] * non_zero_phi[0]) * prec_nz;
      }

      for (int k = 0; k < m; ++k) {
        phi(nzu[k], j0) = non_zero_phi[k];
        phi(j0, nzu[k]) = non_zero_phi[k];
      }
      phi(j0, j0) = diag_phi;

    } else { // m == 0
      xi = R::rgamma(gamma_shape, 1.0 / (B_jj / 2.0));
      diag_phi = xi;
      phi(j0, j0) = diag_phi;
    }

    // update inv_phi (always runs, uses the just-updated phi column and xi)
    inv_phi(j0, j0) = 1.0 / xi;

    arma::vec phi_col(m_full);
    for (int k = 0; k < m_full; ++k) phi_col[k] = phi(not_j[k], j0);

    arma::vec pp = badgers_blas::matvec(P_11, phi_col);
    arma::mat outer_pp = badgers_blas::outer(pp);
    arma::mat inv_phi_11(m_full, m_full);
    for (int a = 0; a < m_full; ++a)
      for (int b = 0; b < m_full; ++b)
        inv_phi_11(a, b) = P_11(a, b) + inv_phi(j0, j0) * outer_pp(a, b);

    arma::vec ip11_phicol = badgers_blas::matvec(inv_phi_11, phi_col);
    arma::vec inv_phi_12_new(m_full);
    for (int k = 0; k < m_full; ++k) inv_phi_12_new[k] = -ip11_phicol[k] / phi(j0, j0);

    for (int a = 0; a < m_full; ++a)
      for (int b = 0; b < m_full; ++b)
        inv_phi(not_j[a], not_j[b]) = inv_phi_11(a, b);
    for (int k = 0; k < m_full; ++k) {
      inv_phi(not_j[k], j0) = inv_phi_12_new[k];
      inv_phi(j0, not_j[k]) = inv_phi_12_new[k];
    }
  }

  return List::create(_["phi"] = wrap(phi), _["inv_phi"] = wrap(inv_phi));
}

// Shared implementations for blocks touching the (typically very large)
// Wtrain/Wtune/Wtest replicate arrays. The main MCMC loop converts these to
// arma::cube ONCE outside the n_iter loop and calls these directly, instead
// of re-converting on every iteration. The Rcpp::export'd wrappers below
// just convert once and delegate here.
namespace badgers_core {

struct XUpdateResult { NumericMatrix xtrain, xtune, xtest; };

// ---- X update (train/tune/test latent variables, Gibbs) ----
// A_x_inv = gamma + tcrossprod(beta)/sigma2; diag(A_x_inv) += r_i/rho2
// A_x_inv_tune = gamma; diag(A_x_inv_tune) += r_i/rho2 (no beta term)
// x_i = A_x %*% vec + rmnorm_chol(1, 0, A_x_inv_chol, prec_param=TRUE), same
// for tune/test using A_x_tune/A_x_inv_chol_tune (see update_beta_cpp for
// the rmnorm_chol/Accelerate caveat, which applies here too).
inline XUpdateResult update_X_core(const arma::cube &Wtrain, const arma::cube &Wtune,
                                    const arma::cube &Wtest, NumericVector ytrain,
                                    NumericVector beta, NumericVector r_i,
                                    const arma::mat &gamma, double sigma2, double rho2) {
  int p = gamma.n_rows;
  int ntrain = (int) Wtrain.n_cols, ntune = (int) Wtune.n_cols, ntest = (int) Wtest.n_cols;

  arma::vec beta_v = as<arma::vec>(beta);
  arma::mat outer_beta = badgers_blas::outer(beta_v);
  arma::mat A_x_inv(p, p);
  for (int a = 0; a < p; ++a)
    for (int b = 0; b < p; ++b)
      A_x_inv(a, b) = gamma(a, b) + outer_beta(a, b) / sigma2;
  for (int k = 0; k < p; ++k) A_x_inv(k, k) += r_i[k] / rho2;

  arma::mat A_x_inv_chol = arma::chol(A_x_inv);
  arma::mat A_x = badgers_blas::chol2inv(A_x_inv_chol);

  arma::mat A_x_inv_tune = gamma;
  for (int k = 0; k < p; ++k) A_x_inv_tune(k, k) += r_i[k] / rho2;

  arma::mat A_x_inv_chol_tune = arma::chol(A_x_inv_tune);
  arma::mat A_x_tune = badgers_blas::chol2inv(A_x_inv_chol_tune);

  NumericMatrix xtrain(ntrain, p), xtune(ntune, p), xtest(ntest, p);

  for (int ii = 0; ii < ntrain; ++ii) {
    arma::vec w_i_sum(p);
    for (int g = 0; g < p; ++g) {
      long double s = 0.0L; // sum(na.omit(x)): long double accumulator
      for (arma::uword r = 0; r < Wtrain.n_rows; ++r) {
        double v = Wtrain(r, ii, g);
        if (!std::isnan(v)) s += v;
      }
      w_i_sum[g] = (double) s;
    }

    arma::vec vec_(p);
    for (int g = 0; g < p; ++g)
      vec_[g] = beta[g] * ytrain[ii] / sigma2 + w_i_sum[g] / rho2;

    arma::vec mean_part = badgers_blas::matvec(A_x, vec_);

    arma::vec z(p);
    for (int k = 0; k < p; ++k) z[k] = R::norm_rand();
    arma::vec draw = arma::solve(arma::trimatu(A_x_inv_chol), z);

    for (int g = 0; g < p; ++g) xtrain(ii, g) = mean_part[g] + draw[g];
  }

  for (int ii = 0; ii < ntune; ++ii) {
    arma::vec w_i_sum(p);
    for (int g = 0; g < p; ++g) {
      long double s = 0.0L;
      for (arma::uword r = 0; r < Wtune.n_rows; ++r) {
        double v = Wtune(r, ii, g);
        if (!std::isnan(v)) s += v;
      }
      w_i_sum[g] = (double) s;
    }
    arma::vec vec_(p);
    for (int g = 0; g < p; ++g) vec_[g] = w_i_sum[g] / rho2;

    arma::vec mean_part = badgers_blas::matvec(A_x_tune, vec_);
    arma::vec z(p);
    for (int k = 0; k < p; ++k) z[k] = R::norm_rand();
    arma::vec draw = arma::solve(arma::trimatu(A_x_inv_chol_tune), z);

    for (int g = 0; g < p; ++g) xtune(ii, g) = mean_part[g] + draw[g];
  }

  for (int ii = 0; ii < ntest; ++ii) {
    arma::vec w_i_sum(p);
    for (int g = 0; g < p; ++g) {
      long double s = 0.0L;
      for (arma::uword r = 0; r < Wtest.n_rows; ++r) {
        double v = Wtest(r, ii, g);
        if (!std::isnan(v)) s += v;
      }
      w_i_sum[g] = (double) s;
    }
    arma::vec vec_(p);
    for (int g = 0; g < p; ++g) vec_[g] = w_i_sum[g] / rho2;

    arma::vec mean_part = badgers_blas::matvec(A_x_tune, vec_);
    arma::vec z(p);
    for (int k = 0; k < p; ++k) z[k] = R::norm_rand();
    arma::vec draw = arma::solve(arma::trimatu(A_x_inv_chol_tune), z);

    for (int g = 0; g < p; ++g) xtest(ii, g) = mean_part[g] + draw[g];
  }

  return XUpdateResult{xtrain, xtune, xtest};
}

// ---- rho^2 update (Gibbs) ----
// rho2_scale_post = b_rho + sum over j, r=1:r_i[j] of sum((Wtrain[r,,j]-xtrain[,j])^2)/2
// The outer accumulation is R's own interpreted for-loop (plain double); the
// inner sum(...) is a genuine R sum() call (long double accumulator).
// r_i[j] is trusted as-is: Wtrain[1:r_i[j], , j] is assumed NA-free.
inline double update_rho2_core(const arma::cube &Wtrain, const arma::mat &xtrain,
                                NumericVector r_i, double b_rho, double rho2_shape_post) {
  int ntrain = (int) Wtrain.n_cols;
  int p = (int) Wtrain.n_slices;

  double rho2_scale_post = b_rho;
  for (int j = 0; j < p; ++j) {
    int rij = (int) r_i[j];
    for (int r = 0; r < rij; ++r) {
      long double s = 0.0L;
      for (int i = 0; i < ntrain; ++i) {
        double diff = Wtrain(r, i, j) - xtrain(i, j);
        s += diff * diff;
      }
      rho2_scale_post += (double) s / 2.0;
    }
  }

  return 1.0 / R::rgamma(rho2_shape_post, 1.0 / rho2_scale_post);
}

} // namespace badgers_core

// [[Rcpp::export]]
List update_X_cpp(NumericVector Wtrain_, NumericVector Wtune_, NumericVector Wtest_,
                   NumericVector ytrain, NumericVector beta, NumericVector r_i,
                   NumericMatrix gamma_, double sigma2, double rho2) {
  arma::cube Wtrain = as<arma::cube>(Wtrain_);
  arma::cube Wtune  = as<arma::cube>(Wtune_);
  arma::cube Wtest  = as<arma::cube>(Wtest_);
  arma::mat gamma = as<arma::mat>(gamma_);

  badgers_core::XUpdateResult res = badgers_core::update_X_core(Wtrain, Wtune, Wtest, ytrain, beta,
                                                         r_i, gamma, sigma2, rho2);
  return List::create(_["xtrain"] = res.xtrain, _["xtune"] = res.xtune, _["xtest"] = res.xtest);
}

// [[Rcpp::export]]
double update_rho2_cpp(NumericVector Wtrain_, NumericMatrix xtrain,
                        NumericVector r_i, double b_rho, double rho2_shape_post) {
  arma::cube Wtrain = as<arma::cube>(Wtrain_);
  arma::mat xtrain_m = as<arma::mat>(xtrain);
  return badgers_core::update_rho2_core(Wtrain, xtrain_m, r_i, b_rho, rho2_shape_post);
}

// ---- Gamma/iota update call site ----
// S <- crossprod(xtrain); wires it into AdaptiveGLasso_cpp.
// [[Rcpp::export]]
List update_gamma_cpp(NumericMatrix xtrain, double n, int p,
                       NumericMatrix gamma_, NumericMatrix inv_gamma_,
                       double a_iota, NumericVector b_iota_edge, double b_iota_diag,
                       IntegerVector upperind, IntegerVector lowerind,
                       IntegerMatrix ind_noi_all) {
  arma::mat X = as<arma::mat>(xtrain);
  arma::mat S = badgers_blas::crossprod(X);
  return AdaptiveGLasso_cpp(wrap(S), n, p, gamma_, inv_gamma_, a_iota, b_iota_edge,
                               b_iota_diag, upperind, lowerind, ind_noi_all);
}

// ---- adapt MH proposal SD (runs only when i %% k == 0) ----
// s = i %/% k; gamma1 = 1/(s+3)^c1; gamma2 = 10*gamma1
// r_alpha[,s] = apply(acs_alpha[,(s*k-k+1):(s*k)], 1, mean)
// adaptFactor = exp(gamma2*(r_alpha[,s]-r_opt)); alpha_sd_prop *= adaptFactor
//
// R's mean() is the two-pass algorithm (long double), not sum/n. pow(x,0.8)
// is a genuine fractional exponent (unlike the `^2` cases elsewhere, which
// lower to x*x) -- std::pow() matches R's `^` here since both fall through
// to the platform libm pow().
// [[Rcpp::export]]
List adapt_alpha_sd_cpp(NumericMatrix acs_alpha, int i, int k, double c1, double r_opt,
                         NumericVector alpha_sd_prop) {
  int p = acs_alpha.nrow();
  int s = i / k; // matches R's i %/% k

  double gamma1 = 1.0 / std::pow((double)(s + 3), c1);
  double gamma2 = 10.0 * gamma1;

  int col_lo = s * k - k; // 0-based first column of the block (R: s*k-k+1, 1-based)
  int col_hi = s * k - 1; // 0-based last column of the block  (R: s*k, 1-based)
  int kk = col_hi - col_lo + 1;

  NumericVector r_alpha_col(p);
  NumericVector new_sd(p);
  for (int j = 0; j < p; ++j) {
    long double sum1 = 0.0L;
    for (int c = col_lo; c <= col_hi; ++c) sum1 += acs_alpha(j, c);
    long double mean_val = sum1 / kk;
    if (R_FINITE((double) mean_val)) {
      long double t = 0.0L;
      for (int c = col_lo; c <= col_hi; ++c) t += (acs_alpha(j, c) - mean_val);
      mean_val += t / kk;
    }
    // r_alpha[,s] is stored as plain double in R before adaptFactor reads it
    // back, so round down to double here too rather than keep using the
    // wider mean_val.
    double mean_val_d = (double) mean_val;
    r_alpha_col[j] = mean_val_d;

    double adaptFactor = std::exp(gamma2 * (mean_val_d - r_opt));
    new_sd[j] = alpha_sd_prop[j] * adaptFactor;
  }

  return List::create(_["alpha_sd_prop"] = new_sd, _["r_alpha_col"] = r_alpha_col);
}

// badgers_mcmc_cpp(): full driver, assembling every block above into one loop
// matching utility_v2.R's MCMC_v2() setup, n_iter loop, and post-processing
// (both lite modes).
// Setup/post-processing pieces specific to this function:
//   - upperind/lowerind/ind_noi_all: pure integer arithmetic replicating
//     Reshape()/pracma::triu()'s column-major traversal.
//   - bfs_distances() replaces igraph::distances().
//   - gamma <- solve(inv_gamma) uses badgers_blas::solve_general (dgesv).
//   - rowMeans() (ytune_hat/ytest_hat) is plain sum/n, unlike mean().
// Wtrain/Wtune/Wtest are converted to arma::cube once here, not per
// iteration. R's `print(i)` in the adapt-SD block is kept as a progress line
// every k iterations; output only, no effect on results.
// [[Rcpp::export]]
List badgers_mcmc_cpp(NumericVector ytrain, NumericVector Wtrain_,
                  NumericVector ytune, NumericVector Wtune_,
                  NumericVector ytest, NumericVector Wtest_,
                  NumericMatrix adj_graph,
                  int n_iter, int burnin, int k,
                  double a_sigma, double b_sigma,
                  double a_iota, double b_iota_diag,
                  double psi,
                  double a_rho, double b_rho,
                  double eta, double zeta, double mu, double nu,
                  double c = 2.0, bool lite = true,
                  bool verbose = false) {

  if (!R_FINITE(c) || c <= 0) {
    stop("c must be a finite positive scalar");
  }

  IntegerVector dim_train = Wtrain_.attr("dim");
  IntegerVector dim_tune  = Wtune_.attr("dim");
  IntegerVector dim_test  = Wtest_.attr("dim");
  int ntrain = dim_train[1];
  int ntune  = dim_tune[1];
  int ntest  = dim_test[1];
  int p      = dim_train[2];

  arma::cube Wtrain = as<arma::cube>(Wtrain_);
  arma::cube Wtune  = as<arma::cube>(Wtune_);
  arma::cube Wtest  = as<arma::cube>(Wtest_);

  // ---- storage ----
  NumericMatrix BETA(p, n_iter);
  NumericMatrix Ytune_hat(ntune, n_iter);
  NumericMatrix Ytest_hat(ntest, n_iter);

  NumericVector SIGMA, RHO;
  NumericMatrix TAU, ALPHA;
  arma::cube IOTA_arr, GAMMA_arr, Xtrain_arr, Xtune_arr, Xtest_arr;
  if (!lite) {
    SIGMA = NumericVector(n_iter);
    RHO   = NumericVector(n_iter);
    TAU   = NumericMatrix(p, n_iter);
    ALPHA = NumericMatrix(p, n_iter);
    IOTA_arr   = arma::cube(p, p, n_iter);
    GAMMA_arr  = arma::cube(p, p, n_iter);
    Xtrain_arr = arma::cube(ntrain, p, n_iter);
    Xtune_arr  = arma::cube(ntune,  p, n_iter);
    Xtest_arr  = arma::cube(ntest,  p, n_iter);
  }

  NumericVector beta(p, 0.0);
  NumericVector tau(p, 1.0);
  NumericVector alpha(p, mu);
  double sigma2 = 1.0, rho2 = 1.0;

  // Initialize the latent true X using the mean of the replicates for each subject-gene pair
  auto init_X = [](const arma::cube &W) -> NumericMatrix {
    int Rmax = (int) W.n_rows, n = (int) W.n_cols, pp = (int) W.n_slices;
    NumericMatrix X(n, pp);
    std::vector<double> vals;
    for (int g = 0; g < pp; ++g) {
      for (int i = 0; i < n; ++i) {
        vals.clear();
        for (int r = 0; r < Rmax; ++r) {
          double v = W(r, i, g);
          if (!std::isnan(v)) vals.push_back(v);
        }
        X(i, g) = badgers_stats::mean_twopass(vals);
      }
    }
    return X;
  };
  NumericMatrix xtrain = init_X(Wtrain);
  NumericMatrix xtune  = init_X(Wtune);
  NumericMatrix xtest  = init_X(Wtest);

  arma::mat graph_a = as<arma::mat>(adj_graph);
  for (int i = 0; i < p; ++i) graph_a(i, i) = 1.0;
  arma::mat phi_a(p, p, arma::fill::eye);
  for (int a = 0; a < p; ++a)
    for (int b = 0; b < p; ++b)
      if (a != b && graph_a(a, b) != 0.0) phi_a(a, b) = 0.05; // fixed initial weight
  arma::mat phi_chol;
  if (!arma::chol(phi_chol, phi_a)) {
    // Keep the graph's zero pattern, but use strict diagonal dominance when
    // the historical 0.05 initialization is not positive definite.
    for (int a = 0; a < p; ++a) {
      double degree = 0.0;
      for (int b = 0; b < p; ++b) if (a != b) degree += std::fabs(phi_a(a, b));
      phi_a(a, a) = 1.0 + degree;
    }
    phi_chol = arma::chol(phi_a);
  }
  arma::mat inv_phi_a = badgers_blas::chol2inv(phi_chol);
  NumericMatrix phi = wrap(phi_a);
  NumericMatrix inv_phi = wrap(inv_phi_a);

  arma::mat xtrain_a0 = as<arma::mat>(xtrain);
  arma::mat S0 = badgers_blas::crossprod(xtrain_a0);
  arma::mat inv_gamma_a = S0 / (double) ntrain;
  for (int i = 0; i < p; ++i) inv_gamma_a(i, i) += 0.5;
  arma::mat gamma_a = badgers_blas::solve_general(inv_gamma_a);
  NumericMatrix gamma = wrap(gamma_a);
  NumericMatrix inv_gamma = wrap(inv_gamma_a);

  // adaptive MH settings for alpha
  NumericVector alpha_sd_prop(p, 1.0);
  NumericMatrix acs_alpha(p, n_iter);
  int n_r_alpha_cols = n_iter / k; // floor(n_iter/k)
  NumericMatrix r_alpha(p, n_r_alpha_cols);
  double r_opt = 0.44, c1 = 0.8;

  NumericVector r_i(p);
  for (int g = 0; g < p; ++g) {
    int cnt = 0;
    for (arma::uword r = 0; r < Wtrain.n_rows; ++r) if (!std::isnan(Wtrain(r, 0, g))) cnt++;
    r_i[g] = cnt;
  }

  std::vector<int> upperind_vec, lowerind_vec;
  for (int col = 0; col < p; ++col)
    for (int row = 0; row < p; ++row)
      if (row < col) upperind_vec.push_back(col * p + row + 1);
  for (int col = 0; col < p; ++col)
    for (int row = 0; row < p; ++row)
      if (row < col) lowerind_vec.push_back(row * p + col + 1);
  IntegerVector upperind = wrap(upperind_vec);
  IntegerVector lowerind = wrap(lowerind_vec);

  IntegerMatrix ind_noi_all(p - 1, p);
  for (int j0 = 0; j0 < p; ++j0) {
    int row = 0;
    for (int v = 1; v <= p; ++v) { if (v - 1 == j0) continue; ind_noi_all(row++, j0) = v; }
  }

  if (!R_FINITE(psi) || psi <= 0) {
    stop("psi must be a finite positive scalar");
  }

  // Peterson et al. network-informed rate construction.
  arma::mat adj_gamma_d = as<arma::mat>(adj_graph);
  for (int i = 0; i < p; ++i) adj_gamma_d(i, i) = 0.0;
  arma::Mat<int> adj_int(p, p);
  for (int a = 0; a < p; ++a)
    for (int b = 0; b < p; ++b)
      adj_int(a, b) = (adj_gamma_d(a, b) != 0.0) ? 1 : 0;
  arma::mat dist_graph = badgers_stats::bfs_distances(adj_int);

  double disconnected_rate = std::pow(10.0, -psi);  // b_{iota,jk} = 10^-psi for pairs without a path
  double graph_rate_scale = std::pow(10.0, std::log10(disconnected_rate) + c);
  arma::mat b_iota_graph(p, p);
  b_iota_graph.fill(disconnected_rate);
  for (int a = 0; a < p; ++a) {
    for (int b = 0; b < p; ++b) {
      double d = dist_graph(a, b);
      if (std::isfinite(d) && d > 0) {
        b_iota_graph(a, b) = std::max(disconnected_rate, graph_rate_scale / d);
      }
    }
  }
  NumericVector b_iota_edge(upperind.size());
  for (int t2 = 0; t2 < upperind.size(); ++t2) b_iota_edge[t2] = b_iota_graph.memptr()[upperind[t2] - 1];

  double sigma2_shape_post = a_sigma + (ntrain + p) / 2.0;
  long double sum_r_i_ld = 0.0L;
  for (int g = 0; g < p; ++g) sum_r_i_ld += r_i[g];
  double sum_r_i = (double) sum_r_i_ld;
  double rho2_shape_post = a_rho + ntrain * sum_r_i / 2.0;

  // =========================== MCMC ===========================
  for (int i = 1; i <= n_iter; ++i) {
    Rcpp::checkUserInterrupt();
    int i0 = i - 1;

    // ---- sigma^2 ----
    sigma2 = update_sigma2_cpp(ytrain, xtrain, beta, tau, b_sigma, sigma2_shape_post);
    if (!lite) SIGMA[i0] = sigma2;

    // ---- alpha ----
    {
      List out = update_alpha_cpp(alpha, tau, phi, alpha_sd_prop, nu, mu);
      alpha = out["alpha"];
      IntegerVector accepted = out["accepted"];
      for (int j = 0; j < p; ++j) acs_alpha(j, i0) = accepted[j];
    }
    if (!lite) ALPHA(_, i0) = alpha;

    // ---- beta ----
    beta = update_beta_cpp(xtrain, ytrain, tau, sigma2);
    BETA(_, i0) = beta;

    // ---- tau ----
    tau = update_tau_cpp(beta, alpha, sigma2);
    if (!lite) TAU(_, i0) = tau;

    // ---- Phi ----
    {
      List out = update_phi_cpp(phi, inv_phi, alpha, mu, eta, zeta, nu);
      phi = as<NumericMatrix>(out["phi"]);
      inv_phi = as<NumericMatrix>(out["inv_phi"]);
    }

    // ---- X ----
    {
      arma::mat gamma_arma = as<arma::mat>(gamma);
      badgers_core::XUpdateResult xr = badgers_core::update_X_core(Wtrain, Wtune, Wtest, ytrain, beta,
                                                            r_i, gamma_arma, sigma2, rho2);
      xtrain = xr.xtrain; xtune = xr.xtune; xtest = xr.xtest;
    }
    if (!lite) {
      Xtrain_arr.slice(i0) = as<arma::mat>(xtrain);
      Xtune_arr.slice(i0)  = as<arma::mat>(xtune);
      Xtest_arr.slice(i0)  = as<arma::mat>(xtest);
    }

    // ---- Gamma and iota ----
    NumericVector iota_edge, iota_diag;
    {
      List out = update_gamma_cpp(xtrain, (double) ntrain, p, gamma, inv_gamma,
                                   a_iota, b_iota_edge, b_iota_diag, upperind, lowerind, ind_noi_all);
      gamma = as<NumericMatrix>(out["gamma"]);
      inv_gamma = as<NumericMatrix>(out["inv_gamma"]);
      iota_edge = out["iota_edge"];
      iota_diag = out["iota_diag"];
    }
    if (!lite) {
      GAMMA_arr.slice(i0) = as<arma::mat>(gamma);
      arma::mat iota_matrix(p, p, arma::fill::zeros);
      for (int t2 = 0; t2 < upperind.size(); ++t2) {
        iota_matrix.memptr()[upperind[t2] - 1] = iota_edge[t2];
        iota_matrix.memptr()[lowerind[t2] - 1] = iota_edge[t2];
      }
      for (int d2 = 0; d2 < p; ++d2) iota_matrix(d2, d2) = iota_diag[d2];
      IOTA_arr.slice(i0) = iota_matrix;
    }

    // ---- rho^2 ----
    {
      arma::mat xtrain_arma = as<arma::mat>(xtrain);
      rho2 = badgers_core::update_rho2_core(Wtrain, xtrain_arma, r_i, b_rho, rho2_shape_post);
    }
    if (!lite) RHO[i0] = rho2;

    // ---- adapt MH proposal SD ----
    if (i % k == 0) {
      int s = i / k;
      List out = adapt_alpha_sd_cpp(acs_alpha, i, k, c1, r_opt, alpha_sd_prop);
      alpha_sd_prop = out["alpha_sd_prop"];
      NumericVector r_col = out["r_alpha_col"];
      for (int j = 0; j < p; ++j) r_alpha(j, s - 1) = r_col[j];
      if (verbose) Rcpp::Rcout << i << "\n";
    }

    // ---- save tuning/test posterior ----
    {
      arma::mat xtune_a = as<arma::mat>(xtune);
      arma::mat xtest_a = as<arma::mat>(xtest);
      arma::vec beta_a = as<arma::vec>(beta);
      arma::vec ytune_pred = badgers_blas::matvec(xtune_a, beta_a);
      arma::vec ytest_pred = badgers_blas::matvec(xtest_a, beta_a);
      for (int ii = 0; ii < ntune; ++ii) Ytune_hat(ii, i0) = ytune_pred[ii];
      for (int ii = 0; ii < ntest; ++ii) Ytest_hat(ii, i0) = ytest_pred[ii];
    }
  }

  // ===== calculate error rate and pick gene =====
  int n_keep = n_iter - burnin;

  arma::mat Ytune_hat_a = as<arma::mat>(Ytune_hat);
  arma::mat Ytest_hat_a = as<arma::mat>(Ytest_hat);
  NumericVector ytune_hat(ntune), ytest_hat(ntest);
  for (int ii = 0; ii < ntune; ++ii) {
    long double s = 0.0L;
    for (int it = burnin; it < n_iter; ++it) s += Ytune_hat_a(ii, it);
    ytune_hat[ii] = (double) (s / n_keep);
  }
  for (int ii = 0; ii < ntest; ++ii) {
    long double s = 0.0L;
    for (int it = burnin; it < n_iter; ++it) s += Ytest_hat_a(ii, it);
    ytest_hat[ii] = (double) (s / n_keep);
  }

  std::vector<double> sq_tune(ntune), sq_test(ntest);
  for (int ii = 0; ii < ntune; ++ii) { double d = ytune[ii] - ytune_hat[ii]; sq_tune[ii] = d * d; }
  for (int ii = 0; ii < ntest; ++ii) { double d = ytest[ii] - ytest_hat[ii]; sq_test[ii] = d * d; }
  double MSPE_tune = badgers_stats::mean_twopass(sq_tune);
  double MSPE_test = badgers_stats::mean_twopass(sq_test);

  arma::mat BETA_a = as<arma::mat>(BETA);
  NumericVector beta_mean(p), Beta_95(p), Beta_90(p), Beta_80(p);
  for (int j = 0; j < p; ++j) {
    std::vector<double> row_vals(n_keep);
    for (int it = 0; it < n_keep; ++it) row_vals[it] = BETA_a(j, burnin + it);
    beta_mean[j] = badgers_stats::mean_twopass(row_vals);

    double q025 = badgers_stats::quantile7(row_vals, 0.025);
    double q975 = badgers_stats::quantile7(row_vals, 0.975);
    Beta_95[j] = (q025 <= 0 && q975 >= 0) ? 0.0 : badgers_stats::mean_twopass(row_vals);

    double q05 = badgers_stats::quantile7(row_vals, 0.05);
    double q95 = badgers_stats::quantile7(row_vals, 0.95);
    Beta_90[j] = (q05 <= 0 && q95 >= 0) ? 0.0 : badgers_stats::mean_twopass(row_vals);

    double q10 = badgers_stats::quantile7(row_vals, 0.10);
    double q90 = badgers_stats::quantile7(row_vals, 0.90);
    Beta_80[j] = (q10 <= 0 && q90 >= 0) ? 0.0 : badgers_stats::mean_twopass(row_vals);
  }

  // ===== return results =====
  List result = List::create(
    _["CI95"] = Beta_95,
    _["CI90"] = Beta_90,
    _["CI80"] = Beta_80,
    _["MSPE_test"] = MSPE_test,
    _["MSPE_tune"] = MSPE_tune,
    _["mu"] = mu,
    _["nu"] = nu,
    _["psi"] = psi,
    _["c"] = c,
    _["BETA"] = BETA,
    _["r_alpha"] = r_alpha,
    _["acs_alpha"] = acs_alpha,
    _["alpha_sd_prop"] = alpha_sd_prop,
    _["beta_mean"] = beta_mean
  );
  result["Xtrain"] = lite ? (SEXP) R_NilValue : (SEXP) wrap(Xtrain_arr);
  result["Xtune"]  = lite ? (SEXP) R_NilValue : (SEXP) wrap(Xtune_arr);
  result["Xtest"]  = lite ? (SEXP) R_NilValue : (SEXP) wrap(Xtest_arr);
  result["IOTA"]   = lite ? (SEXP) R_NilValue : (SEXP) wrap(IOTA_arr);
  result["SIGMA2"] = lite ? (SEXP) R_NilValue : (SEXP) wrap(SIGMA);
  result["ALPHA"]  = lite ? (SEXP) R_NilValue : (SEXP) wrap(ALPHA);
  result["RHO2"]   = lite ? (SEXP) R_NilValue : (SEXP) wrap(RHO);
  result["TAU"]    = lite ? (SEXP) R_NilValue : (SEXP) wrap(TAU);
  result["GAMMA"]  = lite ? (SEXP) R_NilValue : (SEXP) wrap(GAMMA_arr);

  return result;
}
