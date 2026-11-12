// [[Rcpp::depends(RcppArmadillo)]]

#define ARMA_WARN_LEVEL 1
#include <RcppArmadillo.h>
#include <queue>
using namespace Rcpp;

// Sampler updates use Armadillo for linear algebra and R's RNG for draws.
// Keep draw order stable so numerical changes can be compared under a fixed seed.
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
  arma::vec vec(p);
  for (int k = 0; k < p; ++k) {
    vec[k] = (k == j0) ? (alpha[j0] + alpha_j_prop - two_mu)
                       : (alpha[k] + alpha[k] - two_mu);
  }

  // comp3 = exp((-1/(2*nu)) * (alpha_j_prop - alpha_j) * (phi_j %*% vec))
  // phi_j = phi[j, ] is read directly from phi's column-major storage with
  // stride p (no row copy needed).
  arma::mat phi_view(phi.begin(), p, p, false, true);
  double quad = arma::dot(phi_view.row(j0), vec);
  double comp3 = std::exp((-1.0 / (2.0 * nu)) * (alpha_j_prop - alpha_j) * quad);

  return comp1 * comp2 * comp3;
}

namespace badgers_graph {

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

} // namespace badgers_graph

// Solve a Gaussian conditional mean using its existing precision factor.
inline arma::vec solve_precision(const arma::mat &upper, const arma::vec &rhs) {
  arma::vec intermediate = arma::solve(arma::trimatl(upper.t()), rhs,
                                       arma::solve_opts::no_approx);
  return arma::solve(arma::trimatu(upper), intermediate, arma::solve_opts::no_approx);
}

// Adaptive graphical lasso update for latent-covariate precision Gamma.
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
    if (p == 1) {
      // There is no off-diagonal Gaussian block in the scalar conditional.
      double rate = (S(0, 0) + iota_diag[0]) / 2.0;
      gamma(0, 0) = R::rgamma(n / 2.0 + 1.0, 1.0 / rate);
      inv_gamma(0, 0) = 1.0 / gamma(0, 0);
      continue;
    }

    arma::uvec ind_noi(p - 1);
    for (int k = 0; k < p - 1; ++k) ind_noi[k] = (arma::uword)(ind_noi_all(k, j0) - 1);

    arma::vec tau_col(p - 1);
    for (int k = 0; k < p - 1; ++k) tau_col[k] = tau1(ind_noi[k], j0);

    arma::mat ig_sub = inv_gamma.submat(ind_noi, ind_noi);
    arma::vec ig_col(p - 1);
    for (int k = 0; k < p - 1; ++k) ig_col[k] = inv_gamma(ind_noi[k], j0);
    double ig_jj = inv_gamma(j0, j0);

    arma::mat outer_ig = (ig_col * ig_col.t());
    arma::mat inv_gamma_11 = ig_sub - outer_ig / ig_jj;

    double s_plus_iota = S(j0, j0) + iota_diag[j0];
    arma::mat gamma_i = s_plus_iota * inv_gamma_11;
    gamma_i.diag() += 1.0 / tau_col;
    arma::mat gamma_i_sym = (gamma_i + gamma_i.t()) / 2.0;

    arma::mat R_chol = arma::chol(gamma_i_sym); // upper triangular

    arma::vec Scol(p - 1);
    for (int k = 0; k < p - 1; ++k) Scol[k] = S(ind_noi[k], j0);
    arma::vec fw = arma::solve(arma::trimatl(R_chol.t()), Scol, arma::solve_opts::no_approx);
    arma::vec mu_gamma_i = -arma::solve(arma::trimatu(R_chol), fw, arma::solve_opts::no_approx);

    arma::vec z(p - 1);
    for (int k = 0; k < p - 1; ++k) z[k] = R::norm_rand();
    arma::vec beta_i = mu_gamma_i + arma::solve(arma::trimatu(R_chol), z, arma::solve_opts::no_approx);

    for (int k = 0; k < p - 1; ++k) {
      gamma(ind_noi[k], j0) = beta_i[k];
      gamma(j0, ind_noi[k]) = beta_i[k];
    }

    double gamma_indiv = R::rgamma(n / 2.0 + 1.0, 1.0 / (s_plus_iota / 2.0));

    arma::vec inv_gamma_11_beta = (inv_gamma_11 * beta_i);
    double quad = arma::dot(beta_i, inv_gamma_11_beta);
    gamma(j0, j0) = gamma_indiv + quad;

    arma::mat outer_igb = (inv_gamma_11_beta * inv_gamma_11_beta.t());
    inv_gamma.submat(ind_noi, ind_noi) = inv_gamma_11 + outer_igb / gamma_indiv;
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

// ---- sigma^2 update (inverse-Gamma conditional) ----
// [[Rcpp::export]]
double update_sigma2_cpp(NumericVector ytrain, NumericMatrix xtrain, NumericVector beta,
                          NumericVector tau, double b_sigma, double sigma2_shape_post) {
  arma::mat X = as<arma::mat>(xtrain);
  arma::vec residual = as<arma::vec>(ytrain) - X * as<arma::vec>(beta);
  arma::vec beta_square = arma::square(as<arma::vec>(beta));
  arma::vec tau_square = arma::square(as<arma::vec>(tau));
  double scale = b_sigma + (arma::accu(arma::square(residual)) +
                            arma::accu(beta_square / tau_square)) / 2.0;
  return 1.0 / R::rgamma(sigma2_shape_post, 1.0 / scale);
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

// ---- beta update (Gaussian conditional) ----
// [[Rcpp::export]]
NumericVector update_beta_cpp(NumericMatrix xtrain, NumericVector ytrain,
                               NumericVector tau, double sigma2) {
  arma::mat X = as<arma::mat>(xtrain);
  arma::vec y = as<arma::vec>(ytrain);
  arma::vec scales = as<arma::vec>(tau);
  arma::mat precision = X.t() * X;
  precision.diag() += 1.0 / arma::square(scales);
  arma::mat upper = arma::chol(precision / sigma2);
  arma::vec mean = solve_precision(upper, X.t() * y / sigma2);

  arma::vec z(X.n_cols);
  for (arma::uword j = 0; j < z.n_elem; ++j) z[j] = R::norm_rand();
  arma::vec draw = mean + arma::solve(arma::trimatu(upper), z,
                                     arma::solve_opts::no_approx);
  return NumericVector(draw.begin(), draw.end());
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

// ---- Phi update, restricted to the working graph ----
// Conditional rate: eta * (11^T + epsilon * I) + (alpha-mu)(alpha-mu)^T / nu.
// [[Rcpp::export]]
List update_phi_cpp(NumericMatrix phi_, NumericMatrix inv_phi_,
                     NumericVector alpha, double mu,
                     double eta, double epsilon, double nu) {
  int p = alpha.size();
  arma::mat phi = as<arma::mat>(phi_);
  arma::mat inv_phi = as<arma::mat>(inv_phi_);

  arma::vec alpha_mu = as<arma::vec>(alpha) - mu;
  arma::mat B(p, p, arma::fill::value(eta));
  B.diag().fill(eta * (epsilon + 1.0));
  B += alpha_mu * alpha_mu.t() / nu;

  double gamma_shape = 1.0 + (eta * (1.0 + epsilon) + 1.0) / 2.0;
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

    arma::mat outer_ip12 = (inv_phi_12 * inv_phi_12.t());
    arma::mat P_11 = inv_phi_11_sub - outer_ip12 / inv_phi_22;

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
        arma::mat outer_ipnzj = (ip_nz_j * ip_nz_j.t());
        arma::mat prec_nz = ip_nz_nz - outer_ipnzj / inv_phi_22;

        arma::mat prec_nz_sqrt = arma::chol(prec_nz); // upper triangular

        double sd_mid = 1.0 / std::sqrt(B_jj);
        arma::vec mid(m);
        for (int k = 0; k < m; ++k) mid[k] = R::rnorm(0.0, sd_mid);

        double neg_inv_Bjj = -1.0 / B_jj;
        arma::vec pt_mean_o1(m);
        for (int k = 0; k < m; ++k) pt_mean_o1[k] = neg_inv_Bjj * B(nzu[k], j0);

        arma::vec pt_mean_o2 = arma::solve(arma::trimatl(prec_nz_sqrt.t()), pt_mean_o1, arma::solve_opts::no_approx);
        arma::vec mean_nz = arma::solve(arma::trimatu(prec_nz_sqrt), pt_mean_o2, arma::solve_opts::no_approx);
        arma::vec bs_mid = arma::solve(arma::trimatu(prec_nz_sqrt), mid, arma::solve_opts::no_approx);
        for (int k = 0; k < m; ++k) non_zero_phi[k] = bs_mid[k] + mean_nz[k];

        xi = R::rgamma(gamma_shape, 1.0 / (B_jj / 2.0));

        arma::vec partial = (prec_nz_sqrt * non_zero_phi);
        double cp = arma::dot(partial, partial);
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

    arma::vec pp = (P_11 * phi_col);
    arma::mat outer_pp = (pp * pp.t());
    arma::mat inv_phi_11 = P_11 + inv_phi(j0, j0) * outer_pp;

    arma::vec ip11_phicol = (inv_phi_11 * phi_col);
    arma::vec inv_phi_12_new = -ip11_phicol / phi(j0, j0);

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

// ---- X update (Gaussian conditionals for training and held-out data) ----
inline XUpdateResult update_X_core(const arma::cube &Wtrain, const arma::cube &Wtune,
                                    const arma::cube &Wtest, NumericVector ytrain,
                                    NumericVector beta, NumericVector r_i,
                                    const arma::mat &gamma, double sigma2, double rho2) {
  arma::vec b = as<arma::vec>(beta);
  arma::vec counts = as<arma::vec>(r_i);
  arma::mat train_precision = gamma + b * b.t() / sigma2;
  train_precision.diag() += counts / rho2;
  arma::mat heldout_precision = gamma;
  heldout_precision.diag() += counts / rho2;
  arma::mat train_upper = arma::chol(train_precision);
  arma::mat heldout_upper = arma::chol(heldout_precision);

  auto sample_rows = [&](const arma::cube &W, const arma::mat &upper,
                         bool training) -> NumericMatrix {
    arma::mat X(W.n_cols, W.n_slices);
    for (arma::uword i = 0; i < W.n_cols; ++i) {
      arma::vec rhs(W.n_slices);
      for (arma::uword j = 0; j < W.n_slices; ++j) {
        rhs[j] = arma::accu(W.slice(j).col(i).head((arma::uword) counts[j])) / rho2;
      }
      if (training) rhs += b * ytrain[i] / sigma2;
      arma::vec mean = solve_precision(upper, rhs);
      arma::vec z(W.n_slices);
      for (arma::uword j = 0; j < z.n_elem; ++j) z[j] = R::norm_rand();
      X.row(i) = (mean + arma::solve(arma::trimatu(upper), z,
                                    arma::solve_opts::no_approx)).t();
    }
    return wrap(X);
  };

  NumericMatrix xtrain = sample_rows(Wtrain, train_upper, true);
  NumericMatrix xtune = sample_rows(Wtune, heldout_upper, false);
  NumericMatrix xtest = sample_rows(Wtest, heldout_upper, false);
  return XUpdateResult{xtrain, xtune, xtest};
}

// ---- rho^2 update (inverse-Gamma conditional) ----
inline double update_rho2_core(const arma::cube &Wtrain, const arma::mat &xtrain,
                                NumericVector r_i, double b_rho, double rho2_shape_post) {
  double scale = b_rho;
  for (arma::uword j = 0; j < Wtrain.n_slices; ++j) {
    for (int r = 0; r < r_i[j]; ++r) {
      arma::vec residual = Wtrain.slice(j).row(r).t() - xtrain.col(j);
      scale += arma::dot(residual, residual) / 2.0;
    }
  }
  return 1.0 / R::rgamma(rho2_shape_post, 1.0 / scale);
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
  arma::mat S = (X.t() * X);
  return AdaptiveGLasso_cpp(wrap(S), n, p, gamma_, inv_gamma_, a_iota, b_iota_edge,
                               b_iota_diag, upperind, lowerind, ind_noi_all);
}

// ---- adapt MH proposal SD at the end of an acceptance block ----
// [[Rcpp::export]]
List adapt_alpha_sd_cpp(NumericMatrix acs_alpha, int i, int k, double c1, double r_opt,
                         NumericVector alpha_sd_prop) {
  arma::mat accepted(acs_alpha.begin(), acs_alpha.nrow(), acs_alpha.ncol(), false, true);
  arma::vec rates = arma::mean(accepted.cols(i - k, i - 1), 1);
  double step = 10.0 * (1.0 / std::pow((double) (i / k + 3), c1));
  arma::vec sd = as<arma::vec>(alpha_sd_prop) % arma::exp(step * (rates - r_opt));
  return List::create(_["alpha_sd_prop"] = NumericVector(sd.begin(), sd.end()),
                       _["r_alpha_col"] = NumericVector(rates.begin(), rates.end()));
}

// Fixed-hyperparameter sampler. Posterior summaries are computed in R.
// Replicate arrays are converted once, outside the MCMC loop.
// [[Rcpp::export]]
List badgers_mcmc_cpp(NumericVector ytrain, NumericVector Wtrain_,
                  NumericVector Wtune_, NumericVector Wtest_,
                  NumericMatrix adj_graph,
                  int n_iter, int k,
                  double a_sigma, double b_sigma,
                  double a_iota, double b_iota,
                  double a_rho, double b_rho,
                  double eta, double epsilon, double mu, double nu,
                  double c_graph = 2.0, bool lite = true,
                  bool verbose = false) {

  if (!R_FINITE(c_graph) || c_graph <= 0 || c_graph >= 6) {
    stop("c_graph must be a finite scalar strictly between 0 and 6");
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
        X(i, g) = arma::mean(arma::vec(vals));
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
  arma::mat inv_phi_a = arma::inv_sympd(phi_a);
  NumericMatrix phi = wrap(phi_a);
  NumericMatrix inv_phi = wrap(inv_phi_a);

  arma::mat xtrain_a0 = as<arma::mat>(xtrain);
  arma::mat S0 = (xtrain_a0.t() * xtrain_a0);
  arma::mat inv_gamma_a = S0 / (double) ntrain;
  for (int i = 0; i < p; ++i) inv_gamma_a(i, i) += 0.5;
  arma::mat gamma_a = arma::inv_sympd(inv_gamma_a);
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

  double b_iota_diag = b_iota;

  // Peterson et al. network-informed rate construction.
  arma::mat adj_gamma_d = as<arma::mat>(adj_graph);
  for (int i = 0; i < p; ++i) adj_gamma_d(i, i) = 0.0;
  arma::Mat<int> adj_int(p, p);
  for (int a = 0; a < p; ++a)
    for (int b = 0; b < p; ++b)
      adj_int(a, b) = (adj_gamma_d(a, b) != 0.0) ? 1 : 0;
  arma::mat dist_graph = badgers_graph::bfs_distances(adj_int);

  double disconnected_rate = 1e-6;
  double graph_rate_scale = std::pow(10.0, -6.0 + c_graph);
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
  double sum_r_i = sum(r_i);
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
      List out = update_phi_cpp(phi, inv_phi, alpha, mu, eta, epsilon, nu);
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
      arma::vec ytune_pred = (xtune_a * beta_a);
      arma::vec ytest_pred = (xtest_a * beta_a);
      for (int ii = 0; ii < ntune; ++ii) Ytune_hat(ii, i0) = ytune_pred[ii];
      for (int ii = 0; ii < ntest; ++ii) Ytest_hat(ii, i0) = ytest_pred[ii];
    }
  }

  List result = List::create(
    _["mu"] = mu,
    _["nu"] = nu,
    _["c_graph"] = c_graph,
    _["BETA"] = BETA,
    _["r_alpha"] = r_alpha,
    _["acs_alpha"] = acs_alpha,
    _["alpha_sd_prop"] = alpha_sd_prop
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

  result["Ytune_hat"] = Ytune_hat;
  result["Ytest_hat"] = Ytest_hat;
  return result;
}
