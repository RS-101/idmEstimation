// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <algorithm>
#include <vector>
#include <cmath>
#include <cfloat>
#include <limits>
using namespace Rcpp;

struct ModelData {
  // --- counts (ints) ---
  int J, C, K_tilde, U, N_tilde, N1_obs_of_T_star, M, W,
  N, N_star, M_mark, I, K, I_mark, L;

  // --- scalars (double) ---
  double s_max, R_max, e_star_max;

  // --- vectors ---
  NumericVector s_j, L_c, t_c, e_k, L_u, t_u, t_m_in_N_tilde,
  L_m, R_m_vec, t_m, E_star, t_star_n,
  L_bar, R_bar, s_j_c, s_j_full, Q, Q_i_mark, A_union;

  // counts per unique time (int vectors)
  IntegerVector d_n;
  arma::rowvec c_k;

  // --- 2-col interval matrices ---
  NumericMatrix A_m, A_u, A_c, full_A_m, Q_i;

  // --- Indicator matrices ---
  LogicalMatrix alpha_ij;  // I_mark x (J+C)
  LogicalMatrix beta_im;   // I x M'
  LogicalMatrix I_rho_mn;  // M x N
  LogicalMatrix I_pi_un;   // U x N
  LogicalMatrix I_sigma_cn; // C x N


  // ---------- (18) alpha: I' x (J+C) ----------
  void cal_alpha() {
    int I = Q_i.nrow();
    int JpC = s_j_full.size();

    alpha_ij = Rcpp::LogicalMatrix(I_mark, JpC);
    for (int i = 0; i < I_mark; ++i) {
      double left = (i < I) ? Q_i(i,0) : Q_i_mark[i - I];
      for (int j = 0; j < JpC; ++j)
        alpha_ij(i,j) = (s_j_full[j] <= left);
    }
  }

  // ---------- (19) beta: I x M' ----------
  void cal_beta() {
    int I = Q_i.nrow();
    int M_mark = full_A_m.nrow();
    beta_im = Rcpp::LogicalMatrix(I, M_mark);

    for (int i = 0; i < I; ++i) {
      double Li = Q_i(i,0), Ri = Q_i(i,1);
      for (int m = 0; m < M_mark; ++m) {
        double Lm = full_A_m(m,0), Rm = full_A_m(m,1);
        beta_im(i,m) = (Lm <= Li && Ri <= Rm);
      }
    }
  }

  // ----- Calc Indicator matrices ------
  void cal_indicator_matrices() {
    I_rho_mn = Rcpp::LogicalMatrix(M, N);  // M x N
    I_pi_un = Rcpp::LogicalMatrix(U, N);   // U x N
    I_sigma_cn = Rcpp::LogicalMatrix(C, N); // C x N

    for(int n = 0; n < N; ++n) {
      for(int m = 0; m < M; ++m) {
        I_rho_mn(m,n) = t_m[m] >= t_star_n[n];
      }
      for(int u = 0; u < U; ++u) {
        I_pi_un(u,n) = t_u[u] >= t_star_n[n];
      }
      for(int c = 0; c < C; ++c) {
        I_sigma_cn(c,n) = t_c[c] >= t_star_n[n];
      }
    }
  }

  IntegerVector lambda_M_plus_u;

  // lambda_M_p_u <- lambda_n[t_star_n == t_u[u]]
  void calc_lambda_M_plus_u() {
    lambda_M_plus_u = IntegerVector(U, NA_INTEGER);
    for(int u = 0; u < U; ++u){
      // lambda_M_p_u <- lambda_n[t_star_n == t_u[u]]
      const double val = t_u[u];
      auto it = std::find_if(t_star_n.begin(), t_star_n.end(), [&](double x){ return std::abs(x - val) < 1e-9; });
      int index = (it != t_star_n.end()) ? std::distance(t_star_n.begin(), it) : -1;
      lambda_M_plus_u[u] = index;
    }
    // Rcpp::Rcout << "[em_fit] lambda_M_plus_u: " << lambda_M_plus_u << std::endl;
  }

  void check() {
    // ---------- dimension & consistency checks ----------
    if ((int)s_j.size() != J) stop("length(s_j) != J");
    if ((int)s_j_full.size() != (J + C)) stop("length(s_j_full) != J + C");

    if (Q_i.ncol() != 2) stop("Q_i must have 2 columns");
    if (Q_i.nrow() != I) stop("nrow(Q_i) != I");

    if (A_m.ncol() != 2 || A_u.ncol() != 2 || A_c.ncol() != 2 || full_A_m.ncol() != 2)
      stop("All A_* matrices must have 2 columns");
    if (A_m.nrow() != M) stop("nrow(A_m) != M");
    if (A_u.nrow() != U) stop("nrow(A_u) != U");
    if (A_c.nrow() != C) stop("nrow(A_c) != C");
    if (full_A_m.nrow() != M_mark) stop("nrow(full_A_m) != M_mark");

    // vectors tied to A_*
    if ((int)L_m.size() != M || (int)R_m_vec.size() != M || (int)t_m.size() != M)
      stop("L_m, R_m, and t_m must all have length M");
    if ((int)L_u.size() != U || (int)t_u.size() != U)
      stop("L_u and t_u must have length U");
    if ((int)L_c.size() != C || (int)t_c.size() != C)
      stop("L_c and t_c must have length C");

    // t_star_n / N
    if ((int)t_star_n.size() != N) stop("length(t_star_n) != N");
    if ((int)d_n.size() != N) stop("length(d_n) != N");
    if (N1_obs_of_T_star < 0 || N1_obs_of_T_star > N)
      stop("N1_obs_of_T_star must be in [0, N]");
    if (N_star <= 0) stop("N_star must be > 0");

    // c_k / E_star live on I_mark - I
    if (K != (I_mark - I)) stop("K must equal I_mark - I");
    if ((int)c_k.size() != K) stop("length(c_k) != K");
    if ((int)E_star.size() != K) stop("length(E_star) != K");

    // Q_i_mark used when i in [I, I_mark)
    if ((int)Q_i_mark.size() != (I_mark - I))
      stop("length(Q_i_mark) != I_mark - I");

    // beta_im indexing uses offsets M+l and W+l
    if (M + U > M_mark)
      stop("M + U must be <= M_mark for beta_im(i, M + l)");
    if (W + C > M_mark)
      stop("W + C must be <= M_mark for beta_im(i, W + l)");
  }


  // ctor from named R list and calculates alpha and beta
  ModelData(const List& x) {
    // ints
    J  = as<int>(x["J"]);
    C  = as<int>(x["C"]);
    K_tilde = as<int>(x["K_tilde"]);
    U  = as<int>(x["U"]);
    N_tilde = as<int>(x["N_tilde"]);
    N1_obs_of_T_star = as<int>(x["N1_obs_of_T_star"]);
    M  = as<int>(x["M"]);
    W  = as<int>(x["W"]);
    N  = as<int>(x["N"]);
    N_star = as<int>(x["N_star"]);
    M_mark = as<int>(x["M_mark"]);
    I  = as<int>(x["I"]);
    K  = as<int>(x["K"]);
    I_mark = as<int>(x["I_mark"]);

    // scalars
    s_max = as<double>(x["s_max"]);
    R_max = as<double>(x["R_max"]);
    e_star_max = as<double>(x["e_star_max"]);

    // vectors
    s_j = as<NumericVector>(x["s_j"]);
    L_c = as<NumericVector>(x["L_c"]);
    t_c = as<NumericVector>(x["t_c"]);
    e_k = as<NumericVector>(x["e_k"]);
    L_u = as<NumericVector>(x["L_u"]);
    t_u = as<NumericVector>(x["t_u"]);
    t_m_in_N_tilde = as<NumericVector>(x["t_m_in_N_tilde"]);
    L_m = as<NumericVector>(x["L_m"]);
    R_m_vec = as<NumericVector>(x["R_m"]);
    t_m = as<NumericVector>(x["t_m"]);
    E_star = as<NumericVector>(x["E_star"]);
    t_star_n = as<NumericVector>(x["t_star_n"]);
    L_bar = as<NumericVector>(x["L_bar"]);
    R_bar = as<NumericVector>(x["R_bar"]);
    s_j_c = as<NumericVector>(x["s_j_c"]);
    s_j_full = as<NumericVector>(x["s_j_full"]);
    Q = as<NumericVector>(x["Q"]);
    Q_i_mark = as<NumericVector>(x["Q_i_mark"]);
    A_union = as<NumericVector>(x["A_union"]);

    // integer vectors
    c_k = Rcpp::as<arma::rowvec>(x["c_k"]);   // coerces integer
    d_n = as<IntegerVector>(x["d_n"]);

    // matrices (2 columns) — the important part
    A_m      = as<NumericMatrix>(x["A_m"]);
    A_u      = as<NumericMatrix>(x["A_u"]);
    A_c      = as<NumericMatrix>(x["A_c"]);
    full_A_m = as<NumericMatrix>(x["full_A_m"]);
    Q_i      = as<NumericMatrix>(x["Q_i"]);

    cal_alpha();
    cal_beta();

    calc_lambda_M_plus_u();
    cal_indicator_matrices();
    L = std::min(M, std::min(U, std::min(C, J)));

    // checks
    check();
  }
};

// [[Rcpp::export]]
SEXP make_model_data(List x) {
  return XPtr<ModelData>(new ModelData(x), true);
}


// ---------- workspace + helpers ----------
struct Workspace {
  // parameters updated by EM
  arma::rowvec lambda_n;  // length N
  arma::rowvec z_i;       // length I_mark

  Rcpp::NumericVector t_sorted;   // sorted t_n
  Rcpp::NumericVector log_prefix; // prefix sum of log(1 - lambda_n) over sorted order
  Rcpp::NumericVector zero_prefix;// prefix count of zeros in (1 - lambda_n)

  // E-step pieces
  arma::mat mu_mi;
  arma::mat mu_bar_ji;
  arma::mat eta_ui;
  arma::mat gamma_ci;
};

void setup_prod(const ModelData& md, Workspace& ws) {
  // Precompute sorted order and prefix arrays.
  // - Sort indices by md.t_star_n
  // - t_sorted[i] = sorted t_star_n
  // - log_prefix[i+1] = sum_{k<=i, lambda_k!=1} log(1 - lambda_k)
  // - zero_prefix[i+1] = count_{k<=i} [lambda_k == 1]
  const R_xlen_t N = md.t_star_n.size();
  if (ws.lambda_n.size() != N) {
    stop("lambda_n length (%d) must match t_star_n length (%d).",
         (int)ws.lambda_n.size(), (int)N);
  }

  // build sorted index of t_star_n
  std::vector<int> idx(N);
  for (R_xlen_t i = 0; i < N; ++i) idx[i] = static_cast<int>(i);
  std::sort(idx.begin(), idx.end(),
            [&](int a, int b){ return md.t_star_n[a] < md.t_star_n[b]; });

  ws.t_sorted  = NumericVector(N);
  ws.log_prefix = NumericVector(N + 1);
  ws.zero_prefix = NumericVector(N + 1);

  ws.log_prefix[0]  = 0.0;
  ws.zero_prefix[0] = 0.0;

  // Treat lambda exactly equal to 1 as a zero factor (product becomes 0 if present in interval)
  const double TOL = 0.0; // set to, e.g., 1e-15 if you want tolerance
  for (R_xlen_t i = 0; i < N; ++i) {
    int k = idx[i];
    const double t   = md.t_star_n[k];
    const double lam = ws.lambda_n[k];

    ws.t_sorted[i] = t;

    const bool is_one = std::abs(1.0-lam) <= TOL;
    ws.zero_prefix[i + 1] = ws.zero_prefix[i] + (is_one ? 1.0 : 0.0);

    // Only add to log-prefix if not zero; if zero, keep running sum unchanged.
    if (!is_one) {
      const double one_minus = 1.0 - lam;
      if (one_minus <= 0.0) {
        // Defensive: log undefined. You can choose to stop() or handle differently.
        stop("Encountered 1 - lambda_n <= 0 at sorted index %d (lambda=%.17g).", (int)i, lam);
      }
      ws.log_prefix[i + 1] = ws.log_prefix[i] + std::log(one_minus);
    } else {
      ws.log_prefix[i + 1] = ws.log_prefix[i];
    }
  }
}

double evaluate(const Workspace& ws, double L, double R) {
  // Evaluate product_{t_n* in (L, R)} (1 - lambda_n).
  // We implement (L, R) with L open, R open by using:
  //   i = upper_bound(t_sorted, L)   -> first index with t > L
  //   j = lower_bound(t_sorted, R)   -> first index with t >= R
  // If any lambda==1 in (i, j-1), product is 0.
  // Else product = exp(log_prefix[j] - log_prefix[i]).
  const R_xlen_t N = ws.t_sorted.size();
  if (ws.log_prefix.size() != N + 1 || ws.zero_prefix.size() != N + 1)
    stop("Workspace not initialized. Call setup_prod first.");

  if (R <= L || N == 0) return 1.0;

  // Binary searches on sorted t
  const double* begin = &ws.t_sorted[0];
  const double* end   = begin + N;

  const R_xlen_t i = std::upper_bound(begin, end, L) - begin; // t > L
  const R_xlen_t j = std::lower_bound(begin, end, R) - begin; // t >= R

  if (j <= i) return 1.0; // empty interval

  const double zeros_in_range = ws.zero_prefix[j] - ws.zero_prefix[i];
  if (zeros_in_range > 0.0) return 0.0;

  const double log_prod = ws.log_prefix[j] - ws.log_prefix[i];
  return std::exp(log_prod);
}

// Move to the next/previous representable double (1 ULP)
inline double next_double(double x) {
  return std::nextafter(x, std::numeric_limits<double>::infinity());
}
inline double prev_double(double x) {
  return std::nextafter(x, -std::numeric_limits<double>::infinity());
}

bool is_double_eq(double a, double b) {
  return ((a - b) < DBL_EPSILON) && ((b - a) < DBL_EPSILON);
}

void print_summary(const ModelData& md, const Workspace& ws) {
  // --- Debug: sizes of core scalars/lengths ---------------------------------
  Rcpp::Rcout << "[calc_all] sizes:"
              << " I=" << md.I
              << " I_mark=" << md.I_mark
              << " N=" << md.N
              << " M=" << md.M
              << " J=" << md.J
              << " U=" << md.U
              << " C=" << md.C
              << " W=" << md.W
              << " t_star_n=" << md.t_star_n.size()
              << " N_star=" << md.N_star
              << std::endl;

  Rcpp::Rcout << "[calc_all] vector lengths:"
              << " z=" << ws.z_i.size()
              << " lambda=" << ws.lambda_n.size()
              << std::endl;

  // Input vectors used inside loops
  Rcpp::Rcout << "[calc_all] input vector lengths:"
              << " ws.z_i=" << ws.z_i.size()
              << " ws.lambda_n=" << ws.lambda_n.size()
              << " R_m_vec=" << md.R_m_vec.size()
              << " t_u=" << md.t_u.size()
              << " t_c=" << md.t_c.size()
              << " E_star=" << md.E_star.size()
              << " c_k=" << md.c_k.size()
              << " d_n=" << md.d_n.size()
              << std::endl;

  Rcpp::Rcout << "[calc_all] matrix shapes (rows x cols):"
  // using arma n_rows/n_cols for arma matrices and nrow()/ncol() for Rcpp matrices
  << " ws.mu_mi=" << ws.mu_mi.n_rows << "x" << ws.mu_mi.n_cols
  << " ws.mu_bar_ji=" << ws.mu_bar_ji.n_rows << "x" << ws.mu_bar_ji.n_cols
  << " ws.eta_ui=" << ws.eta_ui.n_rows << "x" << ws.eta_ui.n_cols
  << " ws.gamma_ci=" << ws.gamma_ci.n_rows << "x" << ws.gamma_ci.n_cols
  << " beta_im=" << md.beta_im.nrow()  << "x" << md.beta_im.ncol()
  << " alpha_ij=" << md.alpha_ij.nrow() << "x" << md.alpha_ij.ncol()
  << " Q_i=" << md.Q_i.nrow() << "x" << md.Q_i.ncol()
  << std::endl;

  Rcpp::Rcout << "[calc_all] L (min of M,U,C,J)=" << md.L << std::endl;


  // Loop ranges summary (avoid printing at every iteration)
  Rcpp::Rcout << "[calc_all] loop ranges:"
              << " i_main: 0.." << (md.I > 0 ? md.I - 1 : -1)
              << " i_extra(I_mark): " << md.I << ".." << (md.I_mark > 0 ? md.I_mark - 1 : -1)
              << " l_main: 0.." << (md.L > 0 ? md.L - 1 : -1)
              << " tails: M:" << md.L << ".." << (md.M > 0 ? md.M - 1 : -1)
              << " J:" << md.L << ".." << (md.J > 0 ? md.J - 1 : -1)
              << " U:" << md.L << ".." << (md.U > 0 ? md.U - 1 : -1)
              << " C:" << md.L << ".." << (md.C > 0 ? md.C - 1 : -1)
              << std::endl;

}

void run_em_once(const ModelData& md, Workspace& ws) {
  // --- Resets temps ----------------------------------------------------------
  ws.mu_mi.zeros(md.M, md.I);
  ws.mu_bar_ji.zeros(md.J, md.I_mark);
  ws.eta_ui.zeros(md.U, md.I_mark);
  ws.gamma_ci.zeros(md.C, md.I_mark);

  setup_prod(md,ws);

  // --- Loops over I ----------------------------------------------------------
  for(int i = 0; i < md.I; ++i) {
    for (int l = 0; l < md.L; ++l) {
      ws.mu_mi(l,i) = md.beta_im(i,l) ?
      ws.z_i[i]*evaluate(ws, md.Q_i(i,1), next_double(md.R_m_vec[l])) : 0;

      ws.mu_bar_ji(l,i) = md.alpha_ij(i,l) ? ws.z_i[i] : 0;

      ws.eta_ui(l,i) = md.beta_im(i, md.M + l) ? ws.lambda_n[md.lambda_M_plus_u[l]] *
        evaluate(ws, md.Q_i(i,1), md.t_u[l]) * ws.z_i[i] : 0 ;

      ws.gamma_ci(l,i) = (md.alpha_ij(i, md.J + l) ?
                            ws.z_i[i] : 0) +
                            (md.beta_im(i, md.W + l) ?
                            evaluate(ws, md.Q_i(i,1), next_double(md.t_c[l])) *
                            ws.z_i[i]: 0);
    }

    // tails
    for (int l = md.L; l < md.M; ++l) {
      ws.mu_mi(l,i) = md.beta_im(i,l) ?
      ws.z_i[i] * evaluate(ws, md.Q_i(i,1), next_double(md.R_m_vec[l])) : 0;
    }
    for (int l = md.L; l < md.J; ++l) {
      ws.mu_bar_ji(l,i) = md.alpha_ij(i,l) ?ws.z_i[i] : 0;
    }
    for (int l = md.L; l < md.U; ++l) {
      ws.eta_ui(l,i) = md.beta_im(i, md.M + l) ? ws.lambda_n[md.lambda_M_plus_u[l]] *
        evaluate(ws, md.Q_i(i,1), md.t_u[l]) * ws.z_i[i] : 0 ;
    }
    for (int l = md.L; l < md.C; ++l) {
      ws.gamma_ci(l,i) = (md.alpha_ij(i, md.J + l) ?
                            ws.z_i[i] : 0) +
                            (md.beta_im(i, md.W + l) ?
                            evaluate(ws, md.Q_i(i,1), next_double(md.t_c[l])) *
                            ws.z_i[i]: 0);
    }
  }

  // --- Extra for I_mark ------------------------------------------------------
  for(int i = md.I; i < md.I_mark; ++i) {
    for (int l = 0; l < md.L; ++l) {
      ws.mu_bar_ji(l,i) = md.alpha_ij(i,l) ?ws.z_i[i] : 0;
      ws.eta_ui(l,i) = is_double_eq(md.t_u[l], md.E_star[i-md.I]) ? ws.z_i[i] : 0;
      ws.gamma_ci(l,i) = md.alpha_ij(i, md.J + l) ?ws.z_i[i] : 0;
    }
    for (int l = md.L; l < md.J; ++l) {
      ws.mu_bar_ji(l,i) = md.alpha_ij(i,l) ?ws.z_i[i] : 0;
    }
    for (int l = md.L; l < md.U; ++l) {
      ws.eta_ui(l,i) = is_double_eq(md.t_u[l],md.E_star[i-md.I]) ? ws.z_i[i] : 0;
    }
    for (int l = md.L; l < md.C; ++l) {
      ws.gamma_ci(l,i) = md.alpha_ij(i, md.J + l) ?ws.z_i[i] : 0;
    }
  }

  // --- Normalizes rows -------------------------------------------------------
  ws.mu_mi = arma::normalise(ws.mu_mi, 1, 1);
  ws.mu_bar_ji = arma::normalise(ws.mu_bar_ji, 1, 1);
  ws.eta_ui = arma::normalise(ws.eta_ui, 1, 1);
  ws.gamma_ci = arma::normalise(ws.gamma_ci, 1, 1);

  double sum_rho_n;
  double sum_pi_n;
  double sum_pi_full_n;
  double sum_sigma_n;
  double d_n_part;

  // --- Loops over N ----------------------------------------------------------
  for(int n = 0; n < md.N; ++n) {
    sum_rho_n = 0;
    sum_pi_n = 0;
    sum_pi_full_n = 0;
    sum_sigma_n = 0;
    d_n_part = 0;


    d_n_part = n <= md.N1_obs_of_T_star ? md.d_n[n] : 0;
    for (int l = 0; l < md.L; ++l) {
      // rho_mn
      if(md.I_rho_mn(l,n)) {
        for(int i = 0; i < md.I; ++i) {
          sum_rho_n +=
            md.L_m[l] <= md.Q_i(i,0) && md.Q_i(i,1) < md.t_star_n[n] ?
            ws.mu_mi(l,i) : 0;
        }
      }
      // pi_un
      if(md.I_pi_un(l,n)) {
        for(int i = 0; i < md.I; ++i) {
          sum_pi_n +=
            md.L_u[l] <= md.Q_i(i,0) && md.Q_i(i,1) < md.t_star_n[n] ?
            ws.eta_ui(l,i) : 0;

          sum_pi_full_n += md.t_u[l] == md.t_star_n[n] ? ws.eta_ui(l,i) : 0;
        }
      }
      // sigma_cn
      if(md.I_sigma_cn(l,n)) {
        for(int i = 0; i < md.I; ++i) {
          sum_sigma_n +=
            md.L_c[l] <= md.Q_i(i,0) && md.Q_i(i,1) < md.t_star_n[n] ?
            ws.gamma_ci(l,i) : 0;
        }
      }

    }
    for (int l = md.L; l < md.M; ++l) {
      // rho_mn
      if(md.I_rho_mn(l,n)) {
        for(int i = 0; i < md.I; ++i) {
          sum_rho_n +=
            md.L_m[l] <= md.Q_i(i,0) && md.Q_i(i,1) < md.t_star_n[n] ?
            ws.mu_mi(l,i) : 0;
        }
      }
    }
    for (int l = md.L; l < md.U; ++l) {
      // pi_un
      if(md.I_pi_un(l,n)) {
        for(int i = 0; i < md.I; ++i) {
          sum_pi_n +=
            md.L_u[l] <= md.Q_i(i,0) && md.Q_i(i,1) < md.t_star_n[n] ?
            ws.eta_ui(l,i) : 0;

          sum_pi_full_n += is_double_eq(md.t_u[l], md.t_star_n[n]) ? ws.eta_ui(l,i) : 0;
        }
      }
    }
    for (int l = md.L; l < md.C; ++l) {
      // sigma_cn
      if(md.I_sigma_cn(l,n)) {
        for(int i = 0; i < md.I; ++i) {
          sum_sigma_n +=
            md.L_c[l] <= md.Q_i(i,0) && md.Q_i(i,1) < md.t_star_n[n] ?
            ws.gamma_ci(l,i) : 0;
        }
      }
    }

    // Rcpp::Rcout << "[em_fit] sum sum_pi_full_n " << sum_pi_full_n << std::endl;

    double demon = sum_rho_n + sum_pi_n + sum_sigma_n;
    //Rcpp::Rcout << "[em_fit] demon is " << demon << std::endl;

    ws.lambda_n[n] = (d_n_part + sum_pi_full_n)/demon;

    // if lambda_n[n] is larger than 1, we set it to 0.05
    if(ws.lambda_n[n] > 1.0){

     // Rcpp::Rcout << "[em_fit] Warning: lambda_n[" << n << "] > 1.0" << std::endl;
     // Rcpp::Rcout << "[em_fit] Corresponding to t_star_n[" << n << "] = " << md.t_star_n[n] << std::endl;
      ws.lambda_n[n] = 0.05;
    }

  }

  // --- Calculates new z ------------------------------------------------------
  arma::rowvec base = arma::sum(ws.mu_bar_ji, 0);
  base += arma::sum(ws.eta_ui, 0);
  base += arma::sum(ws.gamma_ci, 0);
  arma::rowvec new_z = base;

  new_z.head(md.I) += arma::sum(ws.mu_mi, 0);
  new_z.subvec(md.I, md.I_mark - 1) += md.c_k;
  new_z /= md.N_star;

  // // normalize z_i to sum to 1
  // new_z /= arma::accu(new_z);

  ws.z_i = new_z;


  //Rcpp::Rcout << "[em_fit] eta is " << ws.eta_ui << std::endl;
}


// [[Rcpp::export]]
Rcpp::List em_fit(SEXP md_ptr,
                  Rcpp::Nullable<Rcpp::NumericVector> z_init = R_NilValue,
                  Rcpp::Nullable<Rcpp::NumericVector> lambda_init = R_NilValue,
                  int max_iter = 100,
                  double tol = 1e-3,
                  bool verbose = true) {

  Rcpp::XPtr<ModelData> p(md_ptr);
  const ModelData& md = *p;

  Workspace ws;

  const arma::uword I_mark = static_cast<arma::uword>(md.I_mark);
  ws.z_i.set_size(I_mark);
  if (z_init.isNotNull()) {
    Rcpp::NumericVector z = z_init.get();
    if (static_cast<arma::uword>(z.size()) != I_mark)
      Rcpp::stop("length(z_init) must equal I_mark");
    for (arma::uword i = 0; i < I_mark; ++i) ws.z_i[i] = z[i];
  } else {
    ws.z_i.fill(1.0 / static_cast<double>(I_mark));
  }

  // make sure that z_i sums to 1
  ws.z_i /= arma::accu(ws.z_i);
  // print initial sum of z_i for debugging
  if (verbose) {
    double z_sum_init = arma::accu(ws.z_i);
    Rcpp::Rcout << "[em_fit] Initial sum(z_i) = " << z_sum_init << std::endl;
  }


  const arma::uword N = static_cast<arma::uword>(md.t_star_n.size());
  ws.lambda_n.set_size(N);
  if (lambda_init.isNotNull()) {
    Rcpp::NumericVector lam = lambda_init.get();
    if (static_cast<arma::uword>(lam.size()) != N)
      Rcpp::stop("length(lambda_init) must equal N (length(t_star_n))");
    for (arma::uword n = 0; n < N; ++n) ws.lambda_n[n] = lam[n];
  } else {
    ws.lambda_n.fill(0.005);
  }




  double dz = 0.0, dl = 0.0;
  for (int iter = 0; iter < max_iter; ++iter) {
    arma::rowvec prev_lambda_n = ws.lambda_n; // copy
    arma::rowvec prev_z_i      = ws.z_i;      // copy

    run_em_once(md, ws);

    dz = 0.0; dl = 0.0;
    for (arma::uword i = 0; i < I_mark; ++i)
      dz = std::max(dz, std::abs(ws.z_i[i] - prev_z_i[i]));
    for (arma::uword n = 0; n < N; ++n)
      dl = std::max(dl, std::abs(ws.lambda_n[n] - prev_lambda_n[n]));

    if (verbose) {
      Rcpp::Rcout << "iter " << (iter + 1)
                  << ": max|Δz|=" << dz
                  << " max|Δλ|=" << dl << "\n";
      // print sum of z_i for debugging

      double z_sum = arma::accu(ws.z_i);
      Rcpp::Rcout << "[em_fit] sum(z_i) = " << z_sum << std::endl;
    }

    if (std::isfinite(dz) && std::isfinite(dl) && std::max(dz, dl) < tol)
      break;



  }

  if (verbose) print_summary(md, ws);



  return Rcpp::List::create(
    Rcpp::_["z_i"]      = ws.z_i,
    Rcpp::_["lambda_n"] = ws.lambda_n,
    Rcpp::_["alpha_ij"] = md.alpha_ij,
    Rcpp::_["beta_im"]  = md.beta_im,
    Rcpp::_["mu_mi"]        = ws.mu_mi,
    Rcpp::_["mu_bar_ji"]    = ws.mu_bar_ji,
    Rcpp::_["eta_ui"]       = ws.eta_ui,
    Rcpp::_["gamma_ci"]     = ws.gamma_ci
  );
}
