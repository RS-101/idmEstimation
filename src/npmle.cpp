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
  int N_D, N_F, N_C, N_E, N_A, N_A_star, N_AB, N_B, N_ABE,
  N_AE_star, N, N_ABEF, N_CE_star, I, I_mark,
  L; // min of N_AB, N_E, N_F for loops

  // --- scalars (double) ---
  double s_max, R_max, e_star_max;

  // --- vectors ---
  NumericVector t_D, L_F, t_F, t_C, L_E, t_E,
  L_AB, R_AB, t_AB, t_CE_star, t_AE_star,
  t_DF;

  // counts per unique time (int vectors)
  IntegerVector r_A;
  arma::rowvec r_C;

  // --- 2-col interval matrices ---
  NumericMatrix LR_AB, LR_E, LR_F, LR_ABEF, Q_j;

  // --- Indicator matrices ---
  LogicalMatrix alpha_ji;  // I_mark x (N_D+N_F), my alpha + beta
  LogicalMatrix gamma_ji;   // I x N_ABEF, my gamma

  // ---------- (18) alpha: I' x (N_D+N_F) ----------
  void cal_alpha() {
    int I = Q_j.nrow();
    int N_DF = t_DF.size();

    alpha_ji = Rcpp::LogicalMatrix(I_mark, N_DF);
    for (int j = 0; j < I_mark; ++j) {
      double left = (j < I) ? Q_j(j,0) : t_CE_star[j - I];
      for (int i = 0; i < N_DF; ++i)
        alpha_ji(j,i) = (t_DF[i] <= left);
    }
  }

  // ---------- (19) beta: I x N_AB' ----------
  void cal_gamma() {
    int I = Q_j.nrow();
    int N_ABEF = LR_ABEF.nrow();
    gamma_ji = Rcpp::LogicalMatrix(I, N_ABEF);

    for (int i = 0; i < I; ++i) {
      double Li = Q_j(i,0), Ri = Q_j(i,1);
      for (int m = 0; m < N_ABEF; ++m) {
        double Lm = LR_ABEF(m,0), Rm = LR_ABEF(m,1);
        gamma_ji(i,m) = (Lm <= Li && Ri <= Rm);
      }
    }
  }

  IntegerVector lambda_M_plus_u;

  // lambda_M_p_u <- lambda_u[t_AE_star == t_E[u]]
  void calc_lambda_M_plus_u() {
    lambda_M_plus_u = IntegerVector(N_E, NA_INTEGER);
    for(int u = 0; u < N_E; ++u){
      // lambda_M_p_u <- lambda_u[t_AE_star == t_E[u]]
      const double val = t_E[u];
      auto it = std::find_if(t_AE_star.begin(), t_AE_star.end(), [&](double x){ return std::abs(x - val) < 1e-9; });
      int index = (it != t_AE_star.end()) ? std::distance(t_AE_star.begin(), it) : -1;
      lambda_M_plus_u[u] = index;
    }
    // Rcpp::Rcout << "[em_fit] lambda_M_plus_u: " << lambda_M_plus_u << std::endl;
  }

  void check() {
    // ---------- dimension & consistency checks ----------
    if ((int)t_D.size() != N_D) stop("length(t_D) != N_D");
    if ((int)t_DF.size() != (N_D + N_F)) stop("length(t_DF) != N_D + N_F");

    if (Q_j.ncol() != 2) stop("Q_j must have 2 columns");
    if (Q_j.nrow() != I) stop("nrow(Q_j) != I");

    if (LR_AB.ncol() != 2 || LR_E.ncol() != 2 || LR_F.ncol() != 2 || LR_ABEF.ncol() != 2)
      stop("All A_* matrices must have 2 columns");
    if (LR_AB.nrow() != N_AB) stop("nrow(LR_AB) != N_AB");
    if (LR_E.nrow() != N_E) stop("nrow(LR_E) != N_E");
    if (LR_F.nrow() != N_F) stop("nrow(LR_F) != N_F");
    if (LR_ABEF.nrow() != N_ABEF) stop("nrow(LR_ABEF) != N_ABEF");

    // vectors tied to A_*
    if ((int)L_AB.size() != N_AB || (int)R_AB.size() != N_AB || (int)t_AB.size() != N_AB)
      stop("L_AB, R_m, and t_AB must all have length N_AB");
    if ((int)L_E.size() != N_E || (int)t_E.size() != N_E)
      stop("L_E and t_E must have length N_E");
    if ((int)L_F.size() != N_F || (int)t_F.size() != N_F)
      stop("L_F and t_F must have length N_F");

    // t_AE_star / N_AE_star
    if ((int)t_AE_star.size() != N_AE_star) stop("length(t_AE_star) != N_AE_star");
    if ((int)r_A.size() != N_AE_star) stop("length(r_A) != N_AE_star");
    if (N_A_star < 0 || N_A_star > N_AE_star)
      stop("N_A_star must be in [0, N_AE_star]");
    if (N <= 0) stop("N must be > 0");

    // r_C / t_CE_star live on I_mark - I
    if (N_CE_star != (I_mark - I)) stop("N_CE_star must equal I_mark - I");
    if ((int)r_C.size() != N_CE_star) stop("length(r_C) != N_CE_star");
    if ((int)t_CE_star.size() != N_CE_star) stop("length(t_CE_star) != N_CE_star");

    // t_CE_star used when i in [I, I_mark)
    if ((int)t_CE_star.size() != (I_mark - I))
      stop("length(t_CE_star) != I_mark - I");

    // gamma_ji indexing uses offsets N_AB+l and N_ABE+l
    if (N_AB + N_E > N_ABEF)
      stop("N_AB + N_E must be <= N_ABEF for gamma_ji(i, N_AB + l)");
    if (N_ABE + N_F > N_ABEF)
      stop("N_ABE + N_F must be <= N_ABEF for gamma_ji(i, N_ABE + l)");
  }


  // ctor from named R list and calculates alpha and beta
  ModelData(const List& x) {
    // ints
    N_D  = as<int>(x["N_D"]);
    N_F  = as<int>(x["N_F"]);
    N_C = as<int>(x["N_C"]);
    N_E  = as<int>(x["N_E"]);
    N_A = as<int>(x["N_A"]);
    N_A_star = as<int>(x["N_A_star"]);
    N_AB  = as<int>(x["N_AB"]);
    N_ABE  = as<int>(x["N_ABE"]);
    N_AE_star  = as<int>(x["N_AE_star"]);
    N = as<int>(x["N"]);
    N_ABEF = as<int>(x["N_ABEF"]);
    I  = as<int>(x["I"]);
    N_CE_star  = as<int>(x["N_CE_star"]);
    I_mark = as<int>(x["I_mark"]);

    // scalars
    s_max = as<double>(x["s_max"]);
    R_max = as<double>(x["R_max"]);
    e_star_max = as<double>(x["e_star_max"]);

    // vectors
    t_D = as<NumericVector>(x["t_D"]);
    L_F = as<NumericVector>(x["L_F"]);
    t_F = as<NumericVector>(x["t_F"]);
    t_C = as<NumericVector>(x["t_C"]);
    L_E = as<NumericVector>(x["L_E"]);
    t_E = as<NumericVector>(x["t_E"]);
    L_AB = as<NumericVector>(x["L_AB"]);
    R_AB = as<NumericVector>(x["R_AB"]);
    t_AB = as<NumericVector>(x["t_AB"]);
    t_CE_star = as<NumericVector>(x["t_CE_star"]);
    t_AE_star = as<NumericVector>(x["t_AE_star"]);
    t_DF = as<NumericVector>(x["t_DF"]);

    // integer vectors
    r_C = Rcpp::as<arma::rowvec>(x["r_C"]);   // coerces integer
    r_A = as<IntegerVector>(x["r_A"]);

    // matrices (2 columns) — the important part
    LR_AB      = as<NumericMatrix>(x["LR_AB"]);
    LR_E      = as<NumericMatrix>(x["LR_E"]);
    LR_F      = as<NumericMatrix>(x["LR_F"]);
    LR_ABEF = as<NumericMatrix>(x["LR_ABEF"]);
    Q_j      = as<NumericMatrix>(x["Q_j"]);

    cal_alpha();
    cal_gamma();

    calc_lambda_M_plus_u();
    L = std::min(N_AB, std::min(N_E, std::min(N_F, N_D)));

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
  arma::rowvec lambda_u;  // length N_AE_star
  arma::rowvec z_j;       // length I_mark

  Rcpp::NumericVector t_sorted;   // sorted t_n
  Rcpp::NumericVector log_prefix; // prefix sum of log(1 - lambda_u) over sorted order
  Rcpp::NumericVector zero_prefix;// prefix count of zeros in (1 - lambda_u)

  // E-step pieces
  arma::mat P123_AB;
  arma::mat P123_D;
  arma::mat P123_E;
  arma::mat P123_F;
};

void setup_prod(const ModelData& md, Workspace& ws) {
  // Precompute sorted order and prefix arrays.
  // - Sort indices by md.t_AE_star
  // - t_sorted[i] = sorted t_AE_star
  // - log_prefix[i+1] = sum_{k<=i, lambda_k!=1} log(1 - lambda_k)
  // - zero_prefix[i+1] = count_{k<=i} [lambda_k == 1]
  const R_xlen_t N_AE_star = md.t_AE_star.size();
  if (ws.lambda_u.size() != N_AE_star) {
    stop("lambda_u length (%d) must match t_AE_star length (%d).",
         (int)ws.lambda_u.size(), (int)N_AE_star);
  }

  // build sorted index of t_AE_star
  std::vector<int> idx(N_AE_star);
  for (R_xlen_t i = 0; i < N_AE_star; ++i) idx[i] = static_cast<int>(i);
  std::sort(idx.begin(), idx.end(),
            [&](int a, int b){ return md.t_AE_star[a] < md.t_AE_star[b]; });

  ws.t_sorted  = NumericVector(N_AE_star);
  ws.log_prefix = NumericVector(N_AE_star + 1);
  ws.zero_prefix = NumericVector(N_AE_star + 1);

  ws.log_prefix[0]  = 0.0;
  ws.zero_prefix[0] = 0.0;

  // Treat lambda exactly equal to 1 as a zero factor (product becomes 0 if present in interval)
  const double TOL = 0.0; // set to, e.g., 1e-15 if you want tolerance
  for (R_xlen_t i = 0; i < N_AE_star; ++i) {
    int k = idx[i];
    const double t   = md.t_AE_star[k];
    const double lam = ws.lambda_u[k];

    ws.t_sorted[i] = t;

    const bool is_one = std::abs(1.0-lam) <= TOL;
    ws.zero_prefix[i + 1] = ws.zero_prefix[i] + (is_one ? 1.0 : 0.0);

    // Only add to log-prefix if not zero; if zero, keep running sum unchanged.
    if (!is_one) {
      const double one_minus = 1.0 - lam;
      if (one_minus <= 0.0) {
        // Defensive: log undefined. You can choose to stop() or handle differently.
        stop("Encountered 1 - lambda_u <= 0 at sorted index %d (lambda=%.17g).", (int)i, lam);
      }
      ws.log_prefix[i + 1] = ws.log_prefix[i] + std::log(one_minus);
    } else {
      ws.log_prefix[i + 1] = ws.log_prefix[i];
    }
  }
}

double evaluate(const Workspace& ws, double L, double R) {
  // Evaluate product_{t_n* in (L, R)} (1 - lambda_u).
  // We implement (L, R) with L open, R open by using:
  //   i = upper_bound(t_sorted, L)   -> first index with t > L
  //   j = lower_bound(t_sorted, R)   -> first index with t >= R
  // If any lambda==1 in (i, j-1), product is 0.
  // Else product = exp(log_prefix[j] - log_prefix[i]).
  const R_xlen_t N_AE_star = ws.t_sorted.size();
  if (ws.log_prefix.size() != N_AE_star + 1 || ws.zero_prefix.size() != N_AE_star + 1)
    stop("Workspace not initialized. Call setup_prod first.");

  if (R <= L || N_AE_star == 0) return 1.0;

  // Binary searches on sorted t
  const double* begin = &ws.t_sorted[0];
  const double* end   = begin + N_AE_star;

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
              << " N_AE_star=" << md.N_AE_star
              << " N_AB=" << md.N_AB
              << " N_D=" << md.N_D
              << " N_E=" << md.N_E
              << " N_F=" << md.N_F
              << " N_ABE=" << md.N_ABE
              << " t_AE_star=" << md.t_AE_star.size()
              << " N=" << md.N
              << std::endl;

  Rcpp::Rcout << "[calc_all] vector lengths:"
              << " z=" << ws.z_j.size()
              << " lambda=" << ws.lambda_u.size()
              << std::endl;

  // Input vectors used inside loops
  Rcpp::Rcout << "[calc_all] input vector lengths:"
              << " ws.z_j=" << ws.z_j.size()
              << " ws.lambda_u=" << ws.lambda_u.size()
              << " R_AB=" << md.R_AB.size()
              << " t_E=" << md.t_E.size()
              << " t_F=" << md.t_F.size()
              << " t_CE_star=" << md.t_CE_star.size()
              << " r_C=" << md.r_C.size()
              << " r_A=" << md.r_A.size()
              << std::endl;

  Rcpp::Rcout << "[calc_all] matrix shapes (rows x cols):"
  // using arma n_rows/n_cols for arma matrices and nrow()/ncol() for Rcpp matrices
  << " ws.P123_AB=" << ws.P123_AB.n_rows << "x" << ws.P123_AB.n_cols
  << " ws.P123_D=" << ws.P123_D.n_rows << "x" << ws.P123_D.n_cols
  << " ws.P123_E=" << ws.P123_E.n_rows << "x" << ws.P123_E.n_cols
  << " ws.P123_F=" << ws.P123_F.n_rows << "x" << ws.P123_F.n_cols
  << " gamma_ji=" << md.gamma_ji.nrow()  << "x" << md.gamma_ji.ncol()
  << " alpha_ji=" << md.alpha_ji.nrow() << "x" << md.alpha_ji.ncol()
  << " Q_j=" << md.Q_j.nrow() << "x" << md.Q_j.ncol()
  << std::endl;

  Rcpp::Rcout << "[calc_all] L (min of N_AB,N_E,N_F,N_D)=" << md.L << std::endl;


  // Loop ranges summary (avoid printing at every iteration)
  Rcpp::Rcout << "[calc_all] loop ranges:"
              << " i_main: 0.." << (md.I > 0 ? md.I - 1 : -1)
              << " i_extra(I_mark): " << md.I << ".." << (md.I_mark > 0 ? md.I_mark - 1 : -1)
              << " l_main: 0.." << (md.L > 0 ? md.L - 1 : -1)
              << " tails: N_AB:" << md.L << ".." << (md.N_AB > 0 ? md.N_AB - 1 : -1)
              << " N_D:" << md.L << ".." << (md.N_D > 0 ? md.N_D - 1 : -1)
              << " N_E:" << md.L << ".." << (md.N_E > 0 ? md.N_E - 1 : -1)
              << " N_F:" << md.L << ".." << (md.N_F > 0 ? md.N_F - 1 : -1)
              << std::endl;

}

void run_em_once(const ModelData& md, Workspace& ws) {
  // --- Resets temps ----------------------------------------------------------
  ws.P123_AB.zeros(md.N_AB, md.I);
  ws.P123_D.zeros(md.N_D, md.I_mark);
  ws.P123_E.zeros(md.N_E, md.I_mark);
  ws.P123_F.zeros(md.N_F, md.I_mark);

  setup_prod(md,ws);

  // --- Loops over I ----------------------------------------------------------
  for(int j = 0; j < md.I; ++j) {
    for (int i = 0; i < md.N_AB; ++i) {
      ws.P123_AB(i,j) = md.gamma_ji(j,i) ?
      ws.z_j[j] * evaluate(ws, md.Q_j(j,1), next_double(md.R_AB[i])) : 0;
    }
    for (int i = 0; i < md.N_D; ++i) {
      ws.P123_D(i,j) = md.alpha_ji(j,i) ?ws.z_j[j] : 0;
    }
    for (int i = 0; i < md.N_E; ++i) {
      ws.P123_E(i,j) = md.gamma_ji(j, md.N_AB + i) ? ws.lambda_u[md.lambda_M_plus_u[i]] *
        evaluate(ws, md.Q_j(j,1), md.t_E[i]) * ws.z_j[j] : 0 ;
    }
    for (int i = 0; i < md.N_F; ++i) {
      ws.P123_F(i,j) = (md.alpha_ji(j, md.N_D + i) ?
                            ws.z_j[j] : 0) +
                            (md.gamma_ji(j, md.N_ABE + i) ?
                            evaluate(ws, md.Q_j(j,1), next_double(md.t_F[i])) *
                            ws.z_j[j]: 0);
    }
  }

  // --- Extra for I_mark ------------------------------------------------------
  for(int j = md.I; j < md.I_mark; ++j) {
    for (int i = 0; i < md.N_D; ++i) {
      ws.P123_D(i,j) = md.alpha_ji(j,i) ?ws.z_j[j] : 0;
    }
    for (int i = 0; i < md.N_E; ++i) {
      ws.P123_E(i,j) = is_double_eq(md.t_E[i],md.t_CE_star[j-md.I]) ? ws.z_j[j] : 0;
    }
    for (int i = 0; i < md.N_F; ++i) {
      ws.P123_F(i,j) = md.alpha_ji(j, md.N_D + i) ?ws.z_j[j] : 0;
    }
  }

  // --- Normalizes rows -------------------------------------------------------
  ws.P123_AB = arma::normalise(ws.P123_AB, 1, 1);
  ws.P123_D = arma::normalise(ws.P123_D, 1, 1);
  ws.P123_E = arma::normalise(ws.P123_E, 1, 1);
  ws.P123_F = arma::normalise(ws.P123_F, 1, 1);

  double sum_rho_n;
  double sum_pi_n;
  double sum_pi_full_n;
  double sum_sigma_n;
  double d_n_part;

  // --- Loops over N_AE_star ----------------------------------------------------------
  for(int n = 0; n < md.N_AE_star; ++n) {
    sum_rho_n = 0;
    sum_pi_n = 0;
    sum_pi_full_n = 0;
    sum_sigma_n = 0;
    d_n_part = 0;

    d_n_part = n <= md.N_A_star ? md.r_A[n] : 0;
    for (int i = 0; i < md.N_AB; ++i) {
      // rho_mn
      if(md.t_AB[i] >= md.t_AE_star[n]) {
        for(int j = 0; j < md.I; ++j) {
          sum_rho_n +=
            md.L_AB[i] <= md.Q_j(j,0) && md.Q_j(j,1) < md.t_AE_star[n] ?
            ws.P123_AB(i,j) : 0;
        }
      }
    }
    for (int i = 0; i < md.N_E; ++i) {
      // pi_un
      if(md.t_E[i] >= md.t_AE_star[n]) {
        for(int j = 0; j < md.I; ++j) {
          sum_pi_n +=
            md.L_E[i] <= md.Q_j(j,0) && md.Q_j(j,1) < md.t_AE_star[n] ?
            ws.P123_E(i,j) : 0;

          sum_pi_full_n += is_double_eq(md.t_E[i], md.t_AE_star[n]) ? ws.P123_E(i,j) : 0;
        }
      }
    }
    for (int i = 0; i < md.N_F; ++i) {
      // sigma_cn
      if(md.t_F[i] >= md.t_AE_star[n]) {
        for(int j = 0; j < md.I; ++j) {
          sum_sigma_n +=
            md.L_F[i] <= md.Q_j(j,0) && md.Q_j(j,1) < md.t_AE_star[n] ?
            ws.P123_F(i,j) : 0;
        }
      }
    }

    // Rcpp::Rcout << "[em_fit] sum sum_pi_full_n " << sum_pi_full_n << std::endl;

    double demon = sum_rho_n + sum_pi_n + sum_sigma_n;
    //Rcpp::Rcout << "[em_fit] demon is " << demon << std::endl;

    ws.lambda_u[n] = (d_n_part + sum_pi_full_n)/demon;

    // if lambda_u[n] is larger than 1, we set it to 0.05
    if(ws.lambda_u[n] > 1.0){

     // Rcpp::Rcout << "[em_fit] Warning: lambda_u[" << n << "] > 1.0" << std::endl;
     // Rcpp::Rcout << "[em_fit] Corresponding to t_AE_star[" << n << "] = " << md.t_AE_star[n] << std::endl;
      ws.lambda_u[n] = 0.05;
    }

  }

  // --- Calculates new z ------------------------------------------------------
  arma::rowvec base = arma::sum(ws.P123_D, 0);
  base += arma::sum(ws.P123_E, 0);
  base += arma::sum(ws.P123_F, 0);
  arma::rowvec new_z = base;

  new_z.head(md.I) += arma::sum(ws.P123_AB, 0);
  if(md.r_C.n_elem > 0) {
    new_z.subvec(md.I, md.I_mark - 1) += md.r_C;
  }
  new_z /= md.N;

  // // normalize z_j to sum to 1
  // new_z /= arma::accu(new_z);

  ws.z_j = new_z;


  //Rcpp::Rcout << "[em_fit] eta is " << ws.P123_E << std::endl;
}

double calculate_likelihood(const ModelData& md, const Workspace& ws) {
  double loglik = 0.0;


  // Recall that gamma and alpha are logical matrices (indicators)

  // A_23{lambda} part
  for(int u = 0; u < md.N_A_star; ++u) {
    loglik += md.r_A[u] * std::log(ws.lambda_u[u]);
  }

  // AB part
  for(int i = 0; i < md.N_AB; ++i) {
    // sum over I_mark
    loglik += std::log(evaluate(ws, md.R_AB[i], md.t_AB[i]));

    double AB_likelihood = 0.0;
    // Only iterate up to md.I since gamma_ji is I x N_ABEF
    for(int j = 0; j < md.I; ++j) {
      AB_likelihood += md.gamma_ji(j,i) ? ws.z_j[j] * evaluate(ws, md.Q_j(j,1), next_double(md.R_AB[i])) : 0 ;
    }
    loglik += std::log(AB_likelihood);
  }

  // C part
  for(int k = md.I; k < md.I_mark; ++k) {
    loglik += md.r_C[k - md.I] * std::log(ws.z_j[k]);
  }

  // D part
  for(int i = 0; i < md.N_D; ++i) {
    double D_likelihood = 0.0;
    // sum over I_mark multiply alpha z_j
    for(int j = 0; j < md.I_mark; ++j) {
      D_likelihood += md.alpha_ji(j,i) ? ws.z_j[j] : 0;
    }
    loglik += std::log(D_likelihood);
  }
  for(int i = 0; i < md.N_E; ++i) {
    double E_likelihood = 0.0;

    for(int j = 0; j < md.I; ++j) {
      E_likelihood += md.gamma_ji(j, md.N_AB + i) ? ws.lambda_u[md.lambda_M_plus_u[i]] * ws.z_j[j] *
      evaluate(ws, md.Q_j(j,1), md.t_E[i]) : 0;
    }
    for(int j = md.I; j < md.I_mark; ++j) {
      E_likelihood += is_double_eq(md.t_E[i],md.t_CE_star[j - md.I]) ? ws.z_j[j] : 0;
    }

    loglik += std::log(E_likelihood);
  }
  for(int i = 0; i < md.N_F; ++i) {
    double F_likelihood = 0.0;
    for(int j = 0; j < md.I_mark; ++j) {
      F_likelihood += (md.alpha_ji(j, md.N_D + i) ?
                          ws.z_j[j] : 0);

      // Only access gamma_ji when j < md.I since gamma_ji is I x N_ABEF
      if (j < md.I) {
        F_likelihood += (md.gamma_ji(j, md.N_ABE + i) ?
                         evaluate(ws, md.Q_j(j,1), next_double(md.t_F[i])) * ws.z_j[j] : 0);
      }
    }
    loglik += std::log(F_likelihood);
  }
  return loglik;
}






// [[Rcpp::export]]
Rcpp::List em_fit(SEXP md_ptr,
                  Rcpp::Nullable<Rcpp::NumericVector> z_init = R_NilValue,
                  Rcpp::Nullable<Rcpp::NumericVector> lambda_init = R_NilValue,
                  int max_iter = 100,
                  double tol = 1e-3,
                  bool verbose = true,
                  bool eval_likelihood = false) {

  Rcpp::XPtr<ModelData> p(md_ptr);
  const ModelData& md = *p;
  bool converged = false;
  Workspace ws;
  std::vector<double> likelihoods;

  const arma::uword I_mark = static_cast<arma::uword>(md.I_mark);
  ws.z_j.set_size(I_mark);
  if (z_init.isNotNull()) {
    Rcpp::NumericVector z = z_init.get();
    if (static_cast<arma::uword>(z.size()) != I_mark)
      Rcpp::stop("length(z_init) must equal I_mark");
    for (arma::uword i = 0; i < I_mark; ++i) ws.z_j[i] = z[i];
  } else {
    ws.z_j.fill(1.0 / static_cast<double>(I_mark));
  }

  // make sure that z_j sums to 1
  ws.z_j /= arma::accu(ws.z_j);
  // print initial sum of z_j for debugging
  if (verbose) {
    double z_sum_init = arma::accu(ws.z_j);
    Rcpp::Rcout << "[em_fit] Initial sum(z_j) = " << z_sum_init << std::endl;
  }


  const arma::uword N_AE_star = static_cast<arma::uword>(md.t_AE_star.size());
  ws.lambda_u.set_size(N_AE_star);
  if (lambda_init.isNotNull()) {
    Rcpp::NumericVector lam = lambda_init.get();
    if (static_cast<arma::uword>(lam.size()) != N_AE_star)
      Rcpp::stop("length(lambda_init) must equal N_AE_star (length(t_AE_star))");
    for (arma::uword n = 0; n < N_AE_star; ++n) ws.lambda_u[n] = lam[n];
  } else {
    ws.lambda_u.fill(0.005);
  }

  double dz = 0.0, dl = 0.0;
  for (int iter = 0; iter < max_iter; ++iter) {
    arma::rowvec prev_lambda_n = ws.lambda_u; // copy
    arma::rowvec prev_z_i      = ws.z_j;      // copy
    
    run_em_once(md, ws);

    if(eval_likelihood) {
      likelihoods.push_back(calculate_likelihood(md, ws));
    }
    

    dz = 0.0; dl = 0.0;
    for (arma::uword i = 0; i < I_mark; ++i)
      dz = std::max(dz, std::abs(ws.z_j[i] - prev_z_i[i]));
    for (arma::uword n = 0; n < N_AE_star; ++n)
      dl = std::max(dl, std::abs(ws.lambda_u[n] - prev_lambda_n[n]));

    if (verbose) {
      if(iter % 10 == 0 || iter == max_iter - 1){
          Rcpp::Rcout << "iter " << (iter + 1)
                      << ": max|Δz|=" << dz
                      << " max|Δλ|=" << dl << "\n";
          // print sum of z_j for debugging
          double z_sum = arma::accu(ws.z_j);
          Rcpp::Rcout << "[em_fit] sum(z_j) = " << z_sum << std::endl;
      }
    }

    if (std::isfinite(dz) && std::isfinite(dl) && std::max(dz, dl) < tol) {
      if(std::max(dz, dl) < tol) {
        Rcpp::Rcout << "----------------------------------------------------\n";
        Rcpp::Rcout << "Convergence achieved: max|Δz| and max|Δλ| < tol (" << tol << "). After " << (iter + 1) << " iterations." << std::endl;
        converged = true;
      } else {
        Rcpp::Rcout << "----------------------------------------------------\n";
        Rcpp::Rcout << "Stopping: max|Δz| or max|Δλ| is NaN or Inf. After " << (iter + 1) << " iterations." << std::endl;
      }
      break;
    }
    if (iter == max_iter - 1) {
      Rcpp::Rcout << "----------------------------------------------------\n";
      Rcpp::Rcout << "Maximum iterations reached (" << max_iter << ")." << std::endl;
    }
  }


  // if (verbose) print_summary(md, ws);



  return Rcpp::List::create(
    Rcpp::_["z_j"]      = ws.z_j,
    Rcpp::_["lambda_u"] = ws.lambda_u,
    Rcpp::_["alpha_ji"] = md.alpha_ji,
    Rcpp::_["gamma_ji"]  = md.gamma_ji,
    Rcpp::_["converged"]    = converged,
    Rcpp::_["likelihoods"]  = likelihoods
  );
}
