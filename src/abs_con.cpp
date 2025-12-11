// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <algorithm>
#include <vector>
#include <cmath>
#include <cfloat>
#include <limits>
using namespace Rcpp;

struct CaseData {
    // ========================================================================
    // REQUIRED FIELDS - Common to all cases
    // ========================================================================

    // T_obs I-spline matrices (cumulative hazards at observation time)
    arma::mat T_obs_i_spline_mat_12;
    arma::mat T_obs_i_spline_mat_13;
    arma::mat T_obs_i_spline_mat_23;

    // Time values
    arma::vec T_obs_values;
    arma::vec V_healthy_values;

    // ========================================================================
    // OPTIONAL FIELDS - Case-specific
    // ========================================================================

    // M-spline matrices at T_obs (for hazard rates - cases 2 & 4)
    arma::mat T_obs_m_spline_mat_13;  // Case 2: direct death from healthy
    arma::mat T_obs_m_spline_mat_23;  // Cases 2 & 4: death from ill

    // Grid for integration from V_healthy to T_obs (cases 1 & 2)
    arma::mat grid_T_obs_i_spline_mat_12;
    arma::mat grid_T_obs_i_spline_mat_13;
    arma::mat grid_T_obs_i_spline_mat_23;
    arma::mat grid_T_obs_m_spline_mat_12;
    double dx_grid_T_obs;

    // Grid for integration from V_healthy to V_ill (cases 3 & 4)
    arma::mat grid_V_ill_i_spline_mat_12;
    arma::mat grid_V_ill_i_spline_mat_13;
    arma::mat grid_V_ill_i_spline_mat_23;
    arma::mat grid_V_ill_m_spline_mat_12;
    double dx_grid_V_ill;
    arma::vec V_ill_values;

    // Default constructor for empty cases
    CaseData() {}

    // Constructor from List (already exists)
    CaseData(const List& x) {
        // T_obs_i_spline_mat_12 is optional
        if (x.containsElementNamed("T_obs_i_spline_mat_12")) {
            T_obs_i_spline_mat_12 = as<arma::mat>(x["T_obs_i_spline_mat_12"]);
        }
        if (x.containsElementNamed("T_obs_i_spline_mat_13")) {
            T_obs_i_spline_mat_13 = as<arma::mat>(x["T_obs_i_spline_mat_13"]);
        }
        if (x.containsElementNamed("T_obs_i_spline_mat_23")) {
            T_obs_i_spline_mat_23 = as<arma::mat>(x["T_obs_i_spline_mat_23"]);
        }

        if (x.containsElementNamed("T_obs_values")) {
            T_obs_values = as<arma::vec>(x["T_obs_values"]);
        }
        
        if (x.containsElementNamed("V_healthy_values")) {
            V_healthy_values = as<arma::vec>(x["V_healthy_values"]);
        }

        // ====================================================================
        // OPTIONAL FIELDS - Cases 2 & 4
        // ====================================================================
        if (x.containsElementNamed("T_obs_m_spline_mat_13")) {
            T_obs_m_spline_mat_13 = as<arma::mat>(x["T_obs_m_spline_mat_13"]);
        }
        if (x.containsElementNamed("T_obs_m_spline_mat_23")) {
            T_obs_m_spline_mat_23 = as<arma::mat>(x["T_obs_m_spline_mat_23"]);
        }

        // ====================================================================
        // OPTIONAL FIELDS - Cases 1 & 2 (integration from V_healthy to T_obs)
        // ====================================================================
        if (x.containsElementNamed("grid_T_obs_i_spline_mat_12")) {
            grid_T_obs_i_spline_mat_12 = as<arma::mat>(x["grid_T_obs_i_spline_mat_12"]);
        }
        if (x.containsElementNamed("grid_T_obs_i_spline_mat_13")) {
            grid_T_obs_i_spline_mat_13 = as<arma::mat>(x["grid_T_obs_i_spline_mat_13"]);
        }
        if (x.containsElementNamed("grid_T_obs_i_spline_mat_23")) {
            grid_T_obs_i_spline_mat_23 = as<arma::mat>(x["grid_T_obs_i_spline_mat_23"]);
        }
        if (x.containsElementNamed("grid_T_obs_m_spline_mat_12")) {
            grid_T_obs_m_spline_mat_12 = as<arma::mat>(x["grid_T_obs_m_spline_mat_12"]);
        }
        if (x.containsElementNamed("dx_grid_T_obs")) {
            dx_grid_T_obs = as<double>(x["dx_grid_T_obs"]);
        }

        // ====================================================================
        // OPTIONAL FIELDS - Cases 3 & 4 (integration from V_healthy to V_ill)
        // ====================================================================
        if (x.containsElementNamed("grid_V_ill_i_spline_mat_12")) {
            grid_V_ill_i_spline_mat_12 = as<arma::mat>(x["grid_V_ill_i_spline_mat_12"]);
        }
        if (x.containsElementNamed("grid_V_ill_i_spline_mat_13")) {
            grid_V_ill_i_spline_mat_13 = as<arma::mat>(x["grid_V_ill_i_spline_mat_13"]);
        }
        if (x.containsElementNamed("grid_V_ill_i_spline_mat_23")) {
            grid_V_ill_i_spline_mat_23 = as<arma::mat>(x["grid_V_ill_i_spline_mat_23"]);
        }
        if (x.containsElementNamed("grid_V_ill_m_spline_mat_12")) {
            grid_V_ill_m_spline_mat_12 = as<arma::mat>(x["grid_V_ill_m_spline_mat_12"]);
        }
        if (x.containsElementNamed("dx_grid_V_ill")) {
            dx_grid_V_ill = as<double>(x["dx_grid_V_ill"]);
        }
        if (x.containsElementNamed("V_ill_values")) {
            V_ill_values = as<arma::vec>(x["V_ill_values"]);
        }
    }
};

// Full data structure
struct FullData {
    CaseData case_A;
    CaseData case_B;
    CaseData case_C;
    CaseData case_D;
    CaseData case_E;
    CaseData case_F;
    bool has_case_A = false;
    bool has_case_B = false;
    bool has_case_C = false;
    bool has_case_D = false;
    bool has_case_E = false;
    bool has_case_F = false;


    // Constructor taking all data at once
    FullData(const List& data) {
        // Initialize case-specific data
        if (data.containsElementNamed("case_A")) {
            case_A = CaseData(as<List>(data["case_A"]));
            has_case_A = true;
        }
        if (data.containsElementNamed("case_B")) {
            case_B = CaseData(as<List>(data["case_B"]));
            has_case_B = true;
        }
        if (data.containsElementNamed("case_C")) {
            case_C = CaseData(as<List>(data["case_C"]));
            has_case_C = true;
        }
        if (data.containsElementNamed("case_D")) {
            case_D = CaseData(as<List>(data["case_D"]));
            has_case_D = true;
        }
        if (data.containsElementNamed("case_E")) {
            case_E = CaseData(as<List>(data["case_E"]));
            has_case_E = true;
        }
        if (data.containsElementNamed("case_F")) {
            case_F = CaseData(as<List>(data["case_F"]));
            has_case_F = true;
        }
    }
};

// Full data structure
struct PenaltyConfig {
   
    arma::mat penalty_matrix_12;
    arma::mat penalty_matrix_13;
    arma::mat penalty_matrix_23;

    // Constructor taking all data at once
    PenaltyConfig(const List& data) {    
        penalty_matrix_12 = as<arma::mat>(data["penalty_matrix_12"]);
        penalty_matrix_13 = as<arma::mat>(data["penalty_matrix_13"]);
        penalty_matrix_23 = as<arma::mat>(data["penalty_matrix_23"]);
    }

};

// [[Rcpp::export]]
SEXP create_penalty_config(List penalty_data) {
    Rcpp::XPtr<PenaltyConfig> p(new PenaltyConfig(penalty_data), true);
    return p;
}

// [[Rcpp::export]]
SEXP create_full_data(List data) {
    Rcpp::XPtr<FullData> p(new FullData(data), true);
    return p;
}
double calc_case_A_log_likelihood(
    const CaseData& md,
    const arma::vec& theta_12,
    const arma::vec& theta_13,
    const arma::vec& theta_23
) {

    // Additional factor: α23(T) * exp(-A23(T))
    arma::vec A23_T_obs = md.T_obs_i_spline_mat_23 * theta_23;
    arma::vec a23_T_obs = md.T_obs_m_spline_mat_23 * theta_23;
    arma::vec exp_neg_A23_T_times_a23 = arma::exp(-A23_T_obs) % a23_T_obs;

  
    arma::vec integrand = arma::exp(
        -(md.grid_V_ill_i_spline_mat_12 * theta_12) -
        (md.grid_V_ill_i_spline_mat_13 * theta_13) +
        (md.grid_V_ill_i_spline_mat_23 * theta_23))
        % (md.grid_V_ill_m_spline_mat_12 * theta_12);

    // Trapezoidal rule
    arma::vec mid = (integrand.subvec(0, integrand.n_elem-2) +
                     integrand.subvec(1, integrand.n_elem-1)) / 2.0;

    arma::vec integral(mid.n_elem + 1);
    integral(0) = 0;
    integral.subvec(1, integral.n_elem-1) = arma::cumsum(mid * md.dx_grid_V_ill);

    // Get indices for V_healthy and V_ill in the grid
    double grid_min = arma::min(md.V_healthy_values);
    arma::uvec V_healthy_indices = arma::conv_to<arma::uvec>::from(
        arma::round((md.V_healthy_values - grid_min) / md.dx_grid_V_ill));
    arma::uvec V_ill_indices = arma::conv_to<arma::uvec>::from(
        arma::round((md.V_ill_values - grid_min) / md.dx_grid_V_ill));

    // Compute the integral from V_healthy to V_ill
    arma::vec integral_value = integral.elem(V_ill_indices) - integral.elem(V_healthy_indices);

    // Combine: exp(-A23(T)) * α23(T) * integral_value
    arma::vec likelihood = exp_neg_A23_T_times_a23 % integral_value;

    double log_likelihood = arma::accu(arma::log(likelihood));

    return log_likelihood;
}

double calc_case_B_log_likelihood(
    const CaseData& md,
    const arma::vec& theta_12,
    const arma::vec& theta_13,
    const arma::vec& theta_23
) {
   
    // For the integral, we need exp(-A23(T))
    arma::vec A23_T_obs = md.T_obs_i_spline_mat_23 * theta_23;
    arma::vec exp_neg_A23_T = arma::exp(-A23_T_obs);

    // Compute integrand over grid from V_healthy to V_ill:
    // exp(-A12(u) - A13(u)) * α12(u) * exp(-A23(T)) * exp(A23(u))
    // = exp(-A12(u) - A13(u)) * α12(u) * exp(A23(u) - A23(T))

    arma::vec integrand = arma::exp(
        -(md.grid_V_ill_i_spline_mat_12 * theta_12) -
        (md.grid_V_ill_i_spline_mat_13 * theta_13) +
        (md.grid_V_ill_i_spline_mat_23 * theta_23))
        % (md.grid_V_ill_m_spline_mat_12 * theta_12);

    // Trapezoidal rule
    arma::vec mid = (integrand.subvec(0, integrand.n_elem-2) +
                     integrand.subvec(1, integrand.n_elem-1)) / 2.0;

    arma::vec integral(mid.n_elem + 1);
    integral(0) = 0;
    integral.subvec(1, integral.n_elem-1) = arma::cumsum(mid * md.dx_grid_V_ill);

    // Get indices for V_healthy and V_ill in the grid
    double grid_min = arma::min(md.V_healthy_values);
    arma::uvec V_healthy_indices = arma::conv_to<arma::uvec>::from(
        arma::round((md.V_healthy_values - grid_min) / md.dx_grid_V_ill));
    arma::uvec V_ill_indices = arma::conv_to<arma::uvec>::from(
        arma::round((md.V_ill_values - grid_min) / md.dx_grid_V_ill));

    // Compute the integral from V_healthy to V_ill
    arma::vec integral_value = integral.elem(V_ill_indices) - integral.elem(V_healthy_indices);

    // Combine: factored_out * exp(-A23(T)) * integral_value
    arma::vec likelihood = exp_neg_A23_T % integral_value;

    double log_likelihood = arma::accu(arma::log(likelihood));

    return log_likelihood;
}


double calc_case_C_log_likelihood(
    const CaseData& md,
    const arma::vec& theta_12,
    const arma::vec& theta_13,
    const arma::vec& theta_23
) {

    // Additional factor: α23(T) * exp(-A23(T))
    arma::vec a13_T_obs = md.T_obs_m_spline_mat_13 * theta_13;
    arma::vec A12_T_obs = md.T_obs_i_spline_mat_12 * theta_12;
    arma::vec A13_T_obs = md.T_obs_i_spline_mat_13 * theta_13;

    arma::vec likelihood = arma::exp(-A12_T_obs-A13_T_obs) % a13_T_obs;

    double log_likelihood = arma::accu(arma::log(likelihood));

    return log_likelihood;
}
double calc_case_D_log_likelihood(
    const CaseData& md,
    const arma::vec& theta_12,
    const arma::vec& theta_13,
    const arma::vec& theta_23
) {
    arma::vec A12_T_obs = md.T_obs_i_spline_mat_12 * theta_12;
    arma::vec A13_T_obs = md.T_obs_i_spline_mat_13 * theta_13;

    arma::vec likelihood = arma::exp(-A12_T_obs-A13_T_obs);

    double log_likelihood = arma::accu(arma::log(likelihood));

    return log_likelihood;
}


double calc_case_E_log_likelihood(
    const CaseData& md,
    const arma::vec& theta_12,
    const arma::vec& theta_13,
    const arma::vec& theta_23
) {
    // Case E: Died (state 3) at T_obs, was healthy at V_healthy
    // Likelihood: 1/exp(-A12(V_0)-A13(V_0))* (exp(-A12(T_obs) - A13(T_obs)) * a13(T_obs) +
    //             exp(-A23(T_obs)) * a23(T_obs) *
    //             integral from V_healthy to T_obs of exp(-A12(t) - A13(t)) * a12(t) * exp(A23(t)) dt)


    // Term 1: Direct transition 1->3
    arma::vec cum_hazard_12 = md.T_obs_i_spline_mat_12 * theta_12;
    arma::vec cum_hazard_13 = md.T_obs_i_spline_mat_13 * theta_13;
    arma::vec hazard_13 = md.T_obs_m_spline_mat_13 * theta_13;

    arma::vec term_1 = arma::exp(-(cum_hazard_12 + cum_hazard_13)) % hazard_13;

    // Term 2: Transition 1->2->3
    arma::vec hazard_23 = md.T_obs_m_spline_mat_23 * theta_23;
    arma::vec cum_hazard_23 = md.T_obs_i_spline_mat_23 * theta_23;
    arma::vec exp_A23_T_obs = arma::exp(-cum_hazard_23) % hazard_23;

    // Compute integrand: exp(-A12(t) - A13(t)) * a12(t) * exp(A23(t))
    arma::vec integrand = arma::exp(
        -(md.grid_T_obs_i_spline_mat_12 * theta_12) -
        (md.grid_T_obs_i_spline_mat_13 * theta_13))
        % (md.grid_T_obs_m_spline_mat_12 * theta_12)
        % arma::exp(md.grid_T_obs_i_spline_mat_23 * theta_23);

    // Trapezoidal rule
    arma::vec mid = (integrand.subvec(0, integrand.n_elem-2) +
                     integrand.subvec(1, integrand.n_elem-1)) / 2.0;

    arma::vec integral(mid.n_elem + 1);
    integral(0) = 0;
    integral.subvec(1, integral.n_elem-1) = arma::cumsum(mid * md.dx_grid_T_obs);

    // Get indices
    double grid_min = arma::min(md.V_healthy_values);
    arma::uvec T_obs_indices = arma::conv_to<arma::uvec>::from(
        arma::round((md.T_obs_values - grid_min) / md.dx_grid_T_obs));
    arma::uvec V_healthy_indices = arma::conv_to<arma::uvec>::from(
        arma::round((md.V_healthy_values - grid_min) / md.dx_grid_T_obs));

    arma::vec term_2 = exp_A23_T_obs %
        (integral.elem(T_obs_indices) - integral.elem(V_healthy_indices));

    double log_likelihood = arma::accu(arma::log((term_1 + term_2)));

    return log_likelihood;
}


double calc_case_F_log_likelihood(
    const CaseData& md,
    const arma::vec& theta_12,
    const arma::vec& theta_13,
    const arma::vec& theta_23
) {
    // Case F: in state 1 or 2 at T_obs, was healthy at V_healthy
    // Likelihood:(exp(-A12(T_obs) - A13(T_obs)) +
    //             exp(-A23(T_obs)) *
    //             integral from V_healthy to T_obs of exp(-A12(t) - A13(t)) * a12(t) * exp(A23(t)) dt)

    // Factored out term: 1/exp(-A12(V_0) - A13(V_0))
    // Calculate term_1 with correct negative sign
    arma::vec cum_hazard_12 = md.T_obs_i_spline_mat_12 * theta_12;
    arma::vec cum_hazard_13 = md.T_obs_i_spline_mat_13 * theta_13;
    arma::vec term_1 = arma::exp(-(cum_hazard_12 + cum_hazard_13));

    // Include A23 term in factored_out
    arma::vec exp_A23_T_obs = arma::exp(-(md.T_obs_i_spline_mat_23 * theta_23));

    // Compute integrand
    arma::vec integrand = arma::exp(
        - (md.grid_T_obs_i_spline_mat_12 * theta_12)
        - (md.grid_T_obs_i_spline_mat_13 * theta_13))
        % (md.grid_T_obs_m_spline_mat_12 * theta_12).eval()
        % exp(md.grid_T_obs_i_spline_mat_23 * theta_23);

    // Compute midpoints for trapezoidal rule
    arma::vec mid = (integrand.subvec(0, integrand.n_elem-2) +
                        integrand.subvec(1, integrand.n_elem-1)) / 2.0;


    // Compute cumulative integral with leading 0
    arma::vec integral(mid.n_elem + 1);
    integral(0) = 0;
    integral.subvec(1, integral.n_elem-1) = arma::cumsum(mid * md.dx_grid_T_obs);

    // Calculate proper grid indices based on actual time values
    double grid_min = arma::min(md.V_healthy_values);
    arma::uvec T_obs_indices = arma::conv_to<arma::uvec>::from(
        arma::round((md.T_obs_values - grid_min) / md.dx_grid_T_obs));

    arma::uvec V_healthy_indices = arma::conv_to<arma::uvec>::from(
        arma::round((md.V_healthy_values - grid_min) / md.dx_grid_T_obs));


    arma::vec term_2 = exp_A23_T_obs %
        (integral.elem(T_obs_indices) -
         integral.elem(V_healthy_indices));

    double log_likelihood = arma::accu(arma::log(term_1 + term_2));


    return log_likelihood;
}

// [[Rcpp::export]]
double calc_log_likelihood(
    SEXP full_data_ptr,
    const arma::vec& theta_12,
    const arma::vec& theta_13,
    const arma::vec& theta_23
) {
    Rcpp::XPtr<FullData> pointer(full_data_ptr);
    const FullData& mddata = *pointer;

    double log_likelihood = 0.0;

    // calculate log-likelihood contributions from each case

    // check if

    if (mddata.has_case_A ) {
        log_likelihood += calc_case_A_log_likelihood(
            mddata.case_A, theta_12, theta_13, theta_23);
    }
    if (mddata.has_case_B ) {
        log_likelihood += calc_case_B_log_likelihood(
            mddata.case_B, theta_12, theta_13, theta_23);
    }
    if (mddata.has_case_C ) {
        log_likelihood += calc_case_C_log_likelihood(
            mddata.case_C, theta_12, theta_13, theta_23);
    }
    if (mddata.has_case_D ) {
        log_likelihood += calc_case_D_log_likelihood(
            mddata.case_D, theta_12, theta_13, theta_23);
    }
    if (mddata.has_case_E ) {
        log_likelihood += calc_case_E_log_likelihood(
            mddata.case_E, theta_12, theta_13, theta_23);
    }
    if (mddata.has_case_F ) {
        log_likelihood += calc_case_F_log_likelihood(
            mddata.case_F, theta_12, theta_13, theta_23);
    }
    return log_likelihood;
}

// [[Rcpp::export]]
double calc_full_penalty(
    SEXP full_data_ptr,
    const arma::vec& theta_12,
    const arma::vec& theta_13,
    const arma::vec& theta_23,
    double kappa_12,
    double kappa_13,
    double kappa_23

) {
    Rcpp::XPtr<PenaltyConfig> pointer(full_data_ptr);

    const PenaltyConfig& penalty_config = *pointer;

    // Penalty terms
    double penalty_12 = kappa_12 * arma::as_scalar(
        theta_12.t() * penalty_config.penalty_matrix_12 * theta_12);
    double penalty_13 = kappa_13 * arma::as_scalar(
        theta_13.t() * penalty_config.penalty_matrix_13 * theta_13);
    double penalty_23 = kappa_23 * arma::as_scalar(
        theta_23.t() * penalty_config.penalty_matrix_23 * theta_23);
    return penalty_12 + penalty_13 + penalty_23;
}

// [[Rcpp::export]]
List calc_penalized_log_likelihood(
    SEXP full_data_ptr,
    SEXP penalty_config_ptr,
    const arma::vec& theta_12,
    const arma::vec& theta_13,
    const arma::vec& theta_23,
    double kappa_12,
    double kappa_13,
    double kappa_23
) {
    double log_likelihood = calc_log_likelihood(
        full_data_ptr, theta_12, theta_13, theta_23);

    double full_penalty = calc_full_penalty(
        penalty_config_ptr, theta_12, theta_13, theta_23,
        kappa_12, kappa_13, kappa_23);

    double penalized_log_likelihood = log_likelihood - full_penalty;

    return List::create(
        Named("penalized_log_likelihood") = penalized_log_likelihood,
        Named("log_likelihood") = log_likelihood,
        Named("penalty") = full_penalty
    );
}
