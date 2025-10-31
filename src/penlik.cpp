// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <RcppNumerical.h>
#include <RcppEigen.h>
#include <algorithm>
#include <vector>
#include <cmath>
#include <cfloat>
#include <limits>
using namespace Rcpp;

// Penalized likelihood calculations for illness-death model
// 
// Four observation scenarios (cases):
// Case 1: Individual is censored (state 1) at T_obs, was healthy at V_healthy
// Case 2: Individual died (state 3) at T_obs, was healthy at V_healthy
// Case 3: Individual is ill (state 2) at T_obs, was healthy at V_healthy  
// Case 4: Individual died (state 3) at T_obs, was ill at V_ill
//
// Transitions:
// 1->2: healthy to ill (hazard alpha_12, cumulative hazard A12)
// 1->3: healthy to death (hazard alpha_13, cumulative hazard A13)
// 2->3: ill to death (hazard alpha_23, cumulative hazard A23)

struct PenlikModelData {
    // ========================================================================
    // REQUIRED FIELDS - Common to all cases
    // ========================================================================
    
    // V_0 splines (for factored out normalization term)
    arma::mat V_0_i_spline_mat_12;
    arma::mat V_0_i_spline_mat_13;
    
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

    // Constructor - only requires common fields, rest are optional
    PenlikModelData(const List& x) {
        // ====================================================================
        // REQUIRED FIELDS
        // ====================================================================
        V_0_i_spline_mat_12 = as<arma::mat>(x["V_0_i_spline_mat_12"]);
        V_0_i_spline_mat_13 = as<arma::mat>(x["V_0_i_spline_mat_13"]);
        
        T_obs_i_spline_mat_12 = as<arma::mat>(x["T_obs_i_spline_mat_12"]);
        T_obs_i_spline_mat_13 = as<arma::mat>(x["T_obs_i_spline_mat_13"]);
        T_obs_i_spline_mat_23 = as<arma::mat>(x["T_obs_i_spline_mat_23"]);
        
        T_obs_values = as<arma::vec>(x["T_obs_values"]);
        V_healthy_values = as<arma::vec>(x["V_healthy_values"]);
        
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

// [[Rcpp::export]]
SEXP create_penlik_model_data(List x) {
    PenlikModelData* md = new PenlikModelData(x);
    Rcpp::XPtr<PenlikModelData> p(md, true);
    return p;
}

// [[Rcpp::export]]
double calc_case_1_log_likelihood(
    SEXP md_ptr_case_1,
    const arma::vec& theta_12,
    const arma::vec& theta_13,
    const arma::vec& theta_23
) {
    Rcpp::XPtr<PenlikModelData> p1(md_ptr_case_1);
    const PenlikModelData& md = *p1;
    // case 1: healthy at T_obs

    // Factored out term: 1/exp(-A12(V_0) - A13(V_0))
    arma::vec A12_V_0 = md.V_0_i_spline_mat_12 * theta_12;
    arma::vec A13_V_0 = md.V_0_i_spline_mat_13 * theta_13;
    arma::vec factored_out = arma::exp(A12_V_0 + A13_V_0);

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

    double log_likelihood = arma::accu(arma::log(factored_out % (term_1 + term_2)));


    return log_likelihood;
}


// [[Rcpp::export]]
double calc_case_2_log_likelihood(
    SEXP md_ptr_case_2,
    const arma::vec& theta_12,
    const arma::vec& theta_13,
    const arma::vec& theta_23
) {
    Rcpp::XPtr<PenlikModelData> p2(md_ptr_case_2);
    const PenlikModelData& md = *p2;
    
    // Case 2: Died (state 3) at T_obs, was healthy at V_healthy
    // Likelihood: 1/exp(-A12(V_0)-A13(V_0))* (exp(-A12(T_obs) - A13(T_obs)) * a13(T_obs) + 
    //             exp(-A23(T_obs)) * a23(T_obs) * 
    //             integral from V_healthy to T_obs of exp(-A12(t) - A13(t)) * a12(t) * exp(A23(t)) dt)
    

    // Factored out term: 1/exp(-A12(V_0) - A13(V_0))
    arma::vec A12_V_0 = md.V_0_i_spline_mat_12 * theta_12;
    arma::vec A13_V_0 = md.V_0_i_spline_mat_13 * theta_13;
    arma::vec factored_out = arma::exp(A12_V_0 + A13_V_0);

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

    double log_likelihood = arma::accu(arma::log(factored_out % (term_1 + term_2)));

    return log_likelihood;
}

// [[Rcpp::export]]
double calc_case_3_log_likelihood(
    SEXP md_ptr_case_3,
    const arma::vec& theta_12,
    const arma::vec& theta_13,
    const arma::vec& theta_23
) {
    Rcpp::XPtr<PenlikModelData> p3(md_ptr_case_3);
    const PenlikModelData& md = *p3;
    
    // Case 3: Subject is healthy at V_k (V_healthy), ill at V_{k+1} (V_ill), and still alive at T
    // From image formula (iii):
    // L = (1/exp(-A01(V_k) - A02(V_k))) * 
    //     ∫_{V_k}^{V_{k+1}} exp(-A01(u) - A02(u)) * α01(u) * exp(-A12(T) - A12(u)) du
    //
    // Translating to our notation (01→12, 02→13, 12→23):
    // L = (1/exp(-A12(V_0) - A13(V_0))) * 
    //     ∫_{V_healthy}^{V_ill} exp(-A12(u) - A13(u)) * α12(u) * exp(-A23(T) - A23(u)) du
    
    // Factored out term: 1/exp(-A12(V_0) - A13(V_0))
    arma::vec A12_V_0 = md.V_0_i_spline_mat_12 * theta_12;
    arma::vec A13_V_0 = md.V_0_i_spline_mat_13 * theta_13;
    arma::vec factored_out = arma::exp(A12_V_0 + A13_V_0);
    
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
    arma::vec likelihood = factored_out % exp_neg_A23_T % integral_value;
    
    double log_likelihood = arma::accu(arma::log(likelihood));
    
    return log_likelihood;
}

// [[Rcpp::export]]
double calc_case_4_log_likelihood(
    SEXP md_ptr_case_4,
    const arma::vec& theta_12,
    const arma::vec& theta_13,
    const arma::vec& theta_23
) {
    Rcpp::XPtr<PenlikModelData> p4(md_ptr_case_4);
    const PenlikModelData& md = *p4;
    
    // Case 4: Subject is healthy at V_k (V_healthy), ill at V_{k+1} (V_ill), and dies at T
    // From image formula (iv):
    // L = (1/exp(-A01(V_k) - A02(V_k))) * 
    //     ∫_{V_k}^{V_{k+1}} exp(-A01(u) - A02(u)) * α01(u) * exp(-A12(T) - A12(u)) * α12(T) du
    //
    // Translating to our notation (01→12, 02→13, 12→23):
    // L = (1/exp(-A12(V_0) - A13(V_0))) * 
    //     ∫_{V_healthy}^{V_ill} exp(-A12(u) - A13(u)) * α12(u) * exp(-A23(T) - A23(u)) * α23(T) du
    
    // Factored out term: 1/exp(-A12(V_0) - A13(V_0))
    arma::vec A12_V_0 = md.V_0_i_spline_mat_12 * theta_12;
    arma::vec A13_V_0 = md.V_0_i_spline_mat_13 * theta_13;
    arma::vec factored_out = arma::exp(A12_V_0 + A13_V_0);
    
    // Additional factor: α23(T) * exp(-A23(T))
    arma::vec A23_T_obs = md.T_obs_i_spline_mat_23 * theta_23;
    arma::vec a23_T_obs = md.T_obs_m_spline_mat_23 * theta_23;
    arma::vec exp_neg_A23_T_times_a23 = arma::exp(-A23_T_obs) % a23_T_obs;
    
    // Compute integrand over grid from V_healthy to V_ill:
    // exp(-A12(u) - A13(u)) * α12(u) * exp(A23(u) - A23(T)) * α23(T)
    // Note: α23(T) is constant w.r.t. integration variable u
    
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
    
    // Combine: factored_out * exp(-A23(T)) * α23(T) * integral_value
    arma::vec likelihood = factored_out % exp_neg_A23_T_times_a23 % integral_value;
    
    double log_likelihood = arma::accu(arma::log(likelihood));
    
    return log_likelihood;
}

// [[Rcpp::export]]
double calc_penlik_log_likelihood(
    SEXP md_ptr_case_1,
    SEXP md_ptr_case_2,
    SEXP md_ptr_case_3,
    SEXP md_ptr_case_4,
    const arma::vec& theta_12,
    const arma::vec& theta_13,
    const arma::vec& theta_23
) {
    Rcpp::XPtr<PenlikModelData> p1(md_ptr_case_1);
    const PenlikModelData& md_case_1 = *p1;

    Rcpp::XPtr<PenlikModelData> p2(md_ptr_case_2);
    const PenlikModelData& md_case_2 = *p2;

    Rcpp::XPtr<PenlikModelData> p3(md_ptr_case_3);
    const PenlikModelData& md_case_3 = *p3;

    Rcpp::XPtr<PenlikModelData> p4(md_ptr_case_4);
    const PenlikModelData& md_case_4 = *p4;

    double log_likelihood = 0.0;

    // calculate log-likelihood contributions from each case
    if (md_case_1.T_obs_values.n_elem > 0) {
        log_likelihood += calc_case_1_log_likelihood(
            md_ptr_case_1, theta_12, theta_13, theta_23);
    }
    if (md_case_2.T_obs_values.n_elem > 0) {
        log_likelihood += calc_case_2_log_likelihood(
            md_ptr_case_2, theta_12, theta_13, theta_23);
    }
    if (md_case_3.T_obs_values.n_elem > 0) {
        log_likelihood += calc_case_3_log_likelihood(
            md_ptr_case_3, theta_12, theta_13, theta_23);
    }
    if (md_case_4.T_obs_values.n_elem > 0) {
        log_likelihood += calc_case_4_log_likelihood(
            md_ptr_case_4, theta_12, theta_13, theta_23);
    }


    return log_likelihood;
}

// Objective function class for LBFGS optimization
class PenlikObjective : public Numer::MFuncGrad
{
private:
    SEXP md_ptr_case_1;
    SEXP md_ptr_case_2;
    SEXP md_ptr_case_3;
    SEXP md_ptr_case_4;
    int n_theta_12;
    int n_theta_13;
    int n_theta_23;
    
public:
    PenlikObjective(
        SEXP md1, SEXP md2, SEXP md3, SEXP md4,
        int n12, int n13, int n23
    ) : md_ptr_case_1(md1), md_ptr_case_2(md2), 
        md_ptr_case_3(md3), md_ptr_case_4(md4),
        n_theta_12(n12), n_theta_13(n13), n_theta_23(n23) {}
    
    // Evaluate negative log-likelihood (we minimize, so negate)
    double f_grad(Numer::Constvec& theta_all, Numer::Refvec grad) override
    {
        // Extract theta vectors for each transition (copy data, don't reference)
        arma::vec theta_12(n_theta_12);
        arma::vec theta_13(n_theta_13);
        arma::vec theta_23(n_theta_23);
        
        for (int i = 0; i < n_theta_12; i++) {
            theta_12(i) = theta_all(i);
        }
        for (int i = 0; i < n_theta_13; i++) {
            theta_13(i) = theta_all(n_theta_12 + i);
        }
        for (int i = 0; i < n_theta_23; i++) {
            theta_23(i) = theta_all(n_theta_12 + n_theta_13 + i);
        }
        
        // Calculate log-likelihood
        double log_lik = calc_penlik_log_likelihood(
            md_ptr_case_1, md_ptr_case_2, md_ptr_case_3, md_ptr_case_4,
            theta_12, theta_13, theta_23
        );
        
        // Compute gradient using finite differences
        const double epsilon = 1e-7;
        int n_params = theta_all.size();
        
        for (int i = 0; i < n_params; i++) {
            Eigen::VectorXd theta_plus = theta_all;
            theta_plus(i) += epsilon;
            
            arma::vec theta_12_plus(n_theta_12);
            arma::vec theta_13_plus(n_theta_13);
            arma::vec theta_23_plus(n_theta_23);
            
            for (int j = 0; j < n_theta_12; j++) {
                theta_12_plus(j) = theta_plus(j);
            }
            for (int j = 0; j < n_theta_13; j++) {
                theta_13_plus(j) = theta_plus(n_theta_12 + j);
            }
            for (int j = 0; j < n_theta_23; j++) {
                theta_23_plus(j) = theta_plus(n_theta_12 + n_theta_13 + j);
            }
            
            double log_lik_plus = calc_penlik_log_likelihood(
                md_ptr_case_1, md_ptr_case_2, md_ptr_case_3, md_ptr_case_4,
                theta_12_plus, theta_13_plus, theta_23_plus
            );
            
            // Gradient is negative because we return negative objective
            grad(i) = -(log_lik_plus - log_lik) / epsilon;
        }
        
        // Return negative because LBFGS minimizes
        return -log_lik;
    }
};

// [[Rcpp::export]]
List optimize_penlik(
    SEXP md_ptr_case_1,
    SEXP md_ptr_case_2,
    SEXP md_ptr_case_3,
    SEXP md_ptr_case_4,
    const arma::vec& theta_12_init,
    const arma::vec& theta_13_init,
    const arma::vec& theta_23_init,
    int max_iter = 1000,
    double epsilon = 1e-6
) {
    // Get dimensions
    int n_theta_12 = theta_12_init.n_elem;
    int n_theta_13 = theta_13_init.n_elem;
    int n_theta_23 = theta_23_init.n_elem;
    int n_params = n_theta_12 + n_theta_13 + n_theta_23;
    
    // Combine initial parameters into single vector
    Eigen::VectorXd theta_init(n_params);
    for (int i = 0; i < n_theta_12; i++) {
        theta_init(i) = theta_12_init(i);
    }
    for (int i = 0; i < n_theta_13; i++) {
        theta_init(n_theta_12 + i) = theta_13_init(i);
    }
    for (int i = 0; i < n_theta_23; i++) {
        theta_init(n_theta_12 + n_theta_13 + i) = theta_23_init(i);
    }
    
    // Create objective function
    PenlikObjective objective(
        md_ptr_case_1, md_ptr_case_2, md_ptr_case_3, md_ptr_case_4,
        n_theta_12, n_theta_13, n_theta_23
    );
    
    // Set up LBFGS parameters
    double fopt;
    int status = Numer::optim_lbfgs(objective, theta_init, fopt, max_iter, epsilon);
    
    // Extract optimized parameters
    arma::vec theta_12_opt(n_theta_12);
    arma::vec theta_13_opt(n_theta_13);
    arma::vec theta_23_opt(n_theta_23);
    
    for (int i = 0; i < n_theta_12; i++) {
        theta_12_opt(i) = theta_init(i);
    }
    for (int i = 0; i < n_theta_13; i++) {
        theta_13_opt(i) = theta_init(n_theta_12 + i);
    }
    for (int i = 0; i < n_theta_23; i++) {
        theta_23_opt(i) = theta_init(n_theta_12 + n_theta_13 + i);
    }
    
    // Calculate final log-likelihood
    double final_log_lik = calc_penlik_log_likelihood(
        md_ptr_case_1, md_ptr_case_2, md_ptr_case_3, md_ptr_case_4,
        theta_12_opt, theta_13_opt, theta_23_opt
    );
    
    // Return results
    return List::create(
        Named("theta_12") = theta_12_opt,
        Named("theta_13") = theta_13_opt,
        Named("theta_23") = theta_23_opt,
        Named("log_likelihood") = final_log_lik,
        Named("convergence") = status,
        Named("iterations") = max_iter
    );
}
