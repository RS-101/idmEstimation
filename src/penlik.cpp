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
    // I spline matrices
    arma::mat V_healthy_i_spline_mat_12;
    arma::mat V_healthy_i_spline_mat_13;

    arma::mat T_obs_i_spline_mat_12;
    arma::mat T_obs_i_spline_mat_13;
    arma::mat T_obs_i_spline_mat_23;

    arma::mat V_0_i_spline_mat_12;
    arma::mat V_0_i_spline_mat_13;

    arma::mat grid_T_obs_i_spline_mat_12;
    arma::mat grid_T_obs_i_spline_mat_13;
    arma::mat grid_T_obs_i_spline_mat_23;

    arma::mat grid_V_ill_i_spline_mat_12;
    arma::mat grid_V_ill_i_spline_mat_13;
    arma::mat grid_V_ill_i_spline_mat_23;


    arma::mat T_obs_m_spline_mat_12;
    arma::mat T_obs_m_spline_mat_13;
    arma::mat T_obs_m_spline_mat_23;

    double dx_grid_T_obs;
    arma::mat grid_T_obs_m_spline_mat_12;
    arma::mat grid_T_obs_m_spline_mat_13;
    arma::mat grid_T_obs_m_spline_mat_23;

    double dx_grid_V_ill;
    arma::mat grid_V_ill_m_spline_mat_12;
    arma::mat grid_V_ill_m_spline_mat_13;
    arma::mat grid_V_ill_m_spline_mat_23;
    
    // Actual time values
    arma::vec T_obs_values;
    arma::vec V_healthy_values;
    arma::vec V_ill_values;  // Added for case 4

    // constructor of list with matrices
    PenlikModelData(const List& x) {
        V_healthy_i_spline_mat_12 = as<arma::mat>(x["V_healthy_i_spline_mat_12"]);
        V_healthy_i_spline_mat_13 = as<arma::mat>(x["V_healthy_i_spline_mat_13"]);

        T_obs_i_spline_mat_12 = as<arma::mat>(x["T_obs_i_spline_mat_12"]);
        T_obs_i_spline_mat_13 = as<arma::mat>(x["T_obs_i_spline_mat_13"]);
        T_obs_i_spline_mat_23 = as<arma::mat>(x["T_obs_i_spline_mat_23"]);
        V_0_i_spline_mat_12 = as<arma::mat>(x["V_0_i_spline_mat_12"]);
        V_0_i_spline_mat_13 = as<arma::mat>(x["V_0_i_spline_mat_13"]);
        grid_T_obs_i_spline_mat_12 = as<arma::mat>(x["grid_T_obs_i_spline_mat_12"]);
        grid_T_obs_i_spline_mat_13 = as<arma::mat>(x["grid_T_obs_i_spline_mat_13"]);
        grid_T_obs_i_spline_mat_23 = as<arma::mat>(x["grid_T_obs_i_spline_mat_23"]);
        grid_V_ill_i_spline_mat_12 = as<arma::mat>(x["grid_V_ill_i_spline_mat_12"]);
        grid_V_ill_i_spline_mat_13 = as<arma::mat>(x["grid_V_ill_i_spline_mat_13"]);
        grid_V_ill_i_spline_mat_23 = as<arma::mat>(x["grid_V_ill_i_spline_mat_23"]);
        T_obs_m_spline_mat_12 = as<arma::mat>(x["T_obs_m_spline_mat_12"]);
        T_obs_m_spline_mat_13 = as<arma::mat>(x["T_obs_m_spline_mat_13"]);
        T_obs_m_spline_mat_23 = as<arma::mat>(x["T_obs_m_spline_mat_23"]);
        dx_grid_T_obs = as<double>(x["dx_grid_T_obs"]);
        dx_grid_V_ill = as<double>(x["dx_grid_V_ill"]);
        grid_T_obs_m_spline_mat_12 = as<arma::mat>(x["grid_T_obs_m_spline_mat_12"]);
        grid_T_obs_m_spline_mat_13 = as<arma::mat>(x["grid_T_obs_m_spline_mat_13"]);
        grid_T_obs_m_spline_mat_23 = as<arma::mat>(x["grid_T_obs_m_spline_mat_23"]);
        grid_V_ill_m_spline_mat_12 = as<arma::mat>(x["grid_V_ill_m_spline_mat_12"]);
        grid_V_ill_m_spline_mat_13 = as<arma::mat>(x["grid_V_ill_m_spline_mat_13"]);
        grid_V_ill_m_spline_mat_23 = as<arma::mat>(x["grid_V_ill_m_spline_mat_23"]);
        
        // Load time values
        T_obs_values = as<arma::vec>(x["T_obs_values"]);
        V_healthy_values = as<arma::vec>(x["V_healthy_values"]);
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

    // R implementation data are all vectors and A12(T_obs) can be
    // computed as md.T_obs_i_spline_mat_12 %*% theta_12
    // We are calcuating the integral
    // int_V_healthy^T_obs exp(-A12(t) - A13(t)) * a12(t) * exp(A23(t)) dt
    // In the R code we use splinefun and trapezoidal rule to compute this integral


    // term1 <- exp(-A12(T_obs) - A13(T_obs))

    // factored_out <- exp(-A23(T_obs))

    // grid <- seq(min(V_healthy), max(T_obs), length.out = 250)

    // integrand <- exp(-A12(grid) - A13(grid)) * a12(grid) * exp(A23(grid))

    // dx <- diff(grid)
    // mid <- (integrand[-1] + integrand[-length(grid)]) / 2
    // integral_at_grid <- c(0, cumsum(dx * mid))

    // integral_fun <- splinefun(grid, integral_at_grid, method = "natural")

    // term2 <- factored_out * (integral_fun(T_obs) - integral_fun(V_healthy))

    // Calculate term_1 with correct negative sign
    arma::vec cum_hazard_12 = md.T_obs_i_spline_mat_12 * theta_12;
    arma::vec cum_hazard_13 = md.T_obs_i_spline_mat_13 * theta_13;
    arma::vec term_1 = arma::exp(-(cum_hazard_12 + cum_hazard_13));

    // Debug prints (only if we have enough elements)
    if (cum_hazard_12.n_elem >= 3) {
        Rcout << "C++ A12_T_obs first 3: " << cum_hazard_12.subvec(0,2).t();
        Rcout << "C++ A13_T_obs first 3: " << cum_hazard_13.subvec(0,2).t();
        Rcout << "C++ term_1 first 3: " << term_1.subvec(0,2).t();
    }


    // Include A23 term in factored_out
    arma::vec factored_out = arma::exp(-(md.T_obs_i_spline_mat_23 * theta_23));

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
    
    // Debug prints (only if we have enough elements)
    if (T_obs_indices.n_elem >= 3) {
        Rcout << "C++ T_obs_indices first 3: " << T_obs_indices.subvec(0,2).t();
        Rcout << "C++ V_healthy_indices first 3: " << V_healthy_indices.subvec(0,2).t();
    }
    
    arma::vec term_2 = factored_out %
        (integral.elem(T_obs_indices) -
         integral.elem(V_healthy_indices));

    double log_likelihood = arma::accu(arma::log(term_1 + term_2));


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
    // Likelihood: exp(-A12(T_obs) - A13(T_obs)) * a13(T_obs) + 
    //             exp(-A23(T_obs)) * a23(T_obs) * 
    //             integral from V_healthy to T_obs of exp(-A12(t) - A13(t)) * a12(t) * exp(A23(t)) dt
    
    // Term 1: Direct transition 1->3
    arma::vec cum_hazard_12 = md.T_obs_i_spline_mat_12 * theta_12;
    arma::vec cum_hazard_13 = md.T_obs_i_spline_mat_13 * theta_13;
    arma::vec hazard_13 = md.T_obs_m_spline_mat_13 * theta_13;
    
    arma::vec term_1 = arma::exp(-(cum_hazard_12 + cum_hazard_13)) % hazard_13;
    
    // Term 2: Transition 1->2->3
    arma::vec hazard_23 = md.T_obs_m_spline_mat_23 * theta_23;
    arma::vec cum_hazard_23 = md.T_obs_i_spline_mat_23 * theta_23;
    arma::vec factored_out = arma::exp(-cum_hazard_23) % hazard_23;
    
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
    
    arma::vec term_2 = factored_out % 
        (integral.elem(T_obs_indices) - integral.elem(V_healthy_indices));
    
    double log_likelihood = arma::accu(arma::log(term_1 + term_2));
    
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
    

}

// [[Rcpp::export]]
double calc_case_4_log_likelihood(
    SEXP md_ptr_case_4,
    const arma::vec& theta_12,
    const arma::vec& theta_13,
    const arma::vec& theta_23
) {


}

// // [[Rcpp::export]]
// double calc_penlik_log_likelihood(
//     SEXP md_ptr_case_1,
//     SEXP md_ptr_case_2,
//     SEXP md_ptr_case_3,
//     SEXP md_ptr_case_4,
//     const arma::vec& theta_12,
//     const arma::vec& theta_13,
//     const arma::vec& theta_23
// ) {
//     Rcpp::XPtr<PenlikModelData> p1(md_ptr_case_1);
//     const PenlikModelData& md_case_1 = *p1;

//     Rcpp::XPtr<PenlikModelData> p2(md_ptr_case_2);
//     const PenlikModelData& md_case_2 = *p2;

//     Rcpp::XPtr<PenlikModelData> p3(md_ptr_case_3);
//     const PenlikModelData& md_case_3 = *p3;

//     Rcpp::XPtr<PenlikModelData> p4(md_ptr_case_4);
//     const PenlikModelData& md_case_4 = *p4;

//     double log_likelihood = 0.0;



//     return log_likelihood;
// }



// void optimize_penlik() {


//     Rcpp::Numerical::optim_lbfgs();
// }
