// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <algorithm>
#include <vector>
#include <cmath>
#include <cfloat>
#include <limits>
using namespace Rcpp;

// Piecewise constant hazard calculations for illness-death model
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

// Case-specific data structure
struct PwcModelDataCaseSpecific {
    arma::mat V_0_values;
    arma::vec V_healthy_values;
    arma::vec V_ill_values; 
    arma::vec T_obs_values;

    arma::vec knots_12;
    arma::vec knots_13;
    arma::vec knots_23;

    PwcModelDataCaseSpecific() {}

    // Constructor from R list
    PwcModelDataCaseSpecific(const List& data) {
        V_0_values = as<arma::mat>(data["V_0_values"]);
        V_healthy_values = as<arma::vec>(data["V_healthy_values"]);
        V_ill_values = as<arma::vec>(data["V_ill_values"]);
        T_obs_values = as<arma::vec>(data["T_obs_values"]);

        knots_12 = as<arma::vec>(data["knots_12"]);
        knots_13 = as<arma::vec>(data["knots_13"]);
        knots_23 = as<arma::vec>(data["knots_23"]);
    }
};


// Full data structure
struct PwcModelData {
    PwcModelDataCaseSpecific case_1;
    PwcModelDataCaseSpecific case_2;
    PwcModelDataCaseSpecific case_3;
    PwcModelDataCaseSpecific case_4;


    // Constructor taking all data at once
    PwcModelData(const List& data) {
        // Initialize case-specific data
        if (data.containsElementNamed("case_1")) {
            case_1 = PwcModelDataCaseSpecific(as<List>(data["case_1"]));
        }
        if (data.containsElementNamed("case_2")) {
            case_2 = PwcModelDataCaseSpecific(as<List>(data["case_2"]));
        }
        if (data.containsElementNamed("case_3")) {
            case_3 = PwcModelDataCaseSpecific(as<List>(data["case_3"]));
        }
        if (data.containsElementNamed("case_4")) {
            case_4 = PwcModelDataCaseSpecific(as<List>(data["case_4"]));
        }
    }
};

// [[Rcpp::export]]
SEXP create_pwc_model_data(List data) {
    PwcModelData* md = new PwcModelData(data);
    Rcpp::XPtr<PwcModelData> p(md, true);
    return p;
}

// Helper function: cumulative hazard at time t for piecewise-constant hazard
double cumhaz_pc(double t, const arma::vec& knots, const arma::vec& lambda) {
    int M = lambda.n_elem;
    double total = 0.0;
    
    for (int m = 0; m < M; ++m) {
        double left = knots(m);
        double right = knots(m + 1);
        double len = std::max(0.0, std::min(t, right) - left);
        total += lambda(m) * len;
    }
    return total;
}

// Helper function: phi(theta, delta)
inline double phi(double theta, double delta) {
    if (std::abs(theta) < 1e-12) {
        return delta;
    } else {
        return (std::exp(theta * delta) - 1.0) / theta;
    }
}

// Build union grid for [a,b] given all knots
std::vector<double> build_grid(double a, double b, 
                                const arma::vec& knots12,
                                const arma::vec& knots13,
                                const arma::vec& knots23) {
    std::vector<double> grid;
    grid.push_back(a);
    grid.push_back(b);
    
    // Add knots from all three transitions
    for (size_t i = 0; i < knots12.n_elem; ++i) {
        if (knots12(i) > a && knots12(i) < b) {
            grid.push_back(knots12(i));
        }
    }
    for (size_t i = 0; i < knots13.n_elem; ++i) {
        if (knots13(i) > a && knots13(i) < b) {
            grid.push_back(knots13(i));
        }
    }
    for (size_t i = 0; i < knots23.n_elem; ++i) {
        if (knots23(i) > a && knots23(i) < b) {
            grid.push_back(knots23(i));
        }
    }
    
    // Sort and remove duplicates
    std::sort(grid.begin(), grid.end());
    grid.erase(std::unique(grid.begin(), grid.end()), grid.end());
    
    return grid;
}

// Generic I_23(a,b) integral
double I23_integral(double a, double b, double T,
                    const arma::vec& knots12, const arma::vec& lambda12,
                    const arma::vec& knots13, const arma::vec& lambda13,
                    const arma::vec& knots23, const arma::vec& lambda23) {
    
    std::vector<double> tau = build_grid(a, b, knots12, knots13, knots23);
    int R = tau.size() - 1;
    double out = 0.0;
    double A23T = cumhaz_pc(T, knots23, lambda23);
    
    for (int r = 0; r < R; ++r) {
        double left = tau[r];
        double right = tau[r + 1];
        double delta = right - left;
        
        // Find which interval left falls into for each transition
        int i12 = 0, i13 = 0, i23 = 0;
        for (size_t i = 0; i < knots12.n_elem - 1; ++i) {
            if (left >= knots12(i)) i12 = i;
        }
        for (size_t i = 0; i < knots13.n_elem - 1; ++i) {
            if (left >= knots13(i)) i13 = i;
        }
        for (size_t i = 0; i < knots23.n_elem - 1; ++i) {
            if (left >= knots23(i)) i23 = i;
        }
        
        double lam12 = lambda12(i12);
        double lam13 = lambda13(i13);
        double lam23 = lambda23(i23);
        
        double H12 = cumhaz_pc(left, knots12, lambda12);
        double H13 = cumhaz_pc(left, knots13, lambda13);
        double H23 = cumhaz_pc(left, knots23, lambda23);
        double H = H12 + H13 - H23;
        
        double theta = -lam12 - lam13 + lam23;
        
        out += lam23 * std::exp(-A23T - H) * phi(theta, delta);
    }
    
    return out;
}

// Generic I_12(a,b) integral
double I12_integral(double a, double b, double T,
                    const arma::vec& knots12, const arma::vec& lambda12,
                    const arma::vec& knots13, const arma::vec& lambda13,
                    const arma::vec& knots23, const arma::vec& lambda23) {
    
    std::vector<double> tau = build_grid(a, b, knots12, knots13, knots23);
    int R = tau.size() - 1;
    double out = 0.0;
    double A23T = cumhaz_pc(T, knots23, lambda23);
    
    for (int r = 0; r < R; ++r) {
        double left = tau[r];
        double right = tau[r + 1];
        double delta = right - left;
        
        // Find which interval left falls into for each transition
        int i12 = 0, i13 = 0, i23 = 0;
        for (size_t i = 0; i < knots12.n_elem - 1; ++i) {
            if (left >= knots12(i)) i12 = i;
        }
        for (size_t i = 0; i < knots13.n_elem - 1; ++i) {
            if (left >= knots13(i)) i13 = i;
        }
        for (size_t i = 0; i < knots23.n_elem - 1; ++i) {
            if (left >= knots23(i)) i23 = i;
        }
        
        double lam12 = lambda12(i12);
        double lam13 = lambda13(i13);
        double lam23 = lambda23(i23);
        
        double H12 = cumhaz_pc(left, knots12, lambda12);
        double H13 = cumhaz_pc(left, knots13, lambda13);
        double H23 = cumhaz_pc(left, knots23, lambda23);
        double H = H12 + H13 - H23;
        
        double theta = -lam12 - lam13 + lam23;
        
        out += lam12 * std::exp(-A23T - H) * phi(theta, delta);
    }
    
    return out;
}

// Case 1: Censored (healthy at T_obs, was healthy at V_healthy)
double calc_case_1_log_likelihood(
    const PwcModelDataCaseSpecific& md,
    const arma::vec& lambda_12,
    const arma::vec& lambda_13,
    const arma::vec& lambda_23
) {
    int n = md.T_obs_values.n_elem;
    double log_lik = 0.0;
    
    for (int i = 0; i < n; ++i) {
        double V_0 = md.V_0_values(i, 0);
        double V_m = md.V_healthy_values(i);
        double T = md.T_obs_values(i);
        
        double A12_V0 = cumhaz_pc(V_0, md.knots_12, lambda_12);
        double A13_V0 = cumhaz_pc(V_0, md.knots_13, lambda_13);
        double A12_T = cumhaz_pc(T, md.knots_12, lambda_12);
        double A13_T = cumhaz_pc(T, md.knots_13, lambda_13);
        
        double term1 = std::exp(-A12_T - A13_T);
        double term2 = I23_integral(V_m, T, T, md.knots_12, lambda_12,
                                    md.knots_13, lambda_13,
                                    md.knots_23, lambda_23);
        
        double L = std::exp(A12_V0 + A13_V0) * (term1 + term2);
        
        if (L > 0 && std::isfinite(L)) {
            log_lik += std::log(L);
        } else {
            log_lik += -1e10;
        }
    }
    
    return log_lik;
}

// Case 2: Died at T_obs, was healthy at V_healthy
double calc_case_2_log_likelihood(
    const PwcModelDataCaseSpecific& md,
    const arma::vec& lambda_12,
    const arma::vec& lambda_13,
    const arma::vec& lambda_23
) {
    int n = md.T_obs_values.n_elem;
    double log_lik = 0.0;
    
    for (int i = 0; i < n; ++i) {
        double V_0 = md.V_0_values(i, 0);
        double V_m = md.V_healthy_values(i);
        double T = md.T_obs_values(i);
        
        double A12_V0 = cumhaz_pc(V_0, md.knots_12, lambda_12);
        double A13_V0 = cumhaz_pc(V_0, md.knots_13, lambda_13);
        double A12_T = cumhaz_pc(T, md.knots_12, lambda_12);
        double A13_T = cumhaz_pc(T, md.knots_13, lambda_13);
        
        // Find alpha_13(T) and alpha_23(T)
        int i13 = 0, i23 = 0;
        for (size_t j = 0; j < md.knots_13.n_elem - 1; ++j) {
            if (T >= md.knots_13(j)) i13 = j;
        }
        for (size_t j = 0; j < md.knots_23.n_elem - 1; ++j) {
            if (T >= md.knots_23(j)) i23 = j;
        }
        double alpha13_T = lambda_13(i13);
        double alpha23_T = lambda_23(i23);
        
        double term1 = std::exp(-A12_T - A13_T) * alpha13_T;
        double term2 = alpha23_T * I12_integral(V_m, T, T, md.knots_12, lambda_12,
                                                md.knots_13, lambda_13,
                                                md.knots_23, lambda_23);
        
        double L = std::exp(A12_V0 + A13_V0) * (term1 + term2);
        
        if (L > 0 && std::isfinite(L)) {
            log_lik += std::log(L);
        } else {
            log_lik += -1e10;
        }
    }
    
    return log_lik;
}

// Case 3: Ill at T_obs, was healthy at V_healthy, became ill at V_ill
double calc_case_3_log_likelihood(
    const PwcModelDataCaseSpecific& md,
    const arma::vec& lambda_12,
    const arma::vec& lambda_13,
    const arma::vec& lambda_23
) {
    int n = md.T_obs_values.n_elem;
    double log_lik = 0.0;
    
    for (int i = 0; i < n; ++i) {
        double V_0 = md.V_0_values(i, 0);
        double V_k = md.V_healthy_values(i);
        double V_kp1 = md.V_ill_values(i);
        double T = md.T_obs_values(i);
        
        double A12_V0 = cumhaz_pc(V_0, md.knots_12, lambda_12);
        double A13_V0 = cumhaz_pc(V_0, md.knots_13, lambda_13);
        double A23_T = cumhaz_pc(T, md.knots_23, lambda_23);
        
        double integral = I12_integral(V_k, V_kp1, T, md.knots_12, lambda_12,
                                       md.knots_13, lambda_13,
                                       md.knots_23, lambda_23);
        
        double L = std::exp(A12_V0 + A13_V0 - A23_T) * integral;
        
        if (L > 0 && std::isfinite(L)) {
            log_lik += std::log(L);
        } else {
            log_lik += -1e10;
        }
    }
    
    return log_lik;
}

// Case 4: Died at T_obs, was healthy at V_healthy, became ill at V_ill
double calc_case_4_log_likelihood(
    const PwcModelDataCaseSpecific& md,
    const arma::vec& lambda_12,
    const arma::vec& lambda_13,
    const arma::vec& lambda_23
) {
    int n = md.T_obs_values.n_elem;
    double log_lik = 0.0;
    
    for (int i = 0; i < n; ++i) {
        double V_0 = md.V_0_values(i, 0);
        double V_k = md.V_healthy_values(i);
        double V_kp1 = md.V_ill_values(i);
        double T = md.T_obs_values(i);
        
        double A12_V0 = cumhaz_pc(V_0, md.knots_12, lambda_12);
        double A13_V0 = cumhaz_pc(V_0, md.knots_13, lambda_13);
        double A23_T = cumhaz_pc(T, md.knots_23, lambda_23);
        
        // Find alpha_23(T)
        int i23 = 0;
        for (size_t j = 0; j < md.knots_23.n_elem - 1; ++j) {
            if (T >= md.knots_23(j)) i23 = j;
        }
        double alpha23_T = lambda_23(i23);
        
        double integral = I12_integral(V_k, V_kp1, T, md.knots_12, lambda_12,
                                       md.knots_13, lambda_13,
                                       md.knots_23, lambda_23);
        
        double L = std::exp(A12_V0 + A13_V0 - A23_T) * alpha23_T * integral;
        
        if (L > 0 && std::isfinite(L)) {
            log_lik += std::log(L);
        } else {
            log_lik += -1e10;
        }
    }
    
    return log_lik;
}

// [[Rcpp::export]]
double calc_log_likelihood_pwc(
    SEXP md_ptr,
    const arma::vec& lambda_12,
    const arma::vec& lambda_13,
    const arma::vec& lambda_23
) {
    Rcpp::XPtr<PwcModelData> pointer(md_ptr);
    const PwcModelData& mddata = *pointer;

    double log_likelihood = 0.0;

    // calculate log-likelihood contributions from each case
    if (mddata.case_1.T_obs_values.n_elem > 0) {
        log_likelihood += calc_case_1_log_likelihood(
            mddata.case_1, lambda_12, lambda_13, lambda_23);
    }
    if (mddata.case_2.T_obs_values.n_elem > 0) {
        log_likelihood += calc_case_2_log_likelihood(
            mddata.case_2, lambda_12, lambda_13, lambda_23);
    }
    if (mddata.case_3.T_obs_values.n_elem > 0) {
        log_likelihood += calc_case_3_log_likelihood(
            mddata.case_3, lambda_12, lambda_13, lambda_23);
    }
    if (mddata.case_4.T_obs_values.n_elem > 0) {
        log_likelihood += calc_case_4_log_likelihood(
            mddata.case_4, lambda_12, lambda_13, lambda_23);
    }
    return log_likelihood;
}
