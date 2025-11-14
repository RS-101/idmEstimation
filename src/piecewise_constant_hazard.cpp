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


struct PwcModelDataCaseSpecific {
    arma::mat V_0_values;
    arma::vec V_healthy_values;
    arma::vec V_ill_values; 
    arma::vec T_obs_values;

    arma::vec knots_12;
    arma::vec knots_13;
    arma::vec knots_23;

    // Constructor from R list
    PwcModelDataCaseSpecific(List data) {
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

double calc_case_1_log_likelihood(
    const PwcModelDataCaseSpecific& md,
    const arma::vec& lambda_12,
    const arma::vec& lambda_13,
    const arma::vec& lambda_23
) {
 
}

double calc_case_2_log_likelihood(
    const PwcModelDataCaseSpecific& md,
    const arma::vec& lambda_12,
    const arma::vec& lambda_13,
    const arma::vec& lambda_23
) {
 
}double calc_case_3_log_likelihood(
    const PwcModelDataCaseSpecific& md,
    const arma::vec& lambda_12,
    const arma::vec& lambda_13,
    const arma::vec& lambda_23
) {
 
}double calc_case_4_log_likelihood(
    const PwcModelDataCaseSpecific& md,
    const arma::vec& lambda_12,
    const arma::vec& lambda_13,
    const arma::vec& lambda_23
) {
 
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
