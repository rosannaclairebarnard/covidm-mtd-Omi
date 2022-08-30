// convenience.h

#ifndef CONVENIENCE_H
#define CONVENIENCE_H

#include <vector>

struct Parameters;
class Reporter;

// TODO perform rigorous error checking on inputs

// Helper functions for user code
std::vector<double> seq(double x0, double x1, double by = 1);

// binomial log density
double binom(double x, double size, double prob);

// negative binomial log density
double nbinom(double x, double mean, double size);

// poisson log density
double poisson(double x, double lambda);

// beta binomial log density
double bbinom(double k, double n, double mode, double conc);

// negative binomial log density with retrospective confirmation
double nbinom_gammaconf(unsigned int x, double mean, double size, double days_ago, double conf_delay_mean, double conf_delay_shape);

// normal log density
double norm(double x, double mean, double sd);

// skew normal log density
double skewnorm(double x, double xi, double omega, double alpha);

// beta density
double dbeta(double x, double alpha, double beta);

// calculate penalty for gamma-prior fit to data
double rate_penalty(std::vector<double>& t, std::vector<double>& x, Reporter& dyn, const char* y_comp, const char* y_comp_minus, 
    double alpha, unsigned int k, double sigma);

// construct a delay distribution following a gamma distribution with mean mu and shape parameter shape.
// [[Rcpp::export]]
std::vector<double> delay_gamma(double mu, double shape, double t_max, double t_step, double mult = 1.);

// construct a delay distribution following a lognormal distribution with true mean mu and coefficient of variation cv.
std::vector<double> delay_lnorm(double mu, double cv, double t_max, double t_step);

// add (convolve) two delay distributions
std::vector<double> delay_convolve(const std::vector<double> a, const std::vector<double> b);

// estimate the basic reproduction number
double estimate_R0(Parameters& P, Reporter& rep, double t, unsigned int p, unsigned int iter);

// estimate the effective reproduction number
double estimate_Rt(Parameters& P, Reporter& rep, double t, unsigned int p, unsigned int iter, double t_seas = -999);

// clamp a number between two limits
double clamp(double x, double x0 = 0.0, double x1 = 1.0);

// smootherstep function
double smootherstep(double x0, double x1, double y0, double y1, double x);

// protective efficacy fold change function
double efficacy_fold(double eff, double fold);

// infection to severe protection
double infection_to_severe(double ei);

#endif