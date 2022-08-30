// convenience.cpp

// Helper functions
#include "convenience.h"
#include "parameters.h"
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include <cmath>
#include <iostream>
#include <algorithm>
using namespace std;

vector<double> seq(double x0, double x1, double by)
{
    int steps = int((x1 - x0) / by + by * 1e-3);
    if (steps >= 0)
    {
        vector<double> ret(steps);
        for (int i = 0; i < steps; ++i)
            ret[i] = x0 + by * i;
        return ret;
    }
    return vector<double>(0);
}

// binomial log density
double binom(double x, double size, double prob)
{
    return gsl_sf_lngamma(size + 1.) - gsl_sf_lngamma(x + 1.) - gsl_sf_lngamma(size - x + 1.) + 
        x * log(prob) + (size - x) * log(1. - prob);
}

// negative binomial log density
double nbinom(double x, double mean, double size)
{
    double k = x;
    double p = size / (size + mean);
    double n = mean * p / (1 - p);

    return gsl_sf_lngamma(n + k) - gsl_sf_lngamma(k + 1) - gsl_sf_lngamma(n) + n * log(p) + k * log(1 - p);
}

// poisson log density
double poisson(double x, double lambda)
{
    return x * log(lambda) - lambda - gsl_sf_lngamma(x + 1);
}

// beta binomial log density
double bbinom(double k, double n, double mode, double conc)
{
    auto lgamma = [](double x) { return gsl_sf_lngamma(x); };
    double a = mode * (conc - 2) + 1;
    double b = (1 - mode) * (conc - 2) + 1;

    return (lgamma(n + 1) + lgamma(k + a) + lgamma(n - k + b) + lgamma(a + b))
        - (lgamma(k + 1) + lgamma(n - k + 1) + lgamma(n + a + b) + lgamma(a) + lgamma(b));
}

// negative binomial log density with retrospective confirmation
double nbinom_gammaconf(unsigned int x, double mean, double size, double days_ago, double conf_delay_mean, double conf_delay_shape)
{
    double conf_delay_scale = conf_delay_mean / conf_delay_shape;
    double prop_confirmed = gsl_cdf_gamma_P(days_ago, conf_delay_shape, conf_delay_scale);
    return nbinom(x, mean * prop_confirmed, size);
}

// normal log density
double norm(double x, double mean, double sd)
{
    return -(x - mean) * (x - mean) / (2 * sd * sd) - log(sd) - 0.918938533204673;
}

// skew normal log density
double skewnorm(double x, double xi, double omega, double alpha)
{
    double X = (x - xi) / omega;
    return -log(omega) + norm(X, 0., 1.) + gsl_sf_log_erfc(-alpha * X / 1.414213562373095);
}

// beta density
double dbeta(double x, double alpha, double beta)
{
    return gsl_ran_beta_pdf(x, alpha, beta);
}

// calculate penalty for gamma-prior fit to data
double rate_penalty(std::vector<double>& t, std::vector<double>& x, Reporter& dyn, const char* y_comp, const char* y_comp_minus, 
    double alpha, unsigned int k, double sigma)
{
    vector<double> posterior_shape(t.size(), 0.0);
    vector<double> posterior_rate(t.size(), 0.0);
    
    for (unsigned int i = 0; i < t.size(); ++i)
    {
        double prior_mean = 0.0;
        if (!y_comp_minus)
            prior_mean = max(dyn(y_comp, t[i], {}, {}), alpha);
        else
            prior_mean = max(dyn(y_comp, t[i], {}, {}) - dyn(y_comp_minus, t[i], {}, {}), alpha);
        double prior_sd = prior_mean * sigma;
        double prior_shape = (prior_mean * prior_mean) / (prior_sd * prior_sd);
        double prior_rate = prior_mean / (prior_sd * prior_sd);
        posterior_shape[i] = prior_shape + x[i];
        posterior_rate[i] = prior_rate + 1.0;
    }
    
    double ll = 0.0;
    if (k > 1)
    {
        int hk = k / 2;
        for (int i = 0; i < t.size(); ++i)
        {
            double roll_shape = 0.0, n = 0.0;
            for (int j = std::max(0, i - hk); j <= std::min((int)t.size() - 1, i + hk); ++j)
            {
                roll_shape += posterior_shape[i];
                ++n;
            }
            roll_shape /= n;
            ll += nbinom(max(x[i], alpha), roll_shape / posterior_rate[i], roll_shape);
        }
    }
    else
    {
        for (unsigned int i = 0; i < t.size(); ++i)
        {
            ll += nbinom(max(x[i], alpha), posterior_shape[i] / posterior_rate[i], posterior_shape[i]);
        }
    }

    return ll;
}


// construct a delay distribution following a gamma distribution with mean mu and shape parameter shape.
std::vector<double> delay_gamma(double mu, double shape, double t_max, double t_step, double mult)
{
    double scale = mu / shape;
    vector<double> height;

    for (double t = 0.0; t < t_max + 0.5 * t_step; t += t_step)
        height.push_back(mult * (gsl_cdf_gamma_P(t + t_step/2, shape, scale) -
            gsl_cdf_gamma_P(max(0.0, t - t_step/2), shape, scale)));
    return height;
}

// construct a delay distribution following a lognormal distribution with true mean mu and coefficient of variation cv.
std::vector<double> delay_lnorm(double mu, double cv, double t_max, double t_step)
{
    double meanlog = log(mu / sqrt(1 + cv * cv));
    double sdlog = sqrt(log(1 + cv * cv));
    vector<double> height;

    for (double t = 0.0; t < t_max + 0.5 * t_step; t += t_step)
        height.push_back(gsl_cdf_lognormal_P(t + t_step/2, meanlog, sdlog) -
            gsl_cdf_lognormal_P(max(0.0, t - t_step/2), meanlog, sdlog));
    return height;
}

// add (convolve) two delay distributions
std::vector<double> delay_convolve(const std::vector<double> a, const std::vector<double> b)
{
    std::vector<double> c(a.size(), 0.0);
    for (unsigned int i = 0; i < c.size(); ++i)
    {
        for (unsigned int j = 0; j <= i; ++j)
        {
            c[i] += a[j] * b[i - j];
        }
    }
    return c;
}


// estimate the basic reproduction number: three strains
double estimate_R0(Parameters& P, Reporter& dyn, double t, unsigned int p, unsigned int iter)
{
    return estimate_Rt(P, dyn, P.time0 + 1, p, iter, t);
}

//   0    1    2     3     4     5    6      7     8       9     10
//  "S", "E", "Ip", "Is", "Ia", "R", "E2", "Ip2", "Is2", "Ia2", "R2",
//  
//  11     12    13     14    15      16    17     18      19   20      21 
// "E3", "Ip3", "Is3", "Ia3", "R3", "Va1", "Va2", "Va3", "Vb1", "Vb2", "Vb3",

// estimate the effective reproduction number: three strains
double estimate_Rt(Parameters& P, Reporter& dyn, double t, unsigned int p, unsigned int iter, double t_seas)
{
    vector<double> E(P.pop[p].size.size(), 1.0);
    vector<double> L(P.pop[p].size.size(), 1.0);
    vector<double> Ef = E;
    vector<double> Lf = L;

    vector<double> SE(P.pop[p].size.size(), 0.0);
    vector<double> SL(P.pop[p].size.size(), 0.0);
    vector<double> SE2(P.pop[p].size.size(), 0.0);
    vector<double> SL2(P.pop[p].size.size(), 0.0);
    vector<double> SE3(P.pop[p].size.size(), 0.0);
    vector<double> SL3(P.pop[p].size.size(), 0.0);

    // Calculate effective susceptible fraction
    for (unsigned int a = 0; a < SE.size(); ++a) {
        SE[a] = (dyn(t, p, a, 0) +                              // E-susceptible to strain 1 in S compartment
            dyn(t, p, a, 5)  * (1.0 - P.pop[p].pi_r[a])  * (1.0 - P.pop[p].pd_ri[a])  +   // E-susceptible to strain 1 in R compartment
            dyn(t, p, a, 10) * (1.0 - P.pop[p].pi_r2[a]) * (1.0 - P.pop[p].pd_r2i[a]) +   // E-susceptible to strain 1 in R2 compartment
            dyn(t, p, a, 15) * (1.0 - P.pop[p].pi_r3[a]) * (1.0 - P.pop[p].pd_r3i[a]) +   // E-susceptible to strain 1 in R3 compartment
            dyn(t, p, a, 16) * (1.0 - P.pop[p].ei_va1[a]) * (1.0 - P.pop[p].ed_va1i[a]) + // E-susceptible to strain 1 in Va1 compartment
            dyn(t, p, a, 17) * (1.0 - P.pop[p].ei_va2[a]) * (1.0 - P.pop[p].ed_va2i[a]) + // E-susceptible to strain 1 in Va2 compartment
            dyn(t, p, a, 18) * (1.0 - P.pop[p].ei_va3[a]) * (1.0 - P.pop[p].ed_va3i[a]) + // E-susceptible to strain 1 in Va3 compartment
            dyn(t, p, a, 19) * (1.0 - P.pop[p].ei_vb1[a]) * (1.0 - P.pop[p].ed_vb1i[a]) + // E-susceptible to strain 1 in Vb1 compartment
            dyn(t, p, a, 20) * (1.0 - P.pop[p].ei_vb2[a]) * (1.0 - P.pop[p].ed_vb2i[a]) + // E-susceptible to strain 1 in Vb2 compartment
            dyn(t, p, a, 21) * (1.0 - P.pop[p].ei_vb3[a]) * (1.0 - P.pop[p].ed_vb3i[a])   // E-susceptible to strain 1 in Vb3 compartment
            ) / P.pop[p].size[a];
        SL[a] = (
            dyn(t, p, a, 5)  * (1.0 - P.pop[p].pi_r[a])  * (P.pop[p].pd_ri[a]  - P.pop[p].pt_ri[a])  +   // L-susceptible to strain 1 in R compartment
            dyn(t, p, a, 10) * (1.0 - P.pop[p].pi_r2[a]) * (P.pop[p].pd_r2i[a] - P.pop[p].pt_r2i[a]) +   // L-susceptible to strain 1 in R2 compartment
            dyn(t, p, a, 15) * (1.0 - P.pop[p].pi_r3[a]) * (P.pop[p].pd_r3i[a] - P.pop[p].pt_r3i[a]) +   // L-susceptible to strain 1 in R3 compartment
            dyn(t, p, a, 16) * (1.0 - P.pop[p].ei_va1[a]) * (P.pop[p].ed_va1i[a] - P.pop[p].et_va1i[a]) +   // L-susceptible to strain 1 in Va1 compartment
            dyn(t, p, a, 17) * (1.0 - P.pop[p].ei_va2[a]) * (P.pop[p].ed_va2i[a] - P.pop[p].et_va2i[a]) +   // L-susceptible to strain 1 in Va2 compartment
            dyn(t, p, a, 18) * (1.0 - P.pop[p].ei_va3[a]) * (P.pop[p].ed_va3i[a] - P.pop[p].et_va3i[a]) +   // L-susceptible to strain 1 in Va3 compartment
            dyn(t, p, a, 19) * (1.0 - P.pop[p].ei_vb1[a]) * (P.pop[p].ed_vb1i[a] - P.pop[p].et_vb1i[a]) +   // L-susceptible to strain 1 in Vb1 compartment
            dyn(t, p, a, 20) * (1.0 - P.pop[p].ei_vb2[a]) * (P.pop[p].ed_vb2i[a] - P.pop[p].et_vb2i[a]) +   // L-susceptible to strain 1 in Vb2 compartment
            dyn(t, p, a, 21) * (1.0 - P.pop[p].ei_vb3[a]) * (P.pop[p].ed_vb3i[a] - P.pop[p].et_vb3i[a])     // L-susceptible to strain 1 in Vb3 compartment
            ) / P.pop[p].size[a];
        SE2[a] = (dyn(t, p, a, 0) +                             // E-susceptible to strain 2 in S compartment
            dyn(t, p, a, 5)  * (1.0 - P.pop[p].pi2_r[a])  * (1.0 - P.pop[p].pd_ri2[a])  +   // E-susceptible to strain 2 in R compartment
            dyn(t, p, a, 10) * (1.0 - P.pop[p].pi2_r2[a]) * (1.0 - P.pop[p].pd_r2i2[a]) +   // E-susceptible to strain 2 in R2 compartment
            dyn(t, p, a, 15) * (1.0 - P.pop[p].pi2_r3[a]) * (1.0 - P.pop[p].pd_r3i2[a]) +   // E-susceptible to strain 2 in R3 compartment
            dyn(t, p, a, 16) * (1.0 - P.pop[p].ei2_va1[a]) * (1.0 - P.pop[p].ed_va1i2[a]) + // E-susceptible to strain 2 in Va1 compartment
            dyn(t, p, a, 17) * (1.0 - P.pop[p].ei2_va2[a]) * (1.0 - P.pop[p].ed_va2i2[a]) + // E-susceptible to strain 2 in Va2 compartment
            dyn(t, p, a, 18) * (1.0 - P.pop[p].ei2_va3[a]) * (1.0 - P.pop[p].ed_va3i2[a]) + // E-susceptible to strain 2 in Va3 compartment
            dyn(t, p, a, 19) * (1.0 - P.pop[p].ei2_vb1[a]) * (1.0 - P.pop[p].ed_vb1i2[a]) + // E-susceptible to strain 2 in Vb1 compartment
            dyn(t, p, a, 20) * (1.0 - P.pop[p].ei2_vb2[a]) * (1.0 - P.pop[p].ed_vb2i2[a]) + // E-susceptible to strain 2 in Vb2 compartment
            dyn(t, p, a, 21) * (1.0 - P.pop[p].ei2_vb3[a]) * (1.0 - P.pop[p].ed_vb3i2[a])   // E-susceptible to strain 2 in Vb3 compartment
            ) / P.pop[p].size[a];
        SL2[a] = (
            dyn(t, p, a, 5)  * (1.0 - P.pop[p].pi2_r[a])  * (P.pop[p].pd_ri2[a]  - P.pop[p].pt_ri2[a])  +   // L-susceptible to strain 2 in R compartment
            dyn(t, p, a, 10) * (1.0 - P.pop[p].pi2_r2[a]) * (P.pop[p].pd_r2i2[a] - P.pop[p].pt_r2i2[a]) +   // L-susceptible to strain 2 in R2 compartment
            dyn(t, p, a, 15) * (1.0 - P.pop[p].pi2_r3[a]) * (P.pop[p].pd_r3i2[a] - P.pop[p].pt_r3i2[a]) +   // L-susceptible to strain 2 in R3 compartment
            dyn(t, p, a, 16) * (1.0 - P.pop[p].ei2_va1[a]) * (P.pop[p].ed_va1i2[a] - P.pop[p].et_va1i2[a]) + // L-susceptible to strain 2 in Va1 compartment
            dyn(t, p, a, 17) * (1.0 - P.pop[p].ei2_va2[a]) * (P.pop[p].ed_va2i2[a] - P.pop[p].et_va2i2[a]) + // L-susceptible to strain 2 in Va2 compartment
            dyn(t, p, a, 18) * (1.0 - P.pop[p].ei2_va3[a]) * (P.pop[p].ed_va3i2[a] - P.pop[p].et_va3i2[a]) + // L-susceptible to strain 2 in Va3 compartment
            dyn(t, p, a, 19) * (1.0 - P.pop[p].ei2_vb1[a]) * (P.pop[p].ed_vb1i2[a] - P.pop[p].et_vb1i2[a]) + // L-susceptible to strain 2 in Vb1 compartment
            dyn(t, p, a, 20) * (1.0 - P.pop[p].ei2_vb2[a]) * (P.pop[p].ed_vb2i2[a] - P.pop[p].et_vb2i2[a]) + // L-susceptible to strain 2 in Vb2 compartment
            dyn(t, p, a, 21) * (1.0 - P.pop[p].ei2_vb3[a]) * (P.pop[p].ed_vb3i2[a] - P.pop[p].et_vb3i2[a])   // L-susceptible to strain 2 in Vb3 compartment
            ) / P.pop[p].size[a];
        SE3[a] = (dyn(t, p, a, 0) +                             // E-susceptible to strain 3 in S compartment
            dyn(t, p, a, 5)  * (1.0 - P.pop[p].pi3_r[a])  * (1.0 - P.pop[p].pd_ri3[a])  +   // E-susceptible to strain 3 in R compartment
            dyn(t, p, a, 10) * (1.0 - P.pop[p].pi3_r2[a]) * (1.0 - P.pop[p].pd_r2i3[a]) +   // E-susceptible to strain 3 in R2 compartment
            dyn(t, p, a, 15) * (1.0 - P.pop[p].pi3_r3[a]) * (1.0 - P.pop[p].pd_r3i3[a]) +   // E-susceptible to strain 3 in R3 compartment
            dyn(t, p, a, 16) * (1.0 - P.pop[p].ei3_va1[a]) * (1.0 - P.pop[p].ed_va1i3[a]) + // E-susceptible to strain 3 in Va1 compartment
            dyn(t, p, a, 17) * (1.0 - P.pop[p].ei3_va2[a]) * (1.0 - P.pop[p].ed_va2i3[a]) + // E-susceptible to strain 3 in Va2 compartment
            dyn(t, p, a, 18) * (1.0 - P.pop[p].ei3_va3[a]) * (1.0 - P.pop[p].ed_va3i3[a]) + // E-susceptible to strain 3 in Va3 compartment
            dyn(t, p, a, 19) * (1.0 - P.pop[p].ei3_vb1[a]) * (1.0 - P.pop[p].ed_vb1i3[a]) + // E-susceptible to strain 3 in Vb1 compartment
            dyn(t, p, a, 20) * (1.0 - P.pop[p].ei3_vb2[a]) * (1.0 - P.pop[p].ed_vb2i3[a]) + // E-susceptible to strain 3 in Vb2 compartment
            dyn(t, p, a, 21) * (1.0 - P.pop[p].ei3_vb3[a]) * (1.0 - P.pop[p].ed_vb3i3[a])   // E-susceptible to strain 3 in Vb3 compartment
            ) / P.pop[p].size[a];
        SL3[a] = (
            dyn(t, p, a, 5)  * (1.0 - P.pop[p].pi3_r[a])  * (P.pop[p].pd_ri3[a]  - P.pop[p].pt_ri3[a])  +   // L-susceptible to strain 3 in R compartment
            dyn(t, p, a, 10) * (1.0 - P.pop[p].pi3_r2[a]) * (P.pop[p].pd_r2i3[a] - P.pop[p].pt_r2i3[a]) +   // L-susceptible to strain 3 in R2 compartment
            dyn(t, p, a, 15) * (1.0 - P.pop[p].pi3_r3[a]) * (P.pop[p].pd_r3i3[a] - P.pop[p].pt_r3i3[a]) +   // L-susceptible to strain 3 in R3 compartment
            dyn(t, p, a, 16) * (1.0 - P.pop[p].ei3_va1[a]) * (P.pop[p].ed_va1i3[a] - P.pop[p].et_va1i3[a]) + // L-susceptible to strain 3 in Va1 compartment
            dyn(t, p, a, 17) * (1.0 - P.pop[p].ei3_va2[a]) * (P.pop[p].ed_va2i3[a] - P.pop[p].et_va2i3[a]) + // L-susceptible to strain 3 in Va2 compartment
            dyn(t, p, a, 18) * (1.0 - P.pop[p].ei3_va3[a]) * (P.pop[p].ed_va3i3[a] - P.pop[p].et_va3i3[a]) + // L-susceptible to strain 3 in Va3 compartment
            dyn(t, p, a, 19) * (1.0 - P.pop[p].ei3_vb1[a]) * (P.pop[p].ed_vb1i3[a] - P.pop[p].et_vb1i3[a]) + // L-susceptible to strain 3 in Vb1 compartment
            dyn(t, p, a, 20) * (1.0 - P.pop[p].ei3_vb2[a]) * (P.pop[p].ed_vb2i3[a] - P.pop[p].et_vb2i3[a]) + // L-susceptible to strain 3 in Vb2 compartment
            dyn(t, p, a, 21) * (1.0 - P.pop[p].ei3_vb3[a]) * (P.pop[p].ed_vb3i3[a] - P.pop[p].et_vb3i3[a])   // L-susceptible to strain 3 in Vb3 compartment
            ) / P.pop[p].size[a];
    }

    double dIp  = P.pop[p].dIp.Mean() * P.time_step;
    double dIs  = P.pop[p].dIs.Mean() * P.time_step;
    double dIa  = P.pop[p].dIa.Mean() * P.time_step;
    double dIp2 = (P.pop[p].dIp2.weights.size() > 1 ? P.pop[p].dIp2.Mean() : P.pop[p].dIp.Mean()) * P.time_step;
    double dIs2 = (P.pop[p].dIs2.weights.size() > 1 ? P.pop[p].dIs2.Mean() : P.pop[p].dIs.Mean()) * P.time_step;
    double dIa2 = (P.pop[p].dIa2.weights.size() > 1 ? P.pop[p].dIa2.Mean() : P.pop[p].dIa.Mean()) * P.time_step;
    double dIp3 = (P.pop[p].dIp3.weights.size() > 1 ? P.pop[p].dIp3.Mean() : P.pop[p].dIp.Mean()) * P.time_step;
    double dIs3 = (P.pop[p].dIs3.weights.size() > 1 ? P.pop[p].dIs3.Mean() : P.pop[p].dIs.Mean()) * P.time_step;
    double dIa3 = (P.pop[p].dIa3.weights.size() > 1 ? P.pop[p].dIa3.Mean() : P.pop[p].dIa.Mean()) * P.time_step;

    double nE  = P.pop[p].size.size();
    double nL  = P.pop[p].size.size();
    double nEf = 0;
    double nLf = 0;
    double Rt = 1;

    double seas = 1.0;
    if (t_seas == -999) t_seas += P.pop[p].season_A[0] * cos(2. * M_PI * (t      - P.pop[p].season_phi[0]) / P.pop[p].season_T[0]);
    else                t_seas += P.pop[p].season_A[0] * cos(2. * M_PI * (t_seas - P.pop[p].season_phi[0]) / P.pop[p].season_T[0]);

    vector<double> f(P.pop[p].size.size(), 0);
    vector<double> f2 = f;
    vector<double> f3 = f;
    for (unsigned int a = 0; a < f2.size(); ++a)
    {
        double e1 = dyn(t, p, a, 1);
        double e2 = dyn(t, p, a, 6);
        double e3 = dyn(t, p, a, 11);
        if (e1 + e2 + e3 > 0) {
            f[a] = e1 / (e1 + e2 + e3);
            f2[a] = e2 / (e1 + e2 + e3);
            f3[a] = e3 / (e1 + e2 + e3);
        }
    }

    for (unsigned int i = 0; i < iter; ++i)
    {
        nEf = nLf = 0;
        for (unsigned int a = 0; a < E.size(); ++a)
        {
            Ef[a] = 0;
            Lf[a] = 0;
            for (unsigned int b = 0; b < E.size(); ++b)
            {
                Ef[a] +=
                    f[a] * SE[a] * P.pop[p].u[a] * seas * (
                        E[b] * P.pop[p].cm(a, b) * (
                            P.pop[p].y[b] * (P.pop[p].fIp[b] * dIp + P.pop[p].fIs[b] * dIs) +
                            (1 - P.pop[p].y[b]) * P.pop[p].fIa[b] * dIa
                        ) +
                        L[b] * P.pop[p].cm(a, b) * (P.pop[p].fIa[b] * dIa)
                    ) +
                    f2[a] * SE2[a] * P.pop[p].u2[a] * seas * (
                        E[b] * P.pop[p].cm(a, b) * (
                            P.pop[p].y2[b] * (P.pop[p].fIp[b] * dIp2 + P.pop[p].fIs[b] * dIs2) +
                            (1 - P.pop[p].y2[b]) * P.pop[p].fIa[b] * dIa2
                        ) +
                        L[b] * P.pop[p].cm(a, b) * (P.pop[p].fIa[b] * dIa2)
                    ) +
                    f3[a] * SE3[a] * P.pop[p].u3[a] * seas * (
                        E[b] * P.pop[p].cm(a, b) * (
                            P.pop[p].y3[b] * (P.pop[p].fIp[b] * dIp3 + P.pop[p].fIs[b] * dIs3) +
                            (1 - P.pop[p].y3[b]) * P.pop[p].fIa[b] * dIa3
                        ) +
                        L[b] * P.pop[p].cm(a, b) * (P.pop[p].fIa[b] * dIa3)
                    );

                Lf[a] +=
                    f[a] * SL[a] * P.pop[p].u[a] * seas * (
                        E[b] * P.pop[p].cm(a, b) * (
                            P.pop[p].y[b] * (P.pop[p].fIp[b] * dIp + P.pop[p].fIs[b] * dIs) +
                            (1 - P.pop[p].y[b]) * P.pop[p].fIa[b] * dIa
                        ) +
                        L[b] * P.pop[p].cm(a, b) * (P.pop[p].fIa[b] * dIa)
                    ) +
                    f2[a] * SL2[a] * P.pop[p].u2[a] * seas * (
                        E[b] * P.pop[p].cm(a, b) * (
                            P.pop[p].y2[b] * (P.pop[p].fIp[b] * dIp2 + P.pop[p].fIs[b] * dIs2) +
                            (1 - P.pop[p].y2[b]) * P.pop[p].fIa[b] * dIa2
                        ) +
                        L[b] * P.pop[p].cm(a, b) * (P.pop[p].fIa[b] * dIa2)
                    ) + 
                    f3[a] * SL3[a] * P.pop[p].u3[a] * seas * (
                        E[b] * P.pop[p].cm(a, b) * (
                            P.pop[p].y3[b] * (P.pop[p].fIp[b] * dIp3 + P.pop[p].fIs[b] * dIs3) +
                            (1 - P.pop[p].y3[b]) * P.pop[p].fIa[b] * dIa3
                        ) +
                        L[b] * P.pop[p].cm(a, b) * (P.pop[p].fIa[b] * dIa3)
                    );
            }
            nEf += Ef[a];
            nLf += Lf[a];
        }

        Rt = (nEf + nLf) / (nE + nL);

        swap(nEf, nE);
        swap(nLf, nL);
        swap(Ef, E);
        swap(Lf, L);
    }

    return Rt;
}

// clamp a number between two limits
double clamp(double x, double x0, double x1)
{
    return min(max(x, x0), x1);
}

// smootherstep function
double smootherstep(double x0, double x1, double y0, double y1, double x)
{
    x = clamp((x - x0) / (x1 - x0));
    return y0 + x * x * x * (x * (x * 6. - 15.) + 10.) * (y1 - y0);
}

// protective efficacy fold change function
double efficacy_fold(double eff, double fold)
{
    static vector<double> kh_neut = { 
        0.001, 0.004, 0.006, 0.008, 0.011, 0.013, 0.015, 0.017, 0.019, 0.022, 0.024, 0.026, 0.029, 0.031, 
        0.034, 0.037, 0.039, 0.042, 0.045, 0.048, 0.051, 0.054, 0.057, 0.06, 0.064, 0.067, 0.071, 0.074, 
        0.078, 0.082, 0.086, 0.09, 0.094, 0.099, 0.103, 0.108, 0.113, 0.118, 0.123, 0.128, 0.134, 0.14, 
        0.146, 0.152, 0.158, 0.165, 0.172, 0.179, 0.187, 0.194, 0.203, 0.211, 0.22, 0.229, 0.239, 0.249, 
        0.259, 0.27, 0.282, 0.294, 0.306, 0.32, 0.334, 0.348, 0.364, 0.38, 0.398, 0.416, 0.436, 0.456, 
        0.479, 0.502, 0.527, 0.554, 0.583, 0.614, 0.647, 0.684, 0.723, 0.766, 0.813, 0.864, 0.921, 0.983, 
        1.053, 1.132, 1.221, 1.322, 1.438, 1.574, 1.733, 1.925, 2.159, 2.453, 2.832, 3.345, 4.079, 5.235, 7.365, 11 };
    
    static vector<double> kh_eff = {
        0.23, 1.34, 2.2, 3.09, 4.45, 5.36, 6.26, 7.16, 8.04, 9.34, 10.19, 11.02, 12.24, 13.04, 14.21, 15.34, 16.08, 
        17.17, 18.22, 19.24, 20.24, 21.21, 22.15, 23.07, 24.26, 25.13, 26.25, 27.07, 28.13, 29.15, 30.15, 31.11, 32.05, 
        33.18, 34.06, 35.12, 36.14, 37.13, 38.09, 39.02, 40.1, 41.13, 42.13, 43.1, 44.03, 45.07, 46.08, 47.05, 48.11, 
        49.01, 50.12, 51.06, 52.08, 53.05, 54.09, 55.08, 56.03, 57.03, 58.07, 59.06, 60.01, 61.06, 62.05, 63, 64.03, 
        65, 66.04, 67.01, 68.04, 69, 70.04, 71.02, 72.01, 73.02, 74.02, 75.02, 76, 77.02, 78.01, 79.02, 80.02, 81.01, 
        82.02, 83, 84, 85.01, 86.01, 87.01, 88, 89.01, 90, 91, 92, 93, 94, 95, 96, 97, 98, 98.77 };
    
    // Translate eff from proportion to percentage
    eff = eff * 100;
    
    //Rcpp::Rcout << "eff is " << eff << "\n";
    
    // Translate efficacy into neutralisation
    auto eff_right_bound_iter = lower_bound(kh_eff.begin(), kh_eff.end(), eff);
    if (eff_right_bound_iter == kh_eff.end() || eff_right_bound_iter == kh_eff.begin()) {
        Rcpp::stop("efficacy fold: eff is out of bounds");
    }
    
    unsigned long eff_right_i = eff_right_bound_iter - kh_eff.begin();
    unsigned long eff_left_i = eff_right_i - 1;
    double eff_betweenness = (eff - kh_eff[eff_left_i]) / (kh_eff[eff_right_i] - kh_eff[eff_left_i]);
    
    double neut = kh_neut[eff_left_i] * (1.0 - eff_betweenness) + kh_neut[eff_right_i] * eff_betweenness;
    
    // Adjust neutralisation by fold
    double new_neut = neut * fold;
    
    // Translate new neutralisation into efficacy
    auto neut_right_bound_iter = lower_bound(kh_neut.begin(), kh_neut.end(), new_neut);
    if (neut_right_bound_iter == kh_neut.end() || neut_right_bound_iter == kh_neut.begin()) {
        Rcpp::stop("efficacy fold: neut is out of bounds");
    }
    
    unsigned long neut_right_i = neut_right_bound_iter - kh_neut.begin();
    unsigned long neut_left_i = neut_right_i - 1;
    double neut_betweenness = (new_neut - kh_neut[neut_left_i]) / (kh_neut[neut_right_i] - kh_neut[neut_left_i]);
    
    double new_eff = kh_eff[neut_left_i] * (1.0 - neut_betweenness) + kh_eff[neut_right_i] * neut_betweenness;
    
    return new_eff / 100;
}

// infection to severe protection
double infection_to_severe(double ei)
{
    static vector<double> mild = {
        0.1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 
        21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 
        41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 
        61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 
        81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100 };
    
    static vector<double> severe = {
        1.09, 8.67, 14.96, 20.12, 24.55, 28.45, 31.96, 35.14, 38.06, 40.75, 43.25, 45.58, 
        47.77, 49.83, 51.77, 53.61, 55.35, 57.01, 58.59, 60.09, 61.53, 62.9, 64.22, 65.48, 
        66.69, 67.85, 68.97, 70.04, 71.08, 72.08, 73.04, 73.97, 74.87, 75.73, 76.57, 77.38, 
        78.17, 78.93, 79.66, 80.37, 81.06, 81.73, 82.38, 83.01, 83.62, 84.21, 84.78, 85.34, 
        85.88, 86.41, 86.92, 87.41, 87.89, 88.36, 88.81, 89.26, 89.69, 90.1, 90.51, 90.9, 
        91.29, 91.66, 92.02, 92.37, 92.71, 93.05, 93.37, 93.68, 93.99, 94.28, 94.57, 94.85, 
        95.12, 95.39, 95.64, 95.89, 96.13, 96.37, 96.59, 96.81, 97.03, 97.23, 97.43, 97.63, 
        97.81, 97.99, 98.17, 98.34, 98.5, 98.66, 98.81, 98.95, 99.09, 99.23, 99.36, 99.48, 
        99.59, 99.71, 99.81, 99.91, 100 };
    
    // Convert ei to percentage
    ei = ei * 100;
    
    // Translate EI into ES
    auto ei_right_bound_iter = lower_bound(mild.begin(), mild.end(), ei);
    if (ei_right_bound_iter == mild.end() || ei_right_bound_iter == mild.begin()) {
        Rcpp::stop("infection to severe: ei is out of bounds");
    }
    
    unsigned long ei_right_i = ei_right_bound_iter - mild.begin();
    unsigned long ei_left_i = ei_right_i - 1;
    double ei_betweenness = (ei - mild[ei_left_i]) / (mild[ei_right_i] - mild[ei_left_i]);
    
    double es = severe[ei_left_i] * (1.0 - ei_betweenness) + severe[ei_right_i] * ei_betweenness;
    
    return es / 100;
}


