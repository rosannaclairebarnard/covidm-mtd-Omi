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
    return -log(size + 1.) - gsl_sf_lnbeta(size - x + 1., x + 1.) + x * log(prob) + (size - x) * log(1. - prob);
}

// negative binomial log density
double nbinom(unsigned int x, double mean, double size)
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
    return estimate_Rt(P, dyn, P.time0, p, iter, t);
}

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
            dyn(t, p, a, 5)  * (1.0 - P.pop[p].pi_r[a]) +       // E-susceptible to strain 1 in R compartment
            dyn(t, p, a, 10) * (1.0 - P.pop[p].pi_r2[a]) +      // E-susceptible to strain 1 in R2 compartment
            dyn(t, p, a, 15) * (1.0 - P.pop[p].pi_r2[a]) +      // E-susceptible to strain 1 in R3 compartment
            dyn(t, p, a, 16) * (1.0 - P.pop[p].ei_va1[a]) * (1.0 - P.pop[p].ed_va1i[a]) + // E-susceptible to strain 1 in Va1 compartment
            dyn(t, p, a, 17) * (1.0 - P.pop[p].ei_va2[a]) * (1.0 - P.pop[p].ed_va2i[a]) + // E-susceptible to strain 1 in Va2 compartment
            dyn(t, p, a, 18) * (1.0 - P.pop[p].ei_vb1[a]) * (1.0 - P.pop[p].ed_vb1i[a]) + // E-susceptible to strain 1 in Vb1 compartment
            dyn(t, p, a, 19) * (1.0 - P.pop[p].ei_vb2[a]) * (1.0 - P.pop[p].ed_vb2i[a])   // E-susceptible to strain 1 in Vb2 compartment
            ) / P.pop[p].size[a];
        SL[a] = (
            dyn(t, p, a, 16) * (1.0 - P.pop[p].ei_va1[a]) * (P.pop[p].ed_va1i[a] - P.pop[p].et_va1i[a]) +   // L-susceptible to strain 1 in Va1 compartment
            dyn(t, p, a, 17) * (1.0 - P.pop[p].ei_va2[a]) * (P.pop[p].ed_va2i[a] - P.pop[p].et_va2i[a]) +   // L-susceptible to strain 1 in Va2 compartment
            dyn(t, p, a, 18) * (1.0 - P.pop[p].ei_vb1[a]) * (P.pop[p].ed_vb1i[a] - P.pop[p].et_vb1i[a]) +   // L-susceptible to strain 1 in Vb1 compartment
            dyn(t, p, a, 19) * (1.0 - P.pop[p].ei_vb2[a]) * (P.pop[p].ed_vb2i[a] - P.pop[p].et_vb2i[a])     // L-susceptible to strain 1 in Vb2 compartment
            ) / P.pop[p].size[a];
        SE2[a] = (dyn(t, p, a, 0) +                             // E-susceptible to strain 2 in S compartment
            dyn(t, p, a, 5)  * (1.0 - P.pop[p].pi2_r[a]) +      // E-susceptible to strain 2 in R compartment
            dyn(t, p, a, 10) * (1.0 - P.pop[p].pi2_r2[a]) +     // E-susceptible to strain 2 in R2 compartment
            dyn(t, p, a, 15) * (1.0 - P.pop[p].pi2_r3[a]) +     // E-susceptible to strain 2 in R3 compartment
            dyn(t, p, a, 16) * (1.0 - P.pop[p].ei2_va1[a]) * (1.0 - P.pop[p].ed_va1i2[a]) + // E-susceptible to strain 2 in Va1 compartment
            dyn(t, p, a, 17) * (1.0 - P.pop[p].ei2_va2[a]) * (1.0 - P.pop[p].ed_va2i2[a]) + // E-susceptible to strain 2 in Va2 compartment
            dyn(t, p, a, 18) * (1.0 - P.pop[p].ei2_vb1[a]) * (1.0 - P.pop[p].ed_vb1i2[a]) + // E-susceptible to strain 2 in Vb1 compartment
            dyn(t, p, a, 19) * (1.0 - P.pop[p].ei2_vb2[a]) * (1.0 - P.pop[p].ed_vb2i2[a])   // E-susceptible to strain 2 in Vb2 compartment
            ) / P.pop[p].size[a];
        SL2[a] = (
            dyn(t, p, a, 16) * (1.0 - P.pop[p].ei2_va1[a]) * (P.pop[p].ed_va1i2[a] - P.pop[p].et_va1i2[a]) + // L-susceptible to strain 2 in Va1 compartment
            dyn(t, p, a, 17) * (1.0 - P.pop[p].ei2_va2[a]) * (P.pop[p].ed_va2i2[a] - P.pop[p].et_va2i2[a]) + // L-susceptible to strain 2 in Va2 compartment
            dyn(t, p, a, 18) * (1.0 - P.pop[p].ei2_vb1[a]) * (P.pop[p].ed_vb1i2[a] - P.pop[p].et_vb1i2[a]) + // L-susceptible to strain 2 in Vb1 compartment
            dyn(t, p, a, 19) * (1.0 - P.pop[p].ei2_vb2[a]) * (P.pop[p].ed_vb2i2[a] - P.pop[p].et_vb2i2[a])   // L-susceptible to strain 2 in Vb2 compartment
            ) / P.pop[p].size[a];
        SE3[a] = (dyn(t, p, a, 0) +                             // E-susceptible to strain 3 in S compartment
            dyn(t, p, a, 5)  * (1.0 - P.pop[p].pi3_r[a]) +      // E-susceptible to strain 3 in R compartment
            dyn(t, p, a, 10) * (1.0 - P.pop[p].pi3_r2[a]) +     // E-susceptible to strain 3 in R2 compartment
            dyn(t, p, a, 15) * (1.0 - P.pop[p].pi3_r3[a]) +     // E-susceptible to strain 3 in R3 compartment
            dyn(t, p, a, 16) * (1.0 - P.pop[p].ei3_va1[a]) * (1.0 - P.pop[p].ed_va1i3[a]) + // E-susceptible to strain 3 in Va1 compartment
            dyn(t, p, a, 17) * (1.0 - P.pop[p].ei3_va2[a]) * (1.0 - P.pop[p].ed_va2i3[a]) + // E-susceptible to strain 3 in Va2 compartment
            dyn(t, p, a, 18) * (1.0 - P.pop[p].ei3_vb1[a]) * (1.0 - P.pop[p].ed_vb1i3[a]) + // E-susceptible to strain 3 in Vb1 compartment
            dyn(t, p, a, 19) * (1.0 - P.pop[p].ei3_vb2[a]) * (1.0 - P.pop[p].ed_vb2i3[a])   // E-susceptible to strain 3 in Vb2 compartment
            ) / P.pop[p].size[a];
        SL3[a] = (
            dyn(t, p, a, 16) * (1.0 - P.pop[p].ei3_va1[a]) * (P.pop[p].ed_va1i3[a] - P.pop[p].et_va1i3[a]) + // L-susceptible to strain 3 in Va1 compartment
            dyn(t, p, a, 17) * (1.0 - P.pop[p].ei3_va2[a]) * (P.pop[p].ed_va2i3[a] - P.pop[p].et_va2i3[a]) + // L-susceptible to strain 3 in Va2 compartment
            dyn(t, p, a, 18) * (1.0 - P.pop[p].ei3_vb1[a]) * (P.pop[p].ed_vb1i3[a] - P.pop[p].et_vb1i3[a]) + // L-susceptible to strain 3 in Vb1 compartment
            dyn(t, p, a, 19) * (1.0 - P.pop[p].ei3_vb2[a]) * (P.pop[p].ed_vb2i3[a] - P.pop[p].et_vb2i3[a])   // L-susceptible to strain 3 in Vb2 compartment
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
