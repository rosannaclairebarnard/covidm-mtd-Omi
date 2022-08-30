// sim_compartment.cpp

#include "sim_compartment.h"
#include "parameters.h"
#include "reporter.h"
#include "randomizer.h"
#include "user_defined.h"
#include <Rcpp.h>

const bool BOUNDS_CHECK = true;

#define CHECK_NEGATIVE(x) \
    if (x < 0) { \
        cout << "Variable " << #x << " is " << x << " at time " << t << ", group " << a << ", population " << p << "\n." << flush; \
        stop_due_to_negative = true; \
    }

#define CHECK_NEGATIVE2(x,b) \
    if (x < 0) { \
        cout << "Variable " << #x << " is " << x << " at time " << t << ", group " << a << " with group " << b << ", population " << p << "\n." << flush; \
        stop_due_to_negative = true; \
    }

#define CHECK_EFFICACY(x) \
    if (x < 0 || x > 1) { \
        cout << "Efficacy " << #x << " is " << x << " at time " << t << ", group " << a << ", population " << p << "\n." << flush; \
        stop_due_to_negative = true; \
    }                                                          \

//
// MODEL DYNAMICS
//

Population::Population(Parameters& P, unsigned int pindex)
 : seed_row(0), seed_row2(0), seed_row3(0), p(pindex)
{
    // Set up built-in compartments
    N           = P.pop[p].size;
    S           = N;
    R           = vector<double>(S.size(), 0.);
    R2          = vector<double>(S.size(), 0.);
    R3          = vector<double>(S.size(), 0.);
    EV          = vector<double>(S.size(), 0.);
    Va1         = vector<Compartment>(S.size());
    Vb1         = vector<Compartment>(S.size());
    Va2         = vector<Compartment>(S.size());
    Vb2         = vector<Compartment>(S.size());
    Va3         = vector<double>(S.size(), 0.);
    Vb3         = vector<double>(S.size(), 0.);

    E           = vector<Compartment>(S.size());
    L           = vector<Compartment>(S.size());
    Ip          = vector<Compartment>(S.size());
    Ia          = vector<Compartment>(S.size());
    Is          = vector<Compartment>(S.size());

    E2          = vector<Compartment>(S.size());
    L2          = vector<Compartment>(S.size());
    Ip2         = vector<Compartment>(S.size());
    Ia2         = vector<Compartment>(S.size());
    Is2         = vector<Compartment>(S.size());

    E3          = vector<Compartment>(S.size());
    L3          = vector<Compartment>(S.size());
    Ip3         = vector<Compartment>(S.size());
    Ia3         = vector<Compartment>(S.size());
    Is3         = vector<Compartment>(S.size());

    to_death    = vector<Compartment>(S.size());
    to_hosp     = vector<Compartment>(S.size());
    hosp        = vector<Compartment>(S.size());
    to_icu      = vector<Compartment>(S.size());
    icu         = vector<Compartment>(S.size());

    to_death_V1 = vector<Compartment>(S.size());
    to_death_V2 = vector<Compartment>(S.size());
    to_death_V3 = vector<Compartment>(S.size());
    to_hosp_V1  = vector<Compartment>(S.size());
    to_hosp_V2  = vector<Compartment>(S.size());
    to_hosp_V3  = vector<Compartment>(S.size());

    // Initial immunity
    for (unsigned int a = 0; a < S.size(); ++a) {
        double imm = 0;
        if (P.deterministic) imm = S[a] * P.pop[p].imm0[a];
        else                 imm = (unsigned int)(S[a] * P.pop[p].imm0[a] + 0.5);
        S[a] -= imm;
        R[a] += imm;
    }

    // Set up user-specified processes
    unsigned int n_pc = 0;
    for (auto& p : P.processes)
        n_pc += p.ids.size();
    pc = vector<vector<Compartment>>(n_pc, vector<Compartment>(S.size()));
    pci = vector<double>(pc.size(), 0.);
    pco = vector<double>(pc.size(), 0.);
}

// Do seeding and calculate contagiousness
void Population::Contagiousness(Parameters& P, Randomizer& Rand, double t, vector<double>& contag, vector<double>& contag2, vector<double>& contag3)
{
    auto add = [&](unsigned int age, double n)
    {
        n = min(n, S[age]);
        S[age] -= n;
        E[age].Add(P, Rand, n, P.pop[p].dE);
    };

    auto add2 = [&](unsigned int age, double n)
    {
        n = min(n, S[age]);
        S[age] -= n;
        if (P.pop[p].dE2.weights.size() <= 1)
            E2[age].Add(P, Rand, n, P.pop[p].dE);
        else
            E2[age].Add(P, Rand, n, P.pop[p].dE2);
    };

    auto add3 = [&](unsigned int age, double n)
    {
        n = min(n, S[age]);
        S[age] -= n;
        if (P.pop[p].dE3.weights.size() <= 1)
            E3[age].Add(P, Rand, n, P.pop[p].dE);
        else
            E3[age].Add(P, Rand, n, P.pop[p].dE3);
    };

    // Do seeding
    while (seed_row < P.pop[p].seed_times.size() && t >= P.pop[p].seed_times[seed_row])
    {
        if (P.deterministic)
        {
            for (unsigned int a = 0; a < S.size(); ++a)
                add(a, P.pop[p].dist_seed_ages.weights[a]);
        }
        else
        {
            Rand.Multinomial(1, P.pop[p].dist_seed_ages.weights, P.pop[p].dist_seed_ages.storage);
            for (unsigned int a = 0; a < S.size(); ++a)
            {
                if (P.pop[p].dist_seed_ages.storage[a] == 1)
                {
                    add(a, 1);
                    break;
                }
            }
        }
        ++seed_row;
    }

    // Do seeding (2nd strain)
    while (seed_row2 < P.pop[p].seed_times2.size() && t >= P.pop[p].seed_times2[seed_row2])
    {
        if (P.deterministic)
        {
            for (unsigned int a = 0; a < S.size(); ++a)
                add2(a, P.pop[p].dist_seed_ages.weights[a]);
        }
        else
        {
            Rand.Multinomial(1, P.pop[p].dist_seed_ages.weights, P.pop[p].dist_seed_ages.storage);
            for (unsigned int a = 0; a < S.size(); ++a)
            {
                if (P.pop[p].dist_seed_ages.storage[a] == 1)
                {
                    add2(a, 1);
                    break;
                }
            }
        }
        ++seed_row2;
    }

    // Do seeding (3rd strain)
    while (seed_row3 < P.pop[p].seed_times3.size() && t >= P.pop[p].seed_times3[seed_row3])
    {
        if (P.deterministic)
        {
            for (unsigned int a = 0; a < S.size(); ++a)
                add3(a, P.pop[p].dist_seed_ages.weights[a]);
        }
        else
        {
            Rand.Multinomial(1, P.pop[p].dist_seed_ages.weights, P.pop[p].dist_seed_ages.storage);
            for (unsigned int a = 0; a < S.size(); ++a)
            {
                if (P.pop[p].dist_seed_ages.storage[a] == 1)
                {
                    add3(a, 1);
                    break;
                }
            }
        }
        ++seed_row3;
    }

    // Calculate contagiousness from this population
    for (unsigned int a = 0; a < contag.size(); ++a)
        contag[a]  = (N[a] == 0) ? 0 : (P.pop[p].fIp[a] * Ip[a].Size()  + P.pop[p].fIa[a] * Ia[a].Size()  + P.pop[p].fIs[a] * Is[a].Size() ) / N[a];
    for (unsigned int a = 0; a < contag2.size(); ++a)
        contag2[a] = (N[a] == 0) ? 0 : (P.pop[p].fIp[a] * Ip2[a].Size() + P.pop[p].fIa[a] * Ia2[a].Size() + P.pop[p].fIs[a] * Is2[a].Size()) / N[a];
    for (unsigned int a = 0; a < contag3.size(); ++a)
        contag3[a] = (N[a] == 0) ? 0 : (P.pop[p].fIp[a] * Ip3[a].Size() + P.pop[p].fIa[a] * Ia3[a].Size() + P.pop[p].fIs[a] * Is3[a].Size()) / N[a];
}

// Execute one time step's events
void Population::Tick(Parameters& P, Randomizer& Rand, double t, vector<double>& infec, vector<double>& infec2, vector<double>& infec3, Reporter& rep)
{
    // Calculate force of infection in this compartment
    lambda.assign(infec.size(), 0.0);
    lambda2.assign(infec.size(), 0.0);
    lambda3.assign(infec.size(), 0.0);
    for (unsigned int a = 0; a < lambda.size(); ++a)
    {
        for (unsigned int b = 0; b < lambda.size(); ++b)
        {
            lambda[a] += P.pop[p].u[a] * P.pop[p].cm(a,b) * infec[b];
            lambda2[a] += P.pop[p].u2[a] * P.pop[p].cm(a,b) * infec2[b];
            lambda3[a] += P.pop[p].u3[a] * P.pop[p].cm(a,b) * infec3[b];
        }
    }

    // Account for seasonality
    if (P.pop[p].season_A[0] != 0)
    {
        double f = 1.0 + P.pop[p].season_A[0] * cos(2. * M_PI * (t - P.pop[p].season_phi[0]) / P.pop[p].season_T[0]);
        for (unsigned int a = 0; a < lambda.size(); ++a)
        {
            lambda[a]  *= f;
            lambda2[a] *= f;
            lambda3[a] *= f;
        }
    }

    // Account for importation
    // NOTE - assuming no importation of second strain OR third strain
    for (unsigned int a = 0; a < lambda.size(); ++a)
        lambda[a] += P.pop[p].omega[a];

    // Account for transmission adjustment for particle filter
    for (unsigned int a = 0; a < lambda.size(); ++a)
    {
        lambda[a]  *= P.adjustment;
        lambda2[a] *= P.adjustment;
        lambda3[a] *= P.adjustment;
    }

    // Helpers
    auto multinomial = [&](double n, vector<double>& p, vector<double>& nd_out, vector<unsigned int>& ni_out) {
        nd_out.resize(p.size(), 0.);
        if (P.deterministic)
        {
            for (unsigned int i = 0; i < p.size(); ++i)
                nd_out[i] = n * p[i];
        }
        else
        {
            ni_out.resize(p.size(), 0);
            Rand.Multinomial(n, p, ni_out);
            for (unsigned int i = 0; i < p.size(); ++i)
                nd_out[i] = ni_out[i];
        }
    };

    auto poisson = [&](double l) {
        if (P.deterministic)
            return l;
        else
            return (double)Rand.Poisson(l);
    };

    auto binomial = [&](double n, double p) {
        if (P.deterministic)
            return n * p; // min(1.0, n * p)
        else
            return (double)Rand.Binomial(n, p);
    };

    auto num = [&](double n) {
        if (P.deterministic)
            return n;
        else
            return round(n);
    };

    // Do state transitions and reporting for each age group
    for (unsigned int a = 0; a < lambda.size(); ++a)
    {
        // 0. Report prevalences
        if (t == (int)t)
        {
            // Built-in states
            rep(t, p, a, 0) = S[a];
            rep(t, p, a, 1) = E[a].Size();
            rep(t, p, a, 2) = Ip[a].Size();
            rep(t, p, a, 3) = Is[a].Size();
            rep(t, p, a, 4) = Ia[a].Size();
            rep(t, p, a, 5) = R[a];
            rep(t, p, a, 6) = E2[a].Size();
            rep(t, p, a, 7) = Ip2[a].Size();
            rep(t, p, a, 8) = Is2[a].Size();
            rep(t, p, a, 9) = Ia2[a].Size();
            rep(t, p, a, 10) = R2[a];
            rep(t, p, a, 11) = E3[a].Size();
            rep(t, p, a, 12) = Ip3[a].Size();
            rep(t, p, a, 13) = Is3[a].Size();
            rep(t, p, a, 14) = Ia3[a].Size();
            rep(t, p, a, 15) = R3[a];
            rep(t, p, a, 16) = Va1[a].Size();
            rep(t, p, a, 17) = Va2[a].Size();
            rep(t, p, a, 18) = Va3[a];
            rep(t, p, a, 19) = Vb1[a].Size();
            rep(t, p, a, 20) = Vb2[a].Size();
            rep(t, p, a, 21) = Vb3[a];
            rep(t, p, a, 22) = lambda[a];
            rep(t, p, a, 23) = lambda2[a];
            rep(t, p, a, 24) = lambda3[a];
            rep(t, p, a, 27) = hosp[a].Size();
            rep(t, p, a, 29) = icu[a].Size();
            rep(t, p, a, 36) = EV[a];

            // User-specified processes
            for (auto& process : P.processes)
            {
                for (unsigned int i = 0; i < process.p_cols.size(); ++i)
                    rep(t, p, a, process.p_cols[i]) = pc[process.p_ids[i]][a].Size();
            }
        }

        // 1. Built-in states

        // Vaccination, booster vaccines & waning of natural immunity and vaccine protection

        // S -> Va1, Va1 -> S
        double to_vacc_a = max(0.0, min(P.pop[p].va1[a] * P.time_step, N[a] * P.vax_limit - EV[a]));    ///
        double nS_Va1 = 0;                                                                              ///
        if (EV[a] < N[a])                                                                               ///
            nS_Va1 = min(S[a], num(to_vacc_a * min(1.0, S[a] / (N[a] - EV[a]))));                       ///
        double nVa1_S   = Va1[a].RemoveProb(P, Rand, 1.0 - exp(-P.pop[p].wva1[a] * P.time_step));
        EV[a] += to_vacc_a;

        // S -> Vb1, Vb1 -> S
        double to_vacc_b = max(0.0, min((P.pop[p].vb1[a] + P.extra_vb1[a]) * P.time_step, N[a] * P.vax_limit - EV[a]));    ///
        double nS_Vb1 = 0;                                                                              ///
        if (EV[a] < N[a])                                                                               ///
            nS_Vb1 = min(S[a] - nS_Va1, num(to_vacc_b * min(1.0, S[a] / (N[a] - EV[a]))));              ///
        double nVb1_S   = Vb1[a].RemoveProb(P, Rand, 1.0 - exp(-P.pop[p].wvb1[a] * P.time_step));
        EV[a] += to_vacc_b;

        // Va1 -> Va2, Va2 -> S
        double nVa1_Va2 = Va1[a].Mature(P.time_step);
        double nVa2_S   = Va2[a].RemoveProb(P, Rand, 1.0 - exp(-P.pop[p].wva2[a] * P.time_step));

        // Va2 -> Vb2 (a booster from 2-dose AZ takes you back to the start of the 2-dose Pfizer/Moderna compartment)
        double nVa2_exit = Va2[a].Mature(P.time_step);
        double nVa2_Vb2  = nVa2_exit * P.pop[p].p_boost_va2[a];

        // Va2 -> Va3, Va2 -> Vb2 (booster), Va3 -> S
        double nVa2_Va3 = nVa2_exit - nVa2_Vb2;
        double nVa3_S   = binomial(Va3[a], 1.0 - exp(-P.pop[p].wva3[a] * P.time_step));

        // Vb1 -> Vb2, Vb2 -> S
        double nVb1_Vb2 = Vb1[a].Mature(P.time_step);
        double nVb2_S   = Vb2[a].RemoveProb(P, Rand, 1.0 - exp(-P.pop[p].wvb2[a] * P.time_step));

        // Vb2 -> Vb2 (booster dose from 2-dose Pfizer/Moderna takes you back to the start of the same compartment)
        double nVb2_exit = Vb2[a].Mature(P.time_step);
        double nVb2_Vb2  = nVb2_exit * P.pop[p].p_boost_vb2[a];

        // Vb2 -> Vb3, Vb3 -> S
        double nVb2_Vb3 = nVb2_exit - nVb2_Vb2;
        double nVb3_S   = binomial(Vb3[a], 1.0 - exp(-P.pop[p].wvb3[a] * P.time_step));

        S[a]   -= nS_Va1 + nS_Vb1;
        Va1[a].Add(P, Rand, nS_Va1, P.pop[p].dVa1);
        Vb1[a].Add(P, Rand, nS_Vb1, P.pop[p].dVb1);
        S[a]   += nVa1_S + nVb1_S;
        Va2[a].Add(P, Rand, nVa1_Va2, P.pop[p].dVa2);
        S[a]   += nVa2_S;
        Vb2[a].Add(P, Rand, nVb1_Vb2 + nVb2_Vb2 + nVa2_Vb2, P.pop[p].dVb2);
        S[a]   += nVb2_S;
        Va3[a] += nVa2_Va3;
        Va3[a] -= nVa3_S;
        S[a]   += nVa3_S;
        Vb3[a] += nVb2_Vb3;
        Vb3[a] -= nVb3_S;
        S[a]   += nVb3_S;

        // R -> S, R2 -> S, R3 -> S
        double nR_S     = binomial(R[a] , 1.0 - exp(-P.pop[p].wn[a]  * P.time_step));
        double nR2_S    = binomial(R2[a], 1.0 - exp(-P.pop[p].wn2[a] * P.time_step));
        double nR3_S    = binomial(R3[a], 1.0 - exp(-P.pop[p].wn3[a] * P.time_step));

        R[a]   -= nR_S;
        S[a]   += nR_S;
        R2[a]  -= nR2_S;
        S[a]   += nR2_S;
        R3[a]  -= nR3_S;
        S[a]   += nR3_S;

        // Infection rates
        // S -> E, S -> E2, S -> E3
        double lambda_123 = lambda[a] + lambda2[a] + lambda3[a];
        double nS_E123    = binomial(S[a], 1.0 - exp(-lambda_123 * P.time_step));
        double nS_E       = binomial(nS_E123, lambda_123 > 0 ? lambda[a] / lambda_123 : 0);
        double nS_E2E3    = nS_E123 - nS_E;
        double lambda_23  = lambda2[a] + lambda3[a];
        double nS_E2      = binomial(nS_E2E3, lambda_23 > 0 ? lambda2[a] / lambda_23 : 0);
        double nS_E3      = nS_E2E3 - nS_E2;
        
        // Reinfection rates
        // R -> E, R -> E2, R -> E3; R -> L, R -> L2, R -> L3; R -> R; R -> R2; R -> R3
        double foi_r    = lambda[a]  * (1 - P.pop[p].pi_r[a]);
        double foi2_r   = lambda2[a] * (1 - P.pop[p].pi2_r[a]);
        double foi3_r   = lambda3[a] * (1 - P.pop[p].pi3_r[a]);

        double nR_REL123 = binomial(R[a], 1.0 - exp(-(foi_r + foi2_r + foi3_r) * P.time_step));
        double nR_REL    = binomial(nR_REL123, foi_r + foi2_r + foi3_r > 0 ? foi_r / (foi_r + foi2_r + foi3_r) : 0);
        double nR_REL23  = nR_REL123 - nR_REL;
        double nR_REL2   = binomial(nR_REL23, foi2_r + foi3_r > 0 ? foi2_r / (foi2_r + foi3_r) : 0);
        double nR_REL3   = nR_REL23 - nR_REL2;

        double nR_R     = binomial(nR_REL, P.pop[p].pt_ri[a]);
        double nR_X     = binomial(nR_REL - nR_R, (1 - P.pop[p].pd_ri[a]) / (1 - P.pop[p].pt_ri[a]));
        double nR_E     = binomial(nR_REL - nR_R, (1 - P.pop[p].pd_ri[a]));
        double nR_L     = nR_REL - nR_R - nR_E;
        double nR_R2    = binomial(nR_REL2, P.pop[p].pt_ri2[a]);
        double nR_X2    = binomial(nR_REL2 - nR_R2, (1 - P.pop[p].pd_ri2[a]) / (1 - P.pop[p].pt_ri2[a]));
        double nR_E2    = binomial(nR_REL2 - nR_R2, (1 - P.pop[p].pd_ri2[a]));
        double nR_L2    = nR_REL2 - nR_R2 - nR_E2;
        double nR_R3    = binomial(nR_REL3, P.pop[p].pt_ri3[a]);
        double nR_X3    = binomial(nR_REL3 - nR_R3, (1 - P.pop[p].pd_ri3[a]) / (1 - P.pop[p].pt_ri3[a]));
        double nR_E3    = binomial(nR_REL3 - nR_R3, (1 - P.pop[p].pd_ri3[a]));
        double nR_L3    = nR_REL3 - nR_R3 - nR_E3;

        // R2 -> E, R2 -> E2, R2 -> E3; R2 -> L, R2 -> L2, R2 -> L3; R2 -> R; R2 -> R2; R2 -> R3
        double foi_r2    = lambda[a]  * (1 - P.pop[p].pi_r2[a]);
        double foi2_r2   = lambda2[a] * (1 - P.pop[p].pi2_r2[a]);
        double foi3_r2   = lambda3[a] * (1 - P.pop[p].pi3_r2[a]);

        double nR2_REL123 = binomial(R2[a], 1.0 - exp(-(foi_r2 + foi2_r2 + foi3_r2) * P.time_step));
        double nR2_REL    = binomial(nR2_REL123, foi_r2 + foi2_r2 + foi3_r2 > 0 ? foi_r2 / (foi_r2 + foi2_r2 + foi3_r2) : 0);
        double nR2_REL23  = nR2_REL123 - nR2_REL;
        double nR2_REL2   = binomial(nR2_REL23, foi2_r2 + foi3_r2 > 0 ? foi2_r2 / (foi2_r2 + foi3_r2) : 0);
        double nR2_REL3   = nR2_REL23 - nR2_REL2;

        double nR2_R     = binomial(nR2_REL, P.pop[p].pt_r2i[a]);
        double nR2_X     = binomial(nR2_REL - nR2_R, (1 - P.pop[p].pd_r2i[a]) / (1 - P.pop[p].pt_r2i[a]));
        double nR2_E     = binomial(nR2_REL - nR2_R, (1 - P.pop[p].pd_r2i[a]));
        double nR2_L     = nR2_REL - nR2_R - nR2_E;
        double nR2_R2    = binomial(nR2_REL2, P.pop[p].pt_r2i2[a]);
        double nR2_X2    = binomial(nR2_REL2 - nR2_R2, (1 - P.pop[p].pd_r2i2[a]) / (1 - P.pop[p].pt_r2i2[a]));
        double nR2_E2    = binomial(nR2_REL2 - nR2_R2, (1 - P.pop[p].pd_r2i2[a]));
        double nR2_L2    = nR2_REL2 - nR2_R2 - nR2_E2;
        double nR2_R3    = binomial(nR2_REL3, P.pop[p].pt_r2i3[a]);
        double nR2_X3    = binomial(nR2_REL3 - nR2_R3, (1 - P.pop[p].pd_r2i3[a]) / (1 - P.pop[p].pt_r2i3[a]));
        double nR2_E3    = binomial(nR2_REL3 - nR2_R3, (1 - P.pop[p].pd_r2i3[a]));
        double nR2_L3    = nR2_REL3 - nR2_R3 - nR2_E3;

        // R3 -> E, R3 -> E2, R3 -> E3; R3 -> L, R3 -> L2, R3 -> L3; R3 -> R; R3 -> R2; R3 -> R3
        double foi_r3    = lambda[a]  * (1 - P.pop[p].pi_r3[a]);
        double foi2_r3   = lambda2[a] * (1 - P.pop[p].pi2_r3[a]);
        double foi3_r3   = lambda3[a] * (1 - P.pop[p].pi3_r3[a]);

        double nR3_REL123 = binomial(R3[a], 1.0 - exp(-(foi_r3 + foi2_r3 + foi3_r3) * P.time_step));
        double nR3_REL    = binomial(nR3_REL123, foi_r3 + foi2_r3 + foi3_r3 > 0 ? foi_r3 / (foi_r3 + foi2_r3 + foi3_r3) : 0);
        double nR3_REL23  = nR3_REL123 - nR3_REL;
        double nR3_REL2   = binomial(nR3_REL23, foi2_r3 + foi3_r3 > 0 ? foi2_r3 / (foi2_r3 + foi3_r3) : 0);
        double nR3_REL3   = nR3_REL23 - nR3_REL2;

        double nR3_R     = binomial(nR3_REL, P.pop[p].pt_r3i[a]);
        double nR3_X     = binomial(nR3_REL - nR3_R, (1 - P.pop[p].pd_r3i[a]) / (1 - P.pop[p].pt_r3i[a]));
        double nR3_E     = binomial(nR3_REL - nR3_R, (1 - P.pop[p].pd_r3i[a]));
        double nR3_L     = nR3_REL - nR3_R - nR3_E;
        double nR3_R2    = binomial(nR3_REL2, P.pop[p].pt_r3i2[a]);
        double nR3_X2    = binomial(nR3_REL2 - nR3_R2, (1 - P.pop[p].pd_r3i2[a]) / (1 - P.pop[p].pt_r3i2[a]));
        double nR3_E2    = binomial(nR3_REL2 - nR3_R2, (1 - P.pop[p].pd_r3i2[a]));
        double nR3_L2    = nR3_REL2 - nR3_R2 - nR3_E2;
        double nR3_R3    = binomial(nR3_REL3, P.pop[p].pt_r3i3[a]);
        double nR3_X3    = binomial(nR3_REL3 - nR3_R3, (1 - P.pop[p].pd_r3i3[a]) / (1 - P.pop[p].pt_r3i3[a]));
        double nR3_E3    = binomial(nR3_REL3 - nR3_R3, (1 - P.pop[p].pd_r3i3[a]));
        double nR3_L3    = nR3_REL3 - nR3_R3 - nR3_E3;

        // Infection rates following vaccination
        // Va1 -> E, Va1 -> E2, Va1 -> E3, Va1 -> L, Va1 -> L2, Va1 -> L3, Va1 -> R, Va1 -> R2, Va1 -> R3
        double foi_va1    = lambda[a]  * (1 - P.pop[p].ei_va1[a]);
        double foi2_va1   = lambda2[a] * (1 - P.pop[p].ei2_va1[a]);
        double foi3_va1   = lambda3[a] * (1 - P.pop[p].ei3_va1[a]);

        double nVa1_REL123 = Va1[a].RemoveProb(P, Rand, 1.0 - exp(-(foi_va1 + foi2_va1 + foi3_va1) * P.time_step));
        double nVa1_REL    = binomial(nVa1_REL123, foi_va1 + foi2_va1 + foi3_va1 > 0 ? foi_va1 / (foi_va1 + foi2_va1 + foi3_va1) : 0);
        double nVa1_REL23  = nVa1_REL123 - nVa1_REL;
        double nVa1_REL2   = binomial(nVa1_REL23, foi2_va1 + foi3_va1 > 0 ? foi2_va1 / (foi2_va1 + foi3_va1) : 0);
        double nVa1_REL3   = nVa1_REL23 - nVa1_REL2;

        double nVa1_R     = binomial(nVa1_REL, P.pop[p].et_va1i[a]);
        double nVa1_X     = binomial(nVa1_REL - nVa1_R, (1 - P.pop[p].ed_va1i[a]) / (1 - P.pop[p].et_va1i[a]));
        double nVa1_E     = binomial(nVa1_REL - nVa1_R, (1 - P.pop[p].ed_va1i[a]));
        double nVa1_L     = nVa1_REL - nVa1_R - nVa1_E;
        double nVa1_R2    = binomial(nVa1_REL2, P.pop[p].et_va1i2[a]);
        double nVa1_X2    = binomial(nVa1_REL2 - nVa1_R2, (1 - P.pop[p].ed_va1i2[a]) / (1 - P.pop[p].et_va1i2[a]));
        double nVa1_E2    = binomial(nVa1_REL2 - nVa1_R2, (1 - P.pop[p].ed_va1i2[a]));
        double nVa1_L2    = nVa1_REL2 - nVa1_R2 - nVa1_E2;
        double nVa1_R3    = binomial(nVa1_REL3, P.pop[p].et_va1i3[a]);
        double nVa1_X3    = binomial(nVa1_REL3 - nVa1_R3, (1 - P.pop[p].ed_va1i3[a]) / (1 - P.pop[p].et_va1i3[a]));
        double nVa1_E3    = binomial(nVa1_REL3 - nVa1_R3, (1 - P.pop[p].ed_va1i3[a]));
        double nVa1_L3    = nVa1_REL3 - nVa1_R3 - nVa1_E3;

        // Va2 -> E, Va2 -> E2, Va2 -> E3, Va2 -> L, Va2 -> L2, Va2 -> L3, Va2 -> R, Va2 -> R2, Va2 -> R3
        double foi_va2    = lambda[a]  * (1 - P.pop[p].ei_va2[a]);
        double foi2_va2   = lambda2[a] * (1 - P.pop[p].ei2_va2[a]);
        double foi3_va2   = lambda3[a] * (1 - P.pop[p].ei3_va2[a]);

        double nVa2_REL123 = Va2[a].RemoveProb(P, Rand, 1.0 - exp(-(foi_va2 + foi2_va2 + foi3_va2) * P.time_step));
        double nVa2_REL    = binomial(nVa2_REL123, foi_va2 + foi2_va2 + foi3_va2 > 0 ? foi_va2 / (foi_va2 + foi2_va2 + foi3_va2) : 0);
        double nVa2_REL23  = nVa2_REL123 - nVa2_REL;
        double nVa2_REL2   = binomial(nVa2_REL23, foi2_va2 + foi3_va2 > 0 ? foi2_va2 / (foi2_va2 + foi3_va2) : 0);
        double nVa2_REL3   = nVa2_REL23 - nVa2_REL2;

        double nVa2_R     = binomial(nVa2_REL, P.pop[p].et_va2i[a]);
        double nVa2_X     = binomial(nVa2_REL - nVa2_R, (1 - P.pop[p].ed_va2i[a]) / (1 - P.pop[p].et_va2i[a]));
        double nVa2_E     = binomial(nVa2_REL - nVa2_R, (1 - P.pop[p].ed_va2i[a]));
        double nVa2_L     = nVa2_REL - nVa2_R - nVa2_E;
        double nVa2_R2    = binomial(nVa2_REL2, P.pop[p].et_va2i2[a]);
        double nVa2_X2    = binomial(nVa2_REL2 - nVa2_R2, (1 - P.pop[p].ed_va2i2[a]) / (1 - P.pop[p].et_va2i2[a]));
        double nVa2_E2    = binomial(nVa2_REL2 - nVa2_R2, (1 - P.pop[p].ed_va2i2[a]));
        double nVa2_L2    = nVa2_REL2 - nVa2_R2 - nVa2_E2;
        double nVa2_R3    = binomial(nVa2_REL3, P.pop[p].et_va2i3[a]);
        double nVa2_X3    = binomial(nVa2_REL3 - nVa2_R3, (1 - P.pop[p].ed_va2i3[a]) / (1 - P.pop[p].et_va2i3[a]));
        double nVa2_E3    = binomial(nVa2_REL3 - nVa2_R3, (1 - P.pop[p].ed_va2i3[a]));
        double nVa2_L3    = nVa2_REL3 - nVa2_R3 - nVa2_E3;

        // Va3 -> E, Va3 -> E2, Va3 -> E3, Va3 -> L, Va3 -> L2, Va3 -> L3, Va3 -> R, Va3 -> R2, Va3 -> R3
        double foi_va3    = lambda[a]  * (1 - P.pop[p].ei_va3[a]);
        double foi2_va3   = lambda2[a] * (1 - P.pop[p].ei2_va3[a]);
        double foi3_va3   = lambda3[a] * (1 - P.pop[p].ei3_va3[a]);

        double nVa3_REL123 = binomial(Va3[a], 1.0 - exp(-(foi_va3 + foi2_va3 + foi3_va3) * P.time_step));
        double nVa3_REL    = binomial(nVa3_REL123, foi_va3 + foi2_va3 + foi3_va3 > 0 ? foi_va3 / (foi_va3 + foi2_va3 + foi3_va3) : 0);
        double nVa3_REL23  = nVa3_REL123 - nVa3_REL;
        double nVa3_REL2   = binomial(nVa3_REL23, foi2_va3 + foi3_va3 > 0 ? foi2_va3 / (foi2_va3 + foi3_va3) : 0);
        double nVa3_REL3   = nVa3_REL23 - nVa3_REL2;

        double nVa3_R     = binomial(nVa3_REL, P.pop[p].et_va3i[a]);
        double nVa3_X     = binomial(nVa3_REL - nVa3_R, (1 - P.pop[p].ed_va3i[a]) / (1 - P.pop[p].et_va3i[a]));
        double nVa3_E     = binomial(nVa3_REL - nVa3_R, (1 - P.pop[p].ed_va3i[a]));
        double nVa3_L     = nVa3_REL - nVa3_R - nVa3_E;
        double nVa3_R2    = binomial(nVa3_REL2, P.pop[p].et_va3i2[a]);
        double nVa3_X2    = binomial(nVa3_REL2 - nVa3_R2, (1 - P.pop[p].ed_va3i2[a]) / (1 - P.pop[p].et_va3i2[a]));
        double nVa3_E2    = binomial(nVa3_REL2 - nVa3_R2, (1 - P.pop[p].ed_va3i2[a]));
        double nVa3_L2    = nVa3_REL2 - nVa3_R2 - nVa3_E2;
        double nVa3_R3    = binomial(nVa3_REL3, P.pop[p].et_va3i3[a]);
        double nVa3_X3    = binomial(nVa3_REL3 - nVa3_R3, (1 - P.pop[p].ed_va3i3[a]) / (1 - P.pop[p].et_va3i3[a]));
        double nVa3_E3    = binomial(nVa3_REL3 - nVa3_R3, (1 - P.pop[p].ed_va3i3[a]));
        double nVa3_L3    = nVa3_REL3 - nVa3_R3 - nVa3_E3;

        // Vb1 -> E, Vb1 -> E2, Vb1 -> E3, Vb1 -> L, Vb1 -> L2, Vb1 -> L3, Vb1 -> R, Vb1 -> R2, Vb1 -> R3
        double foi_vb1    = lambda[a]  * (1 - P.pop[p].ei_vb1[a]);
        double foi2_vb1   = lambda2[a] * (1 - P.pop[p].ei2_vb1[a]);
        double foi3_vb1   = lambda3[a] * (1 - P.pop[p].ei3_vb1[a]);

        double nVb1_REL123 = Vb1[a].RemoveProb(P, Rand, 1.0 - exp(-(foi_vb1 + foi2_vb1 + foi3_vb1) * P.time_step));
        double nVb1_REL    = binomial(nVb1_REL123, foi_vb1 + foi2_vb1 + foi3_vb1 > 0 ? foi_vb1 / (foi_vb1 + foi2_vb1 + foi3_vb1) : 0);
        double nVb1_REL23  = nVb1_REL123 - nVb1_REL;
        double nVb1_REL2   = binomial(nVb1_REL23, foi2_vb1 + foi3_vb1 > 0 ? foi2_vb1 / (foi2_vb1 + foi3_vb1) : 0);
        double nVb1_REL3   = nVb1_REL23 - nVb1_REL2;

        double nVb1_R     = binomial(nVb1_REL, P.pop[p].et_vb1i[a]);
        double nVb1_X     = binomial(nVb1_REL - nVb1_R, (1 - P.pop[p].ed_vb1i[a]) / (1 - P.pop[p].et_vb1i[a]));
        double nVb1_E     = binomial(nVb1_REL - nVb1_R, (1 - P.pop[p].ed_vb1i[a]));
        double nVb1_L     = nVb1_REL - nVb1_R - nVb1_E;
        double nVb1_R2    = binomial(nVb1_REL2, P.pop[p].et_vb1i2[a]);
        double nVb1_X2    = binomial(nVb1_REL2 - nVb1_R2, (1 - P.pop[p].ed_vb1i2[a]) / (1 - P.pop[p].et_vb1i2[a]));
        double nVb1_E2    = binomial(nVb1_REL2 - nVb1_R2, (1 - P.pop[p].ed_vb1i2[a]));
        double nVb1_L2    = nVb1_REL2 - nVb1_R2 - nVb1_E2;
        double nVb1_R3    = binomial(nVb1_REL3, P.pop[p].et_vb1i3[a]);
        double nVb1_X3    = binomial(nVb1_REL3 - nVb1_R3, (1 - P.pop[p].ed_vb1i3[a]) / (1 - P.pop[p].et_vb1i3[a]));
        double nVb1_E3    = binomial(nVb1_REL3 - nVb1_R3, (1 - P.pop[p].ed_vb1i3[a]));
        double nVb1_L3    = nVb1_REL3 - nVb1_R3 - nVb1_E3;

        // Vb2 -> E, Vb2 -> E2, Vb2 -> E3, Vb2 -> L, Vb2 -> L2, Vb2 -> L3, Vb2 -> R, Vb2 -> R2, Vb2 -> R3
        double foi_vb2    = lambda[a]  * (1 - P.pop[p].ei_vb2[a]);
        double foi2_vb2   = lambda2[a] * (1 - P.pop[p].ei2_vb2[a]);
        double foi3_vb2   = lambda3[a] * (1 - P.pop[p].ei3_vb2[a]);

        double nVb2_REL123 = Vb2[a].RemoveProb(P, Rand, 1.0 - exp(-(foi_vb2 + foi2_vb2 + foi3_vb2) * P.time_step));
        double nVb2_REL    = binomial(nVb2_REL123, foi_vb2 + foi2_vb2 + foi3_vb2 > 0 ? foi_vb2 / (foi_vb2 + foi2_vb2 + foi3_vb2) : 0);
        double nVb2_REL23  = nVb2_REL123 - nVb2_REL;
        double nVb2_REL2   = binomial(nVb2_REL23, foi2_vb2 + foi3_vb2 > 0 ? foi2_vb2 / (foi2_vb2 + foi3_vb2) : 0);
        double nVb2_REL3   = nVb2_REL23 - nVb2_REL2;

        double nVb2_R     = binomial(nVb2_REL, P.pop[p].et_vb2i[a]);
        double nVb2_X     = binomial(nVb2_REL - nVb2_R, (1 - P.pop[p].ed_vb2i[a]) / (1 - P.pop[p].et_vb2i[a]));
        double nVb2_E     = binomial(nVb2_REL - nVb2_R, (1 - P.pop[p].ed_vb2i[a]));
        double nVb2_L     = nVb2_REL - nVb2_R - nVb2_E;
        double nVb2_R2    = binomial(nVb2_REL2, P.pop[p].et_vb2i2[a]);
        double nVb2_X2    = binomial(nVb2_REL2 - nVb2_R2, (1 - P.pop[p].ed_vb2i2[a]) / (1 - P.pop[p].et_vb2i2[a]));
        double nVb2_E2    = binomial(nVb2_REL2 - nVb2_R2, (1 - P.pop[p].ed_vb2i2[a]));
        double nVb2_L2    = nVb2_REL2 - nVb2_R2 - nVb2_E2;
        double nVb2_R3    = binomial(nVb2_REL3, P.pop[p].et_vb2i3[a]);
        double nVb2_X3    = binomial(nVb2_REL3 - nVb2_R3, (1 - P.pop[p].ed_vb2i3[a]) / (1 - P.pop[p].et_vb2i3[a]));
        double nVb2_E3    = binomial(nVb2_REL3 - nVb2_R3, (1 - P.pop[p].ed_vb2i3[a]));
        double nVb2_L3    = nVb2_REL3 - nVb2_R3 - nVb2_E3;

        // Vb3 -> E, Vb3 -> E2, Vb3 -> E3, Vb3 -> L, Vb3 -> L2, Vb3 -> L3, Vb3 -> R, Vb3 -> R2, Vb3 -> R3
        double foi_vb3    = lambda[a]  * (1 - P.pop[p].ei_vb3[a]);
        double foi2_vb3   = lambda2[a] * (1 - P.pop[p].ei2_vb3[a]);
        double foi3_vb3   = lambda3[a] * (1 - P.pop[p].ei3_vb3[a]);

        double nVb3_REL123 = binomial(Vb3[a], 1.0 - exp(-(foi_vb3 + foi2_vb3 + foi3_vb3) * P.time_step));
        double nVb3_REL    = binomial(nVb3_REL123, foi_vb3 + foi2_vb3 + foi3_vb3 > 0 ? foi_vb3 / (foi_vb3 + foi2_vb3 + foi3_vb3) : 0);
        double nVb3_REL23  = nVb3_REL123 - nVb3_REL;
        double nVb3_REL2   = binomial(nVb3_REL23, foi2_vb3 + foi3_vb3 > 0 ? foi2_vb3 / (foi2_vb3 + foi3_vb3) : 0);
        double nVb3_REL3   = nVb3_REL23 - nVb3_REL2;

        double nVb3_R     = binomial(nVb3_REL, P.pop[p].et_vb3i[a]);
        double nVb3_X     = binomial(nVb3_REL - nVb3_R, (1 - P.pop[p].ed_vb3i[a]) / (1 - P.pop[p].et_vb3i[a]));
        double nVb3_E     = binomial(nVb3_REL - nVb3_R, (1 - P.pop[p].ed_vb3i[a]));
        double nVb3_L     = nVb3_REL - nVb3_R - nVb3_E;
        double nVb3_R2    = binomial(nVb3_REL2, P.pop[p].et_vb3i2[a]);
        double nVb3_X2    = binomial(nVb3_REL2 - nVb3_R2, (1 - P.pop[p].ed_vb3i2[a]) / (1 - P.pop[p].et_vb3i2[a]));
        double nVb3_E2    = binomial(nVb3_REL2 - nVb3_R2, (1 - P.pop[p].ed_vb3i2[a]));
        double nVb3_L2    = nVb3_REL2 - nVb3_R2 - nVb3_E2;
        double nVb3_R3    = binomial(nVb3_REL3, P.pop[p].et_vb3i3[a]);
        double nVb3_X3    = binomial(nVb3_REL3 - nVb3_R3, (1 - P.pop[p].ed_vb3i3[a]) / (1 - P.pop[p].et_vb3i3[a]));
        double nVb3_E3    = binomial(nVb3_REL3 - nVb3_R3, (1 - P.pop[p].ed_vb3i3[a]));
        double nVb3_L3    = nVb3_REL3 - nVb3_R3 - nVb3_E3;

        // Infection / reinfection / infection following vaccination changes
        S[a]   -= nS_E123;
        R[a]   -= nR_REL123;
        R2[a]  -= nR2_REL123;
        R3[a]  -= nR3_REL123;
        Va3[a] -= nVa3_REL123;
        Vb3[a] -= nVb3_REL123;

        // Burden compartments
        double fatal    = binomial(nS_E,    P.pop[p].ifr1[a]) +
                          binomial(nVa1_X,  P.pop[p].ifr1[a] * (1. - P.pop[p].em_va1d[a])) +
                          binomial(nVa2_X,  P.pop[p].ifr1[a] * (1  - P.pop[p].em_va2d[a])) +
                          binomial(nVa3_X,  P.pop[p].ifr1[a] * (1  - P.pop[p].em_va3d[a])) +
                          binomial(nVb1_X,  P.pop[p].ifr1[a] * (1  - P.pop[p].em_vb1d[a])) +
                          binomial(nVb2_X,  P.pop[p].ifr1[a] * (1  - P.pop[p].em_vb2d[a])) +
                          binomial(nVb3_X,  P.pop[p].ifr1[a] * (1  - P.pop[p].em_vb3d[a])) +
                          binomial(nR_X,    P.pop[p].ifr1[a] * (1  - P.pop[p].pm_rd[a])) +
                          binomial(nR2_X,   P.pop[p].ifr1[a] * (1  - P.pop[p].pm_r2d[a])) +
                          binomial(nR3_X,   P.pop[p].ifr1[a] * (1  - P.pop[p].pm_r3d[a])) +
                          binomial(nS_E2,   P.pop[p].ifr2[a]) +
                          binomial(nVa1_X2, P.pop[p].ifr2[a] * (1. - P.pop[p].em_va1d2[a])) +
                          binomial(nVa2_X2, P.pop[p].ifr2[a] * (1  - P.pop[p].em_va2d2[a])) +
                          binomial(nVa3_X2, P.pop[p].ifr2[a] * (1  - P.pop[p].em_va3d2[a])) +
                          binomial(nVb1_X2, P.pop[p].ifr2[a] * (1  - P.pop[p].em_vb1d2[a])) +
                          binomial(nVb2_X2, P.pop[p].ifr2[a] * (1  - P.pop[p].em_vb2d2[a])) +
                          binomial(nVb3_X2, P.pop[p].ifr2[a] * (1  - P.pop[p].em_vb3d2[a])) +
                          binomial(nR_X2,   P.pop[p].ifr2[a] * (1  - P.pop[p].pm_rd2[a])) +
                          binomial(nR2_X2,  P.pop[p].ifr2[a] * (1  - P.pop[p].pm_r2d2[a])) +
                          binomial(nR3_X2,  P.pop[p].ifr2[a] * (1  - P.pop[p].pm_r3d2[a])) +
                          binomial(nS_E3,   P.pop[p].ifr3[a]) +
                          binomial(nVa1_X3, P.pop[p].ifr3[a] * (1. - P.pop[p].em_va1d3[a])) +
                          binomial(nVa2_X3, P.pop[p].ifr3[a] * (1  - P.pop[p].em_va2d3[a])) +
                          binomial(nVa3_X3, P.pop[p].ifr3[a] * (1  - P.pop[p].em_va3d3[a])) +
                          binomial(nVb1_X3, P.pop[p].ifr3[a] * (1  - P.pop[p].em_vb1d3[a])) +
                          binomial(nVb2_X3, P.pop[p].ifr3[a] * (1  - P.pop[p].em_vb2d3[a])) +
                          binomial(nVb3_X3, P.pop[p].ifr3[a] * (1  - P.pop[p].em_vb3d3[a])) +
                          binomial(nR_X3,   P.pop[p].ifr3[a] * (1  - P.pop[p].pm_rd3[a])) +
                          binomial(nR2_X3,  P.pop[p].ifr3[a] * (1  - P.pop[p].pm_r2d3[a])) +
                          binomial(nR3_X3,  P.pop[p].ifr3[a] * (1  - P.pop[p].pm_r3d3[a]));
        double severe   = binomial(nS_E,    P.pop[p].ihr1[a]) +
                          binomial(nVa1_X,  P.pop[p].ihr1[a] * (1. - P.pop[p].eh_va1d[a])) +
                          binomial(nVa2_X,  P.pop[p].ihr1[a] * (1  - P.pop[p].eh_va2d[a])) +
                          binomial(nVa3_X,  P.pop[p].ihr1[a] * (1  - P.pop[p].eh_va3d[a])) +
                          binomial(nVb1_X,  P.pop[p].ihr1[a] * (1  - P.pop[p].eh_vb1d[a])) +
                          binomial(nVb2_X,  P.pop[p].ihr1[a] * (1  - P.pop[p].eh_vb2d[a])) +
                          binomial(nVb3_X,  P.pop[p].ihr1[a] * (1  - P.pop[p].eh_vb3d[a])) +
                          binomial(nR_X,    P.pop[p].ihr1[a] * (1  - P.pop[p].ph_rd[a])) +
                          binomial(nR2_X,   P.pop[p].ihr1[a] * (1  - P.pop[p].ph_r2d[a])) +
                          binomial(nR3_X,   P.pop[p].ihr1[a] * (1  - P.pop[p].ph_r3d[a])) +
                          binomial(nS_E2,   P.pop[p].ihr2[a]) +
                          binomial(nVa1_X2, P.pop[p].ihr2[a] * (1. - P.pop[p].eh_va1d2[a])) +
                          binomial(nVa2_X2, P.pop[p].ihr2[a] * (1  - P.pop[p].eh_va2d2[a])) +
                          binomial(nVa3_X2, P.pop[p].ihr2[a] * (1  - P.pop[p].eh_va3d2[a])) +
                          binomial(nVb1_X2, P.pop[p].ihr2[a] * (1  - P.pop[p].eh_vb1d2[a])) +
                          binomial(nVb2_X2, P.pop[p].ihr2[a] * (1  - P.pop[p].eh_vb2d2[a])) +
                          binomial(nVb3_X2, P.pop[p].ihr2[a] * (1  - P.pop[p].eh_vb3d2[a])) +
                          binomial(nR_X2,   P.pop[p].ihr2[a] * (1  - P.pop[p].ph_rd2[a])) +
                          binomial(nR2_X2,  P.pop[p].ihr2[a] * (1  - P.pop[p].ph_r2d2[a])) +
                          binomial(nR3_X2,  P.pop[p].ihr2[a] * (1  - P.pop[p].ph_r3d2[a])) +
                          binomial(nS_E3,   P.pop[p].ihr3[a]) +
                          binomial(nVa1_X3, P.pop[p].ihr3[a] * (1. - P.pop[p].eh_va1d3[a])) +
                          binomial(nVa2_X3, P.pop[p].ihr3[a] * (1  - P.pop[p].eh_va2d3[a])) +
                          binomial(nVa3_X3, P.pop[p].ihr3[a] * (1  - P.pop[p].eh_va3d3[a])) +
                          binomial(nVb1_X3, P.pop[p].ihr3[a] * (1  - P.pop[p].eh_vb1d3[a])) +
                          binomial(nVb2_X3, P.pop[p].ihr3[a] * (1  - P.pop[p].eh_vb2d3[a])) +
                          binomial(nVb3_X3, P.pop[p].ihr3[a] * (1  - P.pop[p].eh_vb3d3[a])) +
                          binomial(nR_X3,   P.pop[p].ihr3[a] * (1  - P.pop[p].ph_rd3[a])) +
                          binomial(nR2_X3,  P.pop[p].ihr3[a] * (1  - P.pop[p].ph_r2d3[a])) +
                          binomial(nR3_X3,  P.pop[p].ihr3[a] * (1  - P.pop[p].ph_r3d3[a]));
        double critical = binomial(nS_E,    P.pop[p].iir1[a]) +
                          binomial(nVa1_X,  P.pop[p].iir1[a] * (1. - P.pop[p].eh_va1d[a])) +
                          binomial(nVa2_X,  P.pop[p].iir1[a] * (1  - P.pop[p].eh_va2d[a])) +
                          binomial(nVa3_X,  P.pop[p].iir1[a] * (1  - P.pop[p].eh_va3d[a])) +
                          binomial(nVb1_X,  P.pop[p].iir1[a] * (1  - P.pop[p].eh_vb1d[a])) +
                          binomial(nVb2_X,  P.pop[p].iir1[a] * (1  - P.pop[p].eh_vb2d[a])) +
                          binomial(nVb3_X,  P.pop[p].iir1[a] * (1  - P.pop[p].eh_vb3d[a])) +
                          binomial(nR_X,    P.pop[p].iir1[a] * (1  - P.pop[p].ph_rd[a])) +
                          binomial(nR2_X,   P.pop[p].iir1[a] * (1  - P.pop[p].ph_r2d[a])) +
                          binomial(nR3_X,   P.pop[p].iir1[a] * (1  - P.pop[p].ph_r3d[a])) +
                          binomial(nS_E2,   P.pop[p].iir2[a]) +
                          binomial(nVa1_X2, P.pop[p].iir2[a] * (1. - P.pop[p].eh_va1d2[a])) +
                          binomial(nVa2_X2, P.pop[p].iir2[a] * (1  - P.pop[p].eh_va2d2[a])) +
                          binomial(nVa3_X2, P.pop[p].iir2[a] * (1  - P.pop[p].eh_va3d2[a])) +
                          binomial(nVb1_X2, P.pop[p].iir2[a] * (1  - P.pop[p].eh_vb1d2[a])) +
                          binomial(nVb2_X2, P.pop[p].iir2[a] * (1  - P.pop[p].eh_vb2d2[a])) +
                          binomial(nVb3_X2, P.pop[p].iir2[a] * (1  - P.pop[p].eh_vb3d2[a])) +
                          binomial(nR_X2,   P.pop[p].iir2[a] * (1  - P.pop[p].ph_rd2[a])) +
                          binomial(nR2_X2,  P.pop[p].iir2[a] * (1  - P.pop[p].ph_r2d2[a])) +
                          binomial(nR3_X2,  P.pop[p].iir2[a] * (1  - P.pop[p].ph_r3d2[a])) +
                          binomial(nS_E3,   P.pop[p].iir3[a]) +
                          binomial(nVa1_X3, P.pop[p].iir3[a] * (1. - P.pop[p].eh_va1d3[a])) +
                          binomial(nVa2_X3, P.pop[p].iir3[a] * (1  - P.pop[p].eh_va2d3[a])) +
                          binomial(nVa3_X3, P.pop[p].iir3[a] * (1  - P.pop[p].eh_va3d3[a])) +
                          binomial(nVb1_X3, P.pop[p].iir3[a] * (1  - P.pop[p].eh_vb1d3[a])) +
                          binomial(nVb2_X3, P.pop[p].iir3[a] * (1  - P.pop[p].eh_vb2d3[a])) +
                          binomial(nVb3_X3, P.pop[p].iir3[a] * (1  - P.pop[p].eh_vb3d3[a])) +
                          binomial(nR_X3,   P.pop[p].iir3[a] * (1  - P.pop[p].ph_rd3[a])) +
                          binomial(nR2_X3,  P.pop[p].iir3[a] * (1  - P.pop[p].ph_r2d3[a])) +
                          binomial(nR3_X3,  P.pop[p].iir3[a] * (1  - P.pop[p].ph_r3d3[a]));

        to_death[a].Add(P, Rand, fatal,    P.pop[p].dDeath);
        to_hosp [a].Add(P, Rand, severe,   P.pop[p].dHosp);
        to_icu  [a].Add(P, Rand, critical, P.pop[p].dICU);

        double died = to_death[a].Mature(P.time_step);
        double hosp_admissions = to_hosp[a].Mature(P.time_step);
        double icu_admissions = to_icu[a].Mature(P.time_step);
        if (P.pop[p].lHosp2.empty())
            hosp[a].Add(P, Rand, hosp_admissions, P.pop[p].lHosp);
        else
            hosp[a].Add(P, Rand, hosp_admissions, P.pop[p].lHosp2[a]);
        icu [a].Add(P, Rand, icu_admissions,  P.pop[p].lICU);
        hosp[a].Mature(P.time_step);
        icu [a].Mature(P.time_step);

        // Tracking burdens among vaccinated individuals
        double fatal_V1   = binomial(nVa1_X,  P.pop[p].ifr1[a] * (1. - P.pop[p].em_va1d [a])) +
                            binomial(nVb1_X,  P.pop[p].ifr1[a] * (1. - P.pop[p].em_vb1d [a])) +
                            binomial(nVa1_X2, P.pop[p].ifr2[a] * (1. - P.pop[p].em_va1d2[a])) +
                            binomial(nVb1_X2, P.pop[p].ifr2[a] * (1. - P.pop[p].em_vb1d2[a])) +
                            binomial(nVa1_X3, P.pop[p].ifr3[a] * (1. - P.pop[p].em_va1d3[a])) +
                            binomial(nVb1_X3, P.pop[p].ifr3[a] * (1. - P.pop[p].em_vb1d3[a]));
        double fatal_V2   = binomial(nVa2_X,  P.pop[p].ifr1[a] * (1. - P.pop[p].em_va2d [a])) +
                            binomial(nVb2_X,  P.pop[p].ifr1[a] * (1. - P.pop[p].em_vb2d [a])) +
                            binomial(nVa2_X2, P.pop[p].ifr2[a] * (1. - P.pop[p].em_va2d2[a])) +
                            binomial(nVb2_X2, P.pop[p].ifr2[a] * (1. - P.pop[p].em_vb2d2[a])) +
                            binomial(nVa2_X3, P.pop[p].ifr3[a] * (1. - P.pop[p].em_va2d3[a])) +
                            binomial(nVb2_X3, P.pop[p].ifr3[a] * (1. - P.pop[p].em_vb2d3[a]));
        double fatal_V3   = binomial(nVa3_X,  P.pop[p].ifr1[a] * (1. - P.pop[p].em_va3d [a])) +
                            binomial(nVb3_X,  P.pop[p].ifr1[a] * (1. - P.pop[p].em_vb3d [a])) +
                            binomial(nVa3_X2, P.pop[p].ifr2[a] * (1. - P.pop[p].em_va3d2[a])) +
                            binomial(nVb3_X2, P.pop[p].ifr2[a] * (1. - P.pop[p].em_vb3d2[a])) +
                            binomial(nVa3_X3, P.pop[p].ifr3[a] * (1. - P.pop[p].em_va3d3[a])) +
                            binomial(nVb3_X3, P.pop[p].ifr3[a] * (1. - P.pop[p].em_vb3d3[a]));

        double severe_V1  = binomial(nVa1_X,  P.pop[p].ihr1[a] * (1. - P.pop[p].eh_va1d [a])) +
                            binomial(nVb1_X,  P.pop[p].ihr1[a] * (1. - P.pop[p].eh_vb1d [a])) +
                            binomial(nVa1_X2, P.pop[p].ihr2[a] * (1. - P.pop[p].eh_va1d2[a])) +
                            binomial(nVb1_X2, P.pop[p].ihr2[a] * (1. - P.pop[p].eh_vb1d2[a])) +
                            binomial(nVa1_X3, P.pop[p].ihr3[a] * (1. - P.pop[p].eh_va1d3[a])) +
                            binomial(nVb1_X3, P.pop[p].ihr3[a] * (1. - P.pop[p].eh_vb1d3[a]));
        double severe_V2  = binomial(nVa2_X,  P.pop[p].ihr1[a] * (1. - P.pop[p].eh_va2d [a])) +
                            binomial(nVb2_X,  P.pop[p].ihr1[a] * (1. - P.pop[p].eh_vb2d [a])) +
                            binomial(nVa2_X2, P.pop[p].ihr2[a] * (1. - P.pop[p].eh_va2d2[a])) +
                            binomial(nVb2_X2, P.pop[p].ihr2[a] * (1. - P.pop[p].eh_vb2d2[a])) +
                            binomial(nVa2_X3, P.pop[p].ihr3[a] * (1. - P.pop[p].eh_va2d3[a])) +
                            binomial(nVb2_X3, P.pop[p].ihr3[a] * (1. - P.pop[p].eh_vb2d3[a]));
        double severe_V3  = binomial(nVa3_X,  P.pop[p].ihr1[a] * (1. - P.pop[p].eh_va3d [a])) +
                            binomial(nVb3_X,  P.pop[p].ihr1[a] * (1. - P.pop[p].eh_vb3d [a])) +
                            binomial(nVa3_X2, P.pop[p].ihr2[a] * (1. - P.pop[p].eh_va3d2[a])) +
                            binomial(nVb3_X2, P.pop[p].ihr2[a] * (1. - P.pop[p].eh_vb3d2[a])) +
                            binomial(nVa3_X3, P.pop[p].ihr3[a] * (1. - P.pop[p].eh_va3d3[a])) +
                            binomial(nVb3_X3, P.pop[p].ihr3[a] * (1. - P.pop[p].eh_vb3d3[a]));

        to_death_V1[a].Add(P, Rand, fatal_V1,    P.pop[p].dDeath);
        to_death_V2[a].Add(P, Rand, fatal_V2,    P.pop[p].dDeath);
        to_death_V3[a].Add(P, Rand, fatal_V3,    P.pop[p].dDeath);
        to_hosp_V1 [a].Add(P, Rand, severe_V1,   P.pop[p].dHosp);
        to_hosp_V2 [a].Add(P, Rand, severe_V2,   P.pop[p].dHosp);
        to_hosp_V3 [a].Add(P, Rand, severe_V3,   P.pop[p].dHosp);

        double died_V1 = to_death_V1[a].Mature(P.time_step);
        double died_V2 = to_death_V2[a].Mature(P.time_step);
        double died_V3 = to_death_V3[a].Mature(P.time_step);
        double hosp_admissions_V1 = to_hosp_V1[a].Mature(P.time_step);
        double hosp_admissions_V2 = to_hosp_V2[a].Mature(P.time_step);
        double hosp_admissions_V3 = to_hosp_V3[a].Mature(P.time_step);

        // Add new exposed.
        E[a].Add(P, Rand, nS_E + nR_E + nR2_E + nR3_E + nVa1_E + nVa2_E + nVa3_E + nVb1_E + nVb2_E + nVb3_E,  P.pop[p].dE);
        // If dE2 has been supplied, use that as latent period for strain 2, otherwise use dE.
        if (P.pop[p].dE2.weights.size() <= 1)
            E2[a].Add(P, Rand, nS_E2 + nR_E2 + nR2_E2 + nR3_E2 + nVa1_E2 + nVa2_E2 + nVa3_E2 + nVb1_E2 + nVb2_E2 + nVb3_E2, P.pop[p].dE);
        else
            E2[a].Add(P, Rand, nS_E2 + nR_E2 + nR2_E2 + nR3_E2 + nVa1_E2 + nVa2_E2 + nVa3_E2 + nVb1_E2 + nVb2_E2 + nVb3_E2, P.pop[p].dE2);
        // If dE3 has been supplied, use that as latent period for strain 3, otherwise use dE.
        if (P.pop[p].dE3.weights.size() <= 1)
            E3[a].Add(P, Rand, nS_E3 + nR_E3 + nR2_E3 + nR3_E3 + nVa1_E3 + nVa2_E3 + nVa3_E3 + nVb1_E3 + nVb2_E3 + nVb3_E3, P.pop[p].dE);
        else
            E3[a].Add(P, Rand, nS_E3 + nR_E3 + nR2_E3 + nR3_E3 + nVa1_E3 + nVa2_E3 + nVa3_E3 + nVb1_E3 + nVb2_E3 + nVb3_E3, P.pop[p].dE3);

        // Add new "latent" (exposed to infection, but protected against disease by vaccine)
        L[a].Add(P, Rand, nR_L + nR2_L + nR3_L + nVa1_L + nVa2_L + nVa3_L + nVb1_L + nVb2_L + nVb3_L, P.pop[p].dE);
        // If dE2 has been supplied, use that as latent period for strain 2, otherwise use dE.
        if (P.pop[p].dE2.weights.size() <= 1)
            L2[a].Add(P, Rand, nR_L2 + nR2_L2 + nR3_L2 + nVa1_L2 + nVa2_L2 + nVa3_L2 + nVb1_L2 + nVb2_L2 + nVb3_L2, P.pop[p].dE);
        else
            L2[a].Add(P, Rand, nR_L2 + nR2_L2 + nR3_L2 + nVa1_L2 + nVa2_L2 + nVa3_L2 + nVb1_L2 + nVb2_L2 + nVb3_L2, P.pop[p].dE2);
        // If dE3 has been supplied, use that as latent period for strain 3, otherwise use dE.
        if (P.pop[p].dE3.weights.size() <= 1)
            L3[a].Add(P, Rand, nR_L3 + nR2_L3 + nR3_L3 + nVa1_L3 + nVa2_L3 + nVa3_L3 + nVb1_L3 + nVb2_L3 + nVb3_L3, P.pop[p].dE);
        else
            L3[a].Add(P, Rand, nR_L3 + nR2_L3 + nR3_L3 + nVa1_L3 + nVa2_L3 + nVa3_L3 + nVb1_L3 + nVb2_L3 + nVb3_L3, P.pop[p].dE3);

        // Add new "skipped" (exposed to infection, but protected against transmitting by vaccine; go straight to R)
        R [a] += nR_R  + nR2_R  + nR3_R  + nVa1_R  + nVa2_R  + nVa3_R  + nVb1_R  + nVb2_R  + nVb3_R;
        R2[a] += nR_R2 + nR2_R2 + nR3_R2 + nVa1_R2 + nVa2_R2 + nVa3_R2 + nVb1_R2 + nVb2_R2 + nVb3_R2;
        R3[a] += nR_R3 + nR2_R3 + nR3_R3 + nVa1_R3 + nVa2_R3 + nVa3_R3 + nVb1_R3 + nVb2_R3 + nVb3_R3;

        // Strain 1
        // E -> Ip/Ia, L -> Ia
        double nE_Ipa = E[a].Mature(P.time_step);
        double nE_Ip  = binomial(nE_Ipa, P.pop[p].y[a]);
        double nE_Ia  = nE_Ipa - nE_Ip;
        double nL_Ia  = L[a].Mature(P.time_step);
        Ip[a].Add(P, Rand, nE_Ip, P.pop[p].dIp);
        Ia[a].Add(P, Rand, nE_Ia + nL_Ia, P.pop[p].dIa);

        // Ip -> Is
        double nIp_Is = Ip[a].Mature(P.time_step);
        Is[a].Add(P, Rand, nIp_Is, P.pop[p].dIs);

        // Is -> R
        double nIs_R = Is[a].Mature(P.time_step);
        R[a] += nIs_R;

        // Ia -> R
        double nIa_R = Ia[a].Mature(P.time_step);
        R[a] += nIa_R;

        // Strain 2
        // E2 -> Ip2/Ia2, L2 -> Ia2
        double nE2_Ipa2 = E2[a].Mature(P.time_step);
        double nE2_Ip2  = binomial(nE2_Ipa2, P.pop[p].y2[a]);
        double nE2_Ia2  = nE2_Ipa2 - nE2_Ip2;
        double nL2_Ia2  = L2[a].Mature(P.time_step);
        Ip2[a].Add(P, Rand, nE2_Ip2, (P.pop[p].dIp2.weights.size() > 1 ? P.pop[p].dIp2 : P.pop[p].dIp));
        Ia2[a].Add(P, Rand, nE2_Ia2 + nL2_Ia2, (P.pop[p].dIa2.weights.size() > 1 ? P.pop[p].dIa2 : P.pop[p].dIa));

        // Ip2 -> Is2
        double nIp2_Is2 = Ip2[a].Mature(P.time_step);
        Is2[a].Add(P, Rand, nIp2_Is2, (P.pop[p].dIs2.weights.size() > 1 ? P.pop[p].dIs2 : P.pop[p].dIs));

        // Is2 -> R2
        double nIs2_R2 = Is2[a].Mature(P.time_step);
        R2[a] += nIs2_R2;

        // Ia2 -> R2
        double nIa2_R2 = Ia2[a].Mature(P.time_step);
        R2[a] += nIa2_R2;

        // Strain 3
        // E3 -> Ip3/Ia3, L3 -> Ia3
        double nE3_Ipa3 = E3[a].Mature(P.time_step);
        double nE3_Ip3  = binomial(nE3_Ipa3, P.pop[p].y3[a]);
        double nE3_Ia3  = nE3_Ipa3 - nE3_Ip3;
        double nL3_Ia3  = L3[a].Mature(P.time_step);
        Ip3[a].Add(P, Rand, nE3_Ip3, (P.pop[p].dIp3.weights.size() > 1 ? P.pop[p].dIp3 : P.pop[p].dIp));
        Ia3[a].Add(P, Rand, nE3_Ia3 + nL3_Ia3, (P.pop[p].dIa3.weights.size() > 1 ? P.pop[p].dIa3 : P.pop[p].dIa));

        // Ip3 -> Is3
        double nIp3_Is3 = Ip3[a].Mature(P.time_step);
        Is3[a].Add(P, Rand, nIp3_Is3, (P.pop[p].dIs3.weights.size() > 1 ? P.pop[p].dIs3 : P.pop[p].dIs));

        // Is3 -> R3
        double nIs3_R3 = Is3[a].Mature(P.time_step);
        R3[a] += nIs3_R3;

        // Ia3 -> R3
        double nIa3_R3 = Ia3[a].Mature(P.time_step);
        R3[a] += nIa3_R3;
        
        if (BOUNDS_CHECK) {
            // -- Check for negative flows / variables --
            bool stop_due_to_negative = false;

            for (unsigned int b = 0; b < lambda.size(); ++b)
            {
                CHECK_NEGATIVE(P.pop[p].u[a])
                CHECK_NEGATIVE(P.pop[p].u2[a])
                CHECK_NEGATIVE(P.pop[p].u3[a])
                CHECK_NEGATIVE2(P.pop[p].cm(a,b), b)
                CHECK_NEGATIVE2(infec[b], b)
                CHECK_NEGATIVE2(infec2[b], b)
                CHECK_NEGATIVE2(infec3[b], b)
            }
            CHECK_NEGATIVE(lambda[a])
            CHECK_NEGATIVE(lambda2[a])
            CHECK_NEGATIVE(lambda3[a])

            CHECK_NEGATIVE(foi2_r)
            CHECK_NEGATIVE(foi2_r2)
            CHECK_NEGATIVE(foi2_r3)
            CHECK_NEGATIVE(foi2_va1)
            CHECK_NEGATIVE(foi2_va2)
            CHECK_NEGATIVE(foi2_va3)
            CHECK_NEGATIVE(foi2_vb1)
            CHECK_NEGATIVE(foi2_vb2)
            CHECK_NEGATIVE(foi2_vb3)
            CHECK_NEGATIVE(foi3_r)
            CHECK_NEGATIVE(foi3_r2)
            CHECK_NEGATIVE(foi3_r3)
            CHECK_NEGATIVE(foi3_va1)
            CHECK_NEGATIVE(foi3_va2)
            CHECK_NEGATIVE(foi3_va3)
            CHECK_NEGATIVE(foi3_vb1)
            CHECK_NEGATIVE(foi3_vb2)
            CHECK_NEGATIVE(foi3_vb3)
            CHECK_NEGATIVE(foi_r)
            CHECK_NEGATIVE(foi_r2)
            CHECK_NEGATIVE(foi_r3)
            CHECK_NEGATIVE(foi_va1)
            CHECK_NEGATIVE(foi_va2)
            CHECK_NEGATIVE(foi_va3)
            CHECK_NEGATIVE(foi_vb1)
            CHECK_NEGATIVE(foi_vb2)
            CHECK_NEGATIVE(foi_vb3)

            CHECK_NEGATIVE(N[a])
            CHECK_NEGATIVE(S[a])
            CHECK_NEGATIVE(R[a])
            CHECK_NEGATIVE(R2[a])
            CHECK_NEGATIVE(R3[a])
            CHECK_NEGATIVE(EV[a])
            CHECK_NEGATIVE(Va3[a])
            CHECK_NEGATIVE(Vb3[a])

            CHECK_NEGATIVE(Va1[a].Size())
            CHECK_NEGATIVE(Va2[a].Size())
            CHECK_NEGATIVE(Vb1[a].Size())
            CHECK_NEGATIVE(Vb2[a].Size())
            CHECK_NEGATIVE(E[a].Size())
            CHECK_NEGATIVE(L[a].Size())
            CHECK_NEGATIVE(Ip[a].Size())
            CHECK_NEGATIVE(Ia[a].Size())
            CHECK_NEGATIVE(Is[a].Size())
            CHECK_NEGATIVE(E2[a].Size())
            CHECK_NEGATIVE(L2[a].Size())
            CHECK_NEGATIVE(Ip2[a].Size())
            CHECK_NEGATIVE(Ia2[a].Size())
            CHECK_NEGATIVE(Is2[a].Size())
            CHECK_NEGATIVE(E3[a].Size())
            CHECK_NEGATIVE(L3[a].Size())
            CHECK_NEGATIVE(Ip3[a].Size())
            CHECK_NEGATIVE(Ia3[a].Size())
            CHECK_NEGATIVE(Is3[a].Size())
            CHECK_NEGATIVE(to_death[a].Size())
            CHECK_NEGATIVE(to_hosp[a].Size())
            CHECK_NEGATIVE(hosp[a].Size())
            CHECK_NEGATIVE(to_icu[a].Size())
            CHECK_NEGATIVE(icu[a].Size())
            CHECK_NEGATIVE(to_death_V1[a].Size())
            CHECK_NEGATIVE(to_death_V2[a].Size())
            CHECK_NEGATIVE(to_death_V3[a].Size())
            CHECK_NEGATIVE(to_hosp_V1[a].Size())
            CHECK_NEGATIVE(to_hosp_V2[a].Size())
            CHECK_NEGATIVE(to_hosp_V3[a].Size())

            CHECK_NEGATIVE(nE2_Ia2)
            CHECK_NEGATIVE(nE2_Ip2)
            CHECK_NEGATIVE(nE2_Ipa2)
            CHECK_NEGATIVE(nE3_Ia3)
            CHECK_NEGATIVE(nE3_Ip3)
            CHECK_NEGATIVE(nE3_Ipa3)
            CHECK_NEGATIVE(nE_Ia)
            CHECK_NEGATIVE(nE_Ip)
            CHECK_NEGATIVE(nE_Ipa)
            CHECK_NEGATIVE(nIa2_R2)
            CHECK_NEGATIVE(nIa3_R3)
            CHECK_NEGATIVE(nIa_R)
            CHECK_NEGATIVE(nIp2_Is2)
            CHECK_NEGATIVE(nIp3_Is3)
            CHECK_NEGATIVE(nIp_Is)
            CHECK_NEGATIVE(nIs2_R2)
            CHECK_NEGATIVE(nIs3_R3)
            CHECK_NEGATIVE(nIs_R)
            CHECK_NEGATIVE(nL2_Ia2)
            CHECK_NEGATIVE(nL3_Ia3)
            CHECK_NEGATIVE(nL_Ia)
            CHECK_NEGATIVE(nR2_E)
            CHECK_NEGATIVE(nR2_E2)
            CHECK_NEGATIVE(nR2_E3)
            CHECK_NEGATIVE(nR2_L)
            CHECK_NEGATIVE(nR2_L2)
            CHECK_NEGATIVE(nR2_L3)
            CHECK_NEGATIVE(nR2_R)
            CHECK_NEGATIVE(nR2_R2)
            CHECK_NEGATIVE(nR2_R3)
            CHECK_NEGATIVE(nR2_REL)
            CHECK_NEGATIVE(nR2_REL123)
            CHECK_NEGATIVE(nR2_REL2)
            CHECK_NEGATIVE(nR2_REL23)
            CHECK_NEGATIVE(nR2_REL3)
            CHECK_NEGATIVE(nR2_S)
            CHECK_NEGATIVE(nR2_X)
            CHECK_NEGATIVE(nR2_X2)
            CHECK_NEGATIVE(nR2_X3)
            CHECK_NEGATIVE(nR3_E)
            CHECK_NEGATIVE(nR3_E2)
            CHECK_NEGATIVE(nR3_E3)
            CHECK_NEGATIVE(nR3_L)
            CHECK_NEGATIVE(nR3_L2)
            CHECK_NEGATIVE(nR3_L3)
            CHECK_NEGATIVE(nR3_R)
            CHECK_NEGATIVE(nR3_R2)
            CHECK_NEGATIVE(nR3_R3)
            CHECK_NEGATIVE(nR3_REL)
            CHECK_NEGATIVE(nR3_REL123)
            CHECK_NEGATIVE(nR3_REL2)
            CHECK_NEGATIVE(nR3_REL23)
            CHECK_NEGATIVE(nR3_REL3)
            CHECK_NEGATIVE(nR3_S)
            CHECK_NEGATIVE(nR3_X)
            CHECK_NEGATIVE(nR3_X2)
            CHECK_NEGATIVE(nR3_X3)
            CHECK_NEGATIVE(nR_E)
            CHECK_NEGATIVE(nR_E2)
            CHECK_NEGATIVE(nR_E3)
            CHECK_NEGATIVE(nR_L)
            CHECK_NEGATIVE(nR_L2)
            CHECK_NEGATIVE(nR_L3)
            CHECK_NEGATIVE(nR_R)
            CHECK_NEGATIVE(nR_R2)
            CHECK_NEGATIVE(nR_R3)
            CHECK_NEGATIVE(nR_REL)
            CHECK_NEGATIVE(nR_REL123)
            CHECK_NEGATIVE(nR_REL2)
            CHECK_NEGATIVE(nR_REL23)
            CHECK_NEGATIVE(nR_REL3)
            CHECK_NEGATIVE(nR_S)
            CHECK_NEGATIVE(nR_X)
            CHECK_NEGATIVE(nR_X2)
            CHECK_NEGATIVE(nR_X3)
            CHECK_NEGATIVE(nS_E)
            CHECK_NEGATIVE(nS_E123)
            CHECK_NEGATIVE(nS_E2)
            CHECK_NEGATIVE(nS_E2E3)
            CHECK_NEGATIVE(nS_E3)
            CHECK_NEGATIVE(nS_Va1)
            CHECK_NEGATIVE(nS_Vb1)
            CHECK_NEGATIVE(nVa1_E)
            CHECK_NEGATIVE(nVa1_E2)
            CHECK_NEGATIVE(nVa1_E3)
            CHECK_NEGATIVE(nVa1_L)
            CHECK_NEGATIVE(nVa1_L2)
            CHECK_NEGATIVE(nVa1_L3)
            CHECK_NEGATIVE(nVa1_R)
            CHECK_NEGATIVE(nVa1_R2)
            CHECK_NEGATIVE(nVa1_R3)
            CHECK_NEGATIVE(nVa1_REL)
            CHECK_NEGATIVE(nVa1_REL123)
            CHECK_NEGATIVE(nVa1_REL2)
            CHECK_NEGATIVE(nVa1_REL23)
            CHECK_NEGATIVE(nVa1_REL3)
            CHECK_NEGATIVE(nVa1_S)
            CHECK_NEGATIVE(nVa1_Va2)
            CHECK_NEGATIVE(nVa1_X)
            CHECK_NEGATIVE(nVa1_X2)
            CHECK_NEGATIVE(nVa1_X3)
            CHECK_NEGATIVE(nVa2_E)
            CHECK_NEGATIVE(nVa2_E2)
            CHECK_NEGATIVE(nVa2_E3)
            CHECK_NEGATIVE(nVa2_L)
            CHECK_NEGATIVE(nVa2_L2)
            CHECK_NEGATIVE(nVa2_L3)
            CHECK_NEGATIVE(nVa2_R)
            CHECK_NEGATIVE(nVa2_R2)
            CHECK_NEGATIVE(nVa2_R3)
            CHECK_NEGATIVE(nVa2_REL)
            CHECK_NEGATIVE(nVa2_REL123)
            CHECK_NEGATIVE(nVa2_REL2)
            CHECK_NEGATIVE(nVa2_REL23)
            CHECK_NEGATIVE(nVa2_REL3)
            CHECK_NEGATIVE(nVa2_S)
            CHECK_NEGATIVE(nVa2_Va3)
            CHECK_NEGATIVE(nVa2_Vb2)
            CHECK_NEGATIVE(nVa2_X)
            CHECK_NEGATIVE(nVa2_X2)
            CHECK_NEGATIVE(nVa2_X3)
            CHECK_NEGATIVE(nVa3_E)
            CHECK_NEGATIVE(nVa3_E2)
            CHECK_NEGATIVE(nVa3_E3)
            CHECK_NEGATIVE(nVa3_L)
            CHECK_NEGATIVE(nVa3_L2)
            CHECK_NEGATIVE(nVa3_L3)
            CHECK_NEGATIVE(nVa3_R)
            CHECK_NEGATIVE(nVa3_R2)
            CHECK_NEGATIVE(nVa3_R3)
            CHECK_NEGATIVE(nVa3_REL)
            CHECK_NEGATIVE(nVa3_REL123)
            CHECK_NEGATIVE(nVa3_REL2)
            CHECK_NEGATIVE(nVa3_REL23)
            CHECK_NEGATIVE(nVa3_REL3)
            CHECK_NEGATIVE(nVa3_S)
            CHECK_NEGATIVE(nVa3_X)
            CHECK_NEGATIVE(nVa3_X2)
            CHECK_NEGATIVE(nVa3_X3)
            CHECK_NEGATIVE(nVb1_E)
            CHECK_NEGATIVE(nVb1_E2)
            CHECK_NEGATIVE(nVb1_E3)
            CHECK_NEGATIVE(nVb1_L)
            CHECK_NEGATIVE(nVb1_L2)
            CHECK_NEGATIVE(nVb1_L3)
            CHECK_NEGATIVE(nVb1_R)
            CHECK_NEGATIVE(nVb1_R2)
            CHECK_NEGATIVE(nVb1_R3)
            CHECK_NEGATIVE(nVb1_REL)
            CHECK_NEGATIVE(nVb1_REL123)
            CHECK_NEGATIVE(nVb1_REL2)
            CHECK_NEGATIVE(nVb1_REL23)
            CHECK_NEGATIVE(nVb1_REL3)
            CHECK_NEGATIVE(nVb1_S)
            CHECK_NEGATIVE(nVb1_Vb2)
            CHECK_NEGATIVE(nVb1_X)
            CHECK_NEGATIVE(nVb1_X2)
            CHECK_NEGATIVE(nVb1_X3)
            CHECK_NEGATIVE(nVb2_E)
            CHECK_NEGATIVE(nVb2_E2)
            CHECK_NEGATIVE(nVb2_E3)
            CHECK_NEGATIVE(nVb2_L)
            CHECK_NEGATIVE(nVb2_L2)
            CHECK_NEGATIVE(nVb2_L3)
            CHECK_NEGATIVE(nVb2_R)
            CHECK_NEGATIVE(nVb2_R2)
            CHECK_NEGATIVE(nVb2_R3)
            CHECK_NEGATIVE(nVb2_REL)
            CHECK_NEGATIVE(nVb2_REL123)
            CHECK_NEGATIVE(nVb2_REL2)
            CHECK_NEGATIVE(nVb2_REL23)
            CHECK_NEGATIVE(nVb2_REL3)
            CHECK_NEGATIVE(nVb2_S)
            CHECK_NEGATIVE(nVb2_Vb2)
            CHECK_NEGATIVE(nVb2_Vb3)
            CHECK_NEGATIVE(nVb2_X)
            CHECK_NEGATIVE(nVb2_X2)
            CHECK_NEGATIVE(nVb2_X3)
            CHECK_NEGATIVE(nVb3_E)
            CHECK_NEGATIVE(nVb3_E2)
            CHECK_NEGATIVE(nVb3_E3)
            CHECK_NEGATIVE(nVb3_L)
            CHECK_NEGATIVE(nVb3_L2)
            CHECK_NEGATIVE(nVb3_L3)
            CHECK_NEGATIVE(nVb3_R)
            CHECK_NEGATIVE(nVb3_R2)
            CHECK_NEGATIVE(nVb3_R3)
            CHECK_NEGATIVE(nVb3_REL)
            CHECK_NEGATIVE(nVb3_REL123)
            CHECK_NEGATIVE(nVb3_REL2)
            CHECK_NEGATIVE(nVb3_REL23)
            CHECK_NEGATIVE(nVb3_REL3)
            CHECK_NEGATIVE(nVb3_S)
            CHECK_NEGATIVE(nVb3_X)
            CHECK_NEGATIVE(nVb3_X2)
            CHECK_NEGATIVE(nVb3_X3)

            CHECK_EFFICACY(P.pop[p].pi_r[a])
            CHECK_EFFICACY(P.pop[p].pi_r2[a])
            CHECK_EFFICACY(P.pop[p].pi_r3[a])
            CHECK_EFFICACY(P.pop[p].pi2_r[a])
            CHECK_EFFICACY(P.pop[p].pi2_r2[a])
            CHECK_EFFICACY(P.pop[p].pi2_r3[a])
            CHECK_EFFICACY(P.pop[p].pi3_r[a])
            CHECK_EFFICACY(P.pop[p].pi3_r2[a])
            CHECK_EFFICACY(P.pop[p].pi3_r3[a])
            CHECK_EFFICACY(P.pop[p].ei_va1[a])
            CHECK_EFFICACY(P.pop[p].ei2_va1[a])
            CHECK_EFFICACY(P.pop[p].ei3_va1[a])
            CHECK_EFFICACY(P.pop[p].ed_va1i[a])
            CHECK_EFFICACY(P.pop[p].ed_va1i2[a])
            CHECK_EFFICACY(P.pop[p].ed_va1i3[a])
            CHECK_EFFICACY(P.pop[p].ei_va2[a])
            CHECK_EFFICACY(P.pop[p].ei2_va2[a])
            CHECK_EFFICACY(P.pop[p].ei3_va2[a])
            CHECK_EFFICACY(P.pop[p].ed_va2i[a])
            CHECK_EFFICACY(P.pop[p].ed_va2i2[a])
            CHECK_EFFICACY(P.pop[p].ed_va2i3[a])
            CHECK_EFFICACY(P.pop[p].ei_va3[a])
            CHECK_EFFICACY(P.pop[p].ei2_va3[a])
            CHECK_EFFICACY(P.pop[p].ei3_va3[a])
            CHECK_EFFICACY(P.pop[p].ed_va3i[a])
            CHECK_EFFICACY(P.pop[p].ed_va3i2[a])
            CHECK_EFFICACY(P.pop[p].ed_va3i3[a])
            CHECK_EFFICACY(P.pop[p].ei_vb1[a])
            CHECK_EFFICACY(P.pop[p].ei2_vb1[a])
            CHECK_EFFICACY(P.pop[p].ei3_vb1[a])
            CHECK_EFFICACY(P.pop[p].ed_vb1i[a])
            CHECK_EFFICACY(P.pop[p].ed_vb1i2[a])
            CHECK_EFFICACY(P.pop[p].ed_vb1i3[a])
            CHECK_EFFICACY(P.pop[p].ei_vb2[a])
            CHECK_EFFICACY(P.pop[p].ei2_vb2[a])
            CHECK_EFFICACY(P.pop[p].ei3_vb2[a])
            CHECK_EFFICACY(P.pop[p].ed_vb2i[a])
            CHECK_EFFICACY(P.pop[p].ed_vb2i2[a])
            CHECK_EFFICACY(P.pop[p].ed_vb2i3[a])
            CHECK_EFFICACY(P.pop[p].ei_vb3[a])
            CHECK_EFFICACY(P.pop[p].ei2_vb3[a])
            CHECK_EFFICACY(P.pop[p].ei3_vb3[a])
            CHECK_EFFICACY(P.pop[p].ed_vb3i[a])
            CHECK_EFFICACY(P.pop[p].ed_vb3i2[a])
            CHECK_EFFICACY(P.pop[p].ed_vb3i3[a])
            CHECK_EFFICACY(P.pop[p].pd_ri[a])
            CHECK_EFFICACY(P.pop[p].pd_ri2[a])
            CHECK_EFFICACY(P.pop[p].pd_ri3[a])
            CHECK_EFFICACY(P.pop[p].pd_r2i[a])
            CHECK_EFFICACY(P.pop[p].pd_r2i2[a])
            CHECK_EFFICACY(P.pop[p].pd_r2i3[a])
            CHECK_EFFICACY(P.pop[p].pd_r3i[a])
            CHECK_EFFICACY(P.pop[p].pd_r3i2[a])
            CHECK_EFFICACY(P.pop[p].pd_r3i3[a])
            CHECK_EFFICACY(P.pop[p].ph_rd[a])
            CHECK_EFFICACY(P.pop[p].ph_rd2[a])
            CHECK_EFFICACY(P.pop[p].ph_rd3[a])
            CHECK_EFFICACY(P.pop[p].ph_r2d[a])
            CHECK_EFFICACY(P.pop[p].ph_r2d2[a])
            CHECK_EFFICACY(P.pop[p].ph_r2d3[a])
            CHECK_EFFICACY(P.pop[p].ph_r3d[a])
            CHECK_EFFICACY(P.pop[p].ph_r3d2[a])
            CHECK_EFFICACY(P.pop[p].ph_r3d3[a])
            CHECK_EFFICACY(P.pop[p].pm_rd[a])
            CHECK_EFFICACY(P.pop[p].pm_rd2[a])
            CHECK_EFFICACY(P.pop[p].pm_rd3[a])
            CHECK_EFFICACY(P.pop[p].pm_r2d[a])
            CHECK_EFFICACY(P.pop[p].pm_r2d2[a])
            CHECK_EFFICACY(P.pop[p].pm_r2d3[a])
            CHECK_EFFICACY(P.pop[p].pm_r3d[a])
            CHECK_EFFICACY(P.pop[p].pm_r3d2[a])
            CHECK_EFFICACY(P.pop[p].pm_r3d3[a])
            CHECK_EFFICACY(P.pop[p].pt_ri[a])
            CHECK_EFFICACY(P.pop[p].pt_ri2[a])
            CHECK_EFFICACY(P.pop[p].pt_ri3[a])
            CHECK_EFFICACY(P.pop[p].pt_r2i[a])
            CHECK_EFFICACY(P.pop[p].pt_r2i2[a])
            CHECK_EFFICACY(P.pop[p].pt_r2i3[a])
            CHECK_EFFICACY(P.pop[p].pt_r3i[a])
            CHECK_EFFICACY(P.pop[p].pt_r3i2[a])
            CHECK_EFFICACY(P.pop[p].pt_r3i3[a])
            CHECK_EFFICACY(P.pop[p].eh_va1d[a])
            CHECK_EFFICACY(P.pop[p].eh_va1d2[a])
            CHECK_EFFICACY(P.pop[p].eh_va1d3[a])
            CHECK_EFFICACY(P.pop[p].eh_va2d[a])
            CHECK_EFFICACY(P.pop[p].eh_va2d2[a])
            CHECK_EFFICACY(P.pop[p].eh_va2d3[a])
            CHECK_EFFICACY(P.pop[p].eh_va3d[a])
            CHECK_EFFICACY(P.pop[p].eh_va3d2[a])
            CHECK_EFFICACY(P.pop[p].eh_va3d3[a])
            CHECK_EFFICACY(P.pop[p].eh_vb1d[a])
            CHECK_EFFICACY(P.pop[p].eh_vb1d2[a])
            CHECK_EFFICACY(P.pop[p].eh_vb1d3[a])
            CHECK_EFFICACY(P.pop[p].eh_vb2d[a])
            CHECK_EFFICACY(P.pop[p].eh_vb2d2[a])
            CHECK_EFFICACY(P.pop[p].eh_vb2d3[a])
            CHECK_EFFICACY(P.pop[p].eh_vb3d[a])
            CHECK_EFFICACY(P.pop[p].eh_vb3d2[a])
            CHECK_EFFICACY(P.pop[p].eh_vb3d3[a])
            CHECK_EFFICACY(P.pop[p].em_va1d[a])
            CHECK_EFFICACY(P.pop[p].em_va1d2[a])
            CHECK_EFFICACY(P.pop[p].em_va1d3[a])
            CHECK_EFFICACY(P.pop[p].em_va2d[a])
            CHECK_EFFICACY(P.pop[p].em_va2d2[a])
            CHECK_EFFICACY(P.pop[p].em_va2d3[a])
            CHECK_EFFICACY(P.pop[p].em_va3d[a])
            CHECK_EFFICACY(P.pop[p].em_va3d2[a])
            CHECK_EFFICACY(P.pop[p].em_va3d3[a])
            CHECK_EFFICACY(P.pop[p].em_vb1d[a])
            CHECK_EFFICACY(P.pop[p].em_vb1d2[a])
            CHECK_EFFICACY(P.pop[p].em_vb1d3[a])
            CHECK_EFFICACY(P.pop[p].em_vb2d[a])
            CHECK_EFFICACY(P.pop[p].em_vb2d2[a])
            CHECK_EFFICACY(P.pop[p].em_vb2d3[a])
            CHECK_EFFICACY(P.pop[p].em_vb3d[a])
            CHECK_EFFICACY(P.pop[p].em_vb3d2[a])
            CHECK_EFFICACY(P.pop[p].em_vb3d3[a])
            CHECK_EFFICACY(P.pop[p].et_va1i[a])
            CHECK_EFFICACY(P.pop[p].et_va1i2[a])
            CHECK_EFFICACY(P.pop[p].et_va1i3[a])
            CHECK_EFFICACY(P.pop[p].et_va2i[a])
            CHECK_EFFICACY(P.pop[p].et_va2i2[a])
            CHECK_EFFICACY(P.pop[p].et_va2i3[a])
            CHECK_EFFICACY(P.pop[p].et_va3i[a])
            CHECK_EFFICACY(P.pop[p].et_va3i2[a])
            CHECK_EFFICACY(P.pop[p].et_va3i3[a])
            CHECK_EFFICACY(P.pop[p].et_vb1i[a])
            CHECK_EFFICACY(P.pop[p].et_vb1i2[a])
            CHECK_EFFICACY(P.pop[p].et_vb1i3[a])
            CHECK_EFFICACY(P.pop[p].et_vb2i[a])
            CHECK_EFFICACY(P.pop[p].et_vb2i2[a])
            CHECK_EFFICACY(P.pop[p].et_vb2i3[a])
            CHECK_EFFICACY(P.pop[p].et_vb3i[a])
            CHECK_EFFICACY(P.pop[p].et_vb3i2[a])
            CHECK_EFFICACY(P.pop[p].et_vb3i3[a])

            if (stop_due_to_negative) {
                Rcpp::stop("Negative");
            }
        }

        // 2. User-specified processes
        fill(pco.begin(), pco.end(), -99999.);

        for (auto& process : P.processes)
        {
            // Determine number of individuals entering the process
            double n_entering = 0.;
            switch (process.source_id)
            {
                case src_newE:
                    n_entering = nS_E + nR_E + nR2_E + nR3_E + nVa1_E + nVa2_E + nVa3_E + nVb1_E + nVb2_E + nVb3_E; break;
                case src_newI:
                    n_entering = nE_Ipa + nL_Ia; break;
                case src_newIp:
                    n_entering = nE_Ip; break;
                case src_newIs:
                    n_entering = nIp_Is; break;
                case src_newIa:
                    n_entering = nE_Ia + nL_Ia; break;
                case src_newE2:
                    n_entering = nS_E2 + nR_E2 + nR2_E2 + nR3_E2 + nVa1_E2 + nVa2_E2 + nVa3_E2 + nVb1_E2 + nVb2_E2 + nVb3_E2; break;
                case src_newI2:
                    n_entering = nE2_Ipa2 + nL2_Ia2; break;
                case src_newIp2:
                    n_entering = nE2_Ip2; break;
                case src_newIs2:
                    n_entering = nIp2_Is2; break;
                case src_newIa2:
                    n_entering = nE2_Ia2 + nL2_Ia2; break;
                case src_newE3:
                    n_entering = nS_E3 + nR_E3 + nR2_E3 + nR3_E3 + nVa1_E3 + nVa2_E3 + nVa3_E3 + nVb1_E3 + nVb2_E3 + nVb3_E3; break;
                case src_newI3:
                    n_entering = nE3_Ipa3 + nL3_Ia3; break;
                case src_newIp3:
                    n_entering = nE3_Ip3; break;
                case src_newIs3:
                    n_entering = nIp3_Is3; break;
                case src_newIa3:
                    n_entering = nE3_Ia3 + nL3_Ia3; break;
                case src_newEE2E3:
                    n_entering = nS_E + nR_E + nR2_E + nR3_E + nVa1_E + nVa2_E + nVa3_E + nVb1_E + nVb2_E + nVb3_E + nS_E2 + nR_E2 + nR2_E2 + nR3_E2 + nVa1_E2 + nVa2_E2 + nVa3_E2 + nVb1_E2 + nVb2_E2 + nVb3_E2 + nS_E3 + nR_E3 + nR2_E3 + nR3_E3 + nVa1_E3 + nVa2_E3 + nVa3_E3 + nVb1_E3 + nVb2_E3 + nVb3_E3; break;
                case src_newII2I3:
                    n_entering = nE_Ipa + nE2_Ipa2 + nE3_Ipa3 + nL_Ia + nL2_Ia2 + nL3_Ia3; break;
                case src_newIpIp2Ip3:
                    n_entering = nE_Ip + nE2_Ip2 + nE3_Ip3; break;
                case src_newIsIs2Is3:
                    n_entering = nIp_Is + nIp2_Is2 + nIp3_Is3; break;
                case src_newIaIa2Ia3:
                    n_entering = nE_Ia + nE2_Ia2 + nE3_Ia3 + nL_Ia + nL2_Ia2 + nL3_Ia3; break;
                case src_newVa1:
                    n_entering = nS_Va1; break;
                case src_newVb1:
                    n_entering = nS_Vb1; break;
                case src_HospAdm:
                    n_entering = hosp_admissions; break;
                case src_newEL123:
                    n_entering = nS_E + nR_E + nR2_E + nR3_E + nVa1_E + nVa2_E + nVa3_E + nVb1_E + nVb2_E + nVb3_E + nS_E2 + nR_E2 + nR2_E2 + nR3_E2 + nVa1_E2 + nVa2_E2 + nVa3_E2 + nVb1_E2 + nVb2_E2 + nVb3_E2 + nS_E3 + nR_E3 + nR2_E3 + nR3_E3 + nVa1_E3 + nVa2_E3 + nVa3_E3 + nVb1_E3 + nVb2_E3 + nVb3_E3 +
                        nR_L + nR2_L + nR3_L + nVa1_L + nVa2_L + nVa3_L + nVb1_L + nVb2_L + nVb3_L + nR_L2 + nR2_L2 + nR3_L2 + nVa1_L2 + nVa2_L2 + nVa3_L2 + nVb1_L2 + nVb2_L2 + nVb3_L2 + nR_L3 + nR2_L3 + nR3_L3 + nVa1_L3 + nVa2_L3 + nVa3_L3 + nVb1_L3 + nVb2_L3 + nVb3_L3; break;
                default:
                    n_entering = pco[process.source_id];
                    if (n_entering == -99999.)
                        throw logic_error("Process sourced from unset user process. Have user processes been specified in the right order?");
                    if (n_entering < 0)
                        throw logic_error("n_entering is negative (" + to_string(n_entering) + ") for process sourced from " + 
                            process.source_name + " and leading to " + process.names[0] + ". Time " + to_string(t) + ", group " + to_string(a) + ".");
                    break;
            }

            multinomial(n_entering, process.prob[a], nd_out, ni_out);

            // Seed and mature this process's compartments
            unsigned int c = 0;
            for (unsigned int compartment_id : process.ids)
            {
                if (compartment_id != Null)
                {
                    pc[compartment_id][a].Add(P, Rand, nd_out[c], process.delays[c]);
                    pci[compartment_id] = nd_out[c];
                    pco[compartment_id] = pc[compartment_id][a].Mature(P.time_step);
                }
                ++c;
            }
        }

        // 3. Report incidence / outcidence

        // Built-in compartments
        rep(t, p, a, 25) += died;
        rep(t, p, a, 26) += hosp_admissions;
        rep(t, p, a, 28) += icu_admissions;
        rep(t, p, a, 30) += died_V1;
        rep(t, p, a, 31) += died_V2;
        rep(t, p, a, 32) += died_V3;
        rep(t, p, a, 33) += hosp_admissions_V1;
        rep(t, p, a, 34) += hosp_admissions_V2;
        rep(t, p, a, 35) += hosp_admissions_V3;
        rep(t, p, a, 37) += nS_E  + nVa1_E  + nVa2_E  + nVa3_E  + nVb1_E  + nVb2_E  + nVb3_E  + nR_E  + nR2_E  + nR3_E  + 
                            nS_E2 + nVa1_E2 + nVa2_E2 + nVa3_E2 + nVb1_E2 + nVb2_E2 + nVb3_E2 + nR_E2 + nR2_E2 + nR3_E2 + 
                            nS_E3 + nVa1_E3 + nVa2_E3 + nVa3_E3 + nVb1_E3 + nVb2_E3 + nVb3_E3 + nR_E3 + nR2_E3 + nR3_E3;
        rep(t, p, a, 38) += severe;
        


        // User-specified processes
        for (auto& process : P.processes)
        {
            for (unsigned int i = 0; i < process.i_cols.size(); ++i)
                rep(t, p, a, process.i_cols[i]) += pci[process.i_ids[i]];
            for (unsigned int i = 0; i < process.o_cols.size(); ++i)
                rep(t, p, a, process.o_cols[i]) += pco[process.o_ids[i]];
        }
    }

    // Births, deaths, aging
    double Ntot = accumulate(N.begin(), N.end(), 0.0);
    for (unsigned int a = N.size() - 1; ; --a)
    {
        // Births
        double B = poisson(Ntot * (exp(P.pop[p].B[a] * P.time_step) - 1.));

        // Deaths (outer compartments)
        double death_prob = 1.0 - exp(-P.pop[p].D[a] * P.time_step);
        double DS         = binomial(S[a],  death_prob);
        double DR         = binomial(R[a],  death_prob);
        double DR2        = binomial(R2[a], death_prob);
        double DR3        = binomial(R3[a], death_prob);
        double DVa3       = binomial(Va3[a],  death_prob);
        double DVb3       = binomial(Vb3[a],  death_prob);

        // Changes
        N[a]   += B;
        S[a]   += B;
        S[a]   -= DS;
        R[a]   -= DR;
        R2[a]  -= DR2;
        R3[a]  -= DR3;
        Va3[a] -= DVa3;
        Vb3[a] -= DVb3;

        // Deaths (inner compartments)
        double DE   = E[a]  .RemoveProb(P, Rand, death_prob);
        double DL   = L[a]  .RemoveProb(P, Rand, death_prob);
        double DIp  = Ip[a] .RemoveProb(P, Rand, death_prob);
        double DIa  = Ia[a] .RemoveProb(P, Rand, death_prob);
        double DIs  = Is[a] .RemoveProb(P, Rand, death_prob);
        double DE2  = E2[a] .RemoveProb(P, Rand, death_prob);
        double DL2  = L2[a] .RemoveProb(P, Rand, death_prob);
        double DIp2 = Ip2[a].RemoveProb(P, Rand, death_prob);
        double DIa2 = Ia2[a].RemoveProb(P, Rand, death_prob);
        double DIs2 = Is2[a].RemoveProb(P, Rand, death_prob);
        double DE3  = E3[a] .RemoveProb(P, Rand, death_prob);
        double DL3  = L3[a] .RemoveProb(P, Rand, death_prob);
        double DIp3 = Ip3[a].RemoveProb(P, Rand, death_prob);
        double DIa3 = Ia3[a].RemoveProb(P, Rand, death_prob);
        double DIs3 = Is3[a].RemoveProb(P, Rand, death_prob);
        double DVa1 = Va1[a].RemoveProb(P, Rand, death_prob);
        double DVb1 = Vb1[a].RemoveProb(P, Rand, death_prob);
        double DVa2 = Va2[a].RemoveProb(P, Rand, death_prob);
        double DVb2 = Vb2[a].RemoveProb(P, Rand, death_prob);


        N[a]  -= DS + DR + DR2 + DR3 + DVa1 + DVa2 + DVa3 + DVb1 + DVb2 + DVb3 + DE + DL + DIp + DIa + DIs + DE2 + DL2 + DIp2 + DIa2 + DIs2 + DE3 + DL3 + DIp3 + DIa3 + DIs3;

        // Agings
        if (a != lambda.size() - 1)
        {
            double age_prob = 1.0 - exp(-P.pop[p].A[a] * P.time_step);
            double AS   = binomial(S[a],  age_prob);
            double AR   = binomial(R[a],  age_prob);
            double AR2  = binomial(R2[a], age_prob);
            double AR3  = binomial(R3[a], age_prob);
            double AVa3 = binomial(Va3[a],  age_prob);
            double AVb3 = binomial(Vb3[a],  age_prob);

            S[a]       -= AS;
            S[a + 1]   += AS;
            R[a]       -= AR;
            R[a + 1]   += AR;
            R2[a]      -= AR2;
            R2[a + 1]  += AR2;
            R3[a]      -= AR3;
            R3[a + 1]  += AR3;
            Va3[a]     -= AVa3;
            Va3[a + 1] += AVa3;
            Vb3[a]     -= AVb3;
            Vb3[a + 1] += AVb3;

            double AE   = E[a]  .MoveProb(E   [a + 1], P, Rand, age_prob);
            double AL   = L[a]  .MoveProb(L   [a + 1], P, Rand, age_prob);
            double AIp  = Ip[a] .MoveProb(Ip  [a + 1], P, Rand, age_prob);
            double AIa  = Ia[a] .MoveProb(Ia  [a + 1], P, Rand, age_prob);
            double AIs  = Is[a] .MoveProb(Is  [a + 1], P, Rand, age_prob);
            double AE2  = E2[a] .MoveProb(E2  [a + 1], P, Rand, age_prob);
            double AL2  = L2[a] .MoveProb(L2  [a + 1], P, Rand, age_prob);
            double AIp2 = Ip2[a].MoveProb(Ip2 [a + 1], P, Rand, age_prob);
            double AIa2 = Ia2[a].MoveProb(Ia2 [a + 1], P, Rand, age_prob);
            double AIs2 = Is2[a].MoveProb(Is2 [a + 1], P, Rand, age_prob);
            double AE3  = E3[a] .MoveProb(E3  [a + 1], P, Rand, age_prob);
            double AL3  = L3[a] .MoveProb(L3  [a + 1], P, Rand, age_prob);
            double AIp3 = Ip3[a].MoveProb(Ip3 [a + 1], P, Rand, age_prob);
            double AIa3 = Ia3[a].MoveProb(Ia3 [a + 1], P, Rand, age_prob);
            double AIs3 = Is3[a].MoveProb(Is3 [a + 1], P, Rand, age_prob);
            double AVa1 = Va1[a].MoveProb(Va1 [a + 1], P, Rand, age_prob);
            double AVb1 = Vb1[a].MoveProb(Vb1 [a + 1], P, Rand, age_prob);
            double AVa2 = Va2[a].MoveProb(Va2 [a + 1], P, Rand, age_prob);
            double AVb2 = Vb2[a].MoveProb(Vb2 [a + 1], P, Rand, age_prob);


            N[a]      -= AS + AR + AR2 + AR3 + AVa1 + AVa2 + AVa3 + AVb1 + AVb2 + AVb3 + AE + AL + AIp + AIa + AIs + AE2 + AL2 + AIp2 + AIa2 + AIs2 + AE3 + AL3 + AIp3 + AIa3 + AIs3;
            N[a + 1]  += AS + AR + AR2 + AR3 + AVa1 + AVa2 + AVa3 + AVb1 + AVb2 + AVb3 + AE + AL + AIp + AIa + AIs + AE2 + AL2 + AIp2 + AIa2 + AIs2 + AE3 + AL3 + AIp3 + AIa3 + AIs3;
        }

        if (a == 0)
            break;
    }
}

// Print full population details
void Population::DebugPrint() const
{
    auto vecprint = [&](const vector<double>& vec, string name) {
        cout << name;
        for (auto& v : vec)
            cout << " " << v;
        cout << "\n";
    };

    auto comprint = [&](const vector<Compartment>& comp, string name) {
        cout << name;
        for (unsigned int c = 0; c < comp.size(); ++c) {
            cout << "element " << c << "\n";
            comp[c].DebugPrint();
        }
    };

    vecprint(lambda, "lambda");
    vecprint(lambda2, "lambda2");
    vecprint(lambda3, "lambda3");
    vecprint(N, "N");
    vecprint(S, "S");
    vecprint(R, "R");
    comprint(Va1, "Va1");
    comprint(Va2, "Va2");
    vecprint(Va3, "Va3");
    comprint(Vb1, "Vb1");
    comprint(Vb2, "Vb2");
    vecprint(Vb3, "Vb3");
    comprint(E, "E");
    comprint(L, "L");
    comprint(Ip, "Ip");
    comprint(Ia, "Ia");
    comprint(Is, "Is");
    vecprint(R2, "R2");
    comprint(E2, "E2");
    comprint(L2, "L2");
    comprint(Ip2, "Ip2");
    comprint(Ia2, "Ia2");
    comprint(Is2, "Is2");
    vecprint(R3, "R3");
    comprint(E3, "E3");
    comprint(L3, "L3");
    comprint(Ip3, "Ip3");
    comprint(Ia2, "Ia3");
    comprint(Is3, "Is3");
    cout << "seed_row " << seed_row << " p " << p << "\n";
    cout << "seed_row2 " << seed_row2 << " p " << p << "\n";
    for (auto& c : pc)
        comprint(c, "User");

    cout << "\n\n";
}


Metapopulation::Metapopulation(Parameters& P)
{
    P.changes.Capture(P);

    for (unsigned int i = 0; i < P.pop.size(); ++i)
        pops.push_back(Population(P, i));
}

// Execute one time step's events
bool Metapopulation::Tick(Parameters& P, Randomizer& Rand, double t, Reporter& rep)
{
    // Apply any changes to parameters
    P.changes.Apply(P, t);

    unsigned int n_ages = P.pop[0].size.size();

    // Calculate contagiousness from each population
    // NOTE -- 'contag' subscripted first by j, then by a.
    // It's the effective number of infectious individuals FROM subpop j of age a.
    contag.assign(pops.size(), vector<double>(n_ages, 0.0));
    contag2.assign(pops.size(), vector<double>(n_ages, 0.0));
    contag3.assign(pops.size(), vector<double>(n_ages, 0.0));
    for (unsigned int j = 0; j < pops.size(); ++j)
        pops[j].Contagiousness(P, Rand, t, contag[j], contag2[j], contag3[j]);

    // note -- 'infec' subscripted first by i, then by a
    // It's the effective number of infectious individuals who are CURRENTLY IN subpop i of age a.
    infec.assign(pops.size(), vector<double>(n_ages, 0.0));
    infec2.assign(pops.size(), vector<double>(n_ages, 0.0));
    infec3.assign(pops.size(), vector<double>(n_ages, 0.0));
    for (unsigned int i = 0; i < pops.size(); ++i)
        for (unsigned int j = 0; j < pops.size(); ++j)
            for (unsigned int a = 0; a < n_ages; ++a)
            {
                infec[i][a]  += P.travel(j, i) * contag[j][a]  * (j != i ? P.pop[j].tau[a] : 1.0);
                infec2[i][a] += P.travel(j, i) * contag2[j][a] * (j != i ? P.pop[j].tau[a] : 1.0);
                infec3[i][a] += P.travel(j, i) * contag3[j][a] * (j != i ? P.pop[j].tau[a] : 1.0);
            }

    // Update populations
    //#pragma omp parallel for schedule(dynamic) reduction(&&:keep_going)
    for (unsigned int i = 0; i < pops.size(); ++i)
        pops[i].Tick(P, Rand, t, infec[i], infec2[i], infec3[i], rep);

    // Run observer at the last time step of each day.
    if (t + P.time_step == int(t + P.time_step))
        return CppObserver(P, Rand, rep, (int)t, x, *this);

    return true;
}

void Metapopulation::Run(Parameters& P, Randomizer& Rand, Reporter& rep, vector<double> x_fit)
{
    x = x_fit;

    // Run simulation
    unsigned int time_steps = (1 + P.time1 - P.time0) / P.time_step;
    for (unsigned int ts = 0; ts < time_steps; ++ts)
    {
        if (!Tick(P, Rand, P.time0 + ts * P.time_step, rep))
            break;
    }
}

void Metapopulation::RunUntil(Parameters& P, Randomizer& Rand, Reporter& rep, vector<double>& x_fit, double t_from, double t_to,
    bool include_to)
{
    x = x_fit;

    // Run simulation
    unsigned int time_steps = (t_to - t_from + (include_to ? 1 : 0)) / P.time_step;
    for (unsigned int ts = 0; ts < time_steps; ++ts)
    {
        if (!Tick(P, Rand, t_from + ts * P.time_step, rep))
            break;
    }
}
