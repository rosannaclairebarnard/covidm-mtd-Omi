// Rcpp_interface.cpp

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppGSL)]]

#include <vector>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <algorithm>
#include <numeric>
#include <ctime>
#include <limits>
#include <omp.h>
#include <Rcpp.h>
using namespace std;

#include "randomizer.h"
#include "distribution.h"
#include "helper.h"
#include "process_spec.h"
#include "changes.h"
#include "parameters.h"
#include "compartment.h"
#include "reporter.h"
#include "sim_compartment.h"
#include "user_defined.h"
#include "Rcpp_interface.h"
#include "mcmc.h"


Reporter RunSimulation(Parameters& P, Randomizer& Rand, vector<double> x)
{
    Reporter rep(P);

    if (P.model == "SEI3R")
    {
        Metapopulation mp(P);
        mp.Run(P, Rand, rep, x);
    } 
    else
    {
        throw std::logic_error("Unrecognized model type.");
    }

    return rep;
}


// [[Rcpp::export]]
Rcpp::List cm_backend_simulate(Rcpp::List parameters, unsigned int n_run = 1, unsigned long int seed = 0, unsigned int n_threads = 1, string file_out = "")
{
    // Enable multithreading
    #ifdef _OPENMP
    if (n_threads > 1)
        omp_set_num_threads(n_threads);
    #endif

    // Initialise parameters for this simulation
    Randomizer rand_master(seed);
    Parameters covidm_parameters;
    SetParameters(covidm_parameters, parameters, rand_master);

    if (!file_out.empty())
    {
        // Components unique to each run
        vector<Randomizer> rand_r;
        for (unsigned int r = 0; r < n_run; ++r)
            rand_r.emplace_back(gsl_rng_get(rand_master.GSL_RNG()));

        // Run the simulation
        #pragma omp parallel for if(n_threads > 1) schedule(dynamic)
        for (unsigned int r = 0; r < n_run; ++r)
        {
            Parameters P = covidm_parameters;
            P.FilterForRun(r);
            Reporter rep = RunSimulation(P, rand_r[r]);
            rep.Save(file_out + to_string(r), rand_r[r].Seed());
        }

        return Rcpp::List::create(
            Rcpp::Named("file_out") = file_out
        );
    }
    else
    {
        // Components unique to each run
        vector<Randomizer> rand_r;
        vector<Reporter> rep_r(n_run, Reporter(covidm_parameters));
        for (unsigned int r = 0; r < n_run; ++r)
            rand_r.emplace_back(gsl_rng_get(rand_master.GSL_RNG()));

        // Run the simulation
        #pragma omp parallel for if(n_threads > 1) schedule(dynamic)
        for (unsigned int r = 0; r < n_run; ++r)
        {
            Parameters P = covidm_parameters;
            P.FilterForRun(r);
            rep_r[r] = RunSimulation(P, rand_r[r]);
        }

        // Assemble results
        Rcpp::List dynamics;
        Rcpp::List csvs;

        for (unsigned int r = 0; r < n_run; ++r)
        {
            Reporter& rep = rep_r.front();

            // Create times
            Rcpp::NumericVector t(rep.n_times * rep.n_populations * rep.n_age_groups, 0.);
            for (unsigned int it = 0; it < rep.n_times; ++it)
                for (unsigned int j = 0; j < rep.n_populations * rep.n_age_groups; ++j)
                    t[it * rep.n_populations * rep.n_age_groups + j] = covidm_parameters.time0 + it * covidm_parameters.time_step * covidm_parameters.report_every;

            // Create identifier columns
            Rcpp::DataFrame dynamics_df = Rcpp::DataFrame::create(
                Rcpp::Named("t") = t,
                Rcpp::Named("population") = Rcpp::rep(Rcpp::rep_each(Rcpp::seq(1, rep.n_populations), rep.n_age_groups), rep.n_times),
                Rcpp::Named("group") = Rcpp::rep(Rcpp::seq(1, rep.n_age_groups), rep.n_times * rep.n_populations)
            );

            // Allocate all columns to the dataframe
            for (unsigned int c = 0; c < rep.col_names.size(); ++c)
                dynamics_df.push_back(rep.data[c], rep.col_names[c]);

            // Add observer columns
            for (unsigned int c = 0; c < rep.obs.size(); ++c)
                dynamics_df.push_back(rep.obs[c], "obs" + to_string(c));

            // Set dataframe as a data.table
            Rcpp::Function setDT("setDT"); 
            setDT(dynamics_df);

            dynamics.push_back(dynamics_df);
            csvs.push_back(rep.csv);

            rep_r.erase(rep_r.begin());
        }

        return Rcpp::List::create(
            Rcpp::Named("dynamics") = dynamics,
            Rcpp::Named("csv") = csvs
        );
    }
}

// [[Rcpp::export]]
Rcpp::DataFrame cm_evaluate_distribution(string dist_code, unsigned int steps = 101, double xmin = 0, double xmax = -1)
{
    Distribution dist(dist_code);
    if (xmax < xmin)
    {
        xmin = dist.LowerBound();
        xmax = dist.UpperBound();
    }

    vector<double> x(steps, 0.);
    vector<double> p(steps, 0.);

    for (unsigned int s = 0; s < steps; ++s)
    {
        x[s] = xmin + ((xmax - xmin) / (steps - 1.)) * s;
        p[s] = exp(dist.LogProbability(x[s]));
    }

    Rcpp::DataFrame results = Rcpp::DataFrame::create(
        Rcpp::Named("x") = x,
        Rcpp::Named("p") = p
    );

    return results;
}

// [[Rcpp::export]]
Rcpp::DataFrame cm_backend_mcmc_test(Rcpp::List R_base_parameters, Rcpp::List params_priors, unsigned long int seed, 
    unsigned int burn_in, unsigned int iterations, unsigned int n_threads, bool classic_gamma)
{
    // Initialise parameters for this simulation
    // TODO Rand also used for setting parameters -- is it actually used? this may cause issues with seeds for sample fit etc
    Randomizer Rand(seed); // randomizer for fitting; randomizers for model runs are created in the Likelihood class.
    Parameters base_parameters;

    SetParameters(base_parameters, R_base_parameters, Rand);

    vector<string> param_names = Rcpp::as<vector<string>>(params_priors.names());
    vector<Distribution> priors;

    for (unsigned int i = 0; i < params_priors.size(); ++i)
        priors.push_back(Distribution(Rcpp::as<string>(params_priors[i])));

    ///---
    unsigned int n_chains = 2 * params_priors.size();
    bool verbose = true;
    bool reeval_likelihood = false;
    bool in_parallel = n_threads > 1;
    ///---

    Likelihood lik(base_parameters, seed);
    MCMCReporter rep(iterations, n_chains, param_names);

    DEMCMC_Priors(Rand, lik, rep, burn_in, iterations, n_chains, priors, verbose, param_names, 
        reeval_likelihood, in_parallel, n_threads, classic_gamma);

    // Get data.frame as a data.table and return
    Rcpp::DataFrame df = Rcpp::DataFrame::create();
    df.push_back(Rcpp::IntegerVector::import(rep.trial.begin(), rep.trial.end()), "trial");
    df.push_back(Rcpp::NumericVector::import(rep.lp.begin(), rep.lp.end()), "lp");
    df.push_back(Rcpp::IntegerVector::import(rep.chain.begin(), rep.chain.end()), "chain");
    df.push_back(Rcpp::NumericVector::import(rep.ll.begin(), rep.ll.end()), "ll");
    for (unsigned int d = 0; d < rep.theta.size(); ++d)
        df.push_back(Rcpp::NumericVector::import(rep.theta[d].begin(), rep.theta[d].end()), rep.pnames[d]);

    return df;
}

// [[Rcpp::export]]
double cm_backend_loglik(Rcpp::List R_base_parameters, Rcpp::NumericVector theta, unsigned int seed)
{
    // Initialise parameters for this simulation
    // TODO Rand also used for setting parameters -- is it actually used? this may cause issues with seeds for sample fit etc
    Randomizer Rand(seed); // randomizer for fitting; randomizers for model runs are created in the Likelihood class.
    Parameters base_parameters;

    SetParameters(base_parameters, R_base_parameters, Rand);

    Likelihood lik(base_parameters, seed);
    return lik(Rcpp::as<vector<double>>(theta));
}

void get_multi_params(Rcpp::List params_priors, unsigned int n_simulations,
    vector<vector<unsigned int>>& x_source, vector<double>& fixed_params, vector<string>& theta_names, vector<Distribution>& theta_priors)
{
    x_source = vector<vector<unsigned int>>(params_priors.size(), vector<unsigned int>(n_simulations, 0));

    for (unsigned int i = 0; i < params_priors.size(); ++i)
    {
        unsigned int length = Rcpp::as<Rcpp::List>(params_priors[i]).size();

        // Use fixed params
        if (Rf_isNumeric(params_priors[i]))
        {
            for (unsigned int j = 0; j < n_simulations; ++j)
            {
                x_source[i][j] = 10000 + fixed_params.size();
                fixed_params.push_back(Rcpp::as<vector<double>>(params_priors[i])[j % length]);
            }
        }
        // Use fitted params
        else
        {
            // Single fitted param for all simulations
            if (length == 1)
            {
                for (unsigned int j = 0; j < n_simulations; ++j)
                {
                    x_source[i][j] = theta_names.size();
                }
                theta_names.push_back(Rcpp::as<vector<string>>(params_priors.names())[i]);
                theta_priors.push_back(Distribution(Rcpp::as<string>(params_priors[i])));
            }
            // Independent fitted param for each simulation
            else if (length == n_simulations)
            {
                for (unsigned int j = 0; j < n_simulations; ++j)
                {
                    x_source[i][j] = theta_names.size();
                    theta_names.push_back(Rcpp::as<vector<string>>(params_priors.names())[i] + to_string(j));
                    theta_priors.push_back(Distribution(Rcpp::as<vector<string>>(params_priors[i])[j]));
                }
            }
            else
            {
                throw("Fitted params -- need either 1 or same as number of simulations");
            }
        }
    }
}

// [[Rcpp::export]]
Rcpp::DataFrame cm_backend_mcmc_multi(Rcpp::List R_base_parameters_list, Rcpp::List params_priors, unsigned long int seed, 
    unsigned int burn_in, unsigned int iterations, unsigned int n_threads, bool classic_gamma)
{
    // Initialise parameters for this simulation
    // TODO Rand also used for setting parameters -- is it actually used? this may cause issues with seeds for sample fit etc
    Randomizer Rand(seed); // randomizer for fitting; randomizers for model runs are created in the Likelihood class.

    // Set each of base_parameters individually
    vector<Parameters> base_parameters;
    for (unsigned int i = 0; i < R_base_parameters_list.size(); ++i)
    {
        base_parameters.push_back(Parameters());
        SetParameters(base_parameters.back(), Rcpp::as<Rcpp::List>(R_base_parameters_list[i]), Rand);
    }

    // xs = x_source[i][j] says what to use for x[i] in simulation j.
    // if xs >= 10000, use fixed_params[xs - 10000]; if xs < 10000, use theta[xs].
    vector<vector<unsigned int>> x_source;
    vector<double> fixed_params;
    vector<string> theta_names;
    vector<Distribution> theta_priors;

    get_multi_params(params_priors, base_parameters.size(), x_source, fixed_params, theta_names, theta_priors);

    ///---
    unsigned int n_chains = 2 * params_priors.size();
    bool verbose = true;
    bool reeval_likelihood = false;
    bool in_parallel = n_threads > 1;
    ///---

    Likelihood lik(base_parameters, x_source, fixed_params, seed);
    MCMCReporter rep(iterations, n_chains, theta_names);

    DEMCMC_Priors(Rand, lik, rep, burn_in, iterations, n_chains, theta_priors, verbose, theta_names, 
        reeval_likelihood, in_parallel, n_threads, classic_gamma);

    // Get data.frame and return
    Rcpp::DataFrame df = Rcpp::DataFrame::create();
    df.push_back(Rcpp::IntegerVector::import(rep.trial.begin(), rep.trial.end()), "trial");
    df.push_back(Rcpp::NumericVector::import(rep.lp.begin(), rep.lp.end()), "lp");
    df.push_back(Rcpp::IntegerVector::import(rep.chain.begin(), rep.chain.end()), "chain");
    df.push_back(Rcpp::NumericVector::import(rep.ll.begin(), rep.ll.end()), "ll");
    for (unsigned int d = 0; d < rep.theta.size(); ++d)
        df.push_back(Rcpp::NumericVector::import(rep.theta[d].begin(), rep.theta[d].end()), rep.pnames[d]);

    return df;
}

// [[Rcpp::export]]
Rcpp::DataFrame cm_backend_optimize_test(Rcpp::List R_base_parameters, Rcpp::List params_priors, 
    unsigned int maxeval, double ftol_abs, 
    unsigned long int seed, unsigned int n_threads)
{
    // Initialise parameters for this simulation
    // TODO Rand also used for setting parameters -- is it actually used? this may cause issues with seeds for sample fit etc
    Randomizer Rand(seed); // randomizer for fitting; randomizers for model runs are created in the Likelihood class.
    Parameters base_parameters;

    SetParameters(base_parameters, R_base_parameters, Rand);

    vector<string> param_names = Rcpp::as<vector<string>>(params_priors.names());
    vector<Distribution> priors;

    for (unsigned int i = 0; i < params_priors.size(); ++i)
        priors.push_back(Distribution(Rcpp::as<string>(params_priors[i])));

    Likelihood lik(base_parameters, seed);
    MCMCReporter rep(1, 1, param_names);

    // Optimize_Priors(Rand, lik, rep, priors,
    //     maxeval, ftol_abs, true, n_threads > 1, n_threads);

//-----
    (void) maxeval;
    (void) ftol_abs;
    (void) n_threads;
    vector<double> initial(priors.size(), 0.0);
    for (unsigned int i = 0; i < priors.size(); ++i)
        initial[i] = priors[i].RandomInit(Rand);

    NelderMead_Priors(Rand, lik, rep, priors,
        initial, false, false, false);
//-----

    // Get data.frame as a data.table and return
    Rcpp::DataFrame df = Rcpp::DataFrame::create();
    df.push_back(Rcpp::IntegerVector::import(rep.trial.begin(), rep.trial.end()), "trial");
    df.push_back(Rcpp::NumericVector::import(rep.lp.begin(), rep.lp.end()), "lp");
    df.push_back(Rcpp::IntegerVector::import(rep.chain.begin(), rep.chain.end()), "chain");
    df.push_back(Rcpp::NumericVector::import(rep.ll.begin(), rep.ll.end()), "ll");
    for (unsigned int d = 0; d < rep.theta.size(); ++d)
        df.push_back(Rcpp::NumericVector::import(rep.theta[d].begin(), rep.theta[d].end()), rep.pnames[d]);

    return df;
}

// [[Rcpp::export]]
Rcpp::List cm_backend_sample_fit_test(Rcpp::List R_base_parameters, Rcpp::DataFrame posterior, 
    unsigned int n, unsigned long int seed, unsigned int n_threads = 1)
{
    // Enable multithreading
    #ifdef _OPENMP
    if (n_threads > 1)
        omp_set_num_threads(n_threads);
    #endif

    // Initialise parameters for this simulation
    Randomizer Rand(seed); // randomizer for fitting; randomizers for model runs are created in the Likelihood class.
    Parameters base_parameters;

    SetParameters(base_parameters, R_base_parameters, Rand);

    // Pick rows
    vector<unsigned int> rows(n, 0);
    for (unsigned int i = 0; i < rows.size(); ++i)
        rows[i] = Rand.Discrete(posterior.nrows());
    
    // Make reporters
    std::vector<Reporter*> reps(n, 0);
    
    // Get thetas
    vector<vector<double>> thetas;
    for (unsigned int i = 0; i < n; ++i)
    {
        vector<double> theta;
        for (unsigned int j = 4; j < posterior.size(); ++j)
            theta.push_back(Rcpp::as<Rcpp::NumericVector>(posterior[j])[rows[i]]);
        thetas.push_back(theta);
    }

    // Run simulations
    #pragma omp parallel for if(n_threads > 1) schedule(dynamic)
    for (unsigned int i = 0; i < n; ++i)
    {
        // TODO separate fitting and model seeds...
        Randomizer r(seed);
        Parameters P(base_parameters);

        CppChanges(thetas[i], P);

        Reporter rep = RunSimulation(P, r, thetas[i]);
        reps[i] = new Reporter(rep);
    }
    
    // Create list of data tables
    Rcpp::List dynamics;
    for (unsigned int i = 0; i < n; ++i)
    {
        Reporter& rep = *reps[i];

        // Create times
        Rcpp::NumericVector t(rep.n_times * rep.n_populations * rep.n_age_groups, 0.);
        for (unsigned int it = 0; it < rep.n_times; ++it)
            for (unsigned int j = 0; j < rep.n_populations * rep.n_age_groups; ++j)
                t[it * rep.n_populations * rep.n_age_groups + j] = base_parameters.time0 + 
                    it * base_parameters.time_step * base_parameters.report_every;

        // Create identifier columns
        int run = i + 1;
        Rcpp::DataFrame dynamics_df = Rcpp::DataFrame::create(
            Rcpp::Named("run") = Rcpp::rep(run, rep.n_times * rep.n_populations * rep.n_age_groups),
            Rcpp::Named("t") = t,
            Rcpp::Named("population") = Rcpp::rep(Rcpp::rep_each(Rcpp::seq(1, rep.n_populations), rep.n_age_groups), rep.n_times),
            Rcpp::Named("group") = Rcpp::rep(Rcpp::seq(1, rep.n_age_groups), rep.n_times * rep.n_populations)
        );

        // Allocate all columns to the dataframe
        for (unsigned int c = 0; c < rep.col_names.size(); ++c)
            dynamics_df.push_back(rep.data[c], rep.col_names[c]);

        // Add observer columns
        for (unsigned int c = 0; c < rep.obs.size(); ++c)
            dynamics_df.push_back(rep.obs[c], "obs" + to_string(c));

        // Set dataframe as a data.table
        Rcpp::Function setDT("setDT"); 
        setDT(dynamics_df);

        dynamics.push_back(dynamics_df);
        delete reps[i];
    }

    return dynamics;
}

// [[Rcpp::export]]
Rcpp::NumericVector cm_backend_sample_fit_rows(Rcpp::List R_base_parameters, Rcpp::DataFrame posterior, unsigned int n, unsigned long int seed)
{
    // Initialise parameters for this simulation
    Randomizer Rand(seed); // randomizer for fitting; randomizers for model runs are created in the Likelihood class.
    Parameters base_parameters;

    SetParameters(base_parameters, R_base_parameters, Rand);

    Rcpp::NumericVector rows(n, 0);

    for (unsigned int i = 0; i < n; ++i)
    {
        unsigned int row = Rand.Discrete(posterior.nrows());
        rows[i] = row + 1;
    }
    return rows;
}

// [[Rcpp::export]]
Rcpp::List cm_backend_sample_fit_2(Rcpp::List R_base_parameters_list, Rcpp::DataFrame posterior, 
    unsigned int n, unsigned long int seed, unsigned int n_threads = 1)
{
    // Enable multithreading
    #ifdef _OPENMP
    if (n_threads > 1)
        omp_set_num_threads(n_threads);
    #endif

    // Initialise parameters for this simulation
    Randomizer Rand(seed); // randomizer for fitting; randomizers for model runs are created in the Likelihood class.
    vector<Parameters> v_base_parameters;

    for (unsigned int i = 0; i < n; ++i) {
        v_base_parameters.push_back(Parameters());
        SetParameters(v_base_parameters.back(), Rcpp::as<Rcpp::List>(R_base_parameters_list[i]), Rand);
    }

    // Make reporters
    std::vector<Reporter*> reps(n, 0);
    
    // Get thetas
    vector<vector<double>> thetas;
    for (unsigned int i = 0; i < n; ++i)
    {
        vector<double> theta;
        for (unsigned int j = 4; j < posterior.size(); ++j)
            theta.push_back(Rcpp::as<Rcpp::NumericVector>(posterior[j])[i % posterior.nrows()]);
        thetas.push_back(theta);
    }

    // Run simulations
    #pragma omp parallel for if(n_threads > 1) schedule(dynamic)
    for (unsigned int i = 0; i < n; ++i)
    {
        // TODO separate fitting and model seeds...
        Randomizer r(seed);
        Parameters P(v_base_parameters[i]);

        CppChanges(thetas[i], P);

        Reporter rep = RunSimulation(P, r, thetas[i]);
        reps[i] = new Reporter(rep);
    }
    
    // Create list of data tables
    Rcpp::List dynamics;
    for (unsigned int i = 0; i < n; ++i)
    {
        Reporter& rep = *reps[i];

        // Create times
        Rcpp::NumericVector t(rep.n_times * rep.n_populations * rep.n_age_groups, 0.);
        for (unsigned int it = 0; it < rep.n_times; ++it)
            for (unsigned int j = 0; j < rep.n_populations * rep.n_age_groups; ++j)
                t[it * rep.n_populations * rep.n_age_groups + j] = v_base_parameters[i].time0 + 
                    it * v_base_parameters[i].time_step * v_base_parameters[i].report_every;

        // Create identifier columns
        int run = i + 1;
        Rcpp::DataFrame dynamics_df = Rcpp::DataFrame::create(
            Rcpp::Named("run") = Rcpp::rep(run, rep.n_times * rep.n_populations * rep.n_age_groups),
            Rcpp::Named("t") = t,
            Rcpp::Named("population") = Rcpp::rep(Rcpp::rep_each(Rcpp::seq(1, rep.n_populations), rep.n_age_groups), rep.n_times),
            Rcpp::Named("group") = Rcpp::rep(Rcpp::seq(1, rep.n_age_groups), rep.n_times * rep.n_populations)
        );

        // Allocate all columns to the dataframe
        for (unsigned int c = 0; c < rep.col_names.size(); ++c)
            dynamics_df.push_back(rep.data[c], rep.col_names[c]);

        // Add observer columns
        for (unsigned int c = 0; c < rep.obs.size(); ++c)
            dynamics_df.push_back(rep.obs[c], "obs" + to_string(c));

        // Set dataframe as a data.table
        Rcpp::Function setDT("setDT"); 
        setDT(dynamics_df);

        dynamics.push_back(dynamics_df);
        delete reps[i];
    }

    return dynamics;
}



// [[Rcpp::export]]
Rcpp::List cm_backend_sample_fit_multi(Rcpp::List R_base_parameters_list, Rcpp::List params_priors, 
    Rcpp::DataFrame posterior, unsigned int n, unsigned long int seed)
{
    // Initialise parameters for this simulation
    // TODO Rand also used for setting parameters -- is it actually used? this may cause issues with seeds for sample fit etc
    Randomizer Rand(seed); // randomizer for fitting; randomizers for model runs are created in the Likelihood class.

    // Set each of base_parameters individually
    vector<Parameters> base_parameters;
    for (unsigned int i = 0; i < R_base_parameters_list.size(); ++i)
    {
        base_parameters.push_back(Parameters());
        SetParameters(base_parameters.back(), Rcpp::as<Rcpp::List>(R_base_parameters_list[i]), Rand);
    }

    // xs = x_source[i][j] says what to use for x[i] in simulation j.
    // if xs >= 10000, use fixed_params[xs - 10000]; if xs < 10000, use theta[xs].
    vector<vector<unsigned int>> x_source;
    vector<double> fixed_params;
    vector<string> theta_names;
    vector<Distribution> theta_priors;

    get_multi_params(params_priors, base_parameters.size(), x_source, fixed_params, theta_names, theta_priors);

    Rcpp::List dynamics;

    for (unsigned int i = 0; i < n; ++i)
    {
        for (unsigned int j = 0; j < base_parameters.size(); ++j)
        {
            // TODO separate fitting and model seeds...
            Randomizer r(seed);
            Parameters P(base_parameters[j]);

            unsigned int row = Rand.Discrete(posterior.nrows());
            vector<double> theta(x_source.size(), 0);
            for (unsigned int k = 0; k < theta.size(); ++k)
            {
                if (x_source[k][j] >= 10000)
                    theta[k] = fixed_params[x_source[k][j] - 10000];
                else
                    theta[k] = Rcpp::as<Rcpp::NumericVector>(posterior[x_source[k][j] + 4])[row];

            }

            CppChanges(theta, P);

            Reporter rep = RunSimulation(P, r, theta);

            // TODO Dataframe construction -- copied from old reporter.cpp ----
            // TODO Can move some of this outside the loop... and refactor with code above which is similar
        
            // Create times
            Rcpp::NumericVector t(rep.n_times * rep.n_populations * rep.n_age_groups, 0.);
            for (unsigned int it = 0; it < rep.n_times; ++it)
                for (unsigned int j = 0; j < rep.n_populations * rep.n_age_groups; ++j)
                    t[it * rep.n_populations * rep.n_age_groups + j] = P.time0 + it * P.time_step * P.report_every;

            // Create identifier columns
            int run = i + 1;
            Rcpp::DataFrame dynamics_df = Rcpp::DataFrame::create(
                Rcpp::Named("run") = Rcpp::rep(run, rep.n_times * rep.n_populations * rep.n_age_groups),
                Rcpp::Named("t") = t,
                Rcpp::Named("population") = Rcpp::rep(Rcpp::rep_each(Rcpp::seq(1, rep.n_populations), rep.n_age_groups), rep.n_times),
                Rcpp::Named("group") = Rcpp::rep(Rcpp::seq(1, rep.n_age_groups), rep.n_times * rep.n_populations)
            );

            // Allocate all columns to the dataframe
            for (unsigned int c = 0; c < rep.col_names.size(); ++c)
                dynamics_df.push_back(rep.data[c], rep.col_names[c]);

            // Add observer columns
            for (unsigned int c = 0; c < rep.obs.size(); ++c)
                dynamics_df.push_back(rep.obs[c], "obs" + to_string(c));

            // Set dataframe as a data.table
            Rcpp::Function setDT("setDT"); 
            setDT(dynamics_df);

            dynamics.push_back(dynamics_df);
        }
    }

    return dynamics;
}


// [[Rcpp::export]]
Rcpp::List cm_backend_pfilter(Rcpp::List R_base_parameters, Rcpp::DataFrame posterior, 
    unsigned int row, unsigned long int seed, unsigned int n_particles, vector<double> times, 
    double adjust_nu, double adjust_scale, double adjust_b, 
    unsigned int n_threads = 1)
{
    // Enable multithreading
    #ifdef _OPENMP
    if (n_threads > 1)
        omp_set_num_threads(n_threads);
    #endif

    // Initialise parameters for this simulation
    Randomizer Rand(seed);
    Parameters base_parameters;
    SetParameters(base_parameters, R_base_parameters, Rand);
    times.push_back(base_parameters.time1);

    // Create posterior sample
    vector<double> x;
    for (unsigned int i = 4; i < posterior.size(); ++i)
        x.push_back(Rcpp::as<Rcpp::NumericVector>(posterior[i])[row]);

    // Storage for particles and associated properties
    vector<Parameters*> params;
    vector<Randomizer*> randomizers;
    vector<Reporter*> reporters;
    vector<Metapopulation*> particles;
    vector<double> particle_ll(n_particles, 0.0);
    vector<double> particle_ll_norm(n_particles, 0.0);
    vector<double> particle_ll_sum(n_particles, 0.0);
    
    // Initialise particle-specific objects
    for (unsigned int j = 0; j < n_particles; ++j)
    {
        params.push_back(new Parameters(base_parameters));
        randomizers.push_back(new Randomizer(Rand.Discrete(std::numeric_limits<unsigned int>::max())));
        reporters.push_back(new Reporter(*params[j]));
        particles.push_back(new Metapopulation(*params[j]));
        CppChanges(x, *params[j]);
    }
    
    // Temp - return
    Rcpp::NumericVector 
        ret_t(times.size() * n_particles, 0.0),
        ret_j(times.size() * n_particles, 0.0),
        ret_l(times.size() * n_particles, 0.0),
        ret_s(times.size() * n_particles, 0.0),
        ret_a(times.size() * n_particles, 0.0),
        ret_i(times.size() * n_particles, 0.0);

    // Step through time periods
    double t_prev = base_parameters.time0;
    
    for (unsigned int i = 0; i < times.size(); ++i)
    {
        Rcpp::Rcout << "\n\nTo time " << times[i] << "...\n";

        // Propagate
        double particle_ll_max = -std::numeric_limits<double>::infinity();
        
        #pragma omp parallel for if(n_threads > 1) schedule(dynamic)
        for (unsigned int j = 0; j < n_particles; ++j)
        {
            particles[j]->RunUntil(*params[j], *randomizers[j], *reporters[j], x, t_prev, times[i], i == times.size() - 1);
            particle_ll[j] = CppLogLikelihood(x, *params[j], *reporters[j], t_prev, times[i]);
            particle_ll_sum[j] += particle_ll[j];
            if (particle_ll[j] > particle_ll_max)
                particle_ll_max = particle_ll[j];
        }

        // Choose particles for resampling
        std::vector<unsigned int> parents(n_particles, 0);
        for (unsigned int j = 0; j < n_particles; ++j)
            particle_ll_norm[j] = exp(particle_ll[j] - particle_ll_max);
        Rand.Sample(n_particles, particle_ll_norm, parents, true);
        
        // Resample particles
        for (unsigned int j = 0; j < parents.size(); ++j)
        {
            unsigned int k = parents[j];
            
            params.push_back(new Parameters(*params[k]));
            params.back()->changes.Recapture(*params.back());
            reporters.push_back(new Reporter(*reporters[k]));
            randomizers.push_back(new Randomizer(*randomizers[k]));
            particles.push_back(new Metapopulation(*particles[k]));
            particle_ll.push_back(particle_ll[k]);
            particle_ll_sum.push_back(particle_ll_sum[k]);
            
            ret_i[i * n_particles + j] = k;
        }
        for (unsigned int j = 0; j < n_particles; ++j)
        {
            delete params[j];
            delete reporters[j];
            delete randomizers[j];
            delete particles[j];
        }
        params.assign(params.begin() + n_particles, params.end());
        reporters.assign(reporters.begin() + n_particles, reporters.end());
        randomizers.assign(randomizers.begin() + n_particles, randomizers.end());
        particles.assign(particles.begin() + n_particles, particles.end());
        particle_ll.assign(particle_ll.begin() + n_particles, particle_ll.end());
        particle_ll_sum.assign(particle_ll_sum.begin() + n_particles, particle_ll_sum.end());
        
        for (unsigned int j = 0; j < n_particles; j++)
        {
            ret_t[i * n_particles + j] = times[i];
            ret_j[i * n_particles + j] = j;
            ret_l[i * n_particles + j] = particle_ll[j];
            ret_s[i * n_particles + j] = particle_ll_sum[j];
            ret_a[i * n_particles + j] = params[j]->adjustment;
        }
        
        // Mutation step
        for (unsigned int j = 0; j < n_particles; ++j)
        {
            params[j]->adjustment = exp(Rand.TDist(adjust_nu) * adjust_scale + (log(params[j]->adjustment) * adjust_b));
            // params[j]->adjustment = exp(Rand.TDist(adjust_nu, exp(log(params[j]->adjustment) * (1.0 - adjust_b)), adjust_scale));
        }
        
        t_prev = times[i];
    }
    
    Rcpp::DataFrame ret = Rcpp::DataFrame::create(
        Rcpp::Named("t") = ret_t,
        Rcpp::Named("particle") = ret_j,
        Rcpp::Named("ll") = ret_l,
        Rcpp::Named("ll_sum") = ret_s,
        Rcpp::Named("adj") = ret_a,
        Rcpp::Named("parent") = ret_i
    );
    
    // Create list of data tables
    Rcpp::List dynamics;
    for (unsigned int i = 0; i < n_particles; ++i)
    {
        Reporter& rep = *reporters[i];

        // Create times
        Rcpp::NumericVector t(rep.n_times * rep.n_populations * rep.n_age_groups, 0.);
        for (unsigned int it = 0; it < rep.n_times; ++it)
            for (unsigned int j = 0; j < rep.n_populations * rep.n_age_groups; ++j)
                t[it * rep.n_populations * rep.n_age_groups + j] = base_parameters.time0 + 
                    it * base_parameters.time_step * base_parameters.report_every;

        // Create identifier columns
        int run = i + 1;
        Rcpp::DataFrame dynamics_df = Rcpp::DataFrame::create(
            Rcpp::Named("run") = Rcpp::rep(run, rep.n_times * rep.n_populations * rep.n_age_groups),
            Rcpp::Named("t") = t,
            Rcpp::Named("population") = Rcpp::rep(Rcpp::rep_each(Rcpp::seq(1, rep.n_populations), rep.n_age_groups), rep.n_times),
            Rcpp::Named("group") = Rcpp::rep(Rcpp::seq(1, rep.n_age_groups), rep.n_times * rep.n_populations)
        );

        // Allocate all columns to the dataframe
        for (unsigned int c = 0; c < rep.col_names.size(); ++c)
            dynamics_df.push_back(rep.data[c], rep.col_names[c]);

        // Add observer columns
        for (unsigned int c = 0; c < rep.obs.size(); ++c)
            dynamics_df.push_back(rep.obs[c], "obs" + to_string(c));

        // Set dataframe as a data.table
        Rcpp::Function setDT("setDT"); 
        setDT(dynamics_df);

        dynamics.push_back(dynamics_df);
        delete params[i];
        delete reporters[i];
        delete randomizers[i];
        delete particles[i];
    }
    
    Rcpp::Rcout << "NB does not evaluate last point. Fix.\n";

    return Rcpp::List::create(dynamics, ret);
}


// [[Rcpp::export]]
std::vector<unsigned int> test_sampler(unsigned int N, std::vector<double> p, bool stratified, unsigned int seed = 0, unsigned int burn_in = 0)
{
    std::vector<unsigned int> i;
    Randomizer Rand(seed);
    for (unsigned int b = 0; b < burn_in; ++b)
        Rand.Sample(N, p, i, stratified);
    Rand.Sample(N, p, i, stratified);
    return i;
}

// [[Rcpp::export]]
std::vector<unsigned int> test_multinomial(unsigned int N, std::vector<double> p, unsigned int seed = 0)
{
    Randomizer Rand(seed);

    std::vector<unsigned int> draw(p.size(), 0);
    Rand.Multinomial(N, p, draw);
    return draw;
}



// [[Rcpp::export]]
double TestEfficacyFold(double eff, double fold)
{
    return efficacy_fold(eff, fold);
}

// [[Rcpp::export]]
double TestInfectionToSevere(double ei)
{
    return infection_to_severe(ei);
}

// [[Rcpp::export]]
void TestBoost(double ei, double ed_i, double eh_d, double em_d, double frac_boosted, double boost_fold)
{
    auto VEfunc = [](double x, double given) { return (x - given) / (1.0 - given); };
    
    auto boost = [=](double& eff_inf, double eff_dis, double& eff_hosp, double& eff_mort, double frac_boosted = 1.0)
    {
        if (eff_inf == 1.0) return;
        double eff_boosted = efficacy_fold(eff_inf, boost_fold);
        double eff_sev = infection_to_severe(eff_boosted);
        eff_inf = (1.0 - frac_boosted) * eff_inf + frac_boosted * eff_boosted;
        eff_hosp = (1.0 - frac_boosted) * eff_hosp + frac_boosted * VEfunc(eff_sev, 1.0 - (1.0 - eff_boosted) * (1.0 - eff_dis));
        eff_mort = (1.0 - frac_boosted) * eff_mort + frac_boosted * VEfunc(eff_sev, 1.0 - (1.0 - eff_boosted) * (1.0 - eff_dis));
    };
    
    boost(ei, ed_i, eh_d, em_d, frac_boosted);
    Rcpp::Rcout << ei << " " << ed_i << " " << eh_d << " " << em_d << "\n";
}
