// parameters.h

#ifndef PARAMETERS_H
#define PARAMETERS_H

// [[Rcpp::plugins(cpp11)]]
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
#include "observer.h"
#include "user_include.h"

struct Matrix;
struct Change;

void ParamSet(Discrete& variable, Rcpp::RObject& value);
void ParamSet(vector<double>& variable, Rcpp::RObject& value);
void ParamSet(Matrix& variable, Rcpp::RObject& value);
void ParamSet(vector<Matrix>& variable, Rcpp::RObject& value);
void ParamSet(vector<ProcessSpec>& variable, Rcpp::RObject& value);

// Population-level parameters for batch of simulations
struct PopulationParameters
{
public:
    PopulationParameters() : needs_recalc(true) {}
    
    // Output the current vaccine efficacy and natural protection parameters, CSV-formatted, to ostream out; 
    // first_line = true prints header as well.
    void OutputVELine(double t, unsigned int p, bool header_only, ostream& out);

    bool needs_recalc;  // does contact matrix need recalculation?
    Matrix cm;          // contact matrix

    Discrete dE;
    Discrete dIp;       // TODO: any need for these to be age-specific?
    Discrete dIa;
    Discrete dIs;
    Discrete dE2;
    Discrete dIp2;
    Discrete dIa2;
    Discrete dIs2;
    Discrete dE3;
    Discrete dIp3;
    Discrete dIa3;
    Discrete dIs3;

    Discrete dVa1;
    Discrete dVb1;
    Discrete dVa2;
    Discrete dVb2;

    vector<double> size;
    vector<double> imm0;
    vector<Matrix> matrices;
    vector<double> contact;
    vector<double> contact_mult;
    vector<double> contact_lowerto;
    vector<double> u, u2, u3;
    vector<double> fIp;
    vector<double> fIa;
    vector<double> fIs;
    vector<double> y, y2, y3;
    vector<double> omega;
    vector<double> tau;
    vector<double> pi_r, pi_r2, pi_r3, pi2_r, pi2_r2, pi2_r3, pi3_r, pi3_r2, pi3_r3;
    vector<double> wn, wn2, wn3;
    vector<double> va1, wva1, ei_va1, ei2_va1, ei3_va1, ed_va1i, ed_va1i2, ed_va1i3;
    vector<double>      wva2, ei_va2, ei2_va2, ei3_va2, ed_va2i, ed_va2i2, ed_va2i3;
    vector<double>      wva3, ei_va3, ei2_va3, ei3_va3, ed_va3i, ed_va3i2, ed_va3i3;
    vector<double> p_boost_va2;
    vector<double> vb1, wvb1, ei_vb1, ei2_vb1, ei3_vb1, ed_vb1i, ed_vb1i2, ed_vb1i3;
    vector<double>      wvb2, ei_vb2, ei2_vb2, ei3_vb2, ed_vb2i, ed_vb2i2, ed_vb2i3;
    vector<double>      wvb3, ei_vb3, ei2_vb3, ei3_vb3, ed_vb3i, ed_vb3i2, ed_vb3i3;
    vector<double> p_boost_vb2;
    vector<double> A;
    vector<double> B;
    vector<double> D;
    vector<double> season_A;    // note - size 1 vector
    vector<double> season_T;    // note - size 1 vector
    vector<double> season_phi;  // note - size 1 vector

    vector<double> ifr1, ifr2, ifr3;
    vector<double> ihr1, ihr2, ihr3;
    vector<double> iir1, iir2, iir3;
    
    vector<double> pd_ri,  pd_ri2,  pd_ri3;
    vector<double> pd_r2i, pd_r2i2, pd_r2i3;
    vector<double> pd_r3i, pd_r3i2, pd_r3i3;
    vector<double> ph_rd,  ph_rd2,  ph_rd3;
    vector<double> ph_r2d, ph_r2d2, ph_r2d3;
    vector<double> ph_r3d, ph_r3d2, ph_r3d3;
    vector<double> pm_rd,  pm_rd2,  pm_rd3;
    vector<double> pm_r2d, pm_r2d2, pm_r2d3;
    vector<double> pm_r3d, pm_r3d2, pm_r3d3;
    vector<double> pt_ri,  pt_ri2,  pt_ri3;
    vector<double> pt_r2i, pt_r2i2, pt_r2i3;
    vector<double> pt_r3i, pt_r3i2, pt_r3i3;
    
    vector<double> eh_va1d, eh_va1d2, eh_va1d3;
    vector<double> eh_va2d, eh_va2d2, eh_va2d3;
    vector<double> eh_va3d, eh_va3d2, eh_va3d3;
    vector<double> eh_vb1d, eh_vb1d2, eh_vb1d3;
    vector<double> eh_vb2d, eh_vb2d2, eh_vb2d3;
    vector<double> eh_vb3d, eh_vb3d2, eh_vb3d3;
    vector<double> em_va1d, em_va1d2, em_va1d3;
    vector<double> em_va2d, em_va2d2, em_va2d3;
    vector<double> em_va3d, em_va3d2, em_va3d3;
    vector<double> em_vb1d, em_vb1d2, em_vb1d3;
    vector<double> em_vb2d, em_vb2d2, em_vb2d3;
    vector<double> em_vb3d, em_vb3d2, em_vb3d3;
    vector<double> et_va1i, et_va1i2, et_va1i3;
    vector<double> et_va2i, et_va2i2, et_va2i3;
    vector<double> et_va3i, et_va3i2, et_va3i3;
    vector<double> et_vb1i, et_vb1i2, et_vb1i3;
    vector<double> et_vb2i, et_vb2i2, et_vb2i3;
    vector<double> et_vb3i, et_vb3i2, et_vb3i3;

    Discrete dDeath, dHosp, lHosp, dICU, lICU;
    vector<Discrete> lHosp2;

    vector<double> seed_times;
    vector<double> seed_times2;
    vector<double> seed_times3;
    Discrete dist_seed_ages;

    ///Observer observer;
    ///vector<ScheduleEntry> schedule;

    string name;
    vector<string> group_names;

    bool Set(Parameters* parent, string& name, Rcpp::RObject& value);
    bool Set(Parameters* parent, string& name, vector<double>& value);
    
    // For boosters - storage . . .
    vector<double> X_pi_r, X_pi_r2, X_pi_r3;
    vector<double> X_pi2_r, X_pi2_r2, X_pi2_r3;
    vector<double> X_pi3_r, X_pi3_r2, X_pi3_r3;
    
    vector<double> X_ei_va1, X_ei2_va1, X_ei3_va1;
    vector<double> X_ei_va2, X_ei2_va2, X_ei3_va2;
    vector<double> X_ei_va3, X_ei2_va3, X_ei3_va3;
    vector<double> X_ei_vb1, X_ei2_vb1, X_ei3_vb1;
    vector<double> X_ei_vb2, X_ei2_vb2, X_ei3_vb2;
    vector<double> X_ei_vb3, X_ei2_vb3, X_ei3_vb3;
    
    vector<double> X_ph_rd,  X_ph_rd2,  X_ph_rd3;
    vector<double> X_ph_r2d, X_ph_r2d2, X_ph_r2d3;
    vector<double> X_ph_r3d, X_ph_r3d2, X_ph_r3d3;
    vector<double> X_pm_rd,  X_pm_rd2,  X_pm_rd3;
    vector<double> X_pm_r2d, X_pm_r2d2, X_pm_r2d3;
    vector<double> X_pm_r3d, X_pm_r3d2, X_pm_r3d3;

    vector<double> X_eh_va1d, X_eh_va1d2, X_eh_va1d3;
    vector<double> X_eh_va2d, X_eh_va2d2, X_eh_va2d3;
    vector<double> X_eh_va3d, X_eh_va3d2, X_eh_va3d3;
    vector<double> X_eh_vb1d, X_eh_vb1d2, X_eh_vb1d3;
    vector<double> X_eh_vb2d, X_eh_vb2d2, X_eh_vb2d3;
    vector<double> X_eh_vb3d, X_eh_vb3d2, X_eh_vb3d3;
    vector<double> X_em_va1d, X_em_va1d2, X_em_va1d3;
    vector<double> X_em_va2d, X_em_va2d2, X_em_va2d3;
    vector<double> X_em_va3d, X_em_va3d2, X_em_va3d3;
    vector<double> X_em_vb1d, X_em_vb1d2, X_em_vb1d3;
    vector<double> X_em_vb2d, X_em_vb2d2, X_em_vb2d3;
    vector<double> X_em_vb3d, X_em_vb3d2, X_em_vb3d3;

    void Recalculate()
    {
        if (needs_recalc)
        {
            if (matrices.empty())
                throw logic_error("No contact matrices defined.");
            if (matrices.size() != contact.size())
                throw logic_error("Number of contact components not equal to number of matrices.");

            unsigned int ncol = matrices[0].nc;
            unsigned int nrow = matrices[0].x.size() / ncol;

            auto c_mult = [&](unsigned int m) { if (!contact_mult.empty()) return contact_mult[m]; return 1.0; };
            auto c_lowerto = [&](unsigned int m) { if (!contact_lowerto.empty()) return contact_lowerto[m]; return std::numeric_limits<double>::max(); };

            cm = Matrix(0, nrow, ncol);

            for (unsigned int r = 0; r < nrow; ++r)
                for (unsigned int c = 0; c < ncol; ++c)
                    for (unsigned int m = 0; m < matrices.size(); ++m)
                        cm(r, c) += matrices[m](r, c) * min(contact[m] * c_mult(m), c_lowerto(m));

            needs_recalc = false;
        }
    }
};

struct Parameters
{
public:
    void FilterForRun(unsigned int r);

    string model;
    double time_step;
    double time0;
    double time1;
    unsigned int report_every;
    bool fast_multinomial;
    bool deterministic;

    vector<PopulationParameters> pop;

    vector<ProcessSpec> processes;
    ///Observer observer;
    Matrix travel;
    ChangeSet changes;

    double vax_limit;
    double adjustment;
    vector<double> extra_vb1;
};

void SetParameters(Parameters& P, Rcpp::List list, Randomizer& Rand);

#endif
