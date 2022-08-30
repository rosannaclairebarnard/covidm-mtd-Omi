// user_defined.h

// Declarations of C++ functions which require definition by the user

#ifndef USER_DEFINED_H
#define USER_DEFINED_H

#include <vector>
#include <limits>
#include <iostream>
#include "parameters.h"
#include "sim_compartment.h"
using namespace std;

// Declarations of C++ functions which require definition by the user
void CppChanges(const vector<double>& x, Parameters& P);
double CppLogLikelihood(const vector<double>& x, Parameters& P, Reporter& dyn, 
    double t_min = -std::numeric_limits<double>::infinity(), double t_max = std::numeric_limits<double>::infinity());
bool CppObserver(Parameters& P, Randomizer& R, Reporter& dyn, double t, vector<double>& x, Metapopulation& mp);

#endif