// sim_compartment.h

#ifndef SIM_COMPARTMENT_H
#define SIM_COMPARTMENT_H

#include <vector>
using namespace std;
class Compartment;

//
// MODEL DYNAMICS
//

struct Parameters;
class Randomizer;
class Reporter;

// A population of individuals, with SVELI3R dynamics.
class Population
{
public:
    // Construct a population with the specified size by age group; initially all uninfected
    Population(Parameters& P, unsigned int pindex);

    // Do seeding and calculate contagiousness
    void Contagiousness(Parameters& P, Randomizer& Rand, double t, vector<double>& contag, vector<double>& contag2, vector<double>& contag3);

    // Execute one time step's events
    void Tick(Parameters& P, Randomizer& Rand, double t, vector<double>& infec, vector<double>& infec2, vector<double>& infec3, Reporter& rep);

    // Print full population details
    void DebugPrint() const;

//private:
    vector<double> lambda;
    vector<double> lambda2;
    vector<double> lambda3;
    vector<double> N, S, R, R2, R3;              // Total number, susceptible, recovered, recovered 2, recovered 3
    vector<double> EV;                           // Ever vaccinated
    vector<Compartment> Va1, Va2, Vb1, Vb2;      // Total number vaccinated a1, a2, b1, b2
    vector<double> Va3, Vb3;                     // Total number vaccinated a3, b3
    vector<Compartment> E, L, Ip, Ia, Is;        // Strain 1 exposed, latent (exposed to asymptomatic), presymptomatic, asymptomatic, symptomatic
    vector<Compartment> E2, L2, Ip2, Ia2, Is2;   // Strain 2 exposed, latent (exposed to asymptomatic), presymptomatic, asymptomatic, symptomatic
    vector<Compartment> E3, L3, Ip3, Ia3, Is3;   // Strain 3 exposed, latent (exposed to asymptomatic), presymptomatic, asymptomatic, symptomatic
    vector<Compartment> to_death, to_hosp, hosp, to_icu, icu;  // Burden compartments
    vector<Compartment> to_death_V1, to_death_V2, to_death_V3; // Burden compartments to track vaccinated deaths
    vector<Compartment> to_hosp_V1, to_hosp_V2, to_hosp_V3;    // Burden compartments to track vaccinated hospital admissions
    unsigned int seed_row, seed_row2, seed_row3; // Which seed event is next (strain 1, 2, 3)
    unsigned int p;                              // Which population this is
    vector<vector<Compartment>> pc;              // User-specified process compartments, indexed by process id, then group
    vector<unsigned int> ni_out;                 // Temporary storage
    vector<double> nd_out;                       // Temporary storage
    vector<double> pci;                          // Temporary storage
    vector<double> pco;                          // Temporary storage
};

// A metapopulation, containing multiple subpopulations.
class Metapopulation
{
public:
    Metapopulation(Parameters& P);

    // Execute one time step's events
    bool Tick(Parameters& P, Randomizer& Rand, double t, unsigned int ts, Reporter& rep);

    // Run the model
    void Run(Parameters& P, Randomizer& Rand, Reporter& rep, vector<double> x_fit = vector<double>());

//private:
    vector<vector<double>> contag, contag2, contag3;
    vector<vector<double>> infec, infec2, infec3;
    vector<Population> pops;
    vector<double> x;
};

#endif
