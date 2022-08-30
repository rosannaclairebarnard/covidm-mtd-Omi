// parameters.cpp

#include "parameters.h"

// Output the current vaccine efficacy and natural protection parameters, CSV-formatted, to ostream out; 
// first_line = true prints header as well.

#define EFFICACY_PRINT_I_SUB(p_or_e, imm, strain) \
    out << t << ", " << p << ", " << g << ", " \
        << #p_or_e << ", " << "i" << ", " << #imm << ", " << "" << ", " << #strain << ", " \
        << p_or_e ## i ## strain ## _ ## imm [g] << ", " << 1 << ", " << p_or_e ## i ## strain ## _ ## imm [g] << "\n";

#define EFFICACY_PRINT_T_SUB(p_or_e, imm, strain) \
{ \
        double t_given_i = p_or_e ## t ## _ ## imm ## i ## strain [g]; \
        double i_overall = p_or_e ## i ## strain ## _ ## imm [g]; \
        double absolute_value = t_given_i * (1.0 - i_overall) + i_overall; \
        out << t << ", " << p << ", " << g << ", " \
            << #p_or_e << ", " << "t" << ", " << #imm << ", " << "i" << ", " << #strain << ", " \
            << t_given_i << ", " << i_overall << ", " << absolute_value << "\n"; \
    }

#define EFFICACY_PRINT_D_SUB(p_or_e, imm, strain) \
    { \
        double d_given_i = p_or_e ## d ## _ ## imm ## i ## strain [g]; \
        double i_overall = p_or_e ## i ## strain ## _ ## imm [g]; \
        double absolute_value = d_given_i * (1.0 - i_overall) + i_overall; \
        out << t << ", " << p << ", " << g << ", " \
            << #p_or_e << ", " << "d" << ", " << #imm << ", " << "i" << ", " << #strain << ", " \
            << d_given_i << ", " << i_overall << ", " << absolute_value << "\n"; \
    }

#define EFFICACY_PRINT_S_SUB(p_or_e, outcome, imm, strain) \
    { \
        double s_given_d = p_or_e ## outcome ## _ ## imm ## d ## strain [g]; \
        double i_overall = p_or_e ## i ## strain ## _ ## imm [g]; \
        double d_given_i = p_or_e ## d ## _ ## imm ## i ## strain [g]; \
        double d_overall = d_given_i * (1.0 - i_overall) + i_overall; \
        double absolute_value = s_given_d * (1.0 - d_overall) + d_overall; \
        out << t << ", " << p << ", " << g << ", " \
            << #p_or_e << ", " << #outcome << ", " << #imm << ", " << "d" << ", " << #strain << ", " \
            << s_given_d << ", " << d_overall << ", " << absolute_value << "\n"; \
    }

#define EFFICACY_PRINT_H_SUB(p_or_e, imm, strain) EFFICACY_PRINT_S_SUB(p_or_e, h, imm, strain)
#define EFFICACY_PRINT_M_SUB(p_or_e, imm, strain) EFFICACY_PRINT_S_SUB(p_or_e, m, imm, strain)
    
#define EFFICACY_PRINT_I(p_or_e, imm)           \
    EFFICACY_PRINT_I_SUB(p_or_e, imm, )         \
    EFFICACY_PRINT_I_SUB(p_or_e, imm, 2)        \
    EFFICACY_PRINT_I_SUB(p_or_e, imm, 3)        

#define EFFICACY_PRINT_T(p_or_e, imm)           \
    EFFICACY_PRINT_T_SUB(p_or_e, imm, )         \
    EFFICACY_PRINT_T_SUB(p_or_e, imm, 2)        \
    EFFICACY_PRINT_T_SUB(p_or_e, imm, 3)        

#define EFFICACY_PRINT_D(p_or_e, imm)          \
    EFFICACY_PRINT_D_SUB(p_or_e, imm,)         \
    EFFICACY_PRINT_D_SUB(p_or_e, imm,2)        \
    EFFICACY_PRINT_D_SUB(p_or_e, imm,3)        

#define EFFICACY_PRINT_H(p_or_e, imm)          \
    EFFICACY_PRINT_H_SUB(p_or_e, imm,)         \
    EFFICACY_PRINT_H_SUB(p_or_e, imm,2)        \
    EFFICACY_PRINT_H_SUB(p_or_e, imm,3)        

#define EFFICACY_PRINT_M(p_or_e, imm)          \
    EFFICACY_PRINT_M_SUB(p_or_e, imm,)         \
    EFFICACY_PRINT_M_SUB(p_or_e, imm,2)        \
    EFFICACY_PRINT_M_SUB(p_or_e, imm,3) 
        
#define EFFICACY_PRINT_ALL(p_or_e, imm)                        \
    EFFICACY_PRINT_I(p_or_e, imm)                              \
    EFFICACY_PRINT_T(p_or_e, imm)                              \
    EFFICACY_PRINT_D(p_or_e, imm)                              \
    EFFICACY_PRINT_H(p_or_e, imm)                              \
    EFFICACY_PRINT_M(p_or_e, imm)


void PopulationParameters::OutputVELine(double t, unsigned int p, bool header_only, ostream& out)
{
    if (header_only) {
        out << "t, population, group, p_or_e, outcome, imm, given, strain, conditional_val, given_val, absolute_val\n";
        return;
    }
    
    for (unsigned int g = 0; g < size.size(); ++g)
    {
        EFFICACY_PRINT_ALL(e, va1)
        EFFICACY_PRINT_ALL(e, va2)
        EFFICACY_PRINT_ALL(e, va3)
        
        EFFICACY_PRINT_ALL(e, vb1)
        EFFICACY_PRINT_ALL(e, vb2)
        EFFICACY_PRINT_ALL(e, vb3)

        EFFICACY_PRINT_ALL(p, r)
        EFFICACY_PRINT_ALL(p, r2)
        EFFICACY_PRINT_ALL(p, r3)
    }
}

// Helpers to set parameters
#define _CheckSet(variable) else if (name == #variable) { ParamSet(variable, value); }
#define _CheckSetParent(P, variable) else if (P != 0 && name == #variable) { ParamSet(P->variable, value); return false; }

void ParamSet(double& variable, Rcpp::RObject& value)
{
    variable = Rcpp::as<double>(value);
}

void ParamSet(Discrete& variable, Rcpp::RObject& value)
{
    variable = Rcpp::as<vector<double>>(value);
}

void ParamSet(vector<double>& variable, Rcpp::RObject& value)
{
    variable = Rcpp::as<vector<double>>(value);
}

void ParamSet(Matrix& variable, Rcpp::RObject& value)
{
    Rcpp::NumericMatrix rhs = Rcpp::as<Rcpp::NumericMatrix>(value);
    variable = Matrix(0.0, rhs.nrow(), rhs.ncol());
    for (int r = 0; r < rhs.nrow(); ++r)
        for (int c = 0; c < rhs.ncol(); ++c)
            variable(r, c) = rhs(r, c);
}

void ParamSet(vector<Matrix>& variable, Rcpp::RObject& value)
{
    Rcpp::List ml = Rcpp::as<Rcpp::List>(value);
    for (unsigned int i = 0; i < ml.size(); ++i)
    {
        Rcpp::RObject mat = Rcpp::as<Rcpp::RObject>(ml[i]);
        ParamSet(variable[i], mat);
    }
}

void ParamSet(vector<ProcessSpec>& variable, Rcpp::RObject& value)
{
    Rcpp::List pl = Rcpp::as<Rcpp::List>(value);
    unsigned int pc_id = 0;
    vector<string> pc_names;
    for (unsigned int i = 0; i < pl.size(); ++i)
    {
        Rcpp::List pli = Rcpp::as<Rcpp::List>(pl[i]);

        ProcessSpec process;
        process.source_name = Rcpp::as<string>(pli["source"]);
        process.type = Rcpp::as<string>(pli["type"]);

        process.names = Rcpp::as<vector<string>>(pli["names"]);
        process.ids = vector<unsigned int>(process.names.size(), 0);
        for (unsigned int j = 0; j < process.ids.size(); ++j)
        {
            if (process.names[j] == "null")
            {
                process.ids[j] = Null;
            }
            else
            {
                process.ids[j] = pc_id++;
                pc_names.push_back(process.names[j]);
            }
        }
        process.report = Rcpp::as<vector<string>>(pli["report"]);

        Matrix m_prob, m_delays;
        Rcpp::RObject r_prob = Rcpp::as<Rcpp::RObject>(pli["prob"]);
        Rcpp::RObject r_delays = Rcpp::as<Rcpp::RObject>(pli["delays"]);
        ParamSet(m_prob, r_prob);
        ParamSet(m_delays, r_delays);

        for (unsigned int c = 0; c < m_prob.NCol(); ++c)
        {
            process.prob.push_back(vector<double>(m_prob.NRow(), 0.));
            for (unsigned int r = 0; r < m_prob.NRow(); ++r)
                process.prob[c][r] = m_prob(r, c);
        }

        for (unsigned int r = 0; r < m_delays.NRow(); ++r)
        {
            process.delays.push_back(Discrete());
            std::vector<double> uw(m_delays.NCol(), 0.);
            for (unsigned int c = 0; c < m_delays.NCol(); ++c)
                uw[c] = m_delays(r, c);
            process.delays.back() = uw;
        }

        variable.push_back(process);
    }

    // Set source_ids
    for (auto& pr : variable)
    {
        auto sn = std::find(pc_names.begin(), pc_names.end(), pr.source_name);
        if (sn == pc_names.end())
        {
            if (pr.source_name == "newE")
                pr.source_id = src_newE;
            else if (pr.source_name == "newI")
                pr.source_id = src_newI;
            else if (pr.source_name == "newIp")
                pr.source_id = src_newIp;
            else if (pr.source_name == "newIs")
                pr.source_id = src_newIs;
            else if (pr.source_name == "newIa")
                pr.source_id = src_newIa;
            else if (pr.source_name == "newE2")
                pr.source_id = src_newE2;
            else if (pr.source_name == "newI2")
                pr.source_id = src_newI2;
            else if (pr.source_name == "newIp2")
                pr.source_id = src_newIp2;
            else if (pr.source_name == "newIs2")
                pr.source_id = src_newIs2;
            else if (pr.source_name == "newIa2")
                pr.source_id = src_newIa2;
            else if (pr.source_name == "newE3")
                pr.source_id = src_newE3;
            else if (pr.source_name == "newI3")
                pr.source_id = src_newI3;
            else if (pr.source_name == "newIp3")
                pr.source_id = src_newIp3;
            else if (pr.source_name == "newIs3")
                pr.source_id = src_newIs3;
            else if (pr.source_name == "newIa3")
                pr.source_id = src_newIa3;
            else if (pr.source_name == "newEE2E3")
                pr.source_id = src_newEE2E3;
            else if (pr.source_name == "newII2I3")
                pr.source_id = src_newII2I3;
            else if (pr.source_name == "newIpIp2Ip3")
                pr.source_id = src_newIpIp2Ip3;
            else if (pr.source_name == "newIsIs2Is3")
                pr.source_id = src_newIsIs2Is3;
            else if (pr.source_name == "newIaIa2Ia3")
                pr.source_id = src_newIaIa2Ia3;
            else if (pr.source_name == "newVa1")
                pr.source_id = src_newVa1;
            else if (pr.source_name == "newVb1")
                pr.source_id = src_newVb1;
            else if (pr.source_name == "HospAdm")
                pr.source_id = src_HospAdm;
            else if (pr.source_name == "newEL123")
                pr.source_id = src_newEL123;
            else
                throw logic_error("Unrecognized process source name " + pr.source_name);
        }
        else
        {
            pr.source_id = (unsigned int)(sn - pc_names.begin());
        }
    }
}


bool PopulationParameters::Set(Parameters* parent, string& name, Rcpp::RObject& value)
{
    if (name == "contact" || name == "contact_mult" || name == "contact_lowerto" || name == "matrices")
        needs_recalc = true;

    if (false) {}
    _CheckSet(dE)
    _CheckSet(dIp)
    _CheckSet(dIa)
    _CheckSet(dIs)
    _CheckSet(dE2)
    _CheckSet(dIp2)
    _CheckSet(dIa2)
    _CheckSet(dIs2)
    _CheckSet(dE3)
    _CheckSet(dIp3)
    _CheckSet(dIa3)
    _CheckSet(dIs3)
    _CheckSet(dVa1)
    _CheckSet(dVb1)
    _CheckSet(dVa2)
    _CheckSet(dVb2)
    _CheckSet(size)
    _CheckSet(imm0)
    _CheckSet(matrices)
    _CheckSet(contact)
    _CheckSet(contact_mult)
    _CheckSet(contact_lowerto)
    _CheckSet(u)
    _CheckSet(u2)
    _CheckSet(u3)
    _CheckSet(fIp)
    _CheckSet(fIa)
    _CheckSet(fIs)
    _CheckSet(y)
    _CheckSet(y2)
    _CheckSet(y3)
    _CheckSet(omega)
    _CheckSet(tau)
    _CheckSet(pi_r)
    _CheckSet(pi_r2)
    _CheckSet(pi_r3)
    _CheckSet(pi2_r)
    _CheckSet(pi2_r2)
    _CheckSet(pi2_r3)
    _CheckSet(pi3_r)
    _CheckSet(pi3_r2)
    _CheckSet(pi3_r3)
    _CheckSet(wn)
    _CheckSet(wn2)
    _CheckSet(wn3)
    _CheckSet(va1)
    _CheckSet(wva1)
    _CheckSet(ei_va1)
    _CheckSet(ei2_va1)
    _CheckSet(ei3_va1)
    _CheckSet(ed_va1i)
    _CheckSet(ed_va1i2)
    _CheckSet(ed_va1i3)
    _CheckSet(wva2)
    _CheckSet(ei_va2)
    _CheckSet(ei2_va2)
    _CheckSet(ei3_va2)
    _CheckSet(ed_va2i)
    _CheckSet(ed_va2i2)
    _CheckSet(ed_va2i3)
    _CheckSet(wva3)
    _CheckSet(ei_va3)
    _CheckSet(ei2_va3)
    _CheckSet(ei3_va3)
    _CheckSet(ed_va3i)
    _CheckSet(ed_va3i2)
    _CheckSet(ed_va3i3)
    _CheckSet(p_boost_va2)
    _CheckSet(vb1)
    _CheckSet(wvb1)
    _CheckSet(ei_vb1)
    _CheckSet(ei2_vb1)
    _CheckSet(ei3_vb1)
    _CheckSet(ed_vb1i)
    _CheckSet(ed_vb1i2)
    _CheckSet(ed_vb1i3)
    _CheckSet(wvb2)
    _CheckSet(ei_vb2)
    _CheckSet(ei2_vb2)
    _CheckSet(ei3_vb2)
    _CheckSet(ed_vb2i)
    _CheckSet(ed_vb2i2)
    _CheckSet(ed_vb2i3)
    _CheckSet(wvb3)
    _CheckSet(ei_vb3)
    _CheckSet(ei2_vb3)
    _CheckSet(ei3_vb3)
    _CheckSet(ed_vb3i)
    _CheckSet(ed_vb3i2)
    _CheckSet(ed_vb3i3)
    _CheckSet(p_boost_vb2)
    _CheckSet(A)
    _CheckSet(B)
    _CheckSet(D)
    _CheckSet(season_A)
    _CheckSet(season_T)
    _CheckSet(season_phi)
    _CheckSet(ifr1)
    _CheckSet(ihr1)
    _CheckSet(iir1)
    _CheckSet(ifr2)
    _CheckSet(ihr2)
    _CheckSet(iir2)
    _CheckSet(ifr3)
    _CheckSet(ihr3)
    _CheckSet(iir3)
    _CheckSet(pd_ri)
    _CheckSet(pd_ri2)
    _CheckSet(pd_ri3)
    _CheckSet(pd_r2i)
    _CheckSet(pd_r2i2)
    _CheckSet(pd_r2i3)
    _CheckSet(pd_r3i)
    _CheckSet(pd_r3i2)
    _CheckSet(pd_r3i3)
    _CheckSet(ph_rd)
    _CheckSet(ph_rd2)
    _CheckSet(ph_rd3)
    _CheckSet(ph_r2d)
    _CheckSet(ph_r2d2)
    _CheckSet(ph_r2d3)
    _CheckSet(ph_r3d)
    _CheckSet(ph_r3d2)
    _CheckSet(ph_r3d3)
    _CheckSet(pm_rd)
    _CheckSet(pm_rd2)
    _CheckSet(pm_rd3)
    _CheckSet(pm_r2d)
    _CheckSet(pm_r2d2)
    _CheckSet(pm_r2d3)
    _CheckSet(pm_r3d)
    _CheckSet(pm_r3d2)
    _CheckSet(pm_r3d3)
    _CheckSet(pt_ri)
    _CheckSet(pt_ri2)
    _CheckSet(pt_ri3)
    _CheckSet(pt_r2i)
    _CheckSet(pt_r2i2)
    _CheckSet(pt_r2i3)
    _CheckSet(pt_r3i)
    _CheckSet(pt_r3i2)
    _CheckSet(pt_r3i3)
    _CheckSet(eh_va1d)
    _CheckSet(eh_va1d2)
    _CheckSet(eh_va1d3)
    _CheckSet(eh_va2d)
    _CheckSet(eh_va2d2)
    _CheckSet(eh_va2d3)
    _CheckSet(eh_va3d)
    _CheckSet(eh_va3d2)
    _CheckSet(eh_va3d3)
    _CheckSet(eh_vb1d)
    _CheckSet(eh_vb1d2)
    _CheckSet(eh_vb1d3)
    _CheckSet(eh_vb2d)
    _CheckSet(eh_vb2d2)
    _CheckSet(eh_vb2d3)
    _CheckSet(eh_vb3d)
    _CheckSet(eh_vb3d2)
    _CheckSet(eh_vb3d3)
    _CheckSet(em_va1d)
    _CheckSet(em_va1d2)
    _CheckSet(em_va1d3)
    _CheckSet(em_va2d)
    _CheckSet(em_va2d2)
    _CheckSet(em_va2d3)
    _CheckSet(em_va3d)
    _CheckSet(em_va3d2)
    _CheckSet(em_va3d3)
    _CheckSet(em_vb1d)
    _CheckSet(em_vb1d2)
    _CheckSet(em_vb1d3)
    _CheckSet(em_vb2d)
    _CheckSet(em_vb2d2)
    _CheckSet(em_vb2d3)
    _CheckSet(em_vb3d)
    _CheckSet(em_vb3d2)
    _CheckSet(em_vb3d3)
    _CheckSet(et_va1i)
    _CheckSet(et_va1i2)
    _CheckSet(et_va1i3)
    _CheckSet(et_va2i)
    _CheckSet(et_va2i2)
    _CheckSet(et_va2i3)
    _CheckSet(et_va3i)
    _CheckSet(et_va3i2)
    _CheckSet(et_va3i3)
    _CheckSet(et_vb1i)
    _CheckSet(et_vb1i2)
    _CheckSet(et_vb1i3)
    _CheckSet(et_vb2i)
    _CheckSet(et_vb2i2)
    _CheckSet(et_vb2i3)
    _CheckSet(et_vb3i)
    _CheckSet(et_vb3i2)
    _CheckSet(et_vb3i3)
    _CheckSet(dDeath)
    _CheckSet(dHosp)
    _CheckSet(lHosp)
    _CheckSet(dICU)
    _CheckSet(lICU)
    _CheckSetParent(parent, travel)
    else
    {
        throw logic_error("Unrecognised parameter " + name + ".");
    }

    return true;
}



void ParamSet(double& variable, vector<double>& value)
{
    variable = value[0];
}

void ParamSet(Discrete& variable, vector<double>& value)
{
    variable = value;
}

void ParamSet(vector<double>& variable, vector<double>& value)
{
    variable = value;
}

bool PopulationParameters::Set(Parameters* parent, string& name, vector<double>& value)
{
    if (name == "contact" || name == "contact_mult" || name == "contact_lowerto" || name == "matrices")
        needs_recalc = true;

    if (false) {}
    _CheckSet(dE)
    _CheckSet(dIp)
    _CheckSet(dIa)
    _CheckSet(dIs)
    _CheckSet(dE2)
    _CheckSet(dIp2)
    _CheckSet(dIa2)
    _CheckSet(dIs2)
    _CheckSet(dE3)
    _CheckSet(dIp3)
    _CheckSet(dIa3)
    _CheckSet(dIs3)
    _CheckSet(dVa1)
    _CheckSet(dVb1)
    _CheckSet(dVa2)
    _CheckSet(dVb2)
    _CheckSet(size)
    _CheckSet(imm0)
    //_CheckSet(matrices)
    _CheckSet(contact)
    _CheckSet(contact_mult)
    _CheckSet(contact_lowerto)
    _CheckSet(u)
    _CheckSet(u2)
    _CheckSet(u3)
    _CheckSet(fIp)
    _CheckSet(fIa)
    _CheckSet(fIs)
    _CheckSet(y)
    _CheckSet(y2)
    _CheckSet(y3)
    _CheckSet(omega)
    _CheckSet(tau)
    _CheckSet(pi_r)
    _CheckSet(pi_r2)
    _CheckSet(pi_r3)
    _CheckSet(pi2_r)
    _CheckSet(pi2_r2)
    _CheckSet(pi2_r3)
    _CheckSet(pi3_r)
    _CheckSet(pi3_r2)
    _CheckSet(pi3_r3)
    _CheckSet(wn)
    _CheckSet(wn2)
    _CheckSet(wn3)
    _CheckSet(va1)
    _CheckSet(wva1)
    _CheckSet(ei_va1)
    _CheckSet(ei2_va1)
    _CheckSet(ei3_va1)
    _CheckSet(ed_va1i)
    _CheckSet(ed_va1i2)
    _CheckSet(ed_va1i3)
    _CheckSet(wva2)
    _CheckSet(ei_va2)
    _CheckSet(ei2_va2)
    _CheckSet(ei3_va2)
    _CheckSet(ed_va2i)
    _CheckSet(ed_va2i2)
    _CheckSet(ed_va2i3)
    _CheckSet(wva3)
    _CheckSet(ei_va3)
    _CheckSet(ei2_va3)
    _CheckSet(ei3_va3)
    _CheckSet(ed_va3i)
    _CheckSet(ed_va3i2)
    _CheckSet(ed_va3i3)
    _CheckSet(p_boost_va2)
    _CheckSet(vb1)
    _CheckSet(wvb1)
    _CheckSet(ei_vb1)
    _CheckSet(ei2_vb1)
    _CheckSet(ei3_vb1)
    _CheckSet(ed_vb1i)
    _CheckSet(ed_vb1i2)
    _CheckSet(ed_vb1i3)
    _CheckSet(wvb2)
    _CheckSet(ei_vb2)
    _CheckSet(ei2_vb2)
    _CheckSet(ei3_vb2)
    _CheckSet(ed_vb2i)
    _CheckSet(ed_vb2i2)
    _CheckSet(ed_vb2i3)
    _CheckSet(wvb3)
    _CheckSet(ei_vb3)
    _CheckSet(ei2_vb3)
    _CheckSet(ei3_vb3)
    _CheckSet(ed_vb3i)
    _CheckSet(ed_vb3i2)
    _CheckSet(ed_vb3i3)
    _CheckSet(p_boost_vb2)
    _CheckSet(A)
    _CheckSet(B)
    _CheckSet(D)
    _CheckSet(season_A)
    _CheckSet(season_T)
    _CheckSet(season_phi)
    _CheckSet(ifr1)
    _CheckSet(ihr1)
    _CheckSet(iir1)
    _CheckSet(ifr2)
    _CheckSet(ihr2)
    _CheckSet(iir2)
    _CheckSet(ifr3)
    _CheckSet(ihr3)
    _CheckSet(iir3)
    _CheckSet(pd_ri)
    _CheckSet(pd_ri2)
    _CheckSet(pd_ri3)
    _CheckSet(pd_r2i)
    _CheckSet(pd_r2i2)
    _CheckSet(pd_r2i3)
    _CheckSet(pd_r3i)
    _CheckSet(pd_r3i2)
    _CheckSet(pd_r3i3)
    _CheckSet(ph_rd)
    _CheckSet(ph_rd2)
    _CheckSet(ph_rd3)
    _CheckSet(ph_r2d)
    _CheckSet(ph_r2d2)
    _CheckSet(ph_r2d3)
    _CheckSet(ph_r3d)
    _CheckSet(ph_r3d2)
    _CheckSet(ph_r3d3)
    _CheckSet(pm_rd)
    _CheckSet(pm_rd2)
    _CheckSet(pm_rd3)
    _CheckSet(pm_r2d)
    _CheckSet(pm_r2d2)
    _CheckSet(pm_r2d3)
    _CheckSet(pm_r3d)
    _CheckSet(pm_r3d2)
    _CheckSet(pm_r3d3)
    _CheckSet(pt_ri)
    _CheckSet(pt_ri2)
    _CheckSet(pt_ri3)
    _CheckSet(pt_r2i)
    _CheckSet(pt_r2i2)
    _CheckSet(pt_r2i3)
    _CheckSet(pt_r3i)
    _CheckSet(pt_r3i2)
    _CheckSet(pt_r3i3)
    _CheckSet(eh_va1d)
    _CheckSet(eh_va1d2)
    _CheckSet(eh_va1d3)
    _CheckSet(eh_va2d)
    _CheckSet(eh_va2d2)
    _CheckSet(eh_va2d3)
    _CheckSet(eh_va3d)
    _CheckSet(eh_va3d2)
    _CheckSet(eh_va3d3)
    _CheckSet(eh_vb1d)
    _CheckSet(eh_vb1d2)
    _CheckSet(eh_vb1d3)
    _CheckSet(eh_vb2d)
    _CheckSet(eh_vb2d2)
    _CheckSet(eh_vb2d3)
    _CheckSet(eh_vb3d)
    _CheckSet(eh_vb3d2)
    _CheckSet(eh_vb3d3)
    _CheckSet(em_va1d)
    _CheckSet(em_va1d2)
    _CheckSet(em_va1d3)
    _CheckSet(em_va2d)
    _CheckSet(em_va2d2)
    _CheckSet(em_va2d3)
    _CheckSet(em_va3d)
    _CheckSet(em_va3d2)
    _CheckSet(em_va3d3)
    _CheckSet(em_vb1d)
    _CheckSet(em_vb1d2)
    _CheckSet(em_vb1d3)
    _CheckSet(em_vb2d)
    _CheckSet(em_vb2d2)
    _CheckSet(em_vb2d3)
    _CheckSet(em_vb3d)
    _CheckSet(em_vb3d2)
    _CheckSet(em_vb3d3)
    _CheckSet(et_va1i)
    _CheckSet(et_va1i2)
    _CheckSet(et_va1i3)
    _CheckSet(et_va2i)
    _CheckSet(et_va2i2)
    _CheckSet(et_va2i3)
    _CheckSet(et_va3i)
    _CheckSet(et_va3i2)
    _CheckSet(et_va3i3)
    _CheckSet(et_vb1i)
    _CheckSet(et_vb1i2)
    _CheckSet(et_vb1i3)
    _CheckSet(et_vb2i)
    _CheckSet(et_vb2i2)
    _CheckSet(et_vb2i3)
    _CheckSet(et_vb3i)
    _CheckSet(et_vb3i2)
    _CheckSet(et_vb3i3)
    _CheckSet(dDeath)
    _CheckSet(dHosp)
    _CheckSet(lHosp)
    _CheckSet(dICU)
    _CheckSet(lICU)
    //_CheckSetParent(parent, travel)
    else
    {
        throw logic_error("Unrecognised parameter " + name + ".");
    }

    return true;
}

void Parameters::FilterForRun(unsigned int r)
{
    for (auto i = changes.ch.begin(); i != changes.ch.end(); )
    {
        if (i->runs.empty() || find(i->runs.begin(), i->runs.end(), r) != i->runs.end())
            ++i;
        else
            i = changes.ch.erase(i);
    }
}


// Helpers to set parameters
#define ParamAssign(t, v)               if (list.containsElementNamed(#v)) P.v = Rcpp::as<t>(list[#v]);
#define ParamMatrixAssign(v)            if (list.containsElementNamed(#v)) SetMatrix(P.v, Rcpp::as<Rcpp::NumericMatrix>(list[#v]));
#define ParamPopAssign(t, v, i)         if (popi.containsElementNamed(#v)) P.pop[i].v = Rcpp::as<t>(popi[#v]);
#define ParamPopMatrixAssign(v, i)      if (popi.containsElementNamed(#v)) SetMatrix(P.pop[i].v, Rcpp::as<Rcpp::NumericMatrix>(popi[#v]));

void SetMatrix(Matrix& mat, const Rcpp::NumericMatrix& rhs)
{
    mat = Matrix(0.0, rhs.nrow(), rhs.ncol());
    for (int r = 0; r < rhs.nrow(); ++r)
        for (int c = 0; c < rhs.ncol(); ++c)
            mat(r, c) = rhs(r, c);
}

void SetParameters(Parameters& P, Rcpp::List list, Randomizer& Rand)
{
    P.vax_limit = 1.0;
    P.adjustment = 1.0;
    P.extra_vb1 = vector<double>(24, 0.0);

    // TODO clean up and use namespace Rcpp
    // TODO use the ParamSet functions above for all of these parameter assignments
    ParamAssign(string, model);
    ParamAssign(double, time_step);
    ParamAssign(double, time0);
    ParamAssign(double, time1);
    ParamAssign(unsigned int, report_every);
    ParamAssign(bool, fast_multinomial);
    ParamAssign(bool, deterministic);

    if (P.report_every != 1/P.time_step)
        throw("report_every must be the reciprocal of time_step.");

    if (list.containsElementNamed("pop"))
    {
        Rcpp::List populations = Rcpp::as<Rcpp::List>(list["pop"]);
        unsigned int np = populations.size();
        P.pop.assign(np, PopulationParameters());

        for (unsigned int i = 0; i < np; ++i)
        {
            Rcpp::List popi = Rcpp::as<Rcpp::List>(populations[i]);

            ParamPopAssign(vector<double>, dE, i);
            ParamPopAssign(vector<double>, dIp, i);
            ParamPopAssign(vector<double>, dIa, i);
            ParamPopAssign(vector<double>, dIs, i);
            ParamPopAssign(vector<double>, dE2, i);
            ParamPopAssign(vector<double>, dIp2, i);
            ParamPopAssign(vector<double>, dIa2, i);
            ParamPopAssign(vector<double>, dIs2, i);
            ParamPopAssign(vector<double>, dE3, i);
            ParamPopAssign(vector<double>, dIp3, i);
            ParamPopAssign(vector<double>, dIa3, i);
            ParamPopAssign(vector<double>, dIs3, i);
            ParamPopAssign(vector<double>, dVa1, i);
            ParamPopAssign(vector<double>, dVb1, i);
            ParamPopAssign(vector<double>, dVa2, i);
            ParamPopAssign(vector<double>, dVb2, i);

            if (P.fast_multinomial)
            {
                P.pop[i].dE.mn_approx.Set(P, Rand, P.pop[i].dE.weights);
                P.pop[i].dIp.mn_approx.Set(P, Rand, P.pop[i].dIp.weights);
                P.pop[i].dIa.mn_approx.Set(P, Rand, P.pop[i].dIa.weights);
                P.pop[i].dIs.mn_approx.Set(P, Rand, P.pop[i].dIs.weights);
                P.pop[i].dE2.mn_approx.Set(P, Rand, P.pop[i].dE2.weights);
                P.pop[i].dIp2.mn_approx.Set(P, Rand, P.pop[i].dIp2.weights);
                P.pop[i].dIa2.mn_approx.Set(P, Rand, P.pop[i].dIa2.weights);
                P.pop[i].dIs2.mn_approx.Set(P, Rand, P.pop[i].dIs2.weights);
                P.pop[i].dE3.mn_approx.Set(P, Rand, P.pop[i].dE3.weights);
                P.pop[i].dIp3.mn_approx.Set(P, Rand, P.pop[i].dIp3.weights);
                P.pop[i].dIa3.mn_approx.Set(P, Rand, P.pop[i].dIa3.weights);
                P.pop[i].dIs3.mn_approx.Set(P, Rand, P.pop[i].dIs3.weights);
                P.pop[i].dVa1.mn_approx.Set(P, Rand, P.pop[i].dVa1.weights);
                P.pop[i].dVb1.mn_approx.Set(P, Rand, P.pop[i].dVb1.weights);
                P.pop[i].dVa2.mn_approx.Set(P, Rand, P.pop[i].dVa2.weights);
                P.pop[i].dVb2.mn_approx.Set(P, Rand, P.pop[i].dVb2.weights);
            }

            ParamPopAssign(vector<double>, size, i);
            ParamPopAssign(vector<double>, imm0, i);

            if (popi.containsElementNamed("matrices"))
            {
                Rcpp::List ml = Rcpp::as<Rcpp::List>(popi["matrices"]);
                for (unsigned int m = 0; m < ml.size(); ++m)
                {
                    Matrix mat;
                    SetMatrix(mat, Rcpp::as<Rcpp::NumericMatrix>(ml[m]));
                    P.pop[i].matrices.push_back(mat);
                }
            }

            ParamPopAssign(vector<double>, contact, i);
            ParamPopAssign(vector<double>, contact_mult, i);
            ParamPopAssign(vector<double>, contact_lowerto, i);
            ParamPopAssign(vector<double>, u, i);
            ParamPopAssign(vector<double>, u2, i);
            ParamPopAssign(vector<double>, u3, i);
            ParamPopAssign(vector<double>, fIp, i);
            ParamPopAssign(vector<double>, fIa, i);
            ParamPopAssign(vector<double>, fIs, i);
            ParamPopAssign(vector<double>, y, i);
            ParamPopAssign(vector<double>, y2, i);
            ParamPopAssign(vector<double>, y3, i);
            ParamPopAssign(vector<double>, omega, i);
            ParamPopAssign(vector<double>, tau, i);
            ParamPopAssign(vector<double>, pi_r, i);
            ParamPopAssign(vector<double>, pi_r2, i);
            ParamPopAssign(vector<double>, pi_r3, i);
            ParamPopAssign(vector<double>, pi2_r, i);
            ParamPopAssign(vector<double>, pi2_r2, i);
            ParamPopAssign(vector<double>, pi2_r3, i);
            ParamPopAssign(vector<double>, pi3_r, i);
            ParamPopAssign(vector<double>, pi3_r2, i);
            ParamPopAssign(vector<double>, pi3_r3, i);
            ParamPopAssign(vector<double>, wn, i);
            ParamPopAssign(vector<double>, wn2, i);
            ParamPopAssign(vector<double>, wn3, i);
            ParamPopAssign(vector<double>, va1, i);
            ParamPopAssign(vector<double>, wva1, i);
            ParamPopAssign(vector<double>, ei_va1, i);
            ParamPopAssign(vector<double>, ei2_va1, i);
            ParamPopAssign(vector<double>, ei3_va1, i);
            ParamPopAssign(vector<double>, ed_va1i, i);
            ParamPopAssign(vector<double>, ed_va1i2, i);
            ParamPopAssign(vector<double>, ed_va1i3, i);
            ParamPopAssign(vector<double>, wva2, i);
            ParamPopAssign(vector<double>, ei_va2, i);
            ParamPopAssign(vector<double>, ei2_va2, i);
            ParamPopAssign(vector<double>, ei3_va2, i);
            ParamPopAssign(vector<double>, ed_va2i, i);
            ParamPopAssign(vector<double>, ed_va2i2, i);
            ParamPopAssign(vector<double>, ed_va2i3, i);
            ParamPopAssign(vector<double>, wva3, i);
            ParamPopAssign(vector<double>, ei_va3, i);
            ParamPopAssign(vector<double>, ei2_va3, i);
            ParamPopAssign(vector<double>, ei3_va3, i);
            ParamPopAssign(vector<double>, ed_va3i, i);
            ParamPopAssign(vector<double>, ed_va3i2, i);
            ParamPopAssign(vector<double>, ed_va3i3, i);
            ParamPopAssign(vector<double>, p_boost_va2, i);
            ParamPopAssign(vector<double>, vb1, i);
            ParamPopAssign(vector<double>, wvb1, i);
            ParamPopAssign(vector<double>, ei_vb1, i);
            ParamPopAssign(vector<double>, ei2_vb1, i);
            ParamPopAssign(vector<double>, ei3_vb1, i);
            ParamPopAssign(vector<double>, ed_vb1i, i);
            ParamPopAssign(vector<double>, ed_vb1i2, i);
            ParamPopAssign(vector<double>, ed_vb1i3, i);
            ParamPopAssign(vector<double>, wvb2, i);
            ParamPopAssign(vector<double>, ei_vb2, i);
            ParamPopAssign(vector<double>, ei2_vb2, i);
            ParamPopAssign(vector<double>, ei3_vb2, i);
            ParamPopAssign(vector<double>, ed_vb2i, i);
            ParamPopAssign(vector<double>, ed_vb2i2, i);
            ParamPopAssign(vector<double>, ed_vb2i3, i);
            ParamPopAssign(vector<double>, wvb3, i);
            ParamPopAssign(vector<double>, ei_vb3, i);
            ParamPopAssign(vector<double>, ei2_vb3, i);
            ParamPopAssign(vector<double>, ei3_vb3, i);
            ParamPopAssign(vector<double>, ed_vb3i, i);
            ParamPopAssign(vector<double>, ed_vb3i2, i);
            ParamPopAssign(vector<double>, ed_vb3i3, i);
            ParamPopAssign(vector<double>, p_boost_vb2, i);
            ParamPopAssign(vector<double>, A, i);
            ParamPopAssign(vector<double>, B, i);
            ParamPopAssign(vector<double>, D, i);

            ParamPopAssign(vector<double>, season_A, i);
            ParamPopAssign(vector<double>, season_T, i);
            ParamPopAssign(vector<double>, season_phi, i);

            ParamPopAssign(vector<double>, ifr1, i);
            ParamPopAssign(vector<double>, ihr1, i);
            ParamPopAssign(vector<double>, iir1, i);
            ParamPopAssign(vector<double>, ifr2, i);
            ParamPopAssign(vector<double>, ihr2, i);
            ParamPopAssign(vector<double>, iir2, i);
            ParamPopAssign(vector<double>, ifr3, i);
            ParamPopAssign(vector<double>, ihr3, i);
            ParamPopAssign(vector<double>, iir3, i);
            ParamPopAssign(vector<double>, pd_ri, i);
            ParamPopAssign(vector<double>, pd_ri2, i);
            ParamPopAssign(vector<double>, pd_ri3, i);
            ParamPopAssign(vector<double>, pd_r2i, i);
            ParamPopAssign(vector<double>, pd_r2i2, i);
            ParamPopAssign(vector<double>, pd_r2i3, i);
            ParamPopAssign(vector<double>, pd_r3i, i);
            ParamPopAssign(vector<double>, pd_r3i2, i);
            ParamPopAssign(vector<double>, pd_r3i3, i);
            ParamPopAssign(vector<double>, ph_rd, i);
            ParamPopAssign(vector<double>, ph_rd2, i);
            ParamPopAssign(vector<double>, ph_rd3, i);
            ParamPopAssign(vector<double>, ph_r2d, i);
            ParamPopAssign(vector<double>, ph_r2d2, i);
            ParamPopAssign(vector<double>, ph_r2d3, i);
            ParamPopAssign(vector<double>, ph_r3d, i);
            ParamPopAssign(vector<double>, ph_r3d2, i);
            ParamPopAssign(vector<double>, ph_r3d3, i);
            ParamPopAssign(vector<double>, pm_rd, i);
            ParamPopAssign(vector<double>, pm_rd2, i);
            ParamPopAssign(vector<double>, pm_rd3, i);
            ParamPopAssign(vector<double>, pm_r2d, i);
            ParamPopAssign(vector<double>, pm_r2d2, i);
            ParamPopAssign(vector<double>, pm_r2d3, i);
            ParamPopAssign(vector<double>, pm_r3d, i);
            ParamPopAssign(vector<double>, pm_r3d2, i);
            ParamPopAssign(vector<double>, pm_r3d3, i);
            ParamPopAssign(vector<double>, pt_ri, i);
            ParamPopAssign(vector<double>, pt_ri2, i);
            ParamPopAssign(vector<double>, pt_ri3, i);
            ParamPopAssign(vector<double>, pt_r2i, i);
            ParamPopAssign(vector<double>, pt_r2i2, i);
            ParamPopAssign(vector<double>, pt_r2i3, i);
            ParamPopAssign(vector<double>, pt_r3i, i);
            ParamPopAssign(vector<double>, pt_r3i2, i);
            ParamPopAssign(vector<double>, pt_r3i3, i);
            ParamPopAssign(vector<double>, eh_va1d, i);
            ParamPopAssign(vector<double>, eh_va1d2, i);
            ParamPopAssign(vector<double>, eh_va1d3, i);
            ParamPopAssign(vector<double>, eh_va2d, i);
            ParamPopAssign(vector<double>, eh_va2d2, i);
            ParamPopAssign(vector<double>, eh_va2d3, i);
            ParamPopAssign(vector<double>, eh_va3d, i);
            ParamPopAssign(vector<double>, eh_va3d2, i);
            ParamPopAssign(vector<double>, eh_va3d3, i);
            ParamPopAssign(vector<double>, eh_vb1d, i);
            ParamPopAssign(vector<double>, eh_vb1d2, i);
            ParamPopAssign(vector<double>, eh_vb1d3, i);
            ParamPopAssign(vector<double>, eh_vb2d, i);
            ParamPopAssign(vector<double>, eh_vb2d2, i);
            ParamPopAssign(vector<double>, eh_vb2d3, i);
            ParamPopAssign(vector<double>, eh_vb3d, i);
            ParamPopAssign(vector<double>, eh_vb3d2, i);
            ParamPopAssign(vector<double>, eh_vb3d3, i);
            ParamPopAssign(vector<double>, em_va1d, i);
            ParamPopAssign(vector<double>, em_va1d2, i);
            ParamPopAssign(vector<double>, em_va1d3, i);
            ParamPopAssign(vector<double>, em_va2d, i);
            ParamPopAssign(vector<double>, em_va2d2, i);
            ParamPopAssign(vector<double>, em_va2d3, i);
            ParamPopAssign(vector<double>, em_va3d, i);
            ParamPopAssign(vector<double>, em_va3d2, i);
            ParamPopAssign(vector<double>, em_va3d3, i);
            ParamPopAssign(vector<double>, em_vb1d, i);
            ParamPopAssign(vector<double>, em_vb1d2, i);
            ParamPopAssign(vector<double>, em_vb1d3, i);
            ParamPopAssign(vector<double>, em_vb2d, i);
            ParamPopAssign(vector<double>, em_vb2d2, i);
            ParamPopAssign(vector<double>, em_vb2d3, i);
            ParamPopAssign(vector<double>, em_vb3d, i);
            ParamPopAssign(vector<double>, em_vb3d2, i);
            ParamPopAssign(vector<double>, em_vb3d3, i);
            ParamPopAssign(vector<double>, et_va1i, i);
            ParamPopAssign(vector<double>, et_va1i2, i);
            ParamPopAssign(vector<double>, et_va1i3, i);
            ParamPopAssign(vector<double>, et_va2i, i);
            ParamPopAssign(vector<double>, et_va2i2, i);
            ParamPopAssign(vector<double>, et_va2i3, i);
            ParamPopAssign(vector<double>, et_va3i, i);
            ParamPopAssign(vector<double>, et_va3i2, i);
            ParamPopAssign(vector<double>, et_va3i3, i);
            ParamPopAssign(vector<double>, et_vb1i, i);
            ParamPopAssign(vector<double>, et_vb1i2, i);
            ParamPopAssign(vector<double>, et_vb1i3, i);
            ParamPopAssign(vector<double>, et_vb2i, i);
            ParamPopAssign(vector<double>, et_vb2i2, i);
            ParamPopAssign(vector<double>, et_vb2i3, i);
            ParamPopAssign(vector<double>, et_vb3i, i);
            ParamPopAssign(vector<double>, et_vb3i2, i);
            ParamPopAssign(vector<double>, et_vb3i3, i);
            ParamPopAssign(vector<double>, dDeath, i);
            ParamPopAssign(vector<double>, dHosp, i);
            ParamPopAssign(vector<double>, lHosp, i);
            ParamPopAssign(vector<double>, dICU, i);
            ParamPopAssign(vector<double>, lICU, i);

            ParamPopAssign(vector<double>, seed_times, i);
            ParamPopAssign(vector<double>, seed_times2, i);
            ParamPopAssign(vector<double>, seed_times3, i);
            ParamPopAssign(vector<double>, dist_seed_ages, i);

            // Names
            ParamPopAssign(string, name, i);
            ParamPopAssign(vector<string>, group_names, i);

            // Calculate
            P.pop[i].Recalculate();
        }
    }

    // Read in processes
    if (list.containsElementNamed("processes"))
    {
        Rcpp::RObject r_processes = Rcpp::as<Rcpp::RObject>(list["processes"]);
        ParamSet(P.processes, r_processes);
    }

    ParamMatrixAssign(travel);

    // Read in schedule (change set)
    if (list.containsElementNamed("schedule"))
    {
        Rcpp::List schedule = Rcpp::as<Rcpp::List>(list["schedule"]);
        for (unsigned int j = 0; j < schedule.size(); ++j)
        {
            Rcpp::List sched = Rcpp::as<Rcpp::List>(schedule[j]);

            string param_name = Rcpp::as<string>(sched["parameter"]);
            vector<unsigned int> pops = Rcpp::as<vector<unsigned int>>(sched["pops"]);
            vector<unsigned int> runs;
            if (sched.containsElementNamed("runs"))
                runs = Rcpp::as<vector<unsigned int>>(sched["runs"]);

            string mode_name = Rcpp::as<string>(sched["mode"]);
            Change::Mode mode = Change::Assign;
            if (mode_name == "assign")          mode = Change::Assign;
            else if (mode_name == "add")        mode = Change::Add;
            else if (mode_name == "multiply")   mode = Change::Multiply;
            else if (mode_name == "lowerto")    mode = Change::LowerTo;
            else if (mode_name == "raiseto")    mode = Change::RaiseTo;
            else if (mode_name == "bypass")     mode = Change::Bypass;
            else                                throw std::logic_error("Unrecognized change mode " + mode_name + ".");

            vector<double> times = Rcpp::as<vector<double>>(sched["times"]);
            Rcpp::List values_list = Rcpp::as<Rcpp::List>(sched["values"]);
            vector<vector<double>> values;

            for (unsigned int k = 0; k < values_list.size(); ++k)
                values.push_back(Rcpp::as<vector<double>>(values_list[k]));

            P.changes.ch.push_back(Change(P, pops, runs, param_name, mode, times, values));
        }
    }
}
