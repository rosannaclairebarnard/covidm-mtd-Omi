// changes.cpp

#include "changes.h"
#include "parameters.h"
#include "helper.h"
#include <stdexcept>

#define PVector(x)   x
#define PDiscrete(x) x.weights

#define ParamCapture(param, get_vec) \
    else if (param_name == #param) { \
        if (pops.empty()) { \
            for (auto& pp : P.pop) { \
                param_ptr.push_back(&get_vec(pp.param)); \
                if (reset) param_orig.push_back(get_vec(pp.param)); \
            } \
        } else { \
            for (auto& pi : pops) { \
                if (pi >= P.pop.size()) \
                    throw logic_error("Referring to nonexistent population in Change::Capture."); \
                param_ptr.push_back(&get_vec(P.pop[pi].param)); \
                if (reset) param_orig.push_back(get_vec(P.pop[pi].param)); \
            } \
        } \
    }

// Construct a change impacting parameter pname in populations po of parameters P;
// apply value v with mode m at times t
Change::Change(Parameters& P, vector<unsigned int>& po, vector<unsigned int>& ru, string pname,
    Mode m, vector<double>& t, vector<vector<double>>& v)
 : mode(m), times(t), values(v), param_name(pname), current(-1), pops(po), runs(ru)
{
    if (times.size() != values.size())
        throw logic_error("Change: times and values must be same length.");
}

// Capture parameters to change
void Change::Capture(Parameters& P, bool reset)
{
    if (reset)
    {
        current = -1;
        param_ptr.clear();
        param_orig.clear();
    }
    else
    {
        param_ptr.clear();
    }

    if (false) {}
    ParamCapture(dE,              PDiscrete)
    ParamCapture(dIp,             PDiscrete)
    ParamCapture(dIa,             PDiscrete)
    ParamCapture(dIs,             PDiscrete)
    ParamCapture(dE2,             PDiscrete)
    ParamCapture(dIp2,            PDiscrete)
    ParamCapture(dIa2,            PDiscrete)
    ParamCapture(dIs2,            PDiscrete)
    ParamCapture(dE3,             PDiscrete)
    ParamCapture(dIp3,            PDiscrete)
    ParamCapture(dIa3,            PDiscrete)
    ParamCapture(dIs3,            PDiscrete)
    ParamCapture(dVa1,            PDiscrete)
    ParamCapture(dVb1,            PDiscrete)
    ParamCapture(dVa2,            PDiscrete)
    ParamCapture(dVb2,            PDiscrete)
    ParamCapture(contact,         PVector)
    ParamCapture(contact_mult,    PVector)
    ParamCapture(contact_lowerto, PVector)
    ParamCapture(u,               PVector)
    ParamCapture(u2,              PVector)
    ParamCapture(u3,              PVector)
    ParamCapture(fIp,             PVector)
    ParamCapture(fIa,             PVector)
    ParamCapture(fIs,             PVector)
    ParamCapture(y,               PVector)
    ParamCapture(y2,              PVector)
    ParamCapture(y3,              PVector)
    ParamCapture(omega,           PVector)
    ParamCapture(tau,             PVector)
    ParamCapture(pi_r,            PVector)
    ParamCapture(pi_r2,           PVector)
    ParamCapture(pi_r3,           PVector)
    ParamCapture(pi2_r,           PVector)
    ParamCapture(pi2_r2,          PVector)
    ParamCapture(pi2_r3,          PVector)
    ParamCapture(pi3_r,           PVector)
    ParamCapture(pi3_r2,          PVector)
    ParamCapture(pi3_r3,          PVector)
    ParamCapture(wn,              PVector)
    ParamCapture(wn2,             PVector)
    ParamCapture(wn3,             PVector)
    ParamCapture(va1,             PVector)
    ParamCapture(wva1,            PVector)
    ParamCapture(ei_va1,          PVector)
    ParamCapture(ei2_va1,         PVector)
    ParamCapture(ei3_va1,         PVector)
    ParamCapture(ed_va1i,         PVector)
    ParamCapture(ed_va1i2,        PVector)
    ParamCapture(ed_va1i3,        PVector)
    ParamCapture(wva2,            PVector)
    ParamCapture(ei_va2,          PVector)
    ParamCapture(ei2_va2,         PVector)
    ParamCapture(ei3_va2,         PVector)
    ParamCapture(ed_va2i,         PVector)
    ParamCapture(ed_va2i2,        PVector)
    ParamCapture(ed_va2i3,        PVector)
    ParamCapture(wva3,            PVector)
    ParamCapture(ei_va3,          PVector)
    ParamCapture(ei2_va3,         PVector)
    ParamCapture(ei3_va3,         PVector)
    ParamCapture(ed_va3i,         PVector)
    ParamCapture(ed_va3i2,        PVector)
    ParamCapture(ed_va3i3,        PVector)
    ParamCapture(p_boost_va2,     PVector)
    ParamCapture(vb1,             PVector)
    ParamCapture(wvb1,            PVector)
    ParamCapture(ei_vb1,          PVector)
    ParamCapture(ei2_vb1,         PVector)
    ParamCapture(ei3_vb1,         PVector)
    ParamCapture(ed_vb1i,         PVector)
    ParamCapture(ed_vb1i2,        PVector)
    ParamCapture(ed_vb1i3,        PVector)
    ParamCapture(wvb2,            PVector)
    ParamCapture(ei_vb2,          PVector)
    ParamCapture(ei2_vb2,         PVector)
    ParamCapture(ei3_vb2,         PVector)
    ParamCapture(ed_vb2i,         PVector)
    ParamCapture(ed_vb2i2,        PVector)
    ParamCapture(ed_vb2i3,        PVector)
    ParamCapture(wvb3,            PVector)
    ParamCapture(ei_vb3,          PVector)
    ParamCapture(ei2_vb3,         PVector)
    ParamCapture(ei3_vb3,         PVector)
    ParamCapture(ed_vb3i,         PVector)
    ParamCapture(ed_vb3i2,        PVector)
    ParamCapture(ed_vb3i3,        PVector)
    ParamCapture(p_boost_vb2,     PVector)
    ParamCapture(A,               PVector)
    ParamCapture(B,               PVector)
    ParamCapture(D,               PVector)
    ParamCapture(season_A,        PVector)
    ParamCapture(season_T,        PVector)
    ParamCapture(season_phi,      PVector)
    ParamCapture(ifr1,            PVector)
    ParamCapture(ihr1,            PVector)
    ParamCapture(iir1,            PVector)
    ParamCapture(ifr2,            PVector)
    ParamCapture(ihr2,            PVector)
    ParamCapture(iir2,            PVector)
    ParamCapture(ifr3,            PVector)
    ParamCapture(ihr3,            PVector)
    ParamCapture(iir3,            PVector)
    ParamCapture(pd_ri,           PVector)
    ParamCapture(pd_ri2,          PVector)
    ParamCapture(pd_ri3,          PVector)
    ParamCapture(pd_r2i,          PVector)
    ParamCapture(pd_r2i2,         PVector)
    ParamCapture(pd_r2i3,         PVector)
    ParamCapture(pd_r3i,          PVector)
    ParamCapture(pd_r3i2,         PVector)
    ParamCapture(pd_r3i3,         PVector)
    ParamCapture(ph_rd,           PVector)
    ParamCapture(ph_rd2,          PVector)
    ParamCapture(ph_rd3,          PVector)
    ParamCapture(ph_r2d,          PVector)
    ParamCapture(ph_r2d2,         PVector)
    ParamCapture(ph_r2d3,         PVector)
    ParamCapture(ph_r3d,          PVector)
    ParamCapture(ph_r3d2,         PVector)
    ParamCapture(ph_r3d3,         PVector)
    ParamCapture(pm_rd,           PVector)
    ParamCapture(pm_rd2,          PVector)
    ParamCapture(pm_rd3,          PVector)
    ParamCapture(pm_r2d,          PVector)
    ParamCapture(pm_r2d2,         PVector)
    ParamCapture(pm_r2d3,         PVector)
    ParamCapture(pm_r3d,          PVector)
    ParamCapture(pm_r3d2,         PVector)
    ParamCapture(pm_r3d3,         PVector)
    ParamCapture(pt_ri,           PVector)
    ParamCapture(pt_ri2,          PVector)
    ParamCapture(pt_ri3,          PVector)
    ParamCapture(pt_r2i,          PVector)
    ParamCapture(pt_r2i2,         PVector)
    ParamCapture(pt_r2i3,         PVector)
    ParamCapture(pt_r3i,          PVector)
    ParamCapture(pt_r3i2,         PVector)
    ParamCapture(pt_r3i3,         PVector)
    ParamCapture(eh_va1d,         PVector)
    ParamCapture(eh_va1d2,        PVector)
    ParamCapture(eh_va1d3,        PVector)
    ParamCapture(eh_va2d,         PVector)
    ParamCapture(eh_va2d2,        PVector)
    ParamCapture(eh_va2d3,        PVector)
    ParamCapture(eh_va3d,         PVector)
    ParamCapture(eh_va3d2,        PVector)
    ParamCapture(eh_va3d3,        PVector)
    ParamCapture(eh_vb1d,         PVector)
    ParamCapture(eh_vb1d2,        PVector)
    ParamCapture(eh_vb1d3,        PVector)
    ParamCapture(eh_vb2d,         PVector)
    ParamCapture(eh_vb2d2,        PVector)
    ParamCapture(eh_vb2d3,        PVector)
    ParamCapture(eh_vb3d,         PVector)
    ParamCapture(eh_vb3d2,        PVector)
    ParamCapture(eh_vb3d3,        PVector)
    ParamCapture(em_va1d,         PVector)
    ParamCapture(em_va1d2,        PVector)
    ParamCapture(em_va1d3,        PVector)
    ParamCapture(em_va2d,         PVector)
    ParamCapture(em_va2d2,        PVector)
    ParamCapture(em_va2d3,        PVector)
    ParamCapture(em_va3d,         PVector)
    ParamCapture(em_va3d2,        PVector)
    ParamCapture(em_va3d3,        PVector)
    ParamCapture(em_vb1d,         PVector)
    ParamCapture(em_vb1d2,        PVector)
    ParamCapture(em_vb1d3,        PVector)
    ParamCapture(em_vb2d,         PVector)
    ParamCapture(em_vb2d2,        PVector)
    ParamCapture(em_vb2d3,        PVector)
    ParamCapture(em_vb3d,         PVector)
    ParamCapture(em_vb3d2,        PVector)
    ParamCapture(em_vb3d3,        PVector)
    ParamCapture(et_va1i,         PVector)
    ParamCapture(et_va1i2,        PVector)
    ParamCapture(et_va1i3,        PVector)
    ParamCapture(et_va2i,         PVector)
    ParamCapture(et_va2i2,        PVector)
    ParamCapture(et_va2i3,        PVector)
    ParamCapture(et_va3i,         PVector)
    ParamCapture(et_va3i2,        PVector)
    ParamCapture(et_va3i3,        PVector)
    ParamCapture(et_vb1i,         PVector)
    ParamCapture(et_vb1i2,        PVector)
    ParamCapture(et_vb1i3,        PVector)
    ParamCapture(et_vb2i,         PVector)
    ParamCapture(et_vb2i2,        PVector)
    ParamCapture(et_vb2i3,        PVector)
    ParamCapture(et_vb3i,         PVector)
    ParamCapture(et_vb3i2,        PVector)
    ParamCapture(et_vb3i3,        PVector)
    ParamCapture(dDeath,          PDiscrete)
    ParamCapture(dHosp,           PDiscrete)
    ParamCapture(lHosp,           PDiscrete)
    ParamCapture(dICU,            PDiscrete)
    ParamCapture(lICU,            PDiscrete)

    else if (param_name == "travel") {
        param_ptr.push_back(&P.travel.x);
        param_orig.push_back(P.travel.x);
    }
    else
        throw logic_error("Unsupported parameter change for " + param_name);
}

// Return true if parameters will change at time t
bool Change::Update(double t)
{
    bool refresh = false;
    while (current + 1 < (int)times.size() && t >= times[current + 1])
    {
        ++current;
        refresh = true;
    }
    return refresh;
}

// Apply parameter changes
void Change::Apply(double t)
{
    auto op = [&](vector<double>* ptr) {
        if (ptr->size() == values[current].size() || values[current].empty())
        {
            if (!values[current].empty())
            {
                for (unsigned int j = 0; j < ptr->size(); ++j)
                {
                    if (values[current][j] == -999)
                        continue;

                    switch (mode)
                    {
                        case Assign:    (*ptr)[j] = values[current][j]; break;
                        case Add:       (*ptr)[j] += values[current][j]; break;
                        case Multiply:  (*ptr)[j] *= values[current][j]; break;
                        case LowerTo:   (*ptr)[j] = min((*ptr)[j], values[current][j]); break;
                        case RaiseTo:   (*ptr)[j] = max((*ptr)[j], values[current][j]); break;
                        case Bypass:    break;
                    }
                }
            }
        }
        else
        {
            cout << current << "\n";
            cout << ptr->size() << "\n";
            cout << values[current].size() << "\n";
            throw logic_error("Change::Apply: change value has different size from parameter.");
        }
    };

    if (current >= 0)
        for (auto ptr : param_ptr)
            op(ptr);
}

// Reset all linked parameters to their original values
void Change::Reset()
{
    for (unsigned int i = 0; i < param_ptr.size(); ++i)
        *param_ptr[i] = param_orig[i];
}

// Capture parameters to change
void ChangeSet::Capture(Parameters& P)
{
    for (auto& c : ch)
        c.Capture(P, true);
}

// Reset pointers without resetting changes
void ChangeSet::Recapture(Parameters& P)
{
    for (auto& c : ch)
        c.Capture(P, false);
}


// Apply any needed changes
void ChangeSet::Apply(Parameters& P, double t)
{
    bool refresh = false;

    for (auto& c : ch)
    {
        if (c.Update(t))
        {
            refresh = true;
        }
    }

    if (refresh)
    {
        for (auto& c : ch)
            c.Reset();
        for (auto& c : ch)
            c.Apply(t);
        for (auto& pp : P.pop) { // TODO -- slightly wasteful; Change could keep track if linked parameter sets need refreshing; then again in most cases probably needed
            pp.needs_recalc = true;
            pp.Recalculate();
        }
    }
}
