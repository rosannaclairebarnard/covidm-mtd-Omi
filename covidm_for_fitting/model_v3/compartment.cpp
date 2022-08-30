// compartment.cpp

#include "compartment.h"

// Construct the compartment
Compartment::Compartment()
 : type(Undefined), size(0), rate(0) { }

// Add n individuals to the compartment, with maturation times t
// Note that t[0] means individuals who will be in the compartment for 0 time steps
void Compartment::Add(Parameters& P, Randomizer& Rand, double n, Discrete& mat)
{
    if (type == Undefined && mat.weights.size() > 0)
    {
        if (mat.weights[0] >= 0)
            type = DiscreteDist;
        else
            type = Erlang;
    }
    
    switch (type)
    {
    case Undefined:
        break;
    
    case DiscreteDist:
        // Expand compartment span if needed
        if (contents.size() < mat.weights.size())
            contents.resize(mat.weights.size(), 0);
    
        // Seed compartment
        size += n;
        if (P.deterministic)
        {
            for (unsigned int i = 0; i < mat.weights.size(); ++i)
                contents[i] += mat.weights[i] * n;
        }
        else
        {
            if (P.fast_multinomial)
                mat.mn_approx(n, mat.storage);
            else
                Rand.Multinomial(n, mat.weights, mat.storage);
            for (unsigned int i = 0; i < mat.weights.size(); ++i)
                contents[i] += mat.storage[i];
        }
        break;
        
    case Erlang:
        if (mat.weights.size() > 0)
        {
            rate = -mat.weights[0];

            // Expand compartment span if needed
            if (contents.size() < mat.weights.size())
                contents.resize(mat.weights.size(), 0);
    
            // Seed compartment
            size += n;
            contents.back() += n;
        }
        break;
    }
}

// Remove individuals from the compartment with per-individual probability p; return the number of individuals removed
double Compartment::RemoveProb(Parameters& P, Randomizer& Rand, double p)
{
    double total_removed = 0.0;
    
    switch (type)
    {
    case Undefined:
        break;
    
    case DiscreteDist:
    case Erlang:
        // Clamp value of p
        p = min(1.0, max(0.0, p));
    
        size = 0.0;
        if (P.deterministic)
        {
            for (unsigned int i = 0; i < contents.size(); ++i)
            {
                double to_remove = contents[i] * p;
                contents[i] -= to_remove;
                size += contents[i];
                total_removed += to_remove;
            }
        }
        else
        {
            for (unsigned int i = 0; i < contents.size(); ++i)
            {
                double to_remove = Rand.Binomial(contents[i], p);
                contents[i] -= to_remove;
                size += contents[i];
                total_removed += to_remove;
            }
        }
        break;
    }
    
    return total_removed;
}

// Move individuals from this compartment to another with per-individual probability p; return the number of individuals moved
double Compartment::MoveProb(Compartment& dest, Parameters& P, Randomizer& Rand, double p)
{
    double total_moved = 0.0;
    
    if (type != Undefined && type != dest.type)
    {
        throw runtime_error("Cannot move individuals between compartments of differing types.");
    }

    switch (type)
    {
    case Undefined:
        break;
        
    case DiscreteDist:
    case Erlang:
        // Clamp value of p
        p = min(1.0, max(0.0, p));
    
        size = 0.0;
        dest.size = 0.0;
        if (P.deterministic)
        {
            for (unsigned int i = 0; i < contents.size(); ++i)
            {
                double to_move = contents[i] * p;
                contents[i] -= to_move;
                size += contents[i];
                dest.contents[i] += to_move;
                dest.size += dest.contents[i];
                total_moved += to_move;
            }
        }
        else
        {
            for (unsigned int i = 0; i < contents.size(); ++i)
            {
                double to_move = Rand.Binomial(contents[i], p);
                contents[i] -= to_move;
                size += contents[i];
                dest.contents[i] += to_move;
                dest.size += dest.contents[i];
                total_moved += to_move;
            }
        }

        break;
    }

    return total_moved;
}

// Mature the compartment by one time step, returning the number of individuals who have left the compartment.
double Compartment::Mature(double dt)
{
    double m = 0.0;
    
    switch (type)
    {
    case Undefined:
        break;
        
    case DiscreteDist:
        if (!contents.empty())
        {
            m = contents[0];
            size = max(0., size - m);
            contents.erase(contents.begin());
            contents.push_back(0);
        }
        break;
        
    case Erlang:
        if (!contents.empty())
        {
            double frac = (1.0 - exp(-rate * dt));
            m = contents[0] * frac;
            contents[0] -= m;
            size = max(0., size - m);
            
            for (unsigned int i = 1; i < contents.size(); ++i)
            {
                double x = contents[i] * frac;
                contents[i] -= x;
                contents[i - 1] += x;
            }
        }
        break;
    }

    return m;
}

// Return the total number of individuals in the compartment
double Compartment::Size()
{
    return size;
}

// Print contents of compartment
void Compartment::DebugPrint() const
{
    cout << "size " << size << " compartment size " << contents.size();
    double sum = 0;
    for (auto& c : contents)
    {
        cout << " " << c;
        sum += c;
    }
    cout << " sum " << sum << "\n";
}
