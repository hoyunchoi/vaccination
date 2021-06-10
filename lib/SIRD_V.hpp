#pragma once

#include <math.h>

#include <map>
#include <random>
#include <set>
#include <string>
#include <vector>

#include "CSV.hpp"
#include "Epidemics.hpp"
#include "Networks.hpp"
#include "pcg_random.hpp"

/*
    SIRD + V simulation
    At initial state, vaccination_rate portion of population are vaccinated

    S+I -> I+I with rate SI_II
    I -> R with rate I_R * (1-death_rate)
    I -> D with rate I_R * death_rate
*/

namespace SIRD_V {

struct Individual : public Node_Epidemic<unsigned> {
    //* Member variables
    double deathProb{0.0};

    //* Generator
    Individual() {}
    Individual(const unsigned& t_index, const double& t_deathProb) : Node_Epidemic(t_index), deathProb(t_deathProb) {}
    Individual(const unsigned& t_index,    //* Index of individual
               const double& t_deathProb,  //* Death probability of individual. Choosen randomly at the start, stays constant as dynamics evolves.
               const std::string& t_state) //* Current state of the individual
        : Node_Epidemic(t_index, t_state),
          deathProb(t_deathProb) {}
};

struct Generator {
  protected:
    //* Member variables
    //* Network parameters
    unsigned m_networkSize;
    double m_meanDegree;
    std::vector<Individual> m_nodes;
    std::set<unsigned> m_reactingIndex;
    std::set<unsigned> m_vaccinatedIndex;

    //* Spreading parameters
    unsigned m_seed_Size;
    double m_SI_II, m_I_R;

    //* Informations
    double m_time;
    unsigned m_numS, m_numI, m_numR, m_numD, m_numV;
    const std::map<std::string, int> m_state2int = {{"S", 0}, {"I", 1}, {"R", 2}, {"D", 3}, {"V", 3}};

    //* Random variables
    pcg32 m_randomEngine;
    std::uniform_real_distribution<double> m_probabilityDistribution;
    std::uniform_int_distribution<unsigned> m_nodeDistribution;

  public:
    //* Store dynamics
    std::vector<std::vector<unsigned>> dynamics;

    //* Public methods
  public:
    Generator() {}
    Generator(const Network<unsigned>&,
              const std::map<std::string, double>&,
              const std::set<unsigned>&,
              const pcg32&);

    const unsigned syncRun(const double&,
                           const double&);

    //* Protected method
  protected:
    void updateTransitionRate(const unsigned&);
    void updateTransitionRate();
    void syncUpdate(const double&);
};

Generator::Generator(const Network<unsigned>& t_network,
                     const std::map<std::string, double>& t_rates,
                     const std::set<unsigned>& t_vaccinatedIndex,
                     const pcg32& t_randomEngine)
    : m_randomEngine(t_randomEngine),
      m_vaccinatedIndex(t_vaccinatedIndex) {
    /*
        t_network: Network disease will be spread
        t_rates : Rate parameters
        t_vaccinatedIndex : Index of nodes who will be vaccinated at the start of dynamcis
        t_randomEngine: Random engine used for spreading dynamics
    */

    //* Set network
    m_networkSize = t_network.size;
    m_meanDegree = t_network.meanDegree;
    //* Initialize random variables
    m_nodeDistribution.param(std::uniform_int_distribution<unsigned>::param_type(0, m_networkSize - 1));
    m_probabilityDistribution.param(std::uniform_real_distribution<double>::param_type(0.0, 1.0));

    //* Generate individuals with their own death prob
    m_nodes.reserve(m_networkSize);
    for (unsigned index = 0; index < m_networkSize; ++index) {
        const double deathProb = std::pow(m_probabilityDistribution(m_randomEngine), 3);
        Individual node(index, deathProb, "S");
        node.neighbors = t_network.adjacency[index];
        m_nodes.emplace_back(node);
    }

    //* Set rates
    m_SI_II = t_rates.at("SI_II");
    m_I_R = t_rates.at("I_R");

    //* Initialize model
    m_time = 0.0;
    m_numV = m_vaccinatedIndex.size();
    m_numS = m_networkSize - m_numV - 1;
    m_numI = 1;
    m_numR = 0;
    m_numD = 0;

    //* Vaccination
    for (const unsigned& index : m_vaccinatedIndex) {
        m_nodes[index].state = "V";
    }

    //* Seed the model
    m_nodes[0].state = "I";
    m_reactingIndex.emplace(0);
    for (const unsigned& neighbor : m_nodes[0].neighbors) {
        if (m_nodes[neighbor].state != "V") {
            m_reactingIndex.emplace(neighbor);
        }
    }

    //* Update the transition rate
    updateTransitionRate();
}

//* Run t_maxTime with t_deltaT. Return num death
const unsigned Generator::syncRun(const double& t_deltaT,
                                  const double& t_maxTime) {
    /*
        t_deltaT: delta T of discrete time gillespie algorithm
        t_maxTime:  Maximum time where iteration will continue. It can early stop when disease is vanished
    */
    while (m_time <= t_maxTime) {
        syncUpdate(t_deltaT);
        dynamics.emplace_back(std::vector<unsigned>{m_numS, m_numI, m_numR, m_numD, m_numV});
        if (m_reactingIndex.empty()) {
            break;
        }
    }
    return m_numD;
}

void Generator::updateTransitionRate(const unsigned& t_index) {
    const int int_state = m_state2int.at(m_nodes[t_index].state);
    switch (int_state) {
    //* S process
    case 0: {
        unsigned infectiousNeighbor = 0;
        for (const unsigned& neighbor : m_nodes[t_index].neighbors) {
            if (m_nodes[neighbor].state == "I") {
                ++infectiousNeighbor;
            }
        }
        m_nodes[t_index].transitionRate = m_SI_II * infectiousNeighbor;
        break;
    }
    //* I process
    case 1: {
        m_nodes[t_index].transitionRate = m_I_R;
        break;
    }
    //* R, D process
    case 2:
    case 3: {
        m_nodes[t_index].transitionRate = 0.0;
        break;
    }
    //* V process
    default: {
        break;
    }
    }
}

void Generator::updateTransitionRate() {
    for (const unsigned& index : m_reactingIndex) {
        updateTransitionRate(index);
    }
}

void Generator::syncUpdate(const double& t_deltaT) {
    //* Do reactions according to each transition rate and add I into reacting nodes
    std::set<unsigned> newReactingIndex;
    for (const unsigned& index : m_reactingIndex) {
        const int intState = m_state2int.at(m_nodes[index].state);
        const double transitionRate = m_nodes[index].transitionRate;
        const double transitionProb = 1.0 - std::exp(-1.0 * transitionRate * t_deltaT);

        switch (intState) {
        //* S process
        case 0: {
            //* S -> I
            if (m_probabilityDistribution(m_randomEngine) <= transitionProb) {
                m_nodes[index].state = "I";
                --m_numS;
                ++m_numI;
                newReactingIndex.emplace_hint(newReactingIndex.end(), index);
            }
            break;
        }

        //* I process
        case 1: {
            if (m_probabilityDistribution(m_randomEngine) <= transitionProb) {
                --m_numI;
                //* I -> D
                if (m_probabilityDistribution(m_randomEngine) <= m_nodes[index].deathProb) {
                    m_nodes[index].state = "D";
                    ++m_numD;
                }
                //* I -> R
                else {
                    m_nodes[index].state = "R";
                    ++m_numR;
                }
            }
            //* I -> I
            else {
                newReactingIndex.emplace_hint(newReactingIndex.end(), index);
            }
            break;
        }
        default: {
            break;
        }
        }
    }

    //* Add neighbor S of I node into reacting nodes
    m_reactingIndex = newReactingIndex;
    for (const unsigned& index : newReactingIndex) {
        for (const unsigned& neighbor : m_nodes[index].neighbors) {
            if (m_nodes[neighbor].state == "S") {
                m_reactingIndex.emplace(neighbor);
            }
        }
    }

    //* Update time and transition rate
    m_time += t_deltaT;
    updateTransitionRate();
}

} // namespace SIRD_V
