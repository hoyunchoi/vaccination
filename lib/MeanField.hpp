#pragma once

#include <fstream>
#include <limits>
#include <map>
#include <numeric>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

//* Import from library
#include "linearAlgebra.hpp"
#include "stringFormat.hpp"

namespace SIRD_V {

struct AgeGroup {
    //* Member varialbes
  public:
    unsigned index;
    double size;
    double deathProb;

    //* Member methods
    AgeGroup() {}
    AgeGroup(const unsigned& index,
             const double& t_size,
             const double& t_deathProb)
        : index(index),
          size(t_size),
          deathProb(t_deathProb) {}
};

struct MeanField {
    //* Member variables
  private:
    //* Network variables
    std::vector<std::vector<double>> m_contactMatrix;
    std::vector<AgeGroup> m_ageGroups;
    unsigned m_ageGroupNum;

    //* Dynamics variables
    double m_time;
    std::vector<double> m_strategy;
    double m_totVaccine;
    double m_SI_II{0.7};
    double m_I_R{1.0};
    double m_deltaT{0.05};
    const std::map<std::string, int> m_state2int = {{"S", 0},
                                                    {"I", 1},
                                                    {"R", 2},
                                                    {"D", 3},
                                                    {"V", 4}};

    //* Variables for saving
    std::vector<std::vector<double>> m_states;
    // std::map<std::string, double> m_totState;
    std::vector<double> m_totState;
    std::vector<std::vector<double>> m_totStateHistory;

    //* member methods
  public:
    MeanField() {}
    MeanField(const std::vector<double>&,
              const std::vector<std::vector<double>>&,
              const std::vector<double>&,
              const std::vector<double>&,
              const double&);

    void run(const double&);
    void save(const std::string&) const;

  private:
    std::vector<std::vector<double>> dotStates(const std::vector<std::vector<double>>&) const;
    void updateStates();
};

MeanField::MeanField(const std::vector<double>& t_strategy,
                     const std::vector<std::vector<double>>& t_contactMatrix,
                     const std::vector<double>& t_ageDist,
                     const std::vector<double>& t_deathRate,
                     const double& t_seed = 1e-6)
    : m_contactMatrix(t_contactMatrix),
      m_strategy(t_strategy) {

    /*
        t_strategy: strategy to distribute vaccine. sum(t_strategy) = total number of vaccines
        t_contactMatrix: adjacency matrix between modules. t_ageDist.size() x t_ageDist.size()
        t_ageDist: age distribution of modular network
        t_deathRate: death rate of each age groups
        t_seed: Initial seed ratio of each age groups
    */

    //* Construct network to model
    m_ageGroupNum = t_ageDist.size();
    m_ageGroups.reserve(m_ageGroupNum);
    m_states.assign(m_ageGroupNum, std::vector<double>(4, 0.0));
    for (unsigned i = 0; i < m_ageGroupNum; ++i) {
        //* Generate default age group
        const double size = t_ageDist[i];
        AgeGroup ageGroup(i, size, t_deathRate[i]);

        //* Vaccinate and seed according to strategy
        m_states[i][m_state2int.at("I")] = t_seed * size;
        m_states[i][m_state2int.at("S")] = (1.0 - t_seed - m_strategy[i]) * size;
        m_states[i][m_state2int.at("R")] = 0.0;
        m_states[i][m_state2int.at("D")] = 0.0;

        //* Save age groupt to model
        m_ageGroups.emplace_back(ageGroup);
    }

    m_totVaccine = std::accumulate(m_strategy.begin(), m_strategy.end(), 0.0);

    //* Initialize dynamics
    m_time = 0;
    m_totState.assign(4, 0.0);
    m_totState[1] = t_seed * m_ageGroupNum;                 //* I
    m_totState[0] = 1.0 - m_totVaccine - m_totState[1];     //* S
    m_totState[2] = 0.0;                                    //* R
    m_totState[3] = 0.0;                                    //* D
    m_totStateHistory.emplace_back(m_totState);
};

std::vector<std::vector<double>> MeanField::dotStates(const std::vector<std::vector<double>>& t_states) const {
    std::vector<std::vector<double>> dotState;
    for (const AgeGroup& ageGroup : m_ageGroups) {
        const unsigned index = ageGroup.index;

        //* dot S
        double dotS = 0.0;
        for (unsigned neighbor = 0; neighbor < m_ageGroupNum; ++neighbor) {
            dotS -= m_SI_II * m_contactMatrix[index][neighbor] * t_states[index][m_state2int.at("S")] * t_states[neighbor][m_state2int.at("I")];
        }

        //* dot I
        const double dotI = -1.0 * dotS - m_I_R * t_states[index][m_state2int.at("I")];

        //* dot R
        const double dotR = m_I_R * (1.0 - ageGroup.deathProb) * t_states[index][m_state2int.at("I")];

        //* dot D
        const double dotD = m_I_R * ageGroup.deathProb * t_states[index][m_state2int.at("I")];
        dotState.emplace_back(std::vector<double>{dotS, dotI, dotR, dotD});
    }
    return dotState;
}

void MeanField::updateStates() {
    using namespace linearAlgebra;
    const std::vector<std::vector<double>> delta_state1 = dotStates(m_states);
    const std::vector<std::vector<double>> delta_state2 = dotStates(m_states + delta_state1 * m_deltaT / 2.0);
    const std::vector<std::vector<double>> delta_state3 = dotStates(m_states + delta_state2 * m_deltaT / 2.0);
    const std::vector<std::vector<double>> delta_state4 = dotStates(m_states + delta_state3 * m_deltaT);

    m_states += (delta_state1 + 2.0 * delta_state2 + 2.0 * delta_state3 + delta_state4) * m_deltaT / 6.0;
}

void MeanField::run(const double& t_eps = 1e-12) {
    while (true) {
        //* Update every states ratio at each age distribution
        updateStates();

        //* update every states ratio in total network
        for (unsigned state=0; state<m_totState.size(); ++state){
            m_totState[state] = 0.0;
            for (unsigned index = 0; index < m_ageGroupNum; ++index) {
                m_totState[state] += m_states[index][state];
            }
        }
        m_totStateHistory.emplace_back(m_totState);

        //* update time
        m_time += m_deltaT;

        //* Stop calculation when infected nodes are less than t_eps
        if (m_totState[1] < t_eps) {
            break;
        }
    }
}

void MeanField::save(const std::string& t_filePath) const {
    std::ofstream writeFile(t_filePath, std::ios_base::app);
    writeFile.precision(std::numeric_limits<double>::digits10 + 1);
    writeFile << m_totVaccine << "," << m_totState[3] << "\n";
}

} // namespace SIRD_V