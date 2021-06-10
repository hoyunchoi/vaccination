#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <random>
#include <set>
#include <string>
#include <vector>

//* Import from library
#include "CSV.hpp"
#include "Networks.hpp"
#include "pcg_random.hpp"

//* Import from current project
#include "SIRD_V.hpp"

namespace SIRD_V {
struct SA {
    //* Member variables
  private:
    //* Network variables
    Network<unsigned> m_network;

    //* Dynamics variables
    unsigned m_vaccineNum;
    std::set<unsigned> m_vaccinatedIndex;
    std::map<std::string, double> m_rates;
    double m_deltaT{0.1};
    double m_maxTime{10.0};

    //* Random variables for SA
    int m_SA_seed;
    pcg32 m_SA_randomEngine;
    std::uniform_real_distribution<double> m_probDistribution;
    std::uniform_int_distribution<int> m_seedDistribution;
    std::uniform_int_distribution<unsigned> m_indexDistribution;

    //* SA variables
    double m_temperature;
    double m_energy;

    //* Member methods
  public:
    SA() {}
    SA(const Network<unsigned>&, const std::map<std::string, double>&, const unsigned&, const int&, const double&);

    void run(const unsigned&);

  private:
    std::set<unsigned> perturbate();
};

SA::SA(const Network<unsigned>& t_network,
       const std::map<std::string, double>& t_rates,
       const unsigned& t_vaccineNum,
       const int& t_seed,
       const double& t_initialEnergy = std::numeric_limits<double>::max())
    : m_network(t_network),
      m_rates(t_rates),
      m_vaccineNum(t_vaccineNum),
      m_SA_seed(t_seed) {
    /*
        t_network: Network type. Where every disease will be spread
        t_rates: Rates for SIRD_V model
        t_vaccineNum: Number of total vaccine
        t_SA_seed: seed for random engine used in SA algorithm
        t_initialEnergy: Initial energy. If not given, defualt is inf
    */
    //* Generate random engine for SA
    m_probDistribution.param(std::uniform_real_distribution<double>::param_type(0.0, 1.0));
    m_seedDistribution.param(std::uniform_int_distribution<int>::param_type(0, 1000000));
    m_indexDistribution.param(std::uniform_int_distribution<unsigned>::param_type(0, m_network.size - 1));
    m_SA_seed == -1 ? m_SA_randomEngine.seed((std::random_device())()) : m_SA_randomEngine.seed(m_SA_seed);

    //* Randomly choose vaccinated individual
    while (m_vaccinatedIndex.size() == t_vaccineNum) {
        const unsigned vaccinated = m_indexDistribution(m_SA_randomEngine);
        m_vaccinatedIndex.emplace(vaccinated);
    }
}

void SA::run(const unsigned& t_ensembleSize) {
    for (unsigned ensemble = 0; ensemble < t_ensembleSize; ++ensemble) {
        //* Generate random engine for single ensemble
        const int randomEngineSeed = m_seedDistribution(m_SA_randomEngine);
        pcg32 randomEngine(randomEngineSeed);

        //* Perturbation of vaccinated index
        const std::set<unsigned> newVaccinatedIndex = perturbate();

        //* Initialize SIRD_V model and run
        Generator generator(m_network, m_rates, newVaccinatedIndex, randomEngine);
        const unsigned totalDeath = generator.syncRun(m_deltaT, m_maxTime);

        //* Change resulting energy or not
        const double deltaE = (double)totalDeath - m_energy;
        if (m_probDistribution(m_SA_randomEngine) < std::exp(-1.0 * deltaE / m_temperature)){
            m_energy = (double)totalDeath;
            m_vaccinatedIndex = newVaccinatedIndex;
        }
    }
}

std::set<unsigned> SA::perturbate() {
    std::set<unsigned> new_vaccinatedIndex = m_vaccinatedIndex;
    std::uniform_int_distribution<unsigned> vaccineDistribution(0U, m_vaccineNum);

    //* Randomly remove single vaccinated index
    auto it = new_vaccinatedIndex.begin();
    std::advance(it, vaccineDistribution(m_SA_randomEngine));
    new_vaccinatedIndex.erase(it);

    //* Randomly choose new index to vaccinate
    unsigned newVaccinated;
    do {
        newVaccinated = m_indexDistribution(m_SA_randomEngine);
    } while (m_vaccinatedIndex.find(newVaccinated) != m_vaccinatedIndex.end());

    return new_vaccinatedIndex;
}

} // namespace SIRD_V