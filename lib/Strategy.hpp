#pragma once

#include <algorithm>
#include <vector>
#include <string>

namespace SIRD_V {
std::vector<double> getAgeStrategy(const double& t_totVaccineRatio,
                                   const std::vector<double>& t_ageDist) {
    /*
        t_totVaccineRatio: num_vaccine/tot_population. [0,1]
        t_ageDist: age distribution of modular network
    */

    std::vector<double> strategy(t_ageDist.size(), 0.0);

    //* Iterate over inverse of age distribution: older age group first
    double remainingVaccineRatio = t_totVaccineRatio;
    for (unsigned index = strategy.size() - 1; index >= 0; --index) {
        if (t_ageDist[index] <= remainingVaccineRatio) {
            strategy[index] = t_ageDist[index];
            remainingVaccineRatio -= t_ageDist[index];
        } else {
            strategy[index] = remainingVaccineRatio;
            break;
        }
    }

    return strategy;
}

std::vector<double> getRateStrategy(const double& t_totVaccineRatio,
                                    const std::vector<double>& t_ageDist,
                                    const std::vector<std::vector<double>>& t_contactMatrix,
                                    const double& t_vaccineUnit = 1e-4) {
    /*
        t_totVaccineRatio: num_vaccine/tot_population. [0,1]
        t_ageDist: age distribution of modular network
        t_contactMatrix: adjacency matrix between modules. t_ageDist.size() x t_ageDist.size()
        t_vaccineUnit: Minimum ratio of vaccine to give to certain age group. Default: 1e-8
    */

    const unsigned ageGroupNum = t_ageDist.size();

    //* Initialize weight of each age group
    std::vector<double> weight(ageGroupNum, 0.0);
    for (unsigned index = 0; index < ageGroupNum; ++index) {
        for (unsigned neighbor = 0; neighbor < ageGroupNum; ++neighbor) {
            weight[index] += t_contactMatrix[index][neighbor] * t_ageDist[index] * t_ageDist[neighbor];
        }
    }

    //* Iterate until remaing vaccines become 0
    std::vector<double> vaccinatedRatio(ageGroupNum, 0.0);
    double remainingVaccineRatio = t_totVaccineRatio;
    while (remainingVaccineRatio >= 0) {
        //* Get age group with maximum weight of current network
        const unsigned target = std::max_element(weight.begin(), weight.end()) - weight.begin();

        //* Give vaccine to target age group
        vaccinatedRatio[target] += t_vaccineUnit / t_ageDist[target];
        remainingVaccineRatio -= t_vaccineUnit;

        //* Update weights
        for (unsigned index=0; index<ageGroupNum; ++index){
            weight[index] = 0.0;
            for (unsigned neighbor=0; neighbor < ageGroupNum; ++neighbor){
                weight[index] += t_contactMatrix[index][target] * (1.0 - vaccinatedRatio[index]) * (1.0 - vaccinatedRatio[neighbor]);
            }
        }
    }

    //* Calculate strategy from vaccinated Ratio
    std::vector<double> strategy;
    strategy.reserve(ageGroupNum);
    for (unsigned index=0; index<ageGroupNum; ++index){
        strategy.emplace_back(vaccinatedRatio[index] * t_ageDist[index]);
    }

    return strategy;
}
} // namespace SIRD_V
