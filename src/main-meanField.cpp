#include <chrono>
#include <iostream>
#include <vector>
#include <string>

//* Import from library
#include "CSV.hpp"

//* Import from current project
#include "MeanField.hpp"
#include "Strategy.hpp"

int main(int argc, char* argv[]) {
    //* Read network
    std::vector<std::vector<double>> contactMatrix;
    std::vector<double> ageDist, deathRate;
    CSV::read("data/us_contact_matrix.txt", contactMatrix);
    CSV::read("data/us_age_distribution.txt", ageDist);
    CSV::read("data/us_mortality.txt", deathRate);

    //* Generate strategy
    const double totVaccineRatio = std::stod(argv[1]);
    const std::string strategyName = argv[2];
    std::vector<double> strategy;
    if (strategyName == "age"){
        strategy = SIRD_V::getAgeStrategy(totVaccineRatio, ageDist);
    }
    else if (strategyName == "rate"){
        strategy = SIRD_V::getRateStrategy(totVaccineRatio, ageDist, contactMatrix);
    }

    const auto start = std::chrono::system_clock::now();
    //* Generate mean field object
    SIRD_V::MeanField meanField(strategy, contactMatrix, ageDist, deathRate);
    meanField.run();

    //* Save the data
    const std::string dataDir = "data/meanField/deathRate/";
    CSV::generateDirectory(dataDir);
    meanField.save(dataDir + strategyName + ".csv");

    //* End of calculation and saving
    const auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> second = end - start;
    std::cout << std::setprecision(6) << second.count() << " seconds\n";

    return 0;
}