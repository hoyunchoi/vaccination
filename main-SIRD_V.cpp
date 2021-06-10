#include <chrono>
#include <iostream>
#include <map>
#include <random>
#include <set>
#include <string>

//* Import from library
#include "CSV.hpp"
#include "Networks.hpp"
#include "pcg_random.hpp"

//* Import from current project
#include "SIRD_V.hpp"
#include "fileName.hpp"

int main(int argc, char* argv[]) {

    //* Random variables: network generation
    const int networkSeed = 0;
    pcg32 networkEngine;
    networkSeed == -1 ? networkEngine.seed((std::random_device())()) : networkEngine.seed(networkSeed);

    //* Random variables: dynamics generation
    const int randomSeed = 0;
    pcg32 randomEngine;
    randomSeed == -1 ? randomEngine.seed((std::random_device())()) : randomEngine.seed(randomSeed);

    //* Network variables
    const std::string networkType = "ER";
    const unsigned networkSize = 100U;
    const double meanDegree = 5.0;
    const unsigned long long linkSize = networkSize * meanDegree / 2;
    const Network<unsigned> network = ER::generate(networkSize, linkSize, networkEngine);

    //* Dynamics variables
    std::set<unsigned> vaccinatedIndex = {1, 2, 3, 4};
    std::map<std::string, double> rates;
    rates["SI_II"] = 1.0;
    rates["I_R"] = 1.0;
    const double deltaT = 0.1;
    const double maxTime = 10.0;

    //* Save variables
    namespace fs = std::filesystem;
    const std::string dataDirectory = "data/" + SIRD_V::getNetworkPrefix(networkType, networkSize, meanDegree, networkSeed) + "/";
    CSV::generateDirectory(dataDirectory);

    auto start = std::chrono::system_clock::now();
    //* --------------------- Initialize and Run -----------------------
    SIRD_V::Generator generator(network, rates, vaccinatedIndex, randomEngine);
    const unsigned totalDeath = generator.syncRun(deltaT, maxTime);
    //* -----------------------------------------------------------------
    std::chrono::duration<double> sec = std::chrono::system_clock::now() - start;
    std::cout << std::setprecision(10) << sec.count() << "seconds\n";

    return 0;
}