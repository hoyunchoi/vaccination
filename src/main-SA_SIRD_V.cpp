#include <chrono>
#include <iostream>

//* Import from library
#include "Networks.hpp"

//* Import from current project
#include "SA.hpp"


int main(int argc, char* argv[]){
    //* Random variables
    const int randomSeed = 0;

    //* Load network
    Network<unsigned> network;
    network.loadAdjacency("temp.csv");

    //* Dynamics variables
    const unsigned vaccineNum = 10;
    std::map<std::string, double> rates;
    rates["SI_II"] = 1.0;
    rates["I_R"] = 1.0;

    //* SA variables
    const unsigned ensembleSize = 100;
    const double initialTemperature = 2.0;


    auto start = std::chrono::system_clock::now();
    //* --------------------- Initialize and Run -----------------------
    SIRD_V::SA sa(network, rates, vaccineNum, randomSeed, initialTemperature);
    sa.run(ensembleSize);
    //* -----------------------------------------------------------------
    std::chrono::duration<double> sec = std::chrono::system_clock::now() - start;
    std::cout << std::setprecision(10) << sec.count() << "seconds\n";


    return 0;

}