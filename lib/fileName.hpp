#pragma once

#include <string>

#include "stringFormat.hpp"

namespace SIRD_V {
const std::string getNetworkPrefix(const std::string& t_type,
                                   const unsigned& t_networkSize,
                                   const double& t_meanDegree,
                                   const int& t_networkSeed) {
    std::string prefix = t_type + "_N" + std::to_string(t_networkSize) + "M" + to_stringWithPrecision(t_meanDegree, 1);
    return t_networkSeed == -1 ? prefix : prefix + "-" + std::to_string(t_networkSeed);
}

// const std::string getDynamicsPrefix() {
// }

} // namespace SIRD_V