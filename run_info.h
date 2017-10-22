#ifndef RUN_INFO_H
#define RUN_INFO_H

#include "configuration.h"
#include <complex>

struct runInfo
{
    double attractorDist;
    double r;
    double a;
    double b;
    AshStateType fixedPoint;
    std::array<std::complex<double>, 3> eigvals;
    bool isFixedPointComputationValid;
    bool isAttractorComputationValid;
    double largestLyapunovExponent;
};

void writeAttractorDataToFile(const int& i, const int& j, const runInfo& curResult, const Configuration& config, const std::vector<AshStateType>& trajectory);

void writeRunInfo(std::ofstream& outFile, const runInfo& curResult, const Configuration& config);

#endif // RUN_INFO_H
