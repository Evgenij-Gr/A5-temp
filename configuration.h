#ifndef CONFIGURATION_H
#define CONFIGURATION_H

#include <iostream>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include "ash5.h"
#include <map>

struct Configuration
{
public:
    double rInit, aInit, bInit;
    double rtol, atol, eventTol, approximateReturnTime, fixPtTol;
    int fixPtMaxIter;
    AshStateType initialFixedPoint;
    std::string v1_name, v2_name;
    double v1_min, v1_max, v2_min, v2_max;
    int v2_N, v1_N, v2_att_stride, v1_att_stride;
    int mapIterSkip, mapIterLast;
    struct UsageFlags
    {
        bool isEigvalsNeeded;
        bool isLyapNeeded;
        bool isDistNeeded;
        bool isPlotNeeded;
        UsageFlags(): isEigvalsNeeded(false), isLyapNeeded(false), isDistNeeded(false), isPlotNeeded(false)
        {
            ;
        }
    };
    UsageFlags flags;
    double crossingDirection;
    bool isConfigurationConsistent;
    double jacobiNewtonEps, jacobiFixPtTypeEps;
    double skipTime;
    double fixPtEpsNewtonStep;
    AshStateType flowInitPoint;
    double flowMaxTime;
    double timeStep;
    /////////////////////////
    Configuration(int argc, char* argv[]);
    std::map<std::string, double> getParameterValues(int i, int j) const;
private:
    bool checkConfigFileForConsistency(boost::property_tree::ptree config) const;
};

#endif // CONFIGURATION_H
