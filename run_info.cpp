#include "run_info.h"
#include <iostream>
#include <iomanip>

void writeAttractorDataToFile(const int& i, const int& j, const runInfo& curResult, const Configuration& config, const std::vector<AshStateType>& trajectory)
{
    std::stringstream fileName;
//     std::cout<<"Saving parameters and attractor data to external file:"<<std::endl;
    fileName<<"att_"<<std::setfill('0')<<std::setw(std::to_string(config.v2_N).length())
        <<i<<"_"<<std::setfill('0')<<std::setw(std::to_string(config.v1_N).length())<<j<<".txt";
    std::ofstream outFile(fileName.str());
    // print header with parameter informantion
    // we use NumPy convention for comments: #
    outFile<<"# r = "<<curResult.r<<std::endl;
    outFile<<"# a = "<<curResult.a<<std::endl;
    outFile<<"# b = "<<curResult.b<<std::endl;
    outFile<<"# i = "<<i<<", "<<i<<" / "<<(config.v2_N-1)<<std::endl;
    outFile<<"# j = "<<j<<", "<<j<<" / "<<(config.v1_N-1)<<std::endl;
    outFile<<std::setprecision(15)<<std::fixed;
    outFile<<"# fixed_pt = [ "<<curResult.fixedPoint[0]<<", "
        <<curResult.fixedPoint[1]<<", "
        <<curResult.fixedPoint[2]<<", "
        <<curResult.fixedPoint[3]<<" ]"<<std::endl;
    if (config.flags.isEigvalsNeeded)
    {
        outFile<<"# eigenvalues = ["<<curResult.eigvals[0]
            <<", "<<curResult.eigvals[1]
            <<", "<<curResult.eigvals[2]<<"]"<<std::endl;
        outFile<<"# abs_of_eigenvalues = ["<<std::abs(curResult.eigvals[0])
            <<", "<<std::abs(curResult.eigvals[1])
            <<", "<<std::abs(curResult.eigvals[2])<<"]"<<std::endl;
    }
    if (config.flags.isLyapNeeded)
    {
        ;
    }
    if (config.flags.isDistNeeded)
    {
        outFile<<"# attractor_dist = "<<curResult.attractorDist<<std::endl;
    }
    outFile<<std::setprecision(15)<<std::fixed;
    for (int k = config.mapIterSkip; k <= config.mapIterLast; k++)
    {
        outFile<<std::setw(21)<<std::left<<trajectory[k][0]<<" "
            <<std::setw(21)<<std::left<<trajectory[k][1]<<" "
            <<std::setw(21)<<std::left<<trajectory[k][2]<<" "
            <<std::setw(21)<<trajectory[k][3]<<std::endl;
    }
}

void writeRunInfo(std::ofstream& outFile, const runInfo& curResult, const Configuration& config)
{
    // r a b minDist
    outFile<<curResult.r<<" "<<curResult.a<<" "<<curResult.b;
    if (curResult.isFixedPointComputationValid)
    {
        outFile<<" "<<curResult.fixedPoint[0]<<" "
            <<curResult.fixedPoint[1]<<" "
            <<curResult.fixedPoint[2]<<" "
            <<curResult.fixedPoint[3];
    }
    else
    {
        outFile<<" None None None";
    }
    if (config.flags.isDistNeeded)
    {
        outFile<<" "<<curResult.attractorDist;
    }
    else
    {
        outFile<<" None";
    }
    if (config.flags.isLyapNeeded)
    {
        //outFile<<" "<<curResult.largestLyapunovExponent;
        outFile<<" None";
    }
    else
    {
        outFile<<" None";
    }
    if (config.flags.isEigvalsNeeded)
    {
        outFile<<" "<<curResult.eigvals[0].real()
               <<" "<<curResult.eigvals[0].imag()
               <<" "<<curResult.eigvals[1].real()
               <<" "<<curResult.eigvals[1].imag()
               <<" "<<curResult.eigvals[2].real()
               <<" "<<curResult.eigvals[2].imag();
        outFile<<" "<<std::abs(curResult.eigvals[0])
            <<" "<<std::abs(curResult.eigvals[1])
            <<" "<<std::abs(curResult.eigvals[2]);
    }
    else
    {
        outFile<<" None None None None None None None None None";
    }
    outFile<<std::endl;
}
