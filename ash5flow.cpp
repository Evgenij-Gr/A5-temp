#include <iostream>
#include <iomanip>
#include "configuration.h"

#include <boost/numeric/odeint.hpp>

using namespace boost::numeric::odeint;

#include "dyn_utils.h"

int main(int argc, char* argv[])
{
    std::cout<<std::setprecision(15);
    if (argc < 2)
    {
        std::cout<<"Usage: "<<argv[0]<<" <path to INI file> [options]"<<std::endl;
    }
    else
    {
        const Configuration config(argc, argv);
        if (!config.isConfigurationConsistent)
        {
            std::cout<<"Config file is not consistent, terminating program"<<std::endl;
            return 1;
        }
        double rCurrent = config.rInit;
        double aCurrent = config.aInit;
        double bCurrent = config.bInit;
        Ashwin5osc system(rCurrent, aCurrent, bCurrent);
        auto dopri5dense = make_dense_output(config.atol , config.rtol , runge_kutta_dopri5< AshStateType >() );
        double startTime = 0.0;
        double timeLength = config.flowMaxTime;
        std::vector<AshStateType> vStates;
        std::vector<double> vTimes;
        simpleWriter<AshStateType> observer(vStates, vTimes);
        AshStateType pt = config.flowInitPoint;
        double dt = config.timeStep;
        integrate_const(dopri5dense, system, pt, startTime, startTime+timeLength, dt, observer);
        std::ofstream outFile("traj-out.txt");
        outFile<<std::setprecision(15)<<std::fixed;
        outFile<<"# 0    1     2     3     4    "<<std::endl;
        outFile<<"# t    Pt[0] Pt[1] Pt[2] Pt[3]"<<std::endl;
        for(int i = 0; i < vStates.size(); i++)
        {
            outFile<<std::setw(21)<<std::left<<vTimes[i]
                   <<std::setw(21)<<std::left<<vStates[i][0]
                   <<std::setw(21)<<std::left<<vStates[i][1]
                   <<std::setw(21)<<std::left<<vStates[i][2]
                   <<std::setw(21)<<std::left<<vStates[i][3]
                   <<std::endl;
        }
    }
    return 0;
}
