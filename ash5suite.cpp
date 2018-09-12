#include "ash5.h"
#include "dyn_utils.h"

#define EIGEN_DONT_PARALLELIZE
#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>
#include <chrono>
#include <complex>

//#include <omp.h>

typedef Eigen::Vector3d AshCrossSectionPtType;

AshStateType reprojectFromPoincareSection (const AshCrossSectionPtType &pt)
{
    return {pt[0], pt[1], 1.2, pt[2]};
}

AshCrossSectionPtType projectOnPoincareSection (const AshStateType &pt)
{
    return {pt[0], pt[1], pt[3]};
}

#include "configuration.h"

#include <boost/numeric/odeint.hpp>

using namespace boost::numeric::odeint;

class PoincareMap
{
private:
    double rCurrent, aCurrent, bCurrent;
    double rtol, atol, eventTol, approximateReturnTime;
    double crossingDirection;
    double skipTime;
public:
    PoincareMap(std::map<std::string, double> actualOdeParams, const Configuration& config)
    {
        rCurrent = actualOdeParams["r"];
        aCurrent = actualOdeParams["a"];
        bCurrent = actualOdeParams["b"];
        ////
        rtol = config.rtol;
        atol = config.atol;
        ////
        eventTol = config.eventTol;
        crossingDirection = config.crossingDirection;
        approximateReturnTime = config.approximateReturnTime;
        skipTime = config.skipTime;
    }
    std::vector<AshStateType> computeTrajectory(const AshCrossSectionPtType& startPt, int iterMax=1) const
    {
        // There are some moments that are sketchy in that function
        // For example: what if we had started not on Poincare section?
        // How to check reliably that trajectory is ok
        Ashwin5osc system(rCurrent, aCurrent, bCurrent);
        Ash_Section_Event evt;
        int maxEvtsCapacity = iterMax+1;
        std::vector<AshStateType> trajectoryStorage;
        trajectoryStorage.reserve(maxEvtsCapacity);
        auto pt = reprojectFromPoincareSection(startPt);
        auto dopri5dense = make_dense_output(atol , rtol , runge_kutta_dopri5< AshStateType >() );
        auto ashPoincareObserver = EventObserver<Ash_Section_Event, decltype(dopri5dense), Ashwin5osc, AshStateType>(system, evt, dopri5dense, trajectoryStorage, crossingDirection, eventTol, skipTime);
        int iter = 0;
        double startTime = 0.0;
        while ((trajectoryStorage.size() < maxEvtsCapacity) && (iter < 2*maxEvtsCapacity)) // the second condition is to avoid infinite cycle when we stop returning to cross-section
        {
            // since the system is autonomous we can do this trick
            integrate_adaptive(dopri5dense, system, pt, startTime, startTime+approximateReturnTime, 1e-12, ashPoincareObserver);
            iter++;
            startTime += approximateReturnTime;
        }
        return trajectoryStorage;
    }
    AshCrossSectionPtType poincareMap(const AshCrossSectionPtType& pt)
    {
        return projectOnPoincareSection(computeTrajectory(pt)[1]);
    }
    AshCrossSectionPtType residueMap(const AshCrossSectionPtType& pt)
    {
        return (poincareMap(pt)-pt);
    }
};

typedef AshCrossSectionPtType PtType;
typedef Eigen::Matrix3d MatType;

#include "newt_sol.h"

#include "run_info.h"

int main(int argc, char* argv[])
{
    std::cout<<std::setprecision(15);
    if (argc < 2)
    {
        std::cout<<"Usage: "<<argv[0]<<" <path to INI file> [options]"<<std::endl;
        std::cout<<"Options:\n"<<std::endl;
        std::cout<<"--eigv Analyze eigenvalues of nearest fixed point"<<std::endl;
        std::cout<<"--lyap Compute Largest Lyapunov Exponent"<<std::endl;
        std::cout<<"--dist Compute distance from attractor to fixed point"<<std::endl;
        std::cout<<"--plot Save attractors to text files"<<std::endl;
        std::cout<<"If no optional arguments provided, the program only continues the fixed point and provides only that information"<<std::endl;
    }
    else
    {
        std::chrono::time_point<std::chrono::system_clock> start, end;
        std::chrono::duration<double> elapsed_seconds;
     // PROCESS COMMAND LINE OPTIONS, READ THE INI-FILE AND INITIALIZE STARTING FIXED POINT AND PARAMETERS
        const Configuration config(argc, argv);
        if (!config.isConfigurationConsistent)
        {
            std::cerr<<"Config file is not consistent, terminating program"<<std::endl;
            return 1;
        }
     // reserve space for storing results
        std::vector<std::vector<runInfo>> results;
        results.resize(config.v2_N);
        for (int i = 0; i < config.v2_N; i++)
        {
            results[i].resize(config.v1_N);
        }
     // #MAIN LOOPS
     // #Initial loop for establishing parallelism
     // we fix the column and iterate rows
        std::cout<<"PRE-COMPUTING FIXED POINTS IN THE FIRST COLUMN\n"<<std::endl;
        AshCrossSectionPtType pt = projectOnPoincareSection(config.initialFixedPoint);
        start = std::chrono::system_clock::now();
        for (int i = 0; i<config.v2_N; i++)
        {
        // assign parameters
            auto params = config.getParameterValues(i, 0);//!
            std::cout<<"\n\n"<<i<<"/"<<(config.v2_N-1)<<" r = "<<params["r"]<<" a = "<<params["a"]<<" b = "<<params["b"]<<std::endl;
            bool fixPtSuccessStatus;
            NewtonSolver<PtType, MatType> ns(config.fixPtTol,
                                             config.fixPtMaxIter,
                                             config.jacobiNewtonEps,
                                             config.fixPtEpsNewtonStep);
            PoincareMap pm(params, config);
            auto F = std::bind(&PoincareMap::residueMap, pm, std::placeholders::_1);
            pt = ns.performNewtonMethod(F, pt);
            fixPtSuccessStatus = ns.getSuccessStatus();
//             pt = performNewtonMethodForFixedPointFinding(params, config, pt, fixPt_success_status);
        // #Save info about fixed point in array
//             std::cout<<"After Newton Method: ";
            results[i][0].fixedPoint = reprojectFromPoincareSection(pt);
//             std::cout<<"fixedPoint = \n"<<projectOnPoincareSection(results[i][0].fixedPoint)<<std::endl;
            results[i][0].isFixedPointComputationValid = fixPtSuccessStatus;
            ////////////////////////
//             PoincareMap pmtest(params, config);
//             std::cout<<"Instant Re-check of precision 1: \n"<<pmtest.residueMap(projectOnPoincareSection(results[i][0].fixedPoint)).norm()<<std::endl;
        }
        end = std::chrono::system_clock::now();
        elapsed_seconds = end - start;
        std::cout<<"Pre-calculating fixed points in first column: "<<elapsed_seconds.count()*1000.0<<"ms"<<std::endl;
     // #FULL SWEEP
        std::cout<<"\n\nFULL SWEEP"<<std::endl;
//        Eigen::initParallel();
        for (int j = 0; j < config.v1_N; j++)
        {
//            #pragma omp parallel for num_threads(1)
//            #pragma omp parallel for
            for (int i = 0; i < config.v2_N; i++)
            {
                auto& curResult = results[i][j];
                std::cout<<"\n\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"<<std::endl;
                auto params = config.getParameterValues(i, j);
                curResult.r = params["r"];
                curResult.a = params["a"];
                curResult.b = params["b"];
                std::cout<<i<<"/"<<(config.v2_N-1)<<" "<<j<<"/"<<(config.v1_N-1)<<" r = "<<curResult.r<<" a = "<<curResult.a<<" b = "<<curResult.b<<std::endl;
            // #Recalculate fixed point
                AshCrossSectionPtType fixPt;
                bool fixPtSuccessStatus;
                if (j==0)
                {
                    fixPt = projectOnPoincareSection(results[i][0].fixedPoint);
                    fixPtSuccessStatus = results[i][0].isFixedPointComputationValid;
//                     std::cout<<"Instant Re-check of fixed point: \n"<<fixPt<<std::endl;
//                     PoincareMap pmtest(params, config);
//                     std::cout<<"Instant Re-check of precision 2: \n"<<pmtest.residueMap(fixPt).norm()<<std::endl;
                }
                else
                {
                    if (results[i][j-1].isFixedPointComputationValid)
                    {
                        start = std::chrono::system_clock::now();
                        AshCrossSectionPtType prevApproximation = projectOnPoincareSection(results[i][j-1].fixedPoint);
                        NewtonSolver<PtType, MatType> ns(config.fixPtTol,
                                                         config.fixPtMaxIter,
                                                         config.jacobiNewtonEps,
                                                         config.fixPtEpsNewtonStep);
                        PoincareMap pm(params, config);
                        auto F = std::bind(&PoincareMap::residueMap, pm, std::placeholders::_1);
                        fixPt = ns.performNewtonMethod(F, prevApproximation);
                        fixPtSuccessStatus = ns.getSuccessStatus();
                        end = std::chrono::system_clock::now();
                        elapsed_seconds = end - start;
                        std::cout<<"Calculating fixed point: "<<elapsed_seconds.count()*1000.0<<"ms"<<std::endl;
                        curResult.fixedPoint = reprojectFromPoincareSection(fixPt);
                    }
                    else
                    {
                        fixPtSuccessStatus = false;
                        std::cout<<"ERROR: FIXED POINT CONTINUATION FAILED"<<std::endl;
                    }
                }
                curResult.isFixedPointComputationValid = fixPtSuccessStatus;
            // #Do other stuff if finding fixed point didn't fail
                bool additionalCalculationsNeeded = config.flags.isEigvalsNeeded || config.flags.isLyapNeeded || config.flags.isDistNeeded || config.flags.isPlotNeeded;
                if (curResult.isFixedPointComputationValid && additionalCalculationsNeeded)
                {
                    if (config.flags.isEigvalsNeeded)
                    {
                        double epsOrt = config.jacobiFixPtTypeEps;
                        AshCrossSectionPtType pt = fixPt;
                        PoincareMap pm(params, config);
                        std::cout<<"residue.norm() = "<<pm.residueMap(pt).norm()<<std::endl;
                        Eigen::Matrix3d Jac;
                        auto F = std::bind(&PoincareMap::poincareMap, pm, std::placeholders::_1);
                        NumericalJacobiMatrix<PtType, MatType> jmat(epsOrt);
                        jmat.computeJacobiMatrix(F, pt, Jac);
                        Eigen::EigenSolver<Eigen::Matrix3d> es(Jac);
                        std::cout<<"Jacobi = \n"<<Jac<<std::endl;
                        std::cout<<"Eigenvalues are \n"<<es.eigenvalues()<<std::endl;
                        for (int k =0; k<3; k++)
                        {
                            curResult.eigvals[k] = es.eigenvalues()[k];
                            std::cout<<"eigv["<<k<<"] = "<<es.eigenvalues()[k]<<" (abs = "<<std::abs(es.eigenvalues()[k])<<")"<<std::endl;
                        }
                    }
                    bool isAttractorCalculationNeeded = config.flags.isLyapNeeded
                                                     || config.flags.isDistNeeded
                                                     || config.flags.isPlotNeeded;
                    bool rowStride(true), colStride(true);
                    if (config.v2_att_stride != 0)
                    {
                        rowStride = ((i % (config.v2_att_stride+1)) == 0);
                    }
                    if (config.v1_att_stride != 0)
                    {
                        colStride = ((j % (config.v1_att_stride+1)) == 0);
                    }
                    bool strideCondition = rowStride && colStride;
                    if (isAttractorCalculationNeeded && strideCondition)
                    {
                    // Calculate attractor
                        AshCrossSectionPtType nearPt = fixPt + AshCrossSectionPtType({1.0, 0.0, 0.0}) * 0.001; // MAGIC CONSTANT HERE!!!
                        PoincareMap pm(params, config);
                        start = std::chrono::system_clock::now();
                        auto trajectory = pm.computeTrajectory(nearPt, config.mapIterLast);
                        end = std::chrono::system_clock::now();
                        elapsed_seconds = end - start;
                        std::cout<<"Calculating attractor: "<<elapsed_seconds.count()<<"s"<<std::endl;
                        // wrong here: successful trajectory computation has at least config.mapIterLast+1 elements so it should be <= instead of
                        curResult.isAttractorComputationValid = ((config.mapIterLast + 1) <= trajectory.size());
                        if (curResult.isAttractorComputationValid)
                        {
                            if (config.flags.isDistNeeded)
                            {
                            // #Find distance from attractor to fixed point if needed
                                start = std::chrono::system_clock::now();
                                double min_dist;
                                min_dist = (projectOnPoincareSection(trajectory[config.mapIterSkip])-fixPt).norm();
                                for (int k=config.mapIterSkip+1; k <= config.mapIterLast; k++)
                                {
                                    double cur_dist = (fixPt-projectOnPoincareSection(trajectory[k])).norm();
                                    min_dist = std::min(cur_dist, min_dist);
                                }
                                end = std::chrono::system_clock::now();
                                elapsed_seconds = end - start;
                                std::cout<<"Calculating distance from fixed point to attractor: "<<elapsed_seconds.count()*1000.0<<"ms"<<std::endl;
                                curResult.attractorDist = min_dist;
                            }
                            if (config.flags.isLyapNeeded)
                            {
                                std::cout<<"Computing Largest Lyapunov Exponent: not implemented"<<std::endl;
                            }
                            if (config.flags.isPlotNeeded)
                            {
                            // #Save attractor data to external file (r, a, b, LLE, dist, type of fixed point, eigenvalues) if needed
                                std::cout<<"Saving parameters and attractor data to external file"<<std::endl;
                                writeAttractorDataToFile(i, j, curResult, config, trajectory);
                            }
                        }
                    }
                }
            }
        }
        // #SAVE SWEEP RESULTS TO FILE
        std::ofstream outFile("tst-res.txt");
        outFile<<std::setprecision(15)<<std::fixed;
        outFile<<"# 0 1 2 3        4        5        6        7       8   9          10         11         12         13         14         15           16           17"<<std::endl;
        outFile<<"# r a b fixPt[0] fixPt[1] fixPt[2] fixPt[3] minDist LLE eigv[0].re eigv[0].im eigv[1].re eigv[1].im eigv[2].re eigv[2].im abs(eigv[0]) abs(eigv[1]) abs(eigv[2])"<<std::endl;
        for(int i = 0; i < config.v2_N; i++)
        {
            for (int j=0; j < config.v1_N; j++)
            {
                auto curResult = results[i][j];
                writeRunInfo(outFile, curResult, config);
            }
        }
    }
    std::cout<<"\nTHE END"<<std::endl;
    return 0;
}
