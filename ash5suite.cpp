#include <boost/version.hpp>

#if BOOST_VERSION == 106400
#include <boost/serialization/array_wrapper.hpp>
#endif

#include "ash5.h"
#include "dyn_utils.h"

#define EIGEN_DONT_PARALLELIZE
#define EIGEN_DONT_ALIGN_STATICALLY
//#include <Eigen/Dense>
#include <chrono>
#include <Eigen/Eigenvalues>
#include <complex>
#include <functional>
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

#include "configuration.h"

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

class NumericalJacobiMatrix
{
private:
    double epsOrt;
public:
    NumericalJacobiMatrix(double eO): epsOrt(eO) {}
    void computeJacobiMatrix(const std::function<AshCrossSectionPtType(AshCrossSectionPtType)>& F,
                        const AshCrossSectionPtType& pt, Eigen::Matrix3d& jacobiMatrix)
    {
        Eigen::Matrix3d Jac;
        std::array<AshCrossSectionPtType, 3> orts = {{{1.,0.,0.},{0.,1.,0.},{0.,0.,1.}}};
        AshCrossSectionPtType F0= F(pt);
        for (int i = 0; i < 3; i++)
        {
            Jac.col(i) = (F(pt + epsOrt * orts[i]) - F0)/epsOrt;
        }
        jacobiMatrix = Jac;
    }
};

AshCrossSectionPtType solveSLE(const Eigen::Matrix3d& A, const AshCrossSectionPtType& b)
{
    AshCrossSectionPtType result = A.colPivHouseholderQr().solve(b);
    return result;
}

double l2norm(AshCrossSectionPtType pt)
{
    return pt.norm();
}

typedef AshCrossSectionPtType PtType;
typedef Eigen::Matrix3d MatType;

class NewtonSolver
{
private:
    double fixPtTol;
    double epsOrt;
    bool successStatus;
    int maxIterNum;
public:
    NewtonSolver(const Configuration& config)
    {
        fixPtTol = config.fixPtTol;
        maxIterNum = config.fixPtMaxIter;
        epsOrt = config.jacobiNewtonEps;
    }
    PtType performNewtonMethod(const std::function<PtType(PtType)>& F, const PtType& initPt)
    {
// //     std::cout<<"[Newton] Entry point"<<std::endl;
        PtType pt = initPt;
        PtType residue = F(pt);
// //     std::cout<<"[Newton] Init residue = \n"<<residue<<std::endl;
// //     std::cout<<"[Newton] Init residue norm = "<<residue.norm()<<std::endl;
        int iterNum = 0;
        while ((l2norm(residue) > fixPtTol) && (iterNum < maxIterNum))
        {
// //         std::cout<<"#############################################################"<<std::endl;
// //         std::cout<<"[Newton] Residue before = \n"<<residue<<std::endl;
// //         std::cout<<"[Newton] Residue.norm() before = "<<residue.norm()<<std::endl;
// //         std::cout<<"[Newton] Pt before = \n"<<pt<<std::endl;
            ++iterNum;
            MatType Jac;
            NumericalJacobiMatrix jmat(epsOrt);
// //         std::cout<<"[Newton] Jacobi = \n"<<Jac<<std::endl;
            //     // there is a trick to lose less precision here, see http://www.crbond.com/nonlinear.htm
            jmat.computeJacobiMatrix(F, pt, Jac);
//     // solve Jac * delta = -residue
            PtType delta = solveSLE(Jac, -residue);
            pt = pt + delta;
            residue = F(pt);
// //         std::cout<<"[Newton] Residue after = \n"<<residue<<std::endl;
// //         std::cout<<"[Newton] Residue.norm() after = "<<residue.norm()<<std::endl;
// //         std::cout<<"[Newton] Pt after = \n"<<pt<<std::endl;
        }
        double precision = l2norm(residue);
// //     std::cout<<"[Newton] Final point: \n"<<pt<<std::endl;
// //     std::cout<<"[Newton] Final precision: "<<precision<<std::endl;
        successStatus = (precision < fixPtTol);
        return pt;
    }
    bool getSuccessStatus()
    {
        return successStatus;
    }
};

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
            NewtonSolver ns(config);
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
                        NewtonSolver ns(config);
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
                        NumericalJacobiMatrix jmat(epsOrt);
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
                    bool isAttractorCalculationNeeded = config.flags.isLyapNeeded || config.flags.isDistNeeded || config.flags.isPlotNeeded;
                    if (isAttractorCalculationNeeded)
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
