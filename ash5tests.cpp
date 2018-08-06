#include "misc_defs.h"

int main(int argc, char* argv[])
{
    std::cout<<std::setprecision(8);
    // Test: finding events of harmonic oscillator
    {
        std::cout<<"[HARMONIC OSCILLATOR TESTS]"<<std::endl;
        bool testResult;
        int iterMax = 10;
        double atol = 1e-13;
        double rtol = 1e-13;
        double approximateReturnTime = 6.5;
        double crossingDirection = 1.0;
        int maxEvtsCapacity = iterMax +1;
        double evtTolerance = 1e-13;
        double timeSkip = 1e-3;
        HarmonicOscillator ho(1.0);
        EventCalculator<HarmonicOscillator,
                HarmonicOscillatorEvent, HarmOscType> evcalc(atol, rtol,
                                                             approximateReturnTime,
                                                             crossingDirection,
                                                             iterMax,evtTolerance,
                                                             timeSkip);
        // ####################################################################
        {
            std::cout<<"\n@Starting on point not from a cross-section"<<std::endl;
            HarmonicOscillatorEvent hoEvent;
            HarmOscType startPt = {3., 0.1};
            auto trajectoryStorage = evcalc.findEvents(ho, hoEvent, startPt);
            std::cout<<"Events: "<<std::endl;
            for (int i = 0; i <trajectoryStorage.size(); i++)
            {
                auto evt = trajectoryStorage[i];
                std::cout<<"Event["<<i<<"] = ("<<evt[0]<<", "<<evt[1]<<")"<<std::endl;
            }
            // ////////////////////////////////////////////////////////////////

            {
                std::cout<<"[Test] Must have found at "
                           "least "<<(iterMax+1)<<" events"<<std::endl;
                testResult = (trajectoryStorage.size() >= maxEvtsCapacity );
                std::cout<<((testResult)?("PASSED"):("FAILED"))<<std::endl;
            }

            {
                std::cout<<"[Test] All X coordinates of "
                           "events should be negative"<<std::endl;
                testResult = std::all_of(trajectoryStorage.begin(),
                                         trajectoryStorage.end(),
                                         [](HarmOscType x){return (x[0]<0.0);});
                std::cout<<((testResult)?("PASSED"):("FAILED"))<<std::endl;
            }

            {
                std::cout<<"[Test] Starting point is not "
                           "the first point of trajectory"<<std::endl;
                testResult = (trajectoryStorage[0]!=startPt);
                std::cout<<((testResult)?("PASSED"):("FAILED"))<<std::endl;
            }

            {
                std::cout<<"[Test] All points "
                           "meet the event tolerance"<<std::endl;
                testResult = std::all_of(trajectoryStorage.begin(),
                                         trajectoryStorage.end(),
                                         [hoEvent, evtTolerance](HarmOscType x)
                                {return (std::abs(hoEvent(x))<evtTolerance);});
                std::cout<<((testResult)?("PASSED"):("FAILED"))<<std::endl;
            }
        }
        // ####################################################################
        {
            std::cout<<"\n@Starting on point from a cross-section"
                   " in wrong direction"<<std::endl;
            HarmonicOscillatorEvent hoEvent;
            HarmOscType startPt = {3., 0.0};
            auto trajectoryStorage = evcalc.findEvents(ho, hoEvent, startPt);
            std::cout<<"Events: "<<std::endl;
            for (int i = 0; i <trajectoryStorage.size(); i++)
            {
                auto evt = trajectoryStorage[i];
                std::cout<<"Event["<<i<<"] = ("<<evt[0]<<", "<<evt[1]<<")"<<std::endl;
            }
            // //////////////////////////////////////////////////////////////

            {
                std::cout<<"[Test] Must have found at "
                       "least "<<(iterMax+1)<<" events"<<std::endl;
                testResult = (trajectoryStorage.size() >= maxEvtsCapacity );
                std::cout<<((testResult)?("PASSED"):("FAILED"))<<std::endl;
            }

            {
                std::cout<<"[Test] All X coordinates of "
                       "events should be negative"<<std::endl;
                testResult = std::all_of(trajectoryStorage.begin(),
                                     trajectoryStorage.end(),
                                     [](HarmOscType x){return (x[0]<0.0);});
                std::cout<<"There is a reason why this one doesn't fail in my problem:"
                " I put pt on a cross-section and choose crossing direction"
                " depending on where it moves"<<std::endl;
                std::cout<<((testResult)?("PASSED"):("FAILED"))<<std::endl;
            }

            {
                std::cout<<"[Test] All points "
                           "meet the event tolerance"<<std::endl;
                testResult = std::all_of(trajectoryStorage.begin(),
                                         trajectoryStorage.end(),
                                         [hoEvent, evtTolerance](HarmOscType x)
                                {return (std::abs(hoEvent(x))<evtTolerance);});
                std::cout<<((testResult)?("PASSED"):("FAILED"))<<std::endl;
            }
        }
        // ####################################################################
        {
            std::cout<<"\n@Event that will never happen"<<std::endl;
            HarmonicOscillatorEvent hoEvent(4.0);
            HarmOscType startPt = {3., 0.0};
            auto trajectoryStorage = evcalc.findEvents(ho, hoEvent, startPt);
            {
                std::cout<<"[Test] Must have found no events"<<std::endl;
                testResult = (trajectoryStorage.size() == 0);
                std::cout<<((testResult)?("PASSED"):("FAILED"))<<std::endl;
            }
        }
    }
    // Test: finding zeroes of simple functions
    {
        std::cout<<"\n\n[ARTIFICIAL FUNCTIONS TEST]"<<std::endl;
        // ######################################
        {
            DummyRHS drhs;
            DummyCos dcos;
            double t0 = 0.0;
            const double PI = std::acos(-1.0);
            double dt = 2.0*PI;
            bool testResult;
            int iterMax = 5;
            double atol = 1e-13;
            double rtol = 1e-13;
            double approximateReturnTime = 6.5;
            double evtTolerance = 1e-13;
            double timeSkip = 1e-3;
            {
                std::cout<<"\n@Finding ascending zeroes of cosine function on segment"<<std::endl;
                double crossingDirection = 1.0;
                double trueEventCoordinate = 3.0*PI/2.0;
                EventCalculator<DummyRHS, DummyCos, double> evcalc(atol, rtol,
                                                       approximateReturnTime,
                                                       crossingDirection,
                                                       iterMax,evtTolerance,
                                                       timeSkip);
                auto trajectoryStorage = evcalc.findEvents(drhs, dcos, t0, dt);
                std::cout<<"Real event coordinate: "<<trueEventCoordinate<<std::endl;
                std::cout<<"Events found: "<<trajectoryStorage.size()<<std::endl;
                std::cout<<"Numerical event coordinate: "<<trajectoryStorage[0]<<std::endl;
                {
                    std::cout<<"[TEST] Must have found one crossing"<<std::endl;
                    testResult = (trajectoryStorage.size() == 1);
                    std::cout<<((testResult)?("PASSED"):("FAILED"))<<std::endl;
                }
                {
                    std::cout<<"[TEST] Crossing is close to real point"<<std::endl;
                    testResult = (std::abs(trajectoryStorage[0]-trueEventCoordinate)<evtTolerance);
                    std::cout<<((testResult)?("PASSED"):("FAILED"))<<std::endl;
                }
            }
            {
                std::cout<<"\n@Finding descending zeroes of cosine function on segment"<<std::endl;
                double crossingDirection = -1.0;
                double trueEventCoordinate = PI/2.0;
                EventCalculator<DummyRHS, DummyCos, double> evcalc(atol, rtol,
                                                       approximateReturnTime,
                                                       crossingDirection,
                                                       iterMax,evtTolerance,
                                                       timeSkip);
                auto trajectoryStorage = evcalc.findEvents(drhs, dcos, t0, dt);
                std::cout<<"Real event coordinate: "<<trueEventCoordinate<<std::endl;
                std::cout<<"Events found: "<<trajectoryStorage.size()<<std::endl;
                std::cout<<"Numerical event coordinate: "<<trajectoryStorage[0]<<std::endl;
                {
                    std::cout<<"[TEST] Must have found one crossing"<<std::endl;
                    testResult = (trajectoryStorage.size() == 1);
                    std::cout<<((testResult)?("PASSED"):("FAILED"))<<std::endl;
                }
                {
                    std::cout<<"[TEST] Crossing is close to real point"<<std::endl;
                    testResult = (std::abs(trajectoryStorage[0]-trueEventCoordinate)<evtTolerance);
                    std::cout<<((testResult)?("PASSED"):("FAILED"))<<std::endl;
                }
            }
        }
    }
    // Test: eigenvalues and Jacobi matrix precision and fixed points for Henon map
    {
        std::cout<<"\n\n[HENON MAP TESTS]"<<std::endl;
        {
            bool testResult;
            double a = 2.0, b = 0.5, c = 4.0;
            HenonMap hm(a,b,c);
            double epsOrt = 1e-7;
            double eigvPrecision = 1e-6;
            double jacFrobPrecision = 1e-6;
            Vec3 pt = {0.,0.,0.};
            NumericalJacobiMatrix<Vec3, Mat3> njm(epsOrt);
            auto F = std::bind(&HenonMap::mapping, hm, std::placeholders::_1);
            Mat3 Jac;
            njm.computeJacobiMatrix(F, pt,Jac);
            Eigen::EigenSolver<Mat3> es(Jac);
            // #################################################################
            std::cout<<"[[Eigenvalues and Jacobi]]"<<std::endl;
            {
                std::cout<<"[Test] Eigenvalues are roots of characteristic"
                           " equation with high precision"<<std::endl;
                std::array<std::complex<double>, 3> eigvals;
                for (int i = 0; i<3; i++)
                {
                    eigvals[i] = es.eigenvalues()[i];
                }
                std::cout<<"Eigenvalues are \n";
                for (int i = 0; i<3; i++)
                {
                    auto ev = eigvals[i];
                    std::cout<<"ev["<<i<<"] = "<<ev<<", "
                               "f(ev["<<i<<"]) = "<<
                               hm.valueOfCharacteristicEquation(ev)<<std::endl;
                }
                std::cout<<"Precision is "<<eigvPrecision<<std::endl;
                testResult = std::all_of(eigvals.begin(), eigvals.end(),
                                         [hm, eigvPrecision]
                                         (std::complex<double> x)
                {return (std::abs(hm.valueOfCharacteristicEquation(x))
                         <eigvPrecision);});
                std::cout<<((testResult)?("PASSED"):("FAILED"))<<std::endl;
            }
            // #################################################################

            {
                std::cout<<"[Test] Numerical Jacobi matrix is close to real "
                           "Jacobi matrix"<<std::endl;
                std::cout<<"Numerical Jacobi = \n"<<Jac<<std::endl;
                std::cout<<"Real Jacobi = \n"<<hm.jacobiMatrix(pt)<<std::endl;
                auto errFrobeniusNorm = (Jac-hm.jacobiMatrix(pt)).norm();
                std::cout<<"Frobenius norm of error: "<<errFrobeniusNorm<<std::endl;
                std::cout<<"Jacobi precision is "<<jacFrobPrecision<<std::endl;
                testResult = (errFrobeniusNorm < jacFrobPrecision );
                std::cout<<((testResult)?("PASSED"):("FAILED"))<<std::endl;
            }

            // #################################################################
            std::cout<<"[[Fixed points]]"<<std::endl;
            {
                std::cout<<"[Test] (t, t, t) where t = (A+B+C) - 1 is a fixed "
                           "point"<<std::endl;
                double t = (a+b+c - 1.0);
                Vec3 pt = {t, t, t};
                double residue = hm.residueMap(pt).norm();
                double precision = 1e-14;
                std::cout<<"Residue of HM(t, t, t): "<<residue<<std::endl;
                std::cout<<"Precision is "<<precision<<std::endl;
                testResult = (residue < precision);
                std::cout<<((testResult)?("PASSED"):("FAILED"))<<std::endl;
            }
            {
                std::cout<<"[Test] (0.5*t, 0.5*t, 0.5*t) where t = (A+B+C) - 1 is not a fixed point"<<std::endl;
                double t = 0.5*(a+b+c - 1.0);
                Vec3 pt = {t, t, t};
                double residue = hm.residueMap(pt).norm();
                double precision = 1e-14;
                std::cout<<"Residue of HM(t, t, t): "<<residue<<std::endl;
                std::cout<<"Precision is "<<precision<<std::endl;
                testResult = (residue > precision);
                std::cout<<((testResult)?("PASSED"):("FAILED"))<<std::endl;
            }
            {
                std::cout<<"\n@Start from (0.9*t, 0.9*t,0.9*t) where "
                           "t = (A+B+C) - 1 -- find a fixed point"<<std::endl;
                double t = 0.9*(a+b+c - 1.0);
                Vec3 pt = {t, t, t};
                double fixPtTol = 1e-14;
                NewtonSolver<Vec3, Mat3> ns(fixPtTol, 20, 1e-8, 1.0);
                auto G = std::bind(&HenonMap::residueMap, hm,
                                   std::placeholders::_1);
                pt = ns.performNewtonMethod(G, pt);
                {
                    std::cout<<"[TEST] Numerical fixed point meets precision for"
                               " residue"<<std::endl;
                    double residue = hm.residueMap(pt).norm();
                    std::cout<<"Residue of HM(t, t, t): "<<residue<<std::endl;
                    std::cout<<"Precision is "<<fixPtTol<<std::endl;
                    testResult = (residue < fixPtTol);
                    std::cout<<((testResult)?("PASSED"):("FAILED"))<<std::endl;
                }
                {
                    std::cout<<"[TEST] Numerical fixed point is close to real "
                               "fixed point"<<std::endl;
                    std::cout<<"Numeric fixed point: \n"<<pt<<std::endl;
                    double w = (a+b+c -1.0);
                    Vec3 fpt = {w, w, w};
                    std::cout<<"Real fixed point: \n"<<fpt<<std::endl;
                    double residue = (fpt - pt).norm();
                    std::cout<<"Norm of residue: "<<residue<<std::endl;
                    testResult = (residue < fixPtTol);
                    std::cout<<((testResult)?("PASSED"):("FAILED"))<<std::endl;
                }
                {
                    std::cout<<"[TEST] Solver returned success status"<<std::endl;
                    testResult = ns.getSuccessStatus();
                    std::cout<<((testResult)?("PASSED"):("FAILED"))<<std::endl;
                }
            }

        }
    }
    return 0;
}
