#ifndef MISC_DEFS_H
#define MISC_DEFS_H

#include "newt_sol.h"
#include "ash5.h"
#include "dyn_utils.h"
#include <Eigen/Dense>

// This dummy system is needed for emulating
// real functions and searching for their zeroes
// just corresponds to scalar equation x' = 1
class DummyRHS
{
public:
    void operator()(const double& X, double& dxdt, double t)
    {
        dxdt = 1;
    }
};

// ///////////////////////////////////////////////

#include <boost/version.hpp>

#if BOOST_VERSION == 106400
#include <boost/serialization/array_wrapper.hpp>
#endif

#include <boost/array.hpp>

typedef boost::array<double, 2> HarmOscType;

class HarmonicOscillator
{
private:
    double omega;
public:
    HarmonicOscillator(const double& omega_)
    {
        omega = omega_;
    }
    void operator()(const HarmOscType& X, HarmOscType& dxdt, double t)
    {
        auto x = X[0];
        auto y = X[1];
        dxdt[0] = y;
        dxdt[1] = -omega*omega*x;
    }
};

class HarmonicOscillatorEvent
{
public:
    double operator()(const HarmOscType& X)
    {
        auto y = X[1];
        return y;
    }
    void adjustToSection(HarmOscType& X)
    {
        X[1] = 0.0;
    }
};

// ////////////////////////////////////////////////////

typedef Eigen::Vector3d Vec3;
typedef Eigen::Matrix3d Mat3;

class HenonMap
{
private:
    double A, B, C;
public:
    HenonMap(const double& a, const double& b, const double& c)
    {
        A = a;
        B = b;
        C = c;
    }
    std::vector<Vec3> computeTrajectory(const Vec3& startPt, const int& iterMax=1)
    {
        int maxEvtsCapacity = iterMax+1;
        std::vector<Vec3> trajectoryStorage;
        trajectoryStorage.reserve(maxEvtsCapacity);
        auto pt = startPt;
        for (int i=0; i <= iterMax; i++)
        {
            trajectoryStorage.push_back(pt);
            pt = mapping(pt);
            if (l2norm(pt) > 1000.0)
            {
                break;
            }
        }
        return trajectoryStorage;
    }
    Vec3 mapping(const Vec3& pt)
    {
        auto x = pt[0];
        auto y = pt[1];
        auto z = pt[2];
        return {y, z, B*x + C*y + A*z - y*y};
    }
    Vec3 residueMap(const Vec3& pt)
    {
        return (mapping(pt)-pt);
    }
    Mat3 jacobiMatrix(const Vec3& pt)
    {
        const auto y = pt[1];
        Mat3 realJac;
        //[0 1 0]
        //[0 0 1]
        //[B C-2y A]
        realJac.col(0) = Vec3({0, 0, B});
        realJac.col(1) = Vec3({1, 0, C - 2*y});
        realJac.col(2) = Vec3({0, 1, A});
        return realJac;
    }
    double valueOfCharacteristicEquation(const std::complex<double> &E)
    {
        //Jacobi matrix at origin is
        //[0 1 0]
        //[0 0 1]
        //[B C A]
        //Characteristic equation is
        // -B - C*x - A*x^2 + x^3 = 0
        return std::abs(-B-C*E-A*E*E+E*E*E);
    }
};

#include <boost/numeric/odeint.hpp>
using namespace boost::numeric::odeint;

template<class System, class Event, class TStateType>
class EventCalculator
{
private:
    double atol;
    double rtol;
    double approximateReturnTime;
    double crossingDirection;
    int iterMax;
    double eventTolerance;
    double timeSkip;
public:
    EventCalculator(double atol_, double rtol_, double approximateReturnTime_,
                    double crossingDirection_, int iterMax_,
                    double eventTolerance_, double timeSkip_)
    {
        atol = atol_;
        rtol = rtol_;
        approximateReturnTime = approximateReturnTime_;
        crossingDirection = crossingDirection_;
        iterMax = iterMax_;
        eventTolerance = eventTolerance_;
        timeSkip = timeSkip_;
    }
    std::vector<TStateType> findEvents(System S, Event E, TStateType startPt)
    {
        std::vector<TStateType> trajectoryStorage;
        int maxEvtsCapacity = iterMax + 1;
        trajectoryStorage.reserve(maxEvtsCapacity);
        auto dopri5dense = make_dense_output(atol,
                                             rtol,
                                             runge_kutta_dopri5<TStateType>());
        auto observer = EventObserver<Event,decltype(dopri5dense),
                System,TStateType>(S,E,dopri5dense, trajectoryStorage,
                                   crossingDirection, eventTolerance,
                                   timeSkip);
        int iter = 0;
        double startTime = 0.0;
        auto pt = startPt;
        while ((trajectoryStorage.size() < maxEvtsCapacity) &&
               (iter < 2*maxEvtsCapacity))
        {
            // since the system is autonomous we can do this trick
            integrate_adaptive(dopri5dense, S, pt, startTime,
                               startTime+approximateReturnTime,
                               1e-12, observer);
            iter++;
            startTime += approximateReturnTime;
        }
        return trajectoryStorage;
    }
};

#endif // MISC_DEFS_H
