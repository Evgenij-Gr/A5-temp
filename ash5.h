#ifndef ASH5
#define ASH5

#include <boost/version.hpp>

#if BOOST_VERSION == 106400
#include <boost/serialization/array_wrapper.hpp>
#endif

#include <boost/array.hpp>

typedef boost::array<double, 4> AshStateType;

class Ashwin5osc
{
private:
    double r;
    double a;
    double b;
public:
    Ashwin5osc(double R, double A, double B): r(R), a(A), b(B){}
    void operator()(const AshStateType &X , AshStateType &dxdt , double t );
};

// this Poincare section is taken from the Ashwin-Orosz-Townley-Wordsworth-2007
class Ash_Section_Event
{
public:
    double operator()(const AshStateType& X);
    AshStateType gradient(const AshStateType& X);
    void adjustToSection(AshStateType& X);
};

#endif
