#include "ash5.h"

using std::sin;
using std::cos;

void Ashwin5osc::operator()(const AshStateType &X , AshStateType &dxdt , double t )
{
    double gam1, gam2, gam3, gam4;
    gam1 = X[0];
    gam2 = X[1];
    gam3 = X[2];
    gam4 = X[3];
    dxdt[0] = -0.2*r*sin( -2.0*gam1+b)+0.2*r*sin( 2.0*gam1+b)-0.2*sin(-gam3+gam1+a)+0.2*sin(-gam3+a)+0.2*sin(-gam1+a)-0.2*r*sin( -2.0*gam4+b)+0.2*r*sin( -2.0*gam3+2.0*gam1+b)+0.2*r*sin( -2.0*gam4+2.0*gam1+b)-0.2*r*sin( -2.0*gam2+b)-0.2*sin( -2.0*gam3+b)*r-0.2*sin(-gam4+gam1+a)-0.2*sin( gam1+a)-0.2*sin( gam1+a-gam2)+0.2*sin( a-gam2)+0.2*sin(-gam4+a)+0.2*r*sin( 2.0*gam1-2.0*gam2+b);
    dxdt[1] = 0.2*r*sin( -2.0*gam4+2.0*gam2+b)-0.2*r*sin( -2.0*gam1+b)+0.2*r*sin( -2.0*gam3+2.0*gam2+b)-0.2*sin(-gam1+a+gam2)-0.2*sin( a+gam2)+0.2*sin(-gam3+a)+0.2*sin(-gam1+a)-0.2*sin(-gam4+a+gam2)-0.2*r*sin( -2.0*gam4+b)-0.2*r*sin( -2.0*gam2+b)+0.2*r*sin( 2.0*gam2+b)-0.2*sin( -2.0*gam3+b)*r+0.2*sin( -2.0*gam1+2.0*gam2+b)*r-0.2*sin(-gam3+a+gam2)+0.2*sin( a-gam2)+0.2*sin(-gam4+a);
    dxdt[2] = -0.2*sin( gam3+a-gam2)+0.2*r*sin( 2.0*gam3-2.0*gam4+b)-0.2*r*sin( -2.0*gam1+b)-0.2*sin( gam3-gam4+a)+0.2*r*sin( 2.0*gam3-2.0*gam2+b)+0.2*sin(-gam3+a)+0.2*sin(-gam1+a)-0.2*r*sin( -2.0*gam4+b)+0.2*r*sin( 2.0*gam3-2.0*gam1+b)-0.2*r*sin( -2.0*gam2+b)-0.2*sin( -2.0*gam3+b)*r-0.2*sin( gam3-gam1+a)-0.2*sin( gam3+a)+0.2*sin( a-gam2)+0.2*r*sin( 2.0*gam3+b)+0.2*sin(-gam4+a);
    dxdt[3] = 0.2*r*sin( 2.0*gam4-2.0*gam2+b)+0.2*r*sin( -2.0*gam3+2.0*gam4+b)-0.2*r*sin( -2.0*gam1+b)-0.2*sin( gam4-gam1+a)+0.2*sin(-gam3+a)+0.2*sin(-gam1+a)+0.2*r*sin( 2.0*gam4+b)-0.2*r*sin( -2.0*gam4+b)+0.2*r*sin( 2.0*gam4-2.0*gam1+b)-0.2*sin( gam4+a)-0.2*r*sin( -2.0*gam2+b)-0.2*sin( gam4+a-gam2)-0.2*sin( -2.0*gam3+b)*r-0.2*sin(-gam3+gam4+a)+0.2*sin( a-gam2)+0.2*sin(-gam4+a);
}

// this Poincare section is taken from the Ashwin-Orosz-Townley-Wordsworth-2007
double Ash_Section_Event::operator()(const AshStateType& X)
{
    double gam1, gam2, gam3, gam4;
    gam1 = X[0];
    gam2 = X[1];
    gam3 = X[2];
    gam4 = X[3];
    return (gam3-1.2);
}
AshStateType Ash_Section_Event::gradient(const AshStateType& X)
{
    AshStateType ret;
    ret[0] = 0.0;
    ret[1] = 0.0;
    ret[2] = 1.0;
    ret[3] = 0.0;
    return ret;
}
void Ash_Section_Event::adjustToSection(AshStateType& X)
{
    X[2] = 1.2;
}
