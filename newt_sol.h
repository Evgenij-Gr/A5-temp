#ifndef NEWT_SOL_H
#define NEWT_SOL_H

#include "num_jac.h"

template<class PtType, class MatType>
PtType solveSLE(const MatType& A, const PtType& b)
{
    PtType result = A.colPivHouseholderQr().solve(b);
    return result;
}

template<class PtType>
double l2norm(const PtType& pt)
{
    return pt.norm();
}

template<class PtType, class MatType>
class NewtonSolver
{
private:
    double fixPtTol;
    double epsOrt;
    bool successStatus;
    int maxIterNum;
    double epsStep;
public:
    NewtonSolver(const double &fPT, const int &mIN, const double &jNE, const double &eS)
    {
        fixPtTol = fPT;
        maxIterNum = mIN;
        epsOrt = jNE;
        epsStep = eS;
    }
    PtType performNewtonMethod(const std::function<PtType(PtType)>& F, const PtType& initPt)
    {
        PtType pt = initPt;
        PtType residue = F(pt);
        int iterNum = 0;
        while ((l2norm(residue) > fixPtTol) && (iterNum < maxIterNum))
        {
            ++iterNum;
            MatType Jac;
            NumericalJacobiMatrix<PtType, MatType> jmat(epsOrt);
            jmat.computeJacobiMatrix(F, pt, Jac);
            PtType delta = solveSLE(Jac, PtType(-residue));
            pt = pt + epsStep*delta;
            residue = F(pt);
        }
        successStatus = (l2norm(residue) < fixPtTol);
        return pt;
    }
    bool getSuccessStatus()
    {
        return successStatus;
    }
};

#endif // NEWT_SOL_H
