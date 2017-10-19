#ifndef NUM_JAC_H
#define NUM_JAC_H

#include <functional>

template<class PtType, class MatType>
class NumericalJacobiMatrix
{
private:
    double epsOrt;
public:
    NumericalJacobiMatrix(double eO): epsOrt(eO) {}
    void computeJacobiMatrix(const std::function<PtType(PtType)>& F,
                        const PtType& pt, MatType& jacobiMatrix) const;
};

template<class PtType, class MatType>
void NumericalJacobiMatrix<PtType, MatType>::computeJacobiMatrix(const std::function<PtType(PtType)>& F,
                                                const PtType& pt, MatType& jacobiMatrix) const
{
    MatType Jac;
    std::array<PtType, 3> orts = {{{1.,0.,0.},{0.,1.,0.},{0.,0.,1.}}};
    PtType F0 = F(pt);
    for (int i = 0; i < 3; i++)
    {
        Jac.col(i) = (F(pt + epsOrt * orts[i]) - F0)/epsOrt;
    }
    jacobiMatrix = Jac;
}

#endif // NUM_JAC_H
