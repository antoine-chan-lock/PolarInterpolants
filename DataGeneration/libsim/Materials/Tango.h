#pragma once
#include "CommonIncludes.h"
class Tango
{
public:
    Tango();
    virtual ~Tango();
    virtual void init(double G0, double K0);
    virtual double psi(const Matrix<double, 3, 3> &F) const;
    virtual Matrix<double, 1, 9> dpsidF(const Matrix<double, 3, 3> &F) const;
    virtual Matrix<double, 9, 9> dpsidFF(const Matrix<double, 3, 3> &F) const;

protected:
    double m_G0, m_K0;
};
