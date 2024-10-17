#pragma once
#include "CommonIncludes.h"
#include "Tools.h"
class CenterOfMassEnergy
{
public:
    CenterOfMassEnergy();
    ~CenterOfMassEnergy();

    void setXip(const MatrixXd &Xip);
    double getEnergy(const MatrixXd &X) const;
    VectorXd getGradient(const MatrixXd &X) const;
    SparseMatrix<double> getHessian(const MatrixXd &X) const;
    SparseMatrix<double> getHessian() const;

protected:
    double m_coefCentMass = 1;
    MatrixXd m_C, m_Xip;
    int m_nbDoFs;
    SparseMatrix<double> m_Hess;
};
