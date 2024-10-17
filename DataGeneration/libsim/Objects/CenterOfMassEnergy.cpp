#include "CenterOfMassEnergy.h"

CenterOfMassEnergy::CenterOfMassEnergy()
{
}

CenterOfMassEnergy::~CenterOfMassEnergy()
{
}

void CenterOfMassEnergy::setXip(const MatrixXd &Xip)
{
    m_Xip = Xip;
    int nbDoFs = Xip.rows() * Xip.cols();
    m_nbDoFs = nbDoFs;
    m_C = MatrixXd::Zero(3, nbDoFs);
    m_C(0, seq(0, last, 3)) = MatrixXd::Ones(1, nbDoFs / 3);
    m_C(1, seq(1, last, 3)) = MatrixXd::Ones(1, nbDoFs / 3);
    m_C(2, seq(2, last, 3)) = MatrixXd::Ones(1, nbDoFs / 3);
    m_Hess = (m_coefCentMass * m_C.transpose() * m_C).sparseView();
}

double CenterOfMassEnergy::getEnergy(const MatrixXd &X) const
{
    VectorXd vX = Tools::vectorize(X);
    VectorXd vXrest = Tools::vectorize(m_Xip);
    double res = 0;
    res += 0.5 * (m_C * vX).transpose() * (m_C * vX);
    res += -0.5 * (m_C * vXrest).transpose() * (m_C * vX);
    res += -0.5 * (m_C * vX).transpose() * (m_C * vXrest);
    res += 0.5 * (m_C * vXrest).transpose() * (m_C * vXrest);
    return m_coefCentMass * res;
}

VectorXd CenterOfMassEnergy::getGradient(const MatrixXd &X) const
{
    VectorXd vX = Tools::vectorize(X);
    VectorXd vXrest = Tools::vectorize(m_Xip);
    VectorXd res = VectorXd::Zero(m_nbDoFs);
    res += (m_C * vX).transpose() * m_C;
    res += -(m_C * vXrest).transpose() * m_C;
    return m_coefCentMass * res;
}

SparseMatrix<double> CenterOfMassEnergy::getHessian(const MatrixXd &X) const
{
    return m_Hess;
}

SparseMatrix<double> CenterOfMassEnergy::getHessian() const
{
    return m_Hess;
}