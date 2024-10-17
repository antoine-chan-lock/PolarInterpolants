#pragma once
#include "CommonIncludes.h"
#include "Tango.h"
#include "Tools.h"
class Hexahedron
{
public:
    Hexahedron();
    Hexahedron(int id, const MatrixXi &connectivity, const MatrixXd &Xrest);
    ~Hexahedron();

    // Full Energy and derivatives
    double getEnergy(const MatrixXd &Xcurrent, double s, double c, double a, double k, double d) const;
    Matrix<double, 1, 24> getGradient(const MatrixXd &Xcurrent, double s, double c, double a, double k, double d) const;
    Matrix<double, 24, 24> getHessian(const MatrixXd &Xcurrent, double s, double c, double a, double k, double d) const;
    Matrix<double, 24, 24> getHessianClamped(const MatrixXd &Xcurrent, double s, double c, double a, double k, double d) const;

    // Elastic energy and derivatives
    double getElasticEnergy(const Matrix<double, 3, 8> &Xcurrent) const;
    Matrix<double, 1, 24> getElasticGradient(const Matrix<double, 3, 8> &Xcurrent) const;
    Matrix<double, 24, 24> getElasticHessian(const Matrix<double, 3, 8> &Xcurrent) const;

    // Bending and derivatives
    Matrix<double, 3, 8> bend(const Matrix<double, 3, 8> &X, double s, double c, double a, double k, double d) const;
    Matrix<double, 3 * 8, 3 * 8> dbend(const Matrix<double, 3, 8> &X, double s, double c, double a, double k, double d) const;
    vector<Matrix<double, 3 * 8, 3 * 8>> ddbend(const Matrix<double, 3, 8> &X, double s, double c, double a, double k, double d) const;

    // clamp Hessian
    Matrix<double, 3 * 8, 3 * 8> clamp(const Matrix<double, 3 * 8, 3 * 8> &H) const;
    // Getter
    double getVolume() const;
    Matrix<int, 1, 8 * 3> getIndexDoFs() const;

    int m_id, m_nbNodes, m_nbDofs;
    double m_detJac[8];
    double m_volume;
    Matrix<double, 8, 3> m_dN[8];
    Matrix<double, 8, 3> m_dFdx[8];
    Matrix<double, 9, 24> m_dvFdvx[8];
    Matrix<int, 1, 8> m_nodes;
    Matrix<double, 3, 8> m_Xrest;
    Tango m_material;
};
