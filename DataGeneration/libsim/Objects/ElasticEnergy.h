#pragma once
#include "CommonIncludes.h"
#include "Hexahedron.h"
class ElasticEnergy
{
public:
    ElasticEnergy();
    ElasticEnergy(const MatrixXd &X, const MatrixXi &T);
    ~ElasticEnergy();

    double getEnergy(const MatrixXd &X, double s, double c, double a, double k, double d) const;
    VectorXd getGradient(const MatrixXd &X, double s, double c, double a, double k, double d) const;
    SparseMatrix<double> getHessian(const MatrixXd &X, double s, double c, double a, double k, double d) const;
    SparseMatrix<double> getHessianClamped(const MatrixXd &X, double s, double c, double a, double k, double d) const;
    double getVolume() const;

    vector<Hexahedron> m_velement;
    int m_nbElements, m_nbDoFs;
};
