#pragma once
#include "CommonIncludes.h"
#include "Tools.h"
#include <igl/find.h>
#include "CenterOfMassEnergy.h"
#include "ElasticEnergy.h"

class Solver
{
public:
    Solver();
    Solver(string matdir);
    ~Solver();
    MatrixXd genReductionMatrix() const; // We encode PBC in the reduction matrix

    // Set coarse
    void setSCAKD(double s, double c, double a, double k, double d); // Set the coarse configuration

    // Energy and derivatives
    double getEnergy(const MatrixXd &Xf) const;
    VectorXd getGradient(const MatrixXd &Xf) const;
    SparseMatrix<double> getHessian(const MatrixXd &Xf) const;
    SparseMatrix<double> getHessianClamped(const MatrixXd &Xf) const;

    // Periodicity
    double getEnergyPeriodic(const VectorXd &vXfRED) const;
    VectorXd getGradientPeriodic(const VectorXd &vXfRED) const;
    SparseMatrix<double> getHessianPeriodic(const VectorXd &vXfRED) const;
    SparseMatrix<double> getHessianClampedPeriodic(const VectorXd &vXfRED) const;

    // numerical methods
    VectorXd NewtonRaphson(const VectorXd &vXf, bool &success) const;
    VectorXd NewtonRaphson(const VectorXd &vXf, double &errorConvergence) const;
    VectorXd solveLinearSystem(SparseMatrix<double> &A, VectorXd &b) const;
    VectorXd linsearch(const VectorXd &vX, const VectorXd &vdx, bool &linesearch_success) const;

    // getter
    double getVolume() const;
    double getElasticEnergy(const MatrixXd &Xf) const;
    double getElasticEnergy(const MatrixXd &Xf, double s, double c, double a, double k, double d) const;
    double getMinimalEnergy(double s, double c, double a, double k, double d, bool &success);
    double getRestEnergy() const;
    double getMinimalEnergyDebug(double s, double c, double a, double k, double d, bool &success);

    MatrixXd relax(double s, double c, double a, double k, double d, double &errorConvergence);

protected:
    string m_matdir;
    MatrixXi m_T, m_PBCDoFs;
    MatrixXd m_Xrest, m_redMatrix, m_invRedMatrix;
    SparseMatrix<double> m_SparseRedMat;
    int m_nbNodes, m_nbElements, m_nbDoFs;
    CenterOfMassEnergy m_come;
    ElasticEnergy m_ee;
    double m_s, m_c, m_a, m_k, m_d;
};
