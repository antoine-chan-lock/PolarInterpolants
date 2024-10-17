#include "Solver.h"

Solver::Solver()
{
}

Solver::~Solver()
{
}

Solver::Solver(string matdir)
{
    m_matdir = matdir;
    m_T = Tools::loadCsvInt(matdir + "/T.csv");
    m_Xrest = Tools::loadCsvDouble(matdir + "/X.csv");
    m_nbNodes = m_Xrest.rows();
    m_nbElements = m_T.rows();
    m_nbDoFs = m_nbNodes * 3;
    m_PBCDoFs = Tools::loadCsvInt(matdir + "/PBC.csv");
    m_redMatrix = genReductionMatrix();
    m_SparseRedMat = m_redMatrix.sparseView();
    assert(!((m_redMatrix.transpose() * m_redMatrix).determinant() < 1e-10) && "reduction matrix is not invertible");
    m_invRedMatrix = (m_redMatrix.transpose() * m_redMatrix).inverse() * m_redMatrix.transpose();
    m_ee = ElasticEnergy(m_Xrest, m_T);
    m_come.setXip(m_Xrest);
}

MatrixXd Solver::genReductionMatrix() const
{
    int nbDoFs = m_nbDoFs;
    MatrixXi PBC2 = m_PBCDoFs;

    Matrix<bool, Dynamic, Dynamic> res(nbDoFs, nbDoFs);
    res.setIdentity();

    for (int i = 0; i < PBC2.rows(); i++)
    {
        res(PBC2(i, 0), PBC2(i, 1)) = 1;
        res(PBC2(i, 1), PBC2(i, 0)) = 1;
    }
    // collapse cols
    bool collapsed = false;
    while (!collapsed)
    {
        collapsed = true;
        for (int i = 0; i < res.rows(); i++)
        {
            VectorXi I;
            igl::find(res.row(i), I);
            if (I.rows() > 1)
            {
                collapsed = false;
                for (int j = 0; j < I.size(); j++)
                {
                    res.col(I(0)) = res.col(I(0)) || res.col(I(j));
                }
                for (int j = 1; j < I.size(); j++)
                {
                    res.col(I(j)) = res.col(I(j)) * 0;
                }
            }
        }
    }
    VectorXi I;
    igl::find(res.colwise().sum(), I);
    MatrixXd res2(res.rows(), I.size());
    res2.setZero();
    res2 = (res(all, I)).select(1, res2);

    return res2;
}

void Solver::setSCAKD(double s, double c, double a, double k, double d)
{
    m_s = s;
    m_c = c;
    m_a = a;
    m_k = k;
    m_d = d;

    double ener = m_ee.getEnergy(m_Xrest, m_s, m_c, m_a, m_k, m_d);
    if (isnan(ener))
    {
        cout << "Solver::setSCAKD::nan:: INVERSION AT START" << endl;
        exit(0);
    }
}

double Solver::getEnergy(const MatrixXd &Xf) const
{
    MatrixXd X = m_Xrest + Xf;

    double elasticEnergy = m_ee.getEnergy(X, m_s, m_c, m_a, m_k, m_d);
    double centerMassEnergy = m_come.getEnergy(X);
    return elasticEnergy + centerMassEnergy;
}

double Solver::getElasticEnergy(const MatrixXd &Xf) const
{
    MatrixXd X = m_Xrest + Xf;
    return m_ee.getEnergy(X, m_s, m_c, m_a, m_k, m_d);
}

double Solver::getElasticEnergy(const MatrixXd &Xf, double s, double c, double a, double k, double d) const
{
    MatrixXd X = m_Xrest + Xf;
    return m_ee.getEnergy(X, s, c, a, k, d);
}

VectorXd Solver::getGradient(const MatrixXd &Xf) const
{
    MatrixXd X = m_Xrest + Xf;
    VectorXd elasticGradient = m_ee.getGradient(X, m_s, m_c, m_a, m_k, m_d);
    VectorXd centerMassGradient = m_come.getGradient(X);
    return elasticGradient + centerMassGradient;
}

SparseMatrix<double> Solver::getHessian(const MatrixXd &Xf) const
{
    MatrixXd X = m_Xrest + Xf;
    SparseMatrix<double> COMEHess = m_come.getHessian(X);
    SparseMatrix<double> EEHess = m_ee.getHessian(X, m_s, m_c, m_a, m_k, m_d);
    return EEHess + COMEHess;
}

SparseMatrix<double> Solver::getHessianClamped(const MatrixXd &Xf) const
{

    MatrixXd X = m_Xrest + Xf;
    SparseMatrix<double> COMEHess = m_come.getHessian(X);
    SparseMatrix<double> EEHess = m_ee.getHessianClamped(X, m_s, m_c, m_a, m_k, m_d);

    return EEHess + COMEHess;
}

VectorXd Solver::linsearch(const VectorXd &vX, const VectorXd &vdx, bool &linesearch_success) const
{
    double energy = getEnergyPeriodic(vX);
    VectorXd vXguess = vX - vdx;
    double energyGuess = getEnergyPeriodic(vXguess);
    energyGuess = isnan(energyGuess) ? 1e14 : energyGuess;
    linesearch_success = true;
    VectorXd vdxTMP = vdx;
    int iter = 0;
    double coef = 1;
    while ((energyGuess > energy) and (iter < 100))
    {
        vdxTMP = vdxTMP / 2;
        coef = coef / 2;
        vXguess = vX - vdxTMP;
        energyGuess = getEnergyPeriodic(vXguess);
        energyGuess = isnan(energyGuess) ? 1e14 : energyGuess;
        iter += 1;
    }
    if (iter >= 100)
    {
        linesearch_success = false;
        return vX;
    }
    return vXguess;
}

VectorXd Solver::NewtonRaphson(const VectorXd &vXf, bool &success) const
{

    success = true;
    int maxIter = 1e3, iter = 0;
    double maxError = 1e-2;
    VectorXd vXfRED = m_invRedMatrix * vXf;
    int nbDoFsRED = vXfRED.rows();
    bool linesearch_success = true;
    double error = 1e2;

    while ((iter < maxIter) and (error > maxError) and success)
    {
        double enerBegin = getEnergyPeriodic(vXfRED);
        VectorXd dpsi = getGradientPeriodic(vXfRED);
        error = dpsi.norm();
        SparseMatrix<double> ddpsi = getHessianClampedPeriodic(vXfRED);
        VectorXd vdxRED = solveLinearSystem(ddpsi, dpsi);
        vXfRED = linsearch(vXfRED, vdxRED, success);
        double enerEnd = getEnergyPeriodic(vXfRED);
        if (enerEnd > enerBegin)
        {
            cout << "Solver::NR::didnt manage to lower energy" << endl;
            exit(0);
        }
        iter += 1;
    }
    if (iter >= maxIter)
    {
        success = false;
    }
    return m_redMatrix * vXfRED;
}

VectorXd Solver::NewtonRaphson(const VectorXd &vXf, double &errorConvergence) const
{
    bool success = true;
    int maxIter = 1e3, iter = 0;
    double maxError = 1e-2;
    VectorXd vXfRED = m_invRedMatrix * vXf;
    int nbDoFsRED = vXfRED.rows();
    double error = 1e2;
    while ((iter < maxIter) and (error > maxError) and success)
    {
        double enerBegin = getEnergyPeriodic(vXfRED);
        VectorXd dpsi = getGradientPeriodic(vXfRED);
        error = dpsi.norm();
        SparseMatrix<double> ddpsi = getHessianClampedPeriodic(vXfRED);
        VectorXd vdxRED = solveLinearSystem(ddpsi, dpsi);
        vXfRED = linsearch(vXfRED, vdxRED, success);
        double enerEnd = getEnergyPeriodic(vXfRED);
        if (enerEnd > enerBegin)
        {
            cout << "Solver::NR::didnt manage to lower energy" << endl;
            exit(0);
        }
        iter += 1;
    }
    if (iter >= maxIter)
    {
        success = false;
    }
    errorConvergence = error;
    return m_redMatrix * vXfRED;
}

double Solver::getMinimalEnergy(double s, double c, double a, double k, double d, bool &success)
{
    setSCAKD(s, c, a, k, d);
    VectorXd vXf = VectorXd::Zero(m_nbDoFs);
    vXf = NewtonRaphson(vXf, success);
    return getElasticEnergy(Tools::devectorize(vXf));
}

MatrixXd Solver::relax(double s, double c, double a, double k, double d, double &errorConvergence)
{
    setSCAKD(s, c, a, k, d);
    VectorXd vXf = VectorXd::Zero(m_nbDoFs);
    vXf = NewtonRaphson(vXf, errorConvergence);
    return Tools::devectorize(vXf);
}

double Solver::getRestEnergy() const
{
    return m_ee.getEnergy(m_Xrest, m_s, m_c, m_a, m_k, m_d);
}

double Solver::getMinimalEnergyDebug(double s, double c, double a, double k, double d, bool &success)
{
    setSCAKD(s, c, a, k, d);
    VectorXd vXf = VectorXd::Zero(m_nbDoFs);
    cout << "grad: " << getGradient(Tools::devectorize(vXf)).norm() << endl;
    vXf = NewtonRaphson(vXf, success);
    cout << "gradSol: " << getGradient(Tools::devectorize(vXf)).norm() << endl;
    return getEnergy(Tools::devectorize(vXf));
}

double Solver::getVolume() const
{
    return m_ee.getVolume();
}

VectorXd Solver::solveLinearSystem(SparseMatrix<double> &A, VectorXd &b) const
{
    // VectorXd res;
    // SparseMatrix<double, RowMajor> ARowMajor(A);
    // //1st strategy: BiCGSTAB
    // BiCGSTAB<SparseMatrix<double, RowMajor>, IncompleteLUT<double>> solver;
    // solver.compute(ARowMajor);
    // res = solver.solve(b);
    // double error = (A * res - b).norm();
    // if (isnan(res.norm()) or (error > 1e-6))
    // {
    //     //2nd strategy: SparseQR
    //     Eigen::SparseQR<SparseMatrix<double>, COLAMDOrdering<int>> solver3;
    //     A.makeCompressed();
    //     solver3.compute(A);
    //     if (solver3.info() == Eigen::Success)
    //     {
    //         res = solver3.solve(b);
    //     }
    //     else
    //     {
    //         cout << "Solver::solveLinearSystem::fail" << endl;
    //         exit(0);
    //     }
    // }

    VectorXd res;
    Eigen::SparseQR<SparseMatrix<double>, COLAMDOrdering<int>> solver3;
    A.makeCompressed();
    solver3.compute(A);
    if (solver3.info() == Eigen::Success)
    {
        res = solver3.solve(b);
    }

    if (isnan(res.norm()))
    {
        cout << "Solver::solveLinearSystem::nan" << endl;
        exit(0);
    }

    return res;
}

double Solver::getEnergyPeriodic(const VectorXd &vXfRED) const
{
    return getEnergy(Tools::devectorize(m_redMatrix * vXfRED));
}
VectorXd Solver::getGradientPeriodic(const VectorXd &vXfRED) const
{
    return m_redMatrix.transpose() * getGradient(Tools::devectorize(m_redMatrix * vXfRED));
}
SparseMatrix<double> Solver::getHessianPeriodic(const VectorXd &vXfRED) const
{
    return m_SparseRedMat.transpose() * getHessian(Tools::devectorize(m_redMatrix * vXfRED)) * m_SparseRedMat;
}
SparseMatrix<double> Solver::getHessianClampedPeriodic(const VectorXd &vXfRED) const
{
    return m_SparseRedMat.transpose() * getHessianClamped(Tools::devectorize(m_redMatrix * vXfRED)) * m_SparseRedMat;
}
