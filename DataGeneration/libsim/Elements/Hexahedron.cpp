#include "Hexahedron.h"
#include "bend.inc"
#define THRESHOLD 1e-9

Hexahedron::Hexahedron()
{
}

Hexahedron::Hexahedron(int id, const MatrixXi &connectivity, const MatrixXd &Xrest)
{
    m_id = id;
    m_nodes = connectivity.row(id);
    m_material.init(1.7 * 1e6, 84.43 * 1e6);
    m_Xrest = Xrest(m_nodes, all).transpose();
    m_nbNodes = 8;
    m_nbDofs = m_nbNodes * 3;

    Matrix<double, 8, 3> quadCoord;
    double pos = 1.0 / sqrt(3);
    quadCoord << -pos, -pos, -pos,
        -pos, -pos, pos,
        -pos, pos, -pos,
        -pos, pos, pos,
        pos, -pos, -pos,
        pos, -pos, pos,
        pos, pos, -pos,
        pos, pos, pos;

    m_volume = 0.0;

    for (int i = 0; i < 8; i++)
    {
        m_dN[i](0, 0) = -(1 - quadCoord(i, 1)) * (1 - quadCoord(i, 2));
        m_dN[i](1, 0) = (1 - quadCoord(i, 1)) * (1 - quadCoord(i, 2));
        m_dN[i](2, 0) = (1 + quadCoord(i, 1)) * (1 - quadCoord(i, 2));
        m_dN[i](3, 0) = -(1 + quadCoord(i, 1)) * (1 - quadCoord(i, 2));
        m_dN[i](4, 0) = -(1 - quadCoord(i, 1)) * (1 + quadCoord(i, 2));
        m_dN[i](5, 0) = (1 - quadCoord(i, 1)) * (1 + quadCoord(i, 2));
        m_dN[i](6, 0) = (1 + quadCoord(i, 1)) * (1 + quadCoord(i, 2));
        m_dN[i](7, 0) = -(1 + quadCoord(i, 1)) * (1 + quadCoord(i, 2));

        m_dN[i](0, 1) = (1 - quadCoord(i, 0)) * (-1) * (1 - quadCoord(i, 2));
        m_dN[i](1, 1) = (1 + quadCoord(i, 0)) * (-1) * (1 - quadCoord(i, 2));
        m_dN[i](2, 1) = (1 + quadCoord(i, 0)) * (1) * (1 - quadCoord(i, 2));
        m_dN[i](3, 1) = (1 - quadCoord(i, 0)) * (1) * (1 - quadCoord(i, 2));
        m_dN[i](4, 1) = (1 - quadCoord(i, 0)) * (-1) * (1 + quadCoord(i, 2));
        m_dN[i](5, 1) = (1 + quadCoord(i, 0)) * (-1) * (1 + quadCoord(i, 2));
        m_dN[i](6, 1) = (1 + quadCoord(i, 0)) * (1) * (1 + quadCoord(i, 2));
        m_dN[i](7, 1) = (1 - quadCoord(i, 0)) * (1) * (1 + quadCoord(i, 2));

        m_dN[i](0, 2) = (1 - quadCoord(i, 0)) * (1 - quadCoord(i, 1)) * (-1);
        m_dN[i](1, 2) = (1 + quadCoord(i, 0)) * (1 - quadCoord(i, 1)) * (-1);
        m_dN[i](2, 2) = (1 + quadCoord(i, 0)) * (1 + quadCoord(i, 1)) * (-1);
        m_dN[i](3, 2) = (1 - quadCoord(i, 0)) * (1 + quadCoord(i, 1)) * (-1);
        m_dN[i](4, 2) = (1 - quadCoord(i, 0)) * (1 - quadCoord(i, 1)) * (1);
        m_dN[i](5, 2) = (1 + quadCoord(i, 0)) * (1 - quadCoord(i, 1)) * (1);
        m_dN[i](6, 2) = (1 + quadCoord(i, 0)) * (1 + quadCoord(i, 1)) * (1);
        m_dN[i](7, 2) = (1 - quadCoord(i, 0)) * (1 + quadCoord(i, 1)) * (1);
        m_dN[i] = m_dN[i] / 8;

        m_dFdx[i] = m_dN[i] * (m_Xrest * m_dN[i]).inverse();
        Matrix<double, 3, 8> dFdxT = m_dFdx[i].transpose();

        m_dvFdvx[i].setZero();

        m_dvFdvx[i](0, seq(0, last, 3)) = dFdxT.row(0);
        m_dvFdvx[i](1, seq(1, last, 3)) = dFdxT.row(0);
        m_dvFdvx[i](2, seq(2, last, 3)) = dFdxT.row(0);

        m_dvFdvx[i](3, seq(0, last, 3)) = dFdxT.row(1);
        m_dvFdvx[i](4, seq(1, last, 3)) = dFdxT.row(1);
        m_dvFdvx[i](5, seq(2, last, 3)) = dFdxT.row(1);

        m_dvFdvx[i](6, seq(0, last, 3)) = dFdxT.row(2);
        m_dvFdvx[i](7, seq(1, last, 3)) = dFdxT.row(2);
        m_dvFdvx[i](8, seq(2, last, 3)) = dFdxT.row(2);

        Matrix<double, 3, 3> jac = m_Xrest * m_dN[i];
        m_detJac[i] = jac.determinant();
        m_volume += m_detJac[i];
    }
}

Hexahedron::~Hexahedron()
{
}

Matrix<int, 1, 8 * 3> Hexahedron::getIndexDoFs() const
{
    Matrix<int, 1, 8 * 3> res;
    res << 3 * m_nodes(0), 3 * m_nodes(0) + 1, 3 * m_nodes(0) + 2,
        3 * m_nodes(1), 3 * m_nodes(1) + 1, 3 * m_nodes(1) + 2,
        3 * m_nodes(2), 3 * m_nodes(2) + 1, 3 * m_nodes(2) + 2,
        3 * m_nodes(3), 3 * m_nodes(3) + 1, 3 * m_nodes(3) + 2,
        3 * m_nodes(4), 3 * m_nodes(4) + 1, 3 * m_nodes(4) + 2,
        3 * m_nodes(5), 3 * m_nodes(5) + 1, 3 * m_nodes(5) + 2,
        3 * m_nodes(6), 3 * m_nodes(6) + 1, 3 * m_nodes(6) + 2,
        3 * m_nodes(7), 3 * m_nodes(7) + 1, 3 * m_nodes(7) + 2;
    return res;
}

double Hexahedron::getVolume() const
{
    return m_volume;
}

double Hexahedron::getElasticEnergy(const Matrix<double, 3, 8> &Xcurrent) const
{
    Matrix<double, 3, 3> F;
    double res = 0;
    for (int i = 0; i < 8; i++)
    {
        F = Xcurrent * m_dFdx[i];
        res += m_detJac[i] * m_material.psi(F);
    }
    return res;
}

Matrix<double, 1, 24> Hexahedron::getElasticGradient(const Matrix<double, 3, 8> &Xcurrent) const
{
    Matrix<double, 3, 3> F;
    Matrix<double, 1, 24> res;
    res.setZero();
    for (int i = 0; i < 8; i++)
    {
        F = Xcurrent * m_dFdx[i];
        res += m_detJac[i] * m_material.dpsidF(F) * m_dvFdvx[i];
    }
    return res;
}

Matrix<double, 24, 24> Hexahedron::getElasticHessian(const Matrix<double, 3, 8> &Xcurrent) const
{
    Matrix<double, 3, 3> F;
    Matrix<double, 24, 24> res;
    res.setZero();
    for (int i = 0; i < 8; i++)
    {
        F = Xcurrent * m_dFdx[i];
        res += m_detJac[i] * m_dvFdvx[i].transpose() * m_material.dpsidFF(F) * m_dvFdvx[i];
    }
    return res;
}

double Hexahedron::getEnergy(const MatrixXd &Xcurrent, double s, double c, double a, double k, double d) const
{
    Matrix<double, 3, 8> Xcurrent_elt = Xcurrent(m_nodes, all).transpose();

    Matrix<double, 3, 8> Xbent = bend(Xcurrent_elt, s, c, a, k, d);
    return getElasticEnergy(Xbent);
}
Matrix<double, 1, 24> Hexahedron::getGradient(const MatrixXd &Xcurrent, double s, double c, double a, double k, double d) const
{
    Matrix<double, 3, 8> Xcurrent_elt = Xcurrent(m_nodes, all).transpose();
    Matrix<double, 3, 8> Xbent = bend(Xcurrent_elt, s, c, a, k, d);
    Matrix<double, 8 * 3, 8 * 3> dbenddx = dbend(Xcurrent_elt, s, c, a, k, d);
    return getElasticGradient(Xbent) * dbenddx;
}
Matrix<double, 24, 24> Hexahedron::getHessian(const MatrixXd &Xcurrent, double s, double c, double a, double k, double d) const
{
    Matrix<double, 3, 8> Xcurrent_elt = Xcurrent(m_nodes, all).transpose();
    Matrix<double, 3, 8> Xbent = bend(Xcurrent_elt, s, c, a, k, d);
    Matrix<double, 8 * 3, 8 * 3> dbenddx = dbend(Xcurrent_elt, s, c, a, k, d);
    vector<Matrix<double, 3 * 8, 3 * 8>> ddbenddx = ddbend(Xcurrent_elt, s, c, a, k, d);
    Matrix<double, 24, 24> H1 = dbenddx.transpose() * getElasticHessian(Xbent) * dbenddx;
    Matrix<double, 1, 24> grad = getElasticGradient(Xbent);
    Matrix<double, 24, 24> H2;
    H2.setZero();
    for (int i = 0; i < 24; i++)
    {
        H2 += grad(i) * ddbenddx[i];
    }
    return H1 + H2;
}

Matrix<double, 24, 24> Hexahedron::getHessianClamped(const MatrixXd &Xcurrent, double s, double c, double a, double k, double d) const
{
    Matrix<double, 3, 8> Xcurrent_elt = Xcurrent(m_nodes, all).transpose();
    Matrix<double, 3, 8> Xbent = bend(Xcurrent_elt, s, c, a, k, d);
    Matrix<double, 8 * 3, 8 * 3> dbenddx = dbend(Xcurrent_elt, s, c, a, k, d);
    vector<Matrix<double, 3 * 8, 3 * 8>> ddbenddx = ddbend(Xcurrent_elt, s, c, a, k, d);
    Matrix<double, 24, 24> H1 = dbenddx.transpose() * getElasticHessian(Xbent) * dbenddx;
    Matrix<double, 1, 24> grad = getElasticGradient(Xbent);
    Matrix<double, 24, 24> H2;
    H2.setZero();
    for (int i = 0; i < 24; i++)
    {
        H2 += grad(i) * ddbenddx[i];
    }
    return clamp(H1 + H2);
}

Matrix<double, 3, 8> Hexahedron::bend(const Matrix<double, 3, 8> &X, double s, double c, double a, double k, double d) const
{
    Matrix<double, 3, 8> res = X * 0;
    for (int i = 0; i < 8; i++)
    {
        res(0, i) = symx(a, c, d, k, s, THRESHOLD, X.col(i));
        res(1, i) = symy(a, c, d, s, X.col(i));
        res(2, i) = symz(a, c, d, k, s, THRESHOLD, X.col(i));
    }
    return res;
}

Matrix<double, 3 * 8, 3 * 8> Hexahedron::dbend(const Matrix<double, 3, 8> &X, double s, double c, double a, double k, double d) const
{
    Matrix<double, 3 * 8, 3 * 8> res;
    res.setZero();

    for (int i = 0; i < 8; i++)
    {
        Matrix<double, 3, 3> block;
        VectorXd tmp;
        dsymxdx(a, c, d, k, s, THRESHOLD, X.col(i), tmp);
        block.row(0) = tmp;
        dsymydx(a, c, d, s, tmp);
        block.row(1) = tmp;
        dsymzdx(a, c, d, k, s, THRESHOLD, X.col(i), tmp);
        block.row(2) = tmp;
        res.block<3, 3>(3 * i, 3 * i) = block;
    }
    return res;
}

vector<Matrix<double, 3 * 8, 3 * 8>> Hexahedron::ddbend(const Matrix<double, 3, 8> &X, double s, double c, double a, double k, double d) const
{
    vector<Matrix<double, 3 * 8, 3 * 8>> res(24);
    int cpt = 0;
    for (int i = 0; i < 8; i++)
    {
        MatrixXd tmp;
        MatrixXd block;
        d2symxdx2(a, c, d, k, s, THRESHOLD, X.col(i), block);
        tmp.setZero(24, 24);
        tmp.block<3, 3>(3 * i, 3 * i) = block;
        res[cpt] = tmp;
        cpt++;
        tmp.setZero(24, 24);
        res[cpt] = tmp;
        cpt++;
        d2symzdx2(a, c, d, k, s, THRESHOLD, X.col(i), block);
        tmp.setZero(24, 24);
        tmp.block<3, 3>(3 * i, 3 * i) = block;
        res[cpt] = tmp;
        cpt++;
    }
    return res;
}

Matrix<double, 3 * 8, 3 * 8> Hexahedron::clamp(const Matrix<double, 3 * 8, 3 * 8> &H) const
{
    int n = 3 * 8;
    EigenSolver<MatrixXd> es(H);
    MatrixXd D = es.pseudoEigenvalueMatrix();
    MatrixXd V = es.pseudoEigenvectors();
    for (int i = 0; i < n; i++)
    {
        if (D(i, i) < 0)
            D(i, i) = 0;
    }
    MatrixXd res = V * D * V.inverse();
    return res;
}