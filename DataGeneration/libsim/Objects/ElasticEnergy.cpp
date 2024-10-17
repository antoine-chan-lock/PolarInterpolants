#include "ElasticEnergy.h"

ElasticEnergy::ElasticEnergy()
{
}

ElasticEnergy::ElasticEnergy(const MatrixXd &Xrest, const MatrixXi &T)
{
    vector<Hexahedron> elts;
    for (int i = 0; i < T.rows(); i++)
    {
        elts.push_back(Hexahedron(i, T, Xrest));
    }
    m_nbElements = T.rows();
    m_nbDoFs = Xrest.rows() * Xrest.cols();
    m_velement = elts;
}

ElasticEnergy::~ElasticEnergy()
{
}

double ElasticEnergy::getEnergy(const MatrixXd &X, double s, double c, double a, double k, double d) const
{
    double res = 0;
    for (int e = 0; e < m_velement.size(); e++)
    {
        res += m_velement[e].getEnergy(X, s, c, a, k, d);
    }
    return res;
}

VectorXd ElasticEnergy::getGradient(const MatrixXd &X, double s, double c, double a, double k, double d) const
{

    VectorXd res = VectorXd::Zero(X.rows() * X.cols());
    for (int e = 0; e < m_velement.size(); e++)
    {
        Matrix<double, 1, 8 * 3> localGrad;
        localGrad = m_velement[e].getGradient(X, s, c, a, k, d);
        res(m_velement[e].getIndexDoFs()) += localGrad;
    }
    return res;
}

SparseMatrix<double> ElasticEnergy::getHessian(const MatrixXd &X, double s, double c, double a, double k, double d) const
{
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    tripletList.reserve(24 * 24 * m_nbElements);
    vector<Matrix<double, 24, 24>> v_localHess(m_nbElements);
    // #pragma omp parallel for

    for (int e = 0; e < m_velement.size(); e++)
    {
        v_localHess[e] = m_velement[e].getHessian(X, s, c, a, k, d);
    }

    for (int e = 0; e < m_velement.size(); e++)
    {
        Matrix<double, 24, 24> localHess = v_localHess[e];
        Matrix<int, 1, 24> idx = m_velement[e].getIndexDoFs();
        for (int i = 0; i < 24; i++)
        {
            for (int j = 0; j < 24; j++)
            {
                tripletList.push_back(T(idx(i), idx(j), localHess(i, j)));
            }
        }
    }
    SparseMatrix<double> mat(m_nbDoFs, m_nbDoFs);
    mat.setFromTriplets(tripletList.begin(), tripletList.end());

    return mat;
}

SparseMatrix<double> ElasticEnergy::getHessianClamped(const MatrixXd &X, double s, double c, double a, double k, double d) const
{
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    tripletList.reserve(24 * 24 * m_nbElements);
    vector<Matrix<double, 24, 24>> v_localHess(m_nbElements);
    // #pragma omp parallel for

    for (int e = 0; e < m_velement.size(); e++)
    {
        v_localHess[e] = m_velement[e].getHessianClamped(X, s, c, a, k, d);
    }

    for (int e = 0; e < m_velement.size(); e++)
    {
        Matrix<double, 24, 24> localHess = v_localHess[e];
        Matrix<int, 1, 24> idx = m_velement[e].getIndexDoFs();
        for (int i = 0; i < 24; i++)
        {
            for (int j = 0; j < 24; j++)
            {
                tripletList.push_back(T(idx(i), idx(j), localHess(i, j)));
            }
        }
    }
    SparseMatrix<double> mat(m_nbDoFs, m_nbDoFs);
    mat.setFromTriplets(tripletList.begin(), tripletList.end());

    return mat;
}

double ElasticEnergy::getVolume() const
{
    double res = 0;
    for (int e = 0; e < m_velement.size(); e++)
    {
        res += m_velement[e].getVolume();
    }
    return res;
}
