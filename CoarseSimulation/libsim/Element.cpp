#include "Element.h"

Element::Element()
{
}

Element::Element(int id, const MatrixXi &T, const MatrixXi &TT)
{
    m_id = id;
    m_isIn = 1;
    VectorXi nodes = T(id, all);
    m_neighbors = TT.row(id);
    VectorXi flaps(3);
    m_niexists.setOnes();
    for (int i = 0; i < 3; i++)
    {
        if (TT(id, i) == -1)
        {
            flaps(i) = -1;
            m_niexists(i) = 0;
            m_isIn = 0;
        }
        else
        {
            VectorXi tmp = T(TT(id, i), all);
            for (int j = 0; j < 3; j++)
            {
                if (tmp(j) != nodes(0) && tmp(j) != nodes(1) && tmp(j) != nodes(2))
                {
                    flaps(i) = tmp(j);
                }
            }
        }
    }
    VectorXi tmp(nodes.rows() + flaps.rows());
    tmp << nodes,
        flaps;
    m_T = tmp;
}

Element::~Element()
{
}

Matrix<double, 6, 3> Element::getxflaps(const MatrixXd &x) const
{
    Matrix<double, 6, 3> out;
    out.setZero();
    for (int i = 0; i < 6; i++)
    {
        if (m_T(i) != -1)
        {
            out.row(i) = x.row(m_T(i));
        }
        else
        {
            int id_p1 = (i - 2) % 3;
            int id_p2 = (i - 3) % 3;
            int id_p3 = (i - 1) % 3;
            Vector3d p1 = x.row(m_T(id_p1));
            Vector3d p2 = x.row(m_T(id_p2));
            Vector3d p3 = x.row(m_T(id_p3));
            Vector3d v1 = p2 - p1;
            Vector3d v2 = p3 - p1;
            Vector3d proj = p1 + (v1.dot(v2) / v1.dot(v1)) * v1;
            Vector3d v3 = proj - p3;
            out.row(i) = p3 + 2 * v3;
        }
    }
    return out;
}

Matrix<double, 3, 2> Element::getXrest(const MatrixXd &X) const
{
    Matrix<double, 3, 2> out;
    out.setZero();
    for (int i = 0; i < 3; i++)
    {
        if (m_T(i) != -1)
        {
            out.row(i) = X.row(m_T(i));
        }
        else
        {
            cout << "Element::getXrest::shouldnt be called" << endl;
            exit(0);
        }
    }
    return out;
}

double Element::getArea(const MatrixXd &X) const
{
    Matrix<double, 3, 2> Xrest = getXrest(X);
    Vector2d v1 = Xrest.row(1) - Xrest.row(0);
    Vector2d v2 = Xrest.row(2) - Xrest.row(0);
    return 0.5 * abs(v1(0) * v2(1) - v1(1) * v2(0));
}

Matrix<double, 6 * 3, 6 * 3> Element::clamp(const Matrix<double, 6 * 3, 6 * 3> &H) const
{
    EigenSolver<Matrix<double, 6 * 3, 6 * 3>> es((H + H.transpose()) / 2);

    Matrix<double, 6 * 3, 6 * 3> D = es.pseudoEigenvalueMatrix();
    Matrix<double, 6 * 3, 6 * 3> V = es.pseudoEigenvectors();

    for (int i = 0; i < 6 * 3; i++)
        if (D(i, i) < 0)
            D(i, i) = 0;
    Matrix<double, 6 * 3, 6 * 3> res = V * D * V.inverse();
    if (isnan(res.norm()))
    {
        cout << H << endl;
        cout << "--" << endl;
        cout << D << endl;
        cout << "--" << endl;
        cout << V << endl;
        cout << "--" << endl;
        cout << res << endl;
        cout << "nan alert" << endl;
        exit(0);
    }
    return res;
}

double Element::elasticEnergy(const MatrixXd &x, const MatrixXd &X, Matrix<double, 1, 6 * 3> *grad, Matrix<double, 6 * 3, 6 * 3> *hess, const Material &mat) const
{
    Matrix<double, 6, 6> clampBending;
    clampBending.setIdentity();
    Matrix<double, 1, 6> disp;
    disp.setZero(); 
    if (m_isIn == 0)
    {
        clampBending(3, 3) = 0;
        clampBending(4, 4) = 0;
        clampBending(5, 5) = 0;
    }
    Matrix<double, 6, 3> xflaps = getxflaps(x);
    Matrix<double, 3, 2> Xrest = getXrest(X);
    Matrix<double, 1, 6> S = Strains::symStrain(xflaps, m_niexists, Xrest);
    double area = getArea(X);

    double ener;
    Matrix<double, 1, 6> dener;
    Matrix<double, 6, 6> ddener;
    if (grad or hess)
    {
        ener = mat.elasticEnergy(S * clampBending + disp, &dener, &ddener, m_isIn);
        dener = dener * clampBending;
        ddener = clampBending * ddener * clampBending;
    }
    else
    {
        ener = mat.elasticEnergy(S * clampBending + disp, NULL, NULL, m_isIn);
    }
    Matrix<double, 6, 6 * 3> dSdx;
    if (grad)
    {
        dSdx = Strains::symGrad(xflaps, m_niexists, Xrest);
        (*grad) = dener * dSdx;
        (*grad) *= area;
    }
    vector<Matrix<double, 6 * 3, 6 * 3>> dSidxx;
    if (hess)
    {
        dSidxx = Strains::symHess(xflaps, m_niexists, Xrest);
        (*hess) = dSdx.transpose() * ddener * dSdx;
        for (int i = 0; i < 6; i++)
        {
            (*hess) += dener(0, i) * dSidxx[i];
        }
        (*hess) = clamp((*hess) * area);
    }

    return ener * area;
}