#include "Strains.h"

#include "sym_s.inc"
#include "dsym_sdx.inc"
#include "d2sym_sdx2.inc"

#include "sym_c.inc"
#include "dsym_cdx.inc"
#include "d2sym_cdx2.inc"

#include "sym_a.inc"
#include "dsym_adx.inc"
#include "d2sym_adx2.inc"

#include "sym_k1.inc"
#include "dsym_k1dx.inc"
#include "d2sym_k1dx2.inc"

#include "sym_k2.inc"
#include "dsym_k2dx.inc"
#include "d2sym_k2dx2.inc"

#include "sym_d.inc"
#include "dsym_ddx.inc"
#include "d2sym_ddx2.inc"

#define EPSILON 1e-16

Matrix<double, 1, 6> Strains::symStrain(const Matrix<double, 6, 3> &xflaps, const Vector3d &niexists, const Matrix<double, 3, 2> &Xrest)
{
    VectorXd vx = xflaps.transpose().reshaped();
    VectorXd vXrest = Xrest.transpose().reshaped();

    double s = sym_s(EPSILON, vXrest, vx);
    double c = sym_c(EPSILON, vXrest, vx);
    double a = sym_a(EPSILON, vXrest, vx);
    double k1 = sym_k1(EPSILON, vXrest, niexists, vx);
    double k2 = sym_k2(EPSILON, vXrest, niexists, vx);
    double d = sym_d(EPSILON, vXrest, niexists, vx);

    if ((a < 0) or (a > M_PI))
    {
        cout << "a is not in [0,pi]" << endl;
        cout << a << endl;
        exit(0);
    }
    if ((d < 0) or (d > M_PI))
    {
        cout << "d is not in [0,pi]" << endl;
        cout << d << endl;
        exit(0);
    }

    Matrix<double, 1, 6> strain;
    strain << s, c, a, k1, k2, d;

    return strain;
}

Matrix<double, 6, 6 * 3> Strains::symGrad(const Matrix<double, 6, 3> &xflaps, const Vector3d &niexists, const Matrix<double, 3, 2> &Xrest)
{
    VectorXd vx = xflaps.transpose().reshaped();
    VectorXd vXrest = Xrest.transpose().reshaped();
    VectorXd tmp;
    Matrix<double, 6, 6 * 3> grad;

    dsym_sdx(EPSILON, vXrest, vx, tmp);
    grad.row(0) = tmp;
    dsym_cdx(EPSILON, vXrest, vx, tmp);
    grad.row(1) = tmp;
    dsym_adx(EPSILON, vXrest, vx, tmp);
    grad.row(2) = tmp;
    dsym_k1dx(EPSILON, vXrest, niexists, vx, tmp);
    grad.row(3) = tmp;
    dsym_k2dx(EPSILON, vXrest, niexists, vx, tmp);
    grad.row(4) = tmp;
    dsym_ddx(EPSILON, vXrest, niexists, vx, tmp);
    grad.row(5) = tmp;

    if (grad.hasNaN())
    {
        cout << "grad is NaN" << endl;
        cout << symStrain(xflaps, niexists, Xrest) << endl;
        exit(0);
    }

    return grad;
}

vector<Matrix<double, 6 * 3, 6 * 3>> Strains::symHess(const Matrix<double, 6, 3> &xflaps, const Vector3d &niexists, const Matrix<double, 3, 2> &Xrest)
{
    VectorXd vx = xflaps.transpose().reshaped();
    VectorXd vXrest = Xrest.transpose().reshaped();
    vector<Matrix<double, 6 * 3, 6 * 3>> hess;
    MatrixXd tmp;

    d2sym_sdx2(EPSILON, vXrest, vx, tmp);
    hess.push_back(tmp);
    d2sym_cdx2(EPSILON, vXrest, vx, tmp);
    hess.push_back(tmp);
    d2sym_adx2(EPSILON, vXrest, vx, tmp);
    hess.push_back(tmp);
    d2sym_k1dx2(EPSILON, vXrest, niexists, vx, tmp);
    hess.push_back(tmp);
    d2sym_k2dx2(EPSILON, vXrest, niexists, vx, tmp);
    hess.push_back(tmp);
    d2sym_ddx2(EPSILON, vXrest, niexists, vx, tmp);
    hess.push_back(tmp);

    if (hess[0].hasNaN())
    {
        cout << "hess[0] is NaN" << endl;
        cout << symStrain(xflaps, niexists, Xrest) << endl;
        exit(0);
    }

    return hess;
}
