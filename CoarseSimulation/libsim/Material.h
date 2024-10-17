#pragma once
#include "CommonIncludes.h"
#include "Tools.h"
#define BARRIERRANGE 0.5
class Material
{
public:
    Material();
    ~Material();
    Material(string dirname);

    void setInit(int x) { m_isInit = x; }; // Init energy = l1^2 + l2^2 + k1^2 + k2^2

    // RBFs and derivatives
    double RBF(const Matrix<double, 1, 5> &s, const Matrix<double, 1, 5> &C, double r0) const;
    Matrix<double, 1, 5> dRBF(const Matrix<double, 1, 5> &s, const Matrix<double, 1, 5> &C, double r0) const;
    Matrix<double, 5, 5> ddRBF(const Matrix<double, 1, 5> &s, const Matrix<double, 1, 5> &C, double r0) const;

    // Single curvature stress interpolant and derivatives
    double hiointerp(const Matrix<double, 1, 5> &s, const MatrixXd &C, double r0, const MatrixXd &w) const;
    Matrix<double, 1, 5> dhiointerp(const Matrix<double, 1, 5> &s, const MatrixXd &C, double r0, const MatrixXd &w) const;
    Matrix<double, 5, 5> ddhiointerp(const Matrix<double, 1, 5> &s, const MatrixXd &C, double r0, const MatrixXd &w) const;

    // Polar interpolants and derivatives
    double stressInterp(const Matrix<double, 1, 6> &s) const;
    Matrix<double, 1, 6> dstressInterp(const Matrix<double, 1, 6> &s) const;
    Matrix<double, 6, 6> ddstressInterp(const Matrix<double, 1, 6> &s) const;

    // Barrier functions and derivatives
    double barrierEnergy(const Matrix<double, 1, 6> &s, Matrix<double, 1, 6> &grad, Matrix<double, 6, 6> &hess) const;
    double fbarrier(double d) const;
    double dfbarrier(double d) const;
    double ddfbarrier(double d) const;

    // Init energy and derivatives
    double initEnergy(const Matrix<double, 1, 6> &s, Matrix<double, 1, 6> &grad, Matrix<double, 6, 6> &hess) const;
    // Elastic energy
    double elasticEnergy(const Matrix<double, 1, 6> &s, Matrix<double, 1, 6> *gradient, Matrix<double, 6, 6> *hessian, const int &isIn) const;
    // isotropicSafeEnergy() makes sure that the function stays continuous when the angles a and d are undetermined.
    double isotropicSafeEnergy(const Matrix<double, 1, 6> &s, Matrix<double, 1, 6> *gradient, Matrix<double, 6, 6> *hessian) const;
    // pair__() makes sure that energy(K1) and energy(K2) is pair for stability
    double pairK1(const Matrix<double, 1, 6> &s, Matrix<double, 1, 6> *gradient, Matrix<double, 6, 6> *hessian) const;
    double pairK2(const Matrix<double, 1, 6> &s, Matrix<double, 1, 6> *gradient, Matrix<double, 6, 6> *hessian) const;
    // flipsafe() makes sure that energy(s,c,a) = energy(-c,-s,a+90) for stability
    // s represent elongation, c compression, and a the elongation direction. But this can get flipped
    double flipsafe(const Matrix<double, 1, 6> &s, Matrix<double, 1, 6> *gradient, Matrix<double, 6, 6> *hessian) const;

    // freeze is a smooth continuous function that is 0 when l1=l2, otherwise 1.
    double freeze(double l1, double l2, double h2, double h1, Matrix<double, 1, 2> &gradient, Matrix<double, 2, 2> &hessian) const;
    double freezeIP(double l1, double l2, double h2, Matrix<double, 1, 6> &gradient, Matrix<double, 6, 6> &hessian) const;
    double freezeBEND(double l1, double l2, double h2, Matrix<double, 1, 6> &gradient, Matrix<double, 6, 6> &hessian) const;


    int m_isInit = 0;
    MatrixXd m_C;                                          // RBFs control points
    double m_r0;                                           // RBF radius
    MatrixXd m_w;                                          // RBF weights
    double m_smin, m_smax, m_kmin, m_kmax, m_cmin, m_cmax; // safe strain ranges
    MatrixXd m_range;                                      // safe strain ranges
    Matrix<double, 1, 5> m_vrange;                         // safe strain ranges
    double m_smoothParameter = 50;                         // isotropy blending smoothness
};