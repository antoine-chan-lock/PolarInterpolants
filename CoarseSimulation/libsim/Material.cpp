#include "Material.h"


#define BARRIER

Material::Material()
{
}

Material::~Material()
{
}

Material::Material(string dirname)
{
    MatrixXd tmp = Tools::loadCsvDouble(dirname + "/r0C.csv");
    m_C = tmp(seq(1, last), all);
    m_r0 = tmp(0, 0);
    m_w = Tools::loadCsvDouble(dirname + "/w.csv");
    MatrixXd vrange = Tools::loadCsvDouble(dirname + "/strainRange.csv");
    m_range = vrange;
    m_smin = vrange(0, 0);
    m_smax = vrange(1, 0);
    m_cmin = -0.9; //-1 is double elongation. c isn't clearly sampled, it's a result of a minimal energy search in the Data generation, so we can't give clear bounds
    m_cmax = 0.9; //1 is full compression
    m_kmin = -vrange(1, 3); 
    m_kmax = vrange(1, 3);
    m_vrange = vrange.row(1) - vrange.row(0);
}

double Material::RBF(const Matrix<double, 1, 5> &s, const Matrix<double, 1, 5> &C, double r0) const
{
    Matrix<double, 1, 5> dist = (s - C);
    dist(0, 2) = sin(dist(0, 2));
    dist(0, 4) = sin(dist(0, 4));
    dist = dist.cwiseQuotient(m_vrange);
    double r2 = dist.squaredNorm();
    double r02 = r0 * r0;
    return exp(-r2 / r02);
}

Matrix<double, 1, 5> Material::dRBF(const Matrix<double, 1, 5> &s, const Matrix<double, 1, 5> &C, double r0) const
{
    Matrix<double, 1, 5> dist = (s - C);
    Matrix<double, 1, 5> distMod = dist;
    distMod(0, 2) = sin(distMod(0, 2));
    distMod(0, 4) = sin(distMod(0, 4));
    distMod = distMod.cwiseQuotient(m_vrange);
    double r2 = distMod.squaredNorm();
    double r02 = r0 * r0;

    double f = exp(-r2 / r02);
    double dfdr2 = -f / r02;

    Matrix<double, 1, 5> dr2ds;
    dr2ds << 2 * dist(0, 0) / pow(m_vrange(0, 0), 2),
        2 * dist(0, 1) / pow(m_vrange(0, 1), 2),
        2 * sin(dist(0, 2)) * cos(dist(0, 2)) / pow(m_vrange(0, 2), 2),
        2 * dist(0, 3) / pow(m_vrange(0, 3), 2),
        2 * sin(dist(0, 4)) * cos(dist(0, 4)) / pow(m_vrange(0, 4), 2);
    Matrix<double, 1, 5> dfds = dfdr2 * dr2ds;
    return dfds;
}

Matrix<double, 5, 5> Material::ddRBF(const Matrix<double, 1, 5> &s, const Matrix<double, 1, 5> &C, double r0) const
{
    Matrix<double, 1, 5> dist = (s - C);
    Matrix<double, 1, 5> distMod = dist;
    distMod(0, 2) = sin(distMod(0, 2));
    distMod(0, 4) = sin(distMod(0, 4));
    distMod = distMod.cwiseQuotient(m_vrange);
    double r2 = distMod.squaredNorm();
    double r02 = r0 * r0;

    double f = exp(-r2 / r02);
    double dfdr2 = -f / r02;
    double d2fdr22 = f / (r02 * r02);

    Matrix<double, 1, 5> dr2ds;
    dr2ds << 2 * dist(0, 0) / pow(m_vrange(0, 0), 2),
        2 * dist(0, 1) / pow(m_vrange(0, 1), 2),
        2 * sin(dist(0, 2)) * cos(dist(0, 2)) / pow(m_vrange(0, 2), 2),
        2 * dist(0, 3) / pow(m_vrange(0, 3), 2),
        2 * sin(dist(0, 4)) * cos(dist(0, 4)) / pow(m_vrange(0, 4), 2);
    Matrix<double, 1, 5> dfds = dfdr2 * dr2ds;

    Matrix<double, 5, 5> dr2dss;
    dr2dss.setIdentity();
    dr2dss(2, 2) = cos(2 * dist(0, 2));
    dr2dss(4, 4) = cos(2 * dist(0, 4));
    dr2dss *= 2;
    Matrix<double, 5, 5> tmp = m_vrange.cwiseInverse().cwiseProduct(m_vrange.cwiseInverse()).asDiagonal();
    dr2dss = dr2dss * tmp;

    Matrix<double, 5, 5> d2fds2 = d2fdr22 * dr2ds.transpose() * dr2ds + dfdr2 * dr2dss;

    return d2fds2;
}

double Material::hiointerp(const Matrix<double, 1, 5> &s, const MatrixXd &C, double r0, const MatrixXd &w) const
{
    double res = 0;
    for (int i = 0; i < C.rows(); i++)
    {
        Matrix<double, 1, 5> dist = (s - C.row(i));
        dist(0, 2) = sin(dist(0, 2)) * cos(dist(0, 2));
        dist(0, 4) = sin(dist(0, 4)) * cos(dist(0, 4));
        dist = dist.cwiseQuotient(m_vrange);
        res += RBF(s, C.row(i), r0) * (w.row(i).dot(dist));
    }
    return res;
}

Matrix<double, 1, 5> Material::dhiointerp(const Matrix<double, 1, 5> &s, const MatrixXd &C, double r0, const MatrixXd &w) const
{
    Matrix<double, 1, 5> res;
    res.setZero();
    for (int i = 0; i < C.rows(); i++)
    {
        Matrix<double, 1, 5> dist = (s - C.row(i));
        dist(0, 2) = sin(dist(0, 2)) * cos(dist(0, 2));
        dist(0, 4) = sin(dist(0, 4)) * cos(dist(0, 4));
        dist = dist.cwiseQuotient(m_vrange);
        Matrix<double, 1, 5> tmp2 = w.row(i);
        tmp2(0, 2) *= pow(cos(s(0, 2) - C(i, 2)), 2) - pow(sin(s(0, 2) - C(i, 2)), 2);
        tmp2(0, 4) *= pow(cos(s(0, 4) - C(i, 4)), 2) - pow(sin(s(0, 4) - C(i, 4)), 2);
        tmp2 = tmp2.cwiseQuotient(m_vrange);
        res += RBF(s, C.row(i), r0) * tmp2;
        res += w.row(i).dot(dist) * dRBF(s, C.row(i), r0);
    }
    return res;
}

Matrix<double, 5, 5> Material::ddhiointerp(const Matrix<double, 1, 5> &s, const MatrixXd &C, double r0, const MatrixXd &w) const
{
    Matrix<double, 5, 5> res;
    res.setZero();
    for (int i = 0; i < C.rows(); i++)
    {
        Matrix<double, 1, 5> ds = (s - C.row(i));
        Matrix<double, 1, 5> dsm = ds;
        dsm(0, 2) = sin(dsm(0, 2)) * cos(dsm(0, 2));
        dsm(0, 4) = sin(dsm(0, 4)) * cos(dsm(0, 4));
        dsm = dsm.cwiseQuotient(m_vrange);
        res += (w.row(i).dot(dsm)) * ddRBF(s, C.row(i), r0);
        Matrix<double, 1, 5> d_dot = w.row(i);
        d_dot(2) *= pow(cos(ds(2)), 2) - pow(sin(ds(2)), 2);
        d_dot(4) *= pow(cos(ds(4)), 2) - pow(sin(ds(4)), 2);
        d_dot = d_dot.cwiseQuotient(m_vrange);
        Matrix<double, 5, 5> tmp = d_dot.transpose() * dRBF(s, C.row(i), r0);
        res += tmp + tmp.transpose();
        Matrix<double, 5, 5> dd_dot;
        dd_dot.setZero();
        dd_dot(2, 2) = w(i, 2) * (-4 * cos(ds(2)) * sin(ds(2))) / m_vrange(0, 2); // (-4*cos(ds(2))*sin(ds(2)))
        dd_dot(4, 4) = w(i, 4) * (-4 * cos(ds(4)) * sin(ds(4))) / m_vrange(0, 4); // (-4*cos(ds(4))*sin(ds(4)))
        res += dd_dot.transpose() * RBF(s, C.row(i), r0);
    }
    return res;
}

double Material::stressInterp(const Matrix<double, 1, 6> &s) const
{
    Matrix<double, 1, 5> sk1;
    sk1 << s(0), s(1), s(2), s(3), s(5);
    Matrix<double, 1, 5> sk2;
    sk2 << s(0), s(1), s(2), s(4), s(5) + M_PI / 2.;
    return hiointerp(sk1, m_C, m_r0, m_w) + hiointerp(sk2, m_C, m_r0, m_w);
}

Matrix<double, 1, 6> Material::dstressInterp(const Matrix<double, 1, 6> &s) const
{
    Matrix<double, 1, 5> sk1;
    sk1 << s(0), s(1), s(2), s(3), s(5);
    Matrix<double, 1, 5> sk2;
    sk2 << s(0), s(1), s(2), s(4), s(5) + M_PI / 2.;
    Matrix<double, 1, 5> dsk1 = dhiointerp(sk1, m_C, m_r0, m_w);
    Matrix<double, 1, 5> dsk2 = dhiointerp(sk2, m_C, m_r0, m_w);
    Matrix<double, 1, 6> res;
    res << dsk1(0) + dsk2(0), // 0
        dsk1(1) + dsk2(1),    // 1
        dsk1(2) + dsk2(2),    // 2
        dsk1(3),              // 3
        dsk2(3),              // 4
        dsk1(4) + dsk2(4);    // 5
    return res;
}
Matrix<double, 6, 6> Material::ddstressInterp(const Matrix<double, 1, 6> &s) const
{
    Matrix<double, 1, 5> sk1;
    sk1 << s(0), s(1), s(2), s(3), s(5);
    Matrix<double, 1, 5> sk2;
    sk2 << s(0), s(1), s(2), s(4), s(5) + M_PI / 2.;
    Matrix<double, 5, 5> dsk1 = ddhiointerp(sk1, m_C, m_r0, m_w);
    Matrix<double, 5, 5> dsk2 = ddhiointerp(sk2, m_C, m_r0, m_w);
    Matrix<double, 6, 6> res;
    VectorXi idxk1(5), idxk2(5);
    idxk1 << 0, 1, 2, 3, 5;
    idxk2 << 0, 1, 2, 4, 5;
    res.setZero();
    for (int i = 0; i < 5; i++)
    {
        for (int j = 0; j < 5; j++)
        {
            res(idxk1(i), idxk1(j)) += dsk1(i, j);
            res(idxk2(i), idxk2(j)) += dsk2(i, j);
        }
    }
    return res;
}

double Material::barrierEnergy(const Matrix<double, 1, 6> &s, Matrix<double, 1, 6> &grad, Matrix<double, 6, 6> &hess) const
{
    double res = 0;
    Matrix<double, 8, 1> vd;
    vd << s(0) - m_smin,
        s(0) - m_smax,
        s(1) - m_cmin,
        s(1) - m_cmax,
        s(3) - m_kmin,
        s(3) - m_kmax,
        s(4) - m_kmin,
        s(4) - m_kmax;
    if (s(0) < m_smin)
        res += fbarrier(vd(0));
    if (s(0) > m_smax)
        res += fbarrier(vd(1));

    if (s(1) < m_cmin)
        res += fbarrier(vd(2));
    if (s(1) > m_cmax)
        res += fbarrier(vd(3));

    if (s(3) < m_kmin)
        res += fbarrier(vd(4));
    if (s(3) > m_kmax)
        res += fbarrier(vd(5));

    if (s(4) < m_kmin)
        res += fbarrier(vd(6));
    if (s(4) > m_kmax)
        res += fbarrier(vd(7));

    grad.setZero();
    if (s(0) < m_smin)
        grad(0) += dfbarrier(vd(0));
    if (s(0) > m_smax)
        grad(0) += dfbarrier(vd(1));

    if (s(1) < m_cmin)
        grad(1) += dfbarrier(vd(2));
    if (s(1) > m_cmax)
        grad(1) += dfbarrier(vd(3));

    if (s(3) < m_kmin)
        grad(3) += dfbarrier(vd(4));
    if (s(3) > m_kmax)
        grad(3) += dfbarrier(vd(5));

    if (s(4) < m_kmin)
        grad(4) += dfbarrier(vd(6));
    if (s(4) > m_kmax)
        grad(4) += dfbarrier(vd(7));

    hess.setZero();
    if (s(0) < m_smin)
        hess(0, 0) += ddfbarrier(vd(0));
    if (s(0) > m_smax)
        hess(0, 0) += ddfbarrier(vd(1));

    if (s(1) < m_cmin)
        hess(1, 1) += ddfbarrier(vd(2));
    if (s(1) > m_cmax)
        hess(1, 1) += ddfbarrier(vd(3));

    if (s(3) < m_kmin)
        hess(3, 3) += ddfbarrier(vd(4));
    if (s(3) > m_kmax)
        hess(3, 3) += ddfbarrier(vd(5));

    if (s(4) < m_kmin)
        hess(4, 4) += ddfbarrier(vd(6));
    if (s(4) > m_kmax)
        hess(4, 4) += ddfbarrier(vd(7));

    return res;
}

double Material::fbarrier(double d) const
{
    return 1e6*pow(d,2);
}

double Material::dfbarrier(double d) const
{
    return 2e6*d;
}

double Material::ddfbarrier(double d) const
{
    return 2e6;
}

double Material::initEnergy(const Matrix<double, 1, 6> &s, Matrix<double, 1, 6> &grad, Matrix<double, 6, 6> &hess) const
{
    double superCoef = 1;
    double res = superCoef*s(0) * s(0) + superCoef*s(1) * s(1) + s(3) * s(3) + s(4) * s(4);
    grad << superCoef*2 * s(0), superCoef*2 * s(1), 0, 2 * s(3), 2 * s(4), 0;
    hess.setZero();
    hess(0, 0) = superCoef*2;
    hess(1, 1) = superCoef*2;
    hess(3, 3) = 2;
    hess(4, 4) = 2;

    return res;
}

double Material::elasticEnergy(const Matrix<double, 1, 6> &s, Matrix<double, 1, 6> *gradient, Matrix<double, 6, 6> *hessian,const int &isIn) const
{
    if (m_isInit != 0)
    {
        Matrix<double, 1, 6> gradInit;
        Matrix<double, 6, 6> hessInit;
        double energyInit = initEnergy(s, gradInit, hessInit);


        if (gradient)
        {
            (*gradient) = gradInit ;
        }
        if (hessian)
        {
            (*hessian) = hessInit ;
        }
        return energyInit ;
    }
    else
    {

        return flipsafe(s, gradient,hessian);

    }
}

double Material::isotropicSafeEnergy(const Matrix<double, 1, 6> &s, Matrix<double, 1, 6> *gradient, Matrix<double, 6, 6> *hessian) const
{

    if (m_isInit != 0)
    {
        Matrix<double, 1, 6> gradInit;
        Matrix<double, 6, 6> hessInit;
        double energyInit = initEnergy(s, gradInit, hessInit);
        if (gradient)
        {
            (*gradient) = gradInit;
        }
        if (hessian)
        {
            (*hessian) = hessInit;
        }
        return energyInit;
    }

    Matrix<double, 1, 6> gradBarrier;
    Matrix<double, 6, 6> hessBarrier;
    double energyBarrier = barrierEnergy(s, gradBarrier, hessBarrier);

    Matrix<double, 6, 6> clampIP;
    clampIP.setIdentity();
    clampIP(2, 2) = 0;
    Matrix<double, 6, 6> clampBEND;
    clampBEND.setIdentity();
    clampBEND(5, 5) = 0;
    Matrix<double, 6, 6> clampBOTH = clampIP * clampBEND;


    Matrix<double, 1, 6> dfIP, dufIP, dfBEND, dufBEND, dpsi00, dpsi10, dpsi01, dpsi11;
    Matrix<double, 6, 6> ddfIP, ddufIP, ddfBEND, ddufBEND, ddpsi00, ddpsi10, ddpsi01, ddpsi11;
    double h2 = m_smoothParameter;
    double fIP = freezeIP(s(0), s(1), h2, dfIP, ddfIP);
    double ufIP = 1 - fIP;
    dufIP = -dfIP;
    ddufIP = -ddfIP;

    double fBEND = freezeBEND(s(3), s(4), h2, dfBEND, ddfBEND);
    double ufBEND = 1 - fBEND;
    dufBEND = -dfBEND;
    ddufBEND = -ddfBEND;

    Matrix<double,1,6> s0;
    s0.setZero();

    Matrix<double,1,6> d0;
    d0.setZero();

    double psi11;
    double psi01;
    double psi10;
    double psi00;

    double computationThreshold = 1e-8;

    if ( abs(ufIP * ufBEND)<computationThreshold )
    {
        psi00 = 0;
    }
    else
    {
        psi00 = stressInterp(s * clampBOTH + s0 + d0);
    }

    if ( abs(ufIP * fBEND)<computationThreshold )
    {
        psi01 = 0;
    }
    else
    {
        psi01 = stressInterp(s * clampIP + s0);
    }

    if ( abs(fIP * ufBEND)<computationThreshold )
    {
        psi10 = 0;
    }
    else
    {
        psi10 = stressInterp(s * clampBEND + d0);
    }

    if ( abs(fIP * fBEND)<computationThreshold )
    {
        psi11 = 0;
    }
    else
    {
        psi11 = stressInterp(s);
    }

    double energy = psi00 * ufIP * ufBEND  //(e1)
                    + psi01 * ufIP * fBEND //(e2)
                    + psi10 * fIP * ufBEND //(e3)
                    + psi11 * fIP * fBEND; //(e4)

    if (gradient)
    {
        gradient->setZero();

        if ( abs(ufIP * ufBEND)<computationThreshold )
        {
            dpsi00.setZero();
        }
        else
        {
            dpsi00 = dstressInterp(s * clampBOTH + s0 + d0) * clampBOTH;
        }

        if ( abs(ufIP * fBEND)<computationThreshold )
        {
            dpsi01.setZero();
        }
        else
        {
            dpsi01 = dstressInterp(s * clampIP+s0) * clampIP;
        }

        if ( abs(fIP * ufBEND)<computationThreshold )
        {
            dpsi10.setZero();
        }
        else
        {
            dpsi10 = dstressInterp(s * clampBEND+d0) * clampBEND;
        }

        if ( abs(fIP * fBEND)<computationThreshold )
        {
            dpsi11.setZero();
        }
        else
        {
            dpsi11 = dstressInterp(s);
        }



        (*gradient) = dpsi00 * ufIP * ufBEND + psi00 * dufIP * ufBEND + psi00 * ufIP * dufBEND + dpsi01 * ufIP * fBEND + psi01 * dufIP * fBEND + psi01 * ufIP * dfBEND + dpsi10 * fIP * ufBEND + psi10 * dfIP * ufBEND + psi10 * fIP * dufBEND + dpsi11 * fIP * fBEND + psi11 * dfIP * fBEND + psi11 * fIP * dfBEND;
        (*gradient) += gradBarrier;
    }

    if (hessian)
    {
        hessian->setZero();

        if ( abs(ufIP * ufBEND)<computationThreshold )
        {
            ddpsi00.setZero();
        }
        else
        {
            ddpsi00 = clampBOTH * ddstressInterp(s * clampBOTH + s0 + d0) * clampBOTH;
        }

        if ( abs(ufIP * fBEND)<computationThreshold )
        {
            ddpsi01.setZero();
        }
        else
        {
            ddpsi01 = clampIP * ddstressInterp(s * clampIP + s0) * clampIP;
        }

        if ( abs(fIP * ufBEND)<computationThreshold )
        {
            ddpsi10.setZero();
        }
        else
        {
            ddpsi10 = clampBEND * ddstressInterp(s * clampBEND+d0) * clampBEND;
        }

        if ( abs(fIP * fBEND)<computationThreshold )
        {
            ddpsi11.setZero();
        }
        else
        {
            ddpsi11 = ddstressInterp(s);
        }

        // add transpose() later
        (*hessian) = ddpsi00 * ufIP * ufBEND + dpsi00.transpose() * dufIP * ufBEND + dpsi00.transpose() * ufIP * dufBEND + (dpsi00.transpose() * dufIP * ufBEND).transpose() + psi00 * ddufIP * ufBEND + psi00 * dufIP.transpose() * dufBEND + (dpsi00.transpose() * ufIP * dufBEND).transpose() + (psi00 * dufIP.transpose() * dufBEND).transpose() + psi00 * ufIP * ddufBEND

                     + ddpsi01 * ufIP * fBEND + dpsi01.transpose() * dufIP * fBEND + dpsi01.transpose() * ufIP * dfBEND + (dpsi01.transpose() * dufIP * fBEND).transpose() + psi01 * ddufIP * fBEND + psi01 * dufIP.transpose() * dfBEND + (dpsi01.transpose() * ufIP * dfBEND).transpose() + (psi01 * dufIP.transpose() * dfBEND).transpose() + psi01 * ufIP * ddfBEND

                     + ddpsi10 * fIP * ufBEND + dpsi10.transpose() * dfIP * ufBEND + dpsi10.transpose() * fIP * dufBEND + (dpsi10.transpose() * dfIP * ufBEND).transpose() + psi10 * ddfIP * ufBEND + psi10 * dfIP.transpose() * dufBEND + (dpsi10.transpose() * fIP * dufBEND).transpose() + (psi10 * dfIP.transpose() * dufBEND).transpose() + psi10 * fIP * ddufBEND

                     + ddpsi11 * fIP * fBEND + dpsi11.transpose() * dfIP * fBEND + dpsi11.transpose() * fIP * dfBEND + (dpsi11.transpose() * dfIP * fBEND).transpose() + psi11 * ddfIP * fBEND + psi11 * dfIP.transpose() * dfBEND + (dpsi11.transpose() * fIP * dfBEND).transpose() + (psi11 * dfIP.transpose() * dfBEND).transpose() + psi11 * fIP * ddfBEND;
        (*hessian) += hessBarrier;
    }

    return energy + energyBarrier;
}

double Material::pairK1(const Matrix<double, 1, 6> &s, Matrix<double, 1, 6> *gradient,Matrix<double, 6, 6> *hessian) const
{
    Matrix<double, 1, 6> g1,g2;
    Matrix<double, 6, 6> h1,h2;
    Matrix<double,6,6> mk1; mk1.setIdentity();
    mk1(3,3) = -1;
    double e1 = pairK2(s,&g1,&h1);
    double e2 = pairK2(s*mk1,&g2,&h2);
    g2 = g2*mk1;
    h2 = mk1*h2*mk1;
    double ener = (e1 + e2)/2;
    if (gradient)
    {
        (*gradient) = (g1 + g2)/2;
    }
    if (hessian)
    {
        (*hessian) =  (h1 + h2)/2;
    }
    return ener;
}

double Material::pairK2(const Matrix<double, 1, 6> &s, Matrix<double, 1, 6> *gradient,Matrix<double, 6, 6> *hessian) const
{
    Matrix<double, 1, 6> g1,g2;
    Matrix<double, 6, 6> h1,h2;
    Matrix<double,6,6> mk2; mk2.setIdentity();
    mk2(4,4) = -1;
    double e1 = isotropicSafeEnergy(s,&g1,&h1);
    double e2 = isotropicSafeEnergy(s*mk2,&g2,&h2);
    g2 = g2*mk2;
    h2 = mk2*h2*mk2;
    double ener = (e1 + e2)/2;
    if (gradient)
    {
        (*gradient) = (g1 + g2)/2;
    }
    if (hessian)
    {
        (*hessian) =  (h1 + h2)/2;
    }
    return ener;
}

double Material::flipsafe(const Matrix<double, 1, 6> &s, Matrix<double, 1, 6> *gradient,Matrix<double, 6, 6> *hessian) const
{
    Matrix<double,6,6> flip;
    flip.setIdentity();
    flip(0,0) = 0;
    flip(1,1) = 0;
    flip(1,0) = -1;
    flip(0,1) = -1;

    double ener1,ener2,ener;
    Matrix<double, 1, 6> dener1,dener2;
    Matrix<double, 6, 6> ddener1,ddener2;

    Matrix<double, 1, 6> dis;dis.setZero();
    if (s(2) > (M_PI / 2))
    {
        dis(2) = - M_PI / 2;
    }
    else
    {
        dis(2) =  M_PI / 2;
    }
    //pairK1 
    if (gradient and hessian)
    {
        ener1 = pairK1(s,&dener1,&ddener1);
        ener2 = pairK1(s*flip+dis,&dener2,&ddener2);
        dener2 = dener2*flip;
        ddener2 = flip*ddener2*flip;

        ener = (ener1 + ener2)/2;
        (*gradient) = (dener1 + dener2)/2.;
        (*hessian) =  (ddener1 + ddener2)/2.;
    }
    else
    {
        ener1 = pairK1(s,NULL,NULL);
        ener2 = pairK1(s*flip+dis,NULL,NULL);
        ener = (ener1 + ener2)/2.;
    }
    return ener;
}

double Material::freeze(double l1, double l2, double h2, double h1, Matrix<double, 1, 2> &gradient, Matrix<double, 2, 2> &hessian) const
{
    double delta = l1 - l2;
    double res = 0;
    gradient.setZero();
    hessian.setZero();

    double tmp = exp(-(delta * delta) * h2);
    res = 1 - tmp;
    double dfdd = delta * h2 * tmp * 2.0;
    double dfdd2 = h2 * tmp * 2.0 - (delta * delta) * (h2 * h2) * tmp * 4.0;

    gradient(0, 0) = dfdd;
    gradient(0, 1) = -dfdd;

    hessian(0, 0) = dfdd2;
    hessian(1, 1) = dfdd2;
    hessian(0, 1) = -dfdd2;
    hessian(1, 0) = -dfdd2;

    return res;
}

double Material::freezeIP(double l1, double l2, double h2, Matrix<double, 1, 6> &gradient, Matrix<double, 6, 6> &hessian) const
{
    double f;
    Matrix<double, 1, 2> df;
    Matrix<double, 2, 2> ddf;
    f = freeze(l1, -l2, 10,0, df, ddf);
    gradient.setZero();
    gradient(0, 0) = df(0);
    gradient(0, 1) = -df(1);
    hessian.setZero();
    hessian(0, 0) = ddf(0, 0);
    hessian(1, 1) = ddf(1, 1);
    hessian(0, 1) = -ddf(0, 1);
    hessian(1, 0) = -ddf(1, 0);
    return f;
}

double Material::freezeBEND(double l1, double l2, double h2, Matrix<double, 1, 6> &gradient, Matrix<double, 6, 6> &hessian) const
{
    double f;
    Matrix<double, 1, 2> df;
    Matrix<double, 2, 2> ddf;
    f = freeze(l1, l2, 10,0, df, ddf);
    gradient.setZero();
    gradient(0, 3) = df(0);
    gradient(0, 4) = df(1);
    hessian.setZero();
    hessian(3, 3) = ddf(0, 0);
    hessian(4, 4) = ddf(1, 1);
    hessian(3, 4) = ddf(0, 1);
    hessian(4, 3) = ddf(1, 0);
    return f;
}
