#pragma once
#include "Tools.h"
#include "CommonIncludes.h"
#include "Solver.h"
#include <thread>

class DataGeneration
{
public:
    DataGeneration();
    DataGeneration(const string &matid);
    ~DataGeneration();
    void sampling();
    Matrix<double,1,6> derivativesSafe(double s,double c,double a,double k, double d, double h,double &errorConvergence) const;
    void genDerivatives() const;
    Vector2d findMinEnergy(double s,double &c, double a,double &cmin,double &cmax) const;
    Vector2d orthogonalResponse(double s,double a) const;
    void orthRespProfile(double s,double a) const;

private:
    string m_matdir,m_matid;
    Solver m_slv;
    MatrixXd  m_scakd;
};
