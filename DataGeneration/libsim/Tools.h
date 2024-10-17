#pragma once
#include "CommonIncludes.h"

#define MAXBUFSIZE ((int)1e6)

namespace Tools
{

    MatrixXd loadCsvDouble(const string &path);
    MatrixXi loadCsvInt(const string &path);
    void writeMatrix(const MatrixXd &m, const string &b);
    MatrixXd deleteDuplicateRows(const MatrixXd &m);
    MatrixXd devectorize(const VectorXd &vX);
    VectorXd vectorize(const MatrixXd &X);
    double deg2rad(double deg);
    double rad2deg(double rad);

};
