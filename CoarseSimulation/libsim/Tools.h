#pragma once
#include "CommonIncludes.h"

#define MAXBUFSIZE ((int)1e6)

namespace Tools
{

    MatrixXd loadCsvDouble(const string &path);
    MatrixXi loadCsvInt(const string &path);
    void writeMatrix(const MatrixXd &m, const string &b);
    VectorXi unique(const VectorXi &v1);
    bool rowUnicity(const MatrixXd &m);
    MatrixXd deleteDuplicateRows(const MatrixXd &m);
    void printProgress(double percentage);
    MatrixXd devectorize(const VectorXd &vX);
    VectorXd vectorize(const MatrixXd &X);
    double deg2rad(double deg);
    double rad2deg(double rad);
    int rank(const SparseMatrix<double> &A);
    MatrixXd clamp(const MatrixXd &H);
    double modulo(double x, double y);

};
