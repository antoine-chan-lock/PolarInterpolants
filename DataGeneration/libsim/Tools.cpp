#include "Tools.h"
#include <igl/find.h>

MatrixXd Tools::loadCsvDouble(const string &path)
{
    ifstream indata;
    indata.open(path);
    string line;
    vector<double> values;
    int rows = 0;
    while (getline(indata, line))
    {
        stringstream lineStream(line);
        string cell;
        while (getline(lineStream, cell, ','))
        {
            values.push_back(stod(cell));
        }
        ++rows;
    }
    return Map<const Matrix<double, Dynamic, Dynamic, RowMajor>>(values.data(), rows, values.size() / rows);
}

MatrixXi Tools::loadCsvInt(const string &path)
{
    ifstream indata;
    indata.open(path);
    string line;
    vector<int> values;
    int rows = 0;
    while (getline(indata, line))
    {
        stringstream lineStream(line);
        string cell;
        while (getline(lineStream, cell, ','))
        {
            values.push_back(stod(cell));
        }
        ++rows;
    }
    return Map<const Matrix<int, Dynamic, Dynamic, RowMajor>>(values.data(), rows, values.size() / rows);
}

void Tools::writeMatrix(const MatrixXd &m, const string &b)
{
    const static IOFormat CSVFormat(FullPrecision, DontAlignCols, ", ", "\n");
    ofstream file(b);
    file << m.format(CSVFormat);
}

MatrixXd Tools::deleteDuplicateRows(const MatrixXd &m)
{
    VectorXi idDuplicate = VectorXi::Zero(m.rows());
#pragma omp parallel for
    for (int i = 0; i < m.rows(); i++)
    {
        for (int j = i + 1; j < m.rows(); j++)
        {
            double dist = (m.row(i) - m.row(j)).norm();
            if (dist < 1e-6)
            {
                idDuplicate(j) = 1;
            }
        }
    }
    VectorXi ids;
    igl::find(idDuplicate.array() == 0, ids);
    return m(ids, all);
}

MatrixXd Tools::devectorize(const VectorXd &vX)
{
    int nbCols = 3;
    int nbRows = vX.rows() / nbCols;
    MatrixXd X = MatrixXd::Zero(nbRows, nbCols);
    X.col(0) = vX(seq(0, last, 3));
    X.col(1) = vX(seq(1, last, 3));
    X.col(2) = vX(seq(2, last, 3));
    return X;
}

VectorXd Tools::vectorize(const MatrixXd &X)
{
    int nbDoFs = X.rows() * X.cols();
    VectorXd vx = VectorXd::Zero(nbDoFs);
    vx(seq(0, last, 3)) = X.col(0);
    vx(seq(1, last, 3)) = X.col(1);
    vx(seq(2, last, 3)) = X.col(2);
    return vx;
}

double Tools::deg2rad(double deg)
{
    return deg * M_PI / 180.0;
}

double Tools::rad2deg(double rad)
{
    return rad * 180.0 / M_PI;
}