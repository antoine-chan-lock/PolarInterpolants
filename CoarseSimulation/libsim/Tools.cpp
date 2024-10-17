#include "Tools.h"
#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60
#include <igl/find.h>

MatrixXd Tools::clamp(const MatrixXd &H)
{
    EigenSolver<MatrixXd> es((H + H.transpose()) / 2);

    MatrixXd D = es.pseudoEigenvalueMatrix();
    MatrixXd V = es.pseudoEigenvectors();

    int nbDoFs = H.rows();

    for (int i = 0; i < nbDoFs; i++)
        if (D(i, i) < 0)
            D(i, i) = 0;
    MatrixXd res = V * D * V.inverse();
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

double Tools::modulo(double x, double y)
{
    return x - y * floor(x / y);
}

int Tools::rank(const SparseMatrix<double> &A)
{
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(A);
    double tolerance = 1e-6; // choose your tolerance wisely!
    int rank = 0;
    for (int i = 0; i < svd.singularValues().size(); i++)
    {
        if (svd.singularValues()(i) > tolerance)
            rank++;
    }
    return rank;
}

double Tools::deg2rad(double deg)
{
    return deg * M_PI / 180.0;
}

double Tools::rad2deg(double rad)
{
    return rad * 180.0 / M_PI;
}

MatrixXd Tools::devectorize(const VectorXd &vX)
{
    // 3D only
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

MatrixXd Tools::deleteDuplicateRows(const MatrixXd &m)
{
    VectorXi idDuplicate = VectorXi::Zero(m.rows());
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

bool Tools::rowUnicity(const MatrixXd &m)
{
    for (int i = 0; i < m.rows(); i++)
    {
        for (int j = i + 1; j < m.rows(); j++)
        {
            double dist = (m.row(i) - m.row(j)).norm();
            if (dist < 1e-6)
            {
                cout << "rowUnicity error" << endl;
                return false;
            }
        }
    }
    return true;
}

void Tools::printProgress(double percentage)
{
    int val = (int)(percentage * 100);
    int lpad = (int)(percentage * PBWIDTH);
    int rpad = PBWIDTH - lpad;
    printf("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
    fflush(stdout);
}

VectorXi Tools::unique(const VectorXi &v1)
{
    vector<int> v(v1.data(), v1.data() + v1.rows() * v1.cols());
    auto last = std::unique(v.begin(), v.end());
    v.erase(last, v.end());
    std::sort(v.begin(), v.end());
    last = std::unique(v.begin(), v.end());
    v.erase(last, v.end());
    Eigen::VectorXi v3 = Eigen::Map<Eigen::VectorXi, Eigen::Unaligned>(v.data(), v.size());
    return v3;
}

MatrixXd Tools::loadCsvDouble(const string &path)
{
    ifstream indata;
    indata.open(path);
    // if path not found
    if (!indata)
    {
        cout << "Error: File could not be opened" << endl;
        exit(1);
    }
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
    if (!indata)
    {
        cout << "Error: File could not be opened" << endl;
        exit(1);
    }
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
