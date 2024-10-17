#include "Mesh.h"
#include "Tools.h"
#include <igl/triangle_triangle_adjacency.h>
#include <igl/writeOBJ.h>
#include <Eigen/CholmodSupport>

Mesh::Mesh()
{
}

Mesh::Mesh(string dirname, int matType, int idBC)
{
    m_X = Tools::loadCsvDouble(dirname + "/X.csv");
    m_T = Tools::loadCsvInt(dirname + "/T.csv");
    MatrixXd TTi;
    igl::triangle_triangle_adjacency(m_T, m_TT, TTi);
    m_matType = matType;
    m_idBC = idBC;
    m_material = Material("../Fitting/coefs/mat" + to_string(matType));
}

Mesh::~Mesh()
{
}

void Mesh::initElements()
{
    for (int i = 0; i < m_T.rows(); i++)
    {
        Element e(i, m_T, m_TT);
        m_elements.push_back(e);
    }
}

void Mesh::setMaterialInit(int x)
{
    m_material.setInit(x);
}

double Mesh::psiEdge(const MatrixXd &x, int elt, VectorXd *gradParam, SparseMatrix<double> *hessParam) const
{
    double res = 0;
    VectorXd grad = VectorXd::Zero(x.rows() * x.cols());
    vector<Triplet<double>> coefs;
    coefs.reserve(x.rows() * x.cols() * x.rows() * x.cols() * 3);
    Matrix<double, 1, 6> s;

    vector<Matrix<double, 1, 6>> v_S;
    vector<Matrix<double, 6, 6 * 3>> v_dS;
    vector<vector<Matrix<double, 6 * 3, 6 * 3>>> v_ddS;
    vector<Matrix<double, 6, 6>> v_kill;
    vector<VectorXi> v_dofs;
    // get self in plane strain
    Matrix<double, 6, 3> xflaps = m_elements[elt].getxflaps(x);
    Matrix<double, 3, 2> Xrest = m_elements[elt].getXrest(m_X);
    Matrix<double, 1, 6> Sself = Strains::symStrain(xflaps, m_elements[elt].m_niexists, Xrest);
    Matrix<double, 6, 6 * 3> dSdx = Strains::symGrad(xflaps, m_elements[elt].m_niexists, Xrest);
    vector<Matrix<double, 6 * 3, 6 * 3>> dSidxx = Strains::symHess(xflaps, m_elements[elt].m_niexists, Xrest);
    Matrix<double, 6, 6> kill;
    kill.setIdentity();
    kill(3, 3) = 0;
    kill(4, 4) = 0;
    kill(5, 5) = 0;
    VectorXi dofs(6 * 3), nodes;
    dofs.setOnes();
    dofs *= -1;
    nodes = m_elements[elt].getNodes();
    for (int j = 0; j < nodes.rows(); j++)
    {
        if (nodes(j) != -1)
        {
            dofs(j * 3) = nodes(j) * 3;
            dofs(j * 3 + 1) = nodes(j) * 3 + 1;
            dofs(j * 3 + 2) = nodes(j) * 3 + 2;
        }
    }
    v_dofs.push_back(dofs);
    v_S.push_back(Sself);
    v_dS.push_back(dSdx);
    v_ddS.push_back(dSidxx);
    v_kill.push_back(kill);

    // get neighbors bending strain
    Vector3i neighbors = m_elements[elt].m_neighbors;

    for (int j = 0; j < 3; j++)
    {
        if (neighbors(j) != -1)
        {
            xflaps = m_elements[neighbors(j)].getxflaps(x); 
            Xrest = m_elements[neighbors(j)].getXrest(m_X);
            Sself = Strains::symStrain(xflaps, m_elements[neighbors(j)].m_niexists, Xrest);
            dSdx = Strains::symGrad(xflaps, m_elements[neighbors(j)].m_niexists, Xrest);
            vector<Matrix<double, 6 * 3, 6 * 3>> slt = Strains::symHess(xflaps, m_elements[neighbors(j)].m_niexists, Xrest);
            kill.setZero();
            kill(3, 3) = 1;
            kill(4, 4) = 1;
            kill(5, 5) = 1; 
            nodes = m_elements[neighbors(j)].getNodes();
            dofs.setOnes();
            dofs *= -1;
            for (int k = 0; k < nodes.rows(); k++)
            {
                if (nodes(k) != -1)
                {
                    dofs(k * 3) = nodes(k) * 3;
                    dofs(k * 3 + 1) = nodes(k) * 3 + 1;
                    dofs(k * 3 + 2) = nodes(k) * 3 + 2;
                }
            }
            v_dofs.push_back(dofs);
            v_S.push_back(Sself);
            v_dS.push_back(dSdx);
            v_ddS.push_back(slt);
            v_kill.push_back(kill);
            break;
        }
    }

    s.setZero();
    for (int j = 0; j < v_S.size(); j++)
    {
        s += v_S[j] * v_kill[j];
    }
    Matrix<double, 1, 6> dener;
    Matrix<double, 6, 6> ddener;
    res += m_material.elasticEnergy(s, &dener, &ddener, 1);

    for (int j = 0; j < v_S.size(); j++)
    {
        VectorXd gradj = dener * v_kill[j] * v_dS[j];
        for (int k = 0; k < v_dofs[j].rows(); k++)
        {
            if (v_dofs[j](k) != -1)
            {
                grad(v_dofs[j](k)) += gradj(k);
            }
        }
    }

    Matrix<double, 6 * 3, 6 * 3> hessj = v_dS[0].transpose() * v_kill[0] * ddener * v_kill[0] * v_dS[0];
    VectorXd gradj = dener * v_kill[0];
    for (int kk = 0; kk < 6; kk++)
    {
        hessj += gradj(kk) * v_ddS[0][kk];
    }
    for (int k = 0; k < v_dofs[0].rows(); k++)
    {
        if (v_dofs[0](k) != -1)
        {
            for (int l = 0; l < v_dofs[0].rows(); l++)
            {
                if (v_dofs[0](l) != -1)
                {
                    coefs.push_back(Triplet<double>(v_dofs[0](k), v_dofs[0](l), hessj(k, l)));
                }
            }
        }
    }

    hessj = v_dS[1].transpose() * v_kill[1] * ddener * v_kill[1] * v_dS[1];
    gradj = dener * v_kill[1];
    for (int kk = 0; kk < 6; kk++)
    {
        hessj += gradj(kk) * v_ddS[1][kk];
    }
    for (int k = 0; k < v_dofs[1].rows(); k++)
    {
        if (v_dofs[1](k) != -1)
        {
            for (int l = 0; l < v_dofs[1].rows(); l++)
            {
                if (v_dofs[1](l) != -1)
                {
                    coefs.push_back(Triplet<double>(v_dofs[1](k), v_dofs[1](l), hessj(k, l)));
                }
            }
        }
    }

    hessj = v_dS[0].transpose() * v_kill[0] * ddener * v_kill[1] * v_dS[1];
    for (int k = 0; k < v_dofs[0].rows(); k++)
    {
        if (v_dofs[0](k) != -1)
        {
            for (int l = 0; l < v_dofs[1].rows(); l++)
            {
                if (v_dofs[1](l) != -1)
                {
                    coefs.push_back(Triplet<double>(v_dofs[0](k), v_dofs[1](l), hessj(k, l)));
                }
            }
        }
    }

    hessj = v_dS[1].transpose() * v_kill[1] * ddener * v_kill[0] * v_dS[0];
    for (int k = 0; k < v_dofs[1].rows(); k++)
    {
        if (v_dofs[1](k) != -1)
        {
            for (int l = 0; l < v_dofs[0].rows(); l++)
            {
                if (v_dofs[0](l) != -1)
                {
                    coefs.push_back(Triplet<double>(v_dofs[1](k), v_dofs[0](l), hessj(k, l)));
                }
            }
        }
    }

    vector<int> indexMap;
    int idmax = 0;
    for (int i = 0; i < coefs.size(); i++)
    {
        int idx = coefs[i].row();
        if (idx > idmax)
            idmax = idx;
        if (find(indexMap.begin(), indexMap.end(), idx) == indexMap.end())
        {
            indexMap.push_back(idx);
        }
    }

    vector<int> invMap(idmax + 1);
    for (int i = 0; i < indexMap.size(); i++)
    {
        invMap[indexMap[i]] = i;
    }
    MatrixXd hessDense = MatrixXd::Zero(indexMap.size(), indexMap.size());
    for (int i = 0; i < coefs.size(); i++)
    {
        int idx = coefs[i].row();
        int jdx = coefs[i].col();
        hessDense(invMap[idx], invMap[jdx]) += coefs[i].value();
    }
    hessDense = Tools::clamp(hessDense);
    vector<Triplet<double>> coefs2;
    for (int i = 0; i < hessDense.rows(); i++)
    {
        for (int j = 0; j < hessDense.cols(); j++)
        {
            coefs2.push_back(Triplet<double>(indexMap[i], indexMap[j], hessDense(i, j)));
        }
    }

    SparseMatrix<double> assembled(x.rows() * x.cols(), x.rows() * x.cols());
    assembled.setFromTriplets(coefs2.begin(), coefs2.end());

    (*gradParam) = grad;
    (*hessParam) = assembled;

    return res;
}

double Mesh::elasticEnergy(const MatrixXd &x, VectorXd *grad, SparseMatrix<double> *hess) const
{
    double ener = 0;
    if (grad and hess)
    {
        grad->setZero(x.rows() * x.cols());
        vector<Triplet<double>> coefs;
        coefs.reserve(x.rows() * x.cols() * x.rows() * x.cols() * 3);
        vector<double> energies(m_elements.size());
        vector<Matrix<double, 1, 6 * 3>> gradients(m_elements.size());
        vector<Matrix<double, 6 * 3, 6 * 3>> hessians(m_elements.size());

#pragma omp parallel for
        for (int i = 0; i < m_elements.size(); i++)
        {
            energies[i] = m_elements[i].elasticEnergy(x, m_X, &gradients[i], &hessians[i], m_material);
        }
        for (int i = 0; i < m_elements.size(); i++)
        {
            ener += energies[i];
            // assemble gradient
            VectorXi nodes = m_elements[i].getNodes();
            for (int j = 0; j < nodes.rows(); j++)
            {
                if (nodes(j) != -1)
                {
                    (*grad)(nodes(j) * 3) += gradients[i](j * 3);
                    (*grad)(nodes(j) * 3 + 1) += gradients[i](j * 3 + 1);
                    (*grad)(nodes(j) * 3 + 2) += gradients[i](j * 3 + 2);
                }
            }
            // assemble hessian
            for (int j = 0; j < nodes.rows(); j++)
            {
                if (nodes(j) == -1)
                {
                    continue;
                }
                for (int k = 0; k < nodes.rows(); k++)
                {
                    if (nodes(k) == -1)
                    {
                        continue;
                    }

                    for (int l = 0; l < 3; l++)
                    {
                        for (int m = 0; m < 3; m++)
                        {
                            coefs.push_back(Triplet<double>(nodes(j) * 3 + l, nodes(k) * 3 + m, hessians[i](j * 3 + l, k * 3 + m)));
                        }
                    }
                }
            }
        }

        SparseMatrix<double> assembled(x.rows() * x.cols(), x.rows() * x.cols());
        assembled.setFromTriplets(coefs.begin(), coefs.end());
        *hess = assembled;
    }
    else
    {
        vector<double> energies(m_elements.size());
#pragma omp parallel for
        for (int i = 0; i < m_elements.size(); i++)
            energies[i] = m_elements[i].elasticEnergy(x, m_X, NULL, NULL, m_material);
        for (int i = 0; i < m_elements.size(); i++)
            ener += energies[i];
    }

    if (grad and hess)
    {
        ener += BC(x);
        *grad += dBC(x);
        *hess += ddBC(x);
    }
    else
    {
        ener += BC(x);
    }

    return ener;
}

void Mesh::clampDoFs()
{
    vector<int> clampedDoFs;
    vector<double> clampedDoFsValues;

    double minx = m_X.col(0).minCoeff(), maxx = m_X.col(0).maxCoeff();

    double amplitude = -1;
    if (m_idBC == 1)
    {
        amplitude = (maxx - minx) * 0.125;
    }
    if (m_idBC == 2)
    {
        amplitude = (maxx - minx) * 0.25;
    }
    if (m_idBC == 3)
    {
        amplitude = (maxx - minx) * 0.5;
    }
    cout << "amplitude: " << amplitude << endl;

    for (int i = 0; i < m_X.rows(); i++)
    {

        if (m_idBC == 1)
        {
            if (abs(m_X(i, 0) - minx) < 1e-5)
            {
                clampedDoFs.push_back(i * 3);
                clampedDoFsValues.push_back(m_X(i, 0) - amplitude);
                clampedDoFs.push_back(i * 3 + 1);
                clampedDoFsValues.push_back(m_X(i, 1));
                clampedDoFs.push_back(i * 3 + 2);
                clampedDoFsValues.push_back(0);
            }
            if (abs(m_X(i, 0) - maxx) < 1e-5)
            {
                clampedDoFs.push_back(i * 3);
                clampedDoFsValues.push_back(m_X(i, 0));
                clampedDoFs.push_back(i * 3 + 1);
                clampedDoFsValues.push_back(m_X(i, 1));
                clampedDoFs.push_back(i * 3 + 2);
                clampedDoFsValues.push_back(0);
            }
        }
        else if (m_idBC == 2)
        {
            if (abs(m_X(i, 0) - minx) < 1e-5)
            {
                clampedDoFs.push_back(i * 3);
                clampedDoFsValues.push_back(m_X(i, 0) - amplitude);
                clampedDoFs.push_back(i * 3 + 1);
                clampedDoFsValues.push_back(m_X(i, 1));
                clampedDoFs.push_back(i * 3 + 2);
                clampedDoFsValues.push_back(0);
            }
            if (abs(m_X(i, 0) - maxx) < 1e-5)
            {
                clampedDoFs.push_back(i * 3);
                clampedDoFsValues.push_back(m_X(i, 0));
                clampedDoFs.push_back(i * 3 + 1);
                clampedDoFsValues.push_back(m_X(i, 1));
                clampedDoFs.push_back(i * 3 + 2);
                clampedDoFsValues.push_back(0);
            }
        }
        else if (m_idBC == 3)
        {
            if (abs(m_X(i, 0) - minx) < 1e-5)
            {
                clampedDoFs.push_back(i * 3);
                clampedDoFsValues.push_back(m_X(i, 0) - amplitude);
                clampedDoFs.push_back(i * 3 + 1);
                clampedDoFsValues.push_back(m_X(i, 1));
                clampedDoFs.push_back(i * 3 + 2);
                clampedDoFsValues.push_back(0);
            }
            if (abs(m_X(i, 0) - maxx) < 1e-5)
            {
                clampedDoFs.push_back(i * 3);
                clampedDoFsValues.push_back(m_X(i, 0));
                clampedDoFs.push_back(i * 3 + 1);
                clampedDoFsValues.push_back(m_X(i, 1));
                clampedDoFs.push_back(i * 3 + 2);
                clampedDoFsValues.push_back(0);
            }
        }

        else if (m_idBC == 4)
        {
            if (abs(m_X(i, 0) - minx) < 1e-5)
            {
                clampedDoFs.push_back(i * 3);
                clampedDoFsValues.push_back(m_X(i, 0));
                clampedDoFs.push_back(i * 3 + 1);
                clampedDoFsValues.push_back(m_X(i, 1));
                clampedDoFs.push_back(i * 3 + 2);
                clampedDoFsValues.push_back(0);
            }
            if (abs(m_X(i, 0) - maxx) < 1e-5)
            {
                clampedDoFs.push_back(i * 3);
                clampedDoFsValues.push_back(m_X(i, 0));
                clampedDoFs.push_back(i * 3 + 1);
                clampedDoFsValues.push_back(0);
                clampedDoFs.push_back(i * 3 + 2);
                clampedDoFsValues.push_back(m_X(i, 1));
            }
        }

        else if (m_idBC == 5)
        {
            if (abs(m_X(i, 0) - minx) < 3.)
            {
                clampedDoFs.push_back(i * 3);
                clampedDoFsValues.push_back(0);
                clampedDoFs.push_back(i * 3 + 1);
                clampedDoFsValues.push_back(m_X(i, 1));
                clampedDoFs.push_back(i * 3 + 2);
                clampedDoFsValues.push_back(abs(minx - m_X(i, 0)));
            }
            if (abs(m_X(i, 0) - maxx) < 3.)
            {
                clampedDoFs.push_back(i * 3);
                clampedDoFsValues.push_back(1);
                clampedDoFs.push_back(i * 3 + 1);
                clampedDoFsValues.push_back(m_X(i, 1));
                clampedDoFs.push_back(i * 3 + 2);
                clampedDoFsValues.push_back(abs(maxx - m_X(i, 0)));
            }
        }
        else
        {
            cout << "Dirichlet::wrong idBC" << endl;
            cout << m_idBC << endl;
            exit(1);
        }
    }
    m_clampedDoFs = VectorXi::Map(clampedDoFs.data(), clampedDoFs.size());
    m_clampedDoFsValues = VectorXd::Map(clampedDoFsValues.data(), clampedDoFsValues.size());
}

double Mesh::BC(const MatrixXd &x) const
{
    VectorXd vx = Tools::vectorize(x);
    double res = 0;
    for (int i = 0; i < m_clampedDoFs.rows(); i++)
    {
        res += m_penalty * pow(vx(m_clampedDoFs(i)) - m_clampedDoFsValues(i), 2);
    }
    return res;
}

VectorXd Mesh::dBC(const MatrixXd &x) const
{
    VectorXd vx = Tools::vectorize(x);
    VectorXd res = vx * 0;
    for (int i = 0; i < m_clampedDoFs.rows(); i++)
    {
        res(m_clampedDoFs(i)) = m_penalty * 2. * (vx(m_clampedDoFs(i)) - m_clampedDoFsValues(i));
    }
    return res;
}

SparseMatrix<double> Mesh::ddBC(const MatrixXd &x) const
{
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    tripletList.reserve(m_clampedDoFs.rows());
    for (int i = 0; i < m_clampedDoFs.rows(); i++)
    {
        tripletList.push_back(T(m_clampedDoFs(i), m_clampedDoFs(i), m_penalty * 2));
    }
    int nbDoFs = x.cols() * x.rows();
    SparseMatrix<double> mat(nbDoFs, nbDoFs);
    mat.setFromTriplets(tripletList.begin(), tripletList.end());
    return mat;
}

VectorXd Mesh::newtonRaphson(const VectorXd &vXinit, bool &success) const
{

    success = true;
    int maxIter, iter = 0;
    double maxError = 5e0;
    if (m_material.m_isInit == 1)
    {
        maxIter = 10000;
        maxError = 1e-3;
    }
    else if (m_material.m_isInit == 0)
    {
        maxIter = 40;
        maxError = 1e0;
    }
    else
    {
        cout << "NR::wrong isInit" << endl;
        exit(0);
    }
    bool linsearch_success = true;
    double error = 1e2;
    VectorXd vXsol = vXinit;
    int callsLinSearch = 0;
    m_material.m_smoothParameter = 100;
    cout << "smooth parameter: " << m_material.m_smoothParameter << endl;
    while ((iter < maxIter) and (error > maxError) and success)
    {
        VectorXd grad;
        SparseMatrix<double> hess;
        double energyOld = elasticEnergy(Tools::devectorize(vXsol), &grad, &hess);
        error = grad.norm();
        VectorXd vdx = solveLinearSystem(hess, grad);
        double dispmax = vdx.array().abs().maxCoeff();
        if (dispmax > 1)
        {
            vdx = vdx * 1 / dispmax;
        }

        cout << "iter: " << iter << " error: " << error << " energy: " << energyOld;
        VectorXd vXsoltmp = vXsol;
        vXsol = linsearch(vXsol, vdx, linsearch_success, 20, iter);
        callsLinSearch++;

        if (linsearch_success == false)
        {

            vXsol = linsearch(vXsoltmp, grad, linsearch_success, 40, iter);
            cout << "  GD";
        }
        double energyNew = elasticEnergy(Tools::devectorize(vXsol), NULL, NULL);
        cout << " energyImprovement: " << (energyNew - energyOld);
        iter++;
        cout << endl;
    }
    if (iter == maxIter)
    {
        success = false;
    }
    return vXsol;
}

VectorXd Mesh::linsearch(const VectorXd &vX, const VectorXd &vdx, bool &linsearch_success, int maxIter, int cpt) const
{
    double energy = elasticEnergy(Tools::devectorize(vX), NULL, NULL);
    VectorXd vXguess = vX - vdx;
    double energyGuess = elasticEnergy(Tools::devectorize(vXguess), NULL, NULL);
    energyGuess = isnan(energyGuess) ? 1e20 : energyGuess;
    linsearch_success = true;
    VectorXd vdxTemp = vdx;
    int iter = 0;
    double coef = 1;

    while (energyGuess > energy)
    {
        vdxTemp /= 2.;
        coef /= 2.;
        vXguess = vX - vdxTemp;
        energyGuess = elasticEnergy(Tools::devectorize(vXguess), NULL, NULL);
        energyGuess = isnan(energyGuess) ? 1e20 : energyGuess; // can be nan because the BC penalty is strict and generates big gradients
        iter++;
        if (iter > maxIter)
        {
            linsearch_success = false;
            break;
        }
    }

    cout << "  coef: " << coef;

    return vXguess;
}

VectorXd Mesh::solveLinearSystem(SparseMatrix<double> &A, const VectorXd &b) const
{

    CholmodSupernodalLLT<SparseMatrix<double>, Eigen::Lower> solver;
    solver.analyzePattern(A);
    solver.factorize(A);
    VectorXd res = solver.solve(b);
    return res;

    // VectorXd res;
    // Eigen::SparseQR<SparseMatrix<double>, COLAMDOrdering<int>> solver3;
    // // A.makeCompressed();
    // solver3.compute(A);
    // if (solver3.info() == Eigen::Success)
    // {
    //     res = solver3.solve(b);
    // }
    // else
    // {
    //     cout << "Solver::solveLinearSystem::fail" << endl;
    //     exit(0);
    // }

    return res;
}
