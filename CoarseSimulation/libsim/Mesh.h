#pragma once
#include "CommonIncludes.h"
#include "Element.h"

class Mesh
{
public:
    Mesh();
    Mesh(string dirname, int matType, int idBC);
    ~Mesh();

    // getters
    MatrixXd getRest() const { return m_X; }

    // init
    void initElements();
    void setMaterialInit(int x);

    // material energy and derivatives
    // The Shape operator is Ill-defined for edge elements. In this case, we take the average curavture of the neighboring elements
    // psiEdge is the energy of the edge element
    double psiEdge(const MatrixXd &x, int elt, VectorXd *gradParam, SparseMatrix<double> *hessParam) const;
    // elasticEnergy is the energy non edge elements
    double elasticEnergy(const MatrixXd &x, VectorXd *grad, SparseMatrix<double> *hess) const;

    // boundary conditions
    void clampDoFs(); // Sets the clamped DoFs and their values
    double BC(const MatrixXd &x) const; // Penalty energy
    VectorXd dBC(const MatrixXd &x) const;
    SparseMatrix<double> ddBC(const MatrixXd &x) const;

    // solvers
    VectorXd newtonRaphson(const VectorXd &vXinit, bool &success) const;
    VectorXd linsearch(const VectorXd &vX, const VectorXd &vdx, bool &linsearch_success, int maxIter, int cpt) const;
    VectorXd solveLinearSystem(SparseMatrix<double> &A, const VectorXd &b) const;
    
    vector<Element> m_elements;
    MatrixXd m_X;
    MatrixXi m_T, m_TT;
    int m_matType, m_idBC;
    VectorXi m_clampedDoFs;
    VectorXd m_clampedDoFsValues;
    double m_penalty = 1e8;
    mutable Material m_material;
    VectorXd m_vkill;
};