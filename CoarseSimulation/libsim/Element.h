#pragma once
#include "CommonIncludes.h"
#include "Strains.h"
#include "Material.h"
#include <igl/slice_mask.h>

class Element
{
public:
    Element();
    Element(int id, const MatrixXi &T, const MatrixXi &TT);
    ~Element();

    VectorXi getNodes() const { return m_T; };
    Matrix<double, 6, 3> getxflaps(const MatrixXd &x) const;
    Matrix<double, 3, 2> getXrest(const MatrixXd &X) const;
    double getArea(const MatrixXd &X) const;
    Matrix<double, 6 * 3, 6 * 3> clamp(const Matrix<double, 6 * 3, 6 * 3> &H) const;
    double elasticEnergy(const MatrixXd &x, const MatrixXd &X, Matrix<double, 1, 6 * 3> *grad, Matrix<double, 6 * 3, 6 * 3> *hess, const Material &mat) const;

    // protected:
    int m_id;
    VectorXi m_T;
    VectorXi m_idDofs;
    Vector3d m_niexists;
    double m_isIn;
    VectorXi m_neighbors;
};