#pragma once
#include "CommonIncludes.h"
#include "Tools.h"

namespace Strains
{

    Matrix<double, 1, 6> symStrain(const Matrix<double, 6, 3> &xflaps, const Vector3d &niexists, const Matrix<double, 3, 2> &Xrest);
    Matrix<double, 6, 6 * 3> symGrad(const Matrix<double, 6, 3> &xflaps, const Vector3d &niexists, const Matrix<double, 3, 2> &Xrest);
    vector<Matrix<double, 6 * 3, 6 * 3>> symHess(const Matrix<double, 6, 3> &xflaps, const Vector3d &niexists, const Matrix<double, 3, 2> &Xrest);
}