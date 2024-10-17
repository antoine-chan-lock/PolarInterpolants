double sym_s(double epsilon, const VectorXd& vXrest, const VectorXd& vx) {
    const double x0 = -vXrest[0] + vXrest[2];
    const double x1 = vXrest[0]*vXrest[3] - vXrest[0]*vXrest[5] - vXrest[1]*vXrest[2] + vXrest[1]*vXrest[4] + vXrest[2]*vXrest[5] - vXrest[3]*vXrest[4];
    const double x2 = 1.0/x1;
    const double x3 = x2*(-vx[0] + vx[6]);
    const double x4 = -vXrest[0] + vXrest[4];
    const double x5 = -1/x1;
    const double x6 = x5*(-vx[0] + vx[3]);
    const double x7 = x0*x3 + x4*x6;
    const double x8 = vXrest[1] - vXrest[3];
    const double x9 = vXrest[1] - vXrest[5];
    const double x10 = x3*x8 + x6*x9;
    const double x11 = x2*(-vx[1] + vx[7]);
    const double x12 = x5*(-vx[1] + vx[4]);
    const double x13 = x0*x11 + x12*x4;
    const double x14 = x11*x8 + x12*x9;
    const double x15 = x2*(-vx[2] + vx[8]);
    const double x16 = x5*(-vx[2] + vx[5]);
    const double x17 = x0*x15 + x16*x4;
    const double x18 = x15*x8 + x16*x9;
    const double x19 = pow(x10, 2);
    const double x20 = pow(x14, 2);
    const double x21 = pow(x18, 2);
    const double x22 = (1.0/2.0)*pow(x13, 2) + (1.0/2.0)*pow(x17, 2) + (1.0/2.0)*pow(x7, 2);

    double out;
    out = sqrt((1.0/2.0)*x19 + (1.0/2.0)*x20 + (1.0/2.0)*x21 + x22 + sqrt(epsilon + pow(x10*x7 + x13*x14 + x17*x18, 2) + pow((1.0/2.0)*x19 + (1.0/2.0)*x20 + (1.0/2.0)*x21 - x22, 2))) - 1;
    return out;

}
