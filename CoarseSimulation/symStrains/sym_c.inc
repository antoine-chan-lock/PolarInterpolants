double sym_c(double epsilon, const VectorXd& vXrest, const VectorXd& vx) {
    const double x0 = -vx[0] + vx[6];
    const double x1 = vXrest[0]*vXrest[3] - vXrest[0]*vXrest[5] - vXrest[1]*vXrest[2] + vXrest[1]*vXrest[4] + vXrest[2]*vXrest[5] - vXrest[3]*vXrest[4];
    const double x2 = 1.0/x1;
    const double x3 = x2*(vXrest[1] - vXrest[3]);
    const double x4 = -vx[0] + vx[3];
    const double x5 = -1/x1;
    const double x6 = x5*(vXrest[1] - vXrest[5]);
    const double x7 = x0*x3 + x4*x6;
    const double x8 = pow(x7, 2);
    const double x9 = -vx[1] + vx[7];
    const double x10 = -vx[1] + vx[4];
    const double x11 = x10*x6 + x3*x9;
    const double x12 = pow(x11, 2);
    const double x13 = -vx[2] + vx[8];
    const double x14 = -vx[2] + vx[5];
    const double x15 = x13*x3 + x14*x6;
    const double x16 = pow(x15, 2);
    const double x17 = x2*(-vXrest[0] + vXrest[2]);
    const double x18 = x5*(-vXrest[0] + vXrest[4]);
    const double x19 = x0*x17 + x18*x4;
    const double x20 = x10*x18 + x17*x9;
    const double x21 = x13*x17 + x14*x18;
    const double x22 = (1.0/2.0)*pow(x19, 2) + (1.0/2.0)*pow(x20, 2) + (1.0/2.0)*pow(x21, 2);

    double out;
    out = 1 - sqrt((1.0/2.0)*x12 + (1.0/2.0)*x16 + x22 + (1.0/2.0)*x8 - sqrt(epsilon + pow(x11*x20 + x15*x21 + x19*x7, 2) + pow((1.0/2.0)*x12 + (1.0/2.0)*x16 - x22 + (1.0/2.0)*x8, 2)));
    return out;

}
