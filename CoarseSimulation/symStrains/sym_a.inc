double sym_a(double epsilon, const VectorXd& vXrest, const VectorXd& vx) {
    const double x0 = -vx[0] + vx[6];
    const double x1 = vXrest[0]*vXrest[3] - vXrest[0]*vXrest[5] - vXrest[1]*vXrest[2] + vXrest[1]*vXrest[4] + vXrest[2]*vXrest[5] - vXrest[3]*vXrest[4];
    const double x2 = 1.0/x1;
    const double x3 = x2*(-vXrest[0] + vXrest[2]);
    const double x4 = -vx[0] + vx[3];
    const double x5 = -1/x1;
    const double x6 = x5*(-vXrest[0] + vXrest[4]);
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
    const double x17 = x2*(vXrest[1] - vXrest[3]);
    const double x18 = x5*(vXrest[1] - vXrest[5]);
    const double x19 = x0*x17 + x18*x4;
    const double x20 = pow(x19, 2);
    const double x21 = x10*x18 + x17*x9;
    const double x22 = pow(x21, 2);
    const double x23 = x13*x17 + x14*x18;
    const double x24 = pow(x23, 2);
    const double x25 = x12 + x16 - x20 - x22 - x24 + x8 <= 0;
    const double x26 = x19*x7;
    const double x27 = x11*x21;
    const double x28 = x15*x23;
    const double x29 = x26 + x27 + x28;
    const double x30 = (1.0/2.0)*x8;
    const double x31 = (1.0/2.0)*x12;
    const double x32 = (1.0/2.0)*x16;
    const double x33 = (1.0/2.0)*x20;
    const double x34 = (1.0/2.0)*x22;
    const double x35 = (1.0/2.0)*x24;
    const double x36 = x30 + x31 + x32 - x33 - x34 - x35;
    const double x37 = sqrt(epsilon + pow(x29, 2) + pow(x36, 2));
    const double x38 = x36 + x37;
    const double x39 = ((x25) ? (
        1.0*x26 + 1.0*x27 + 1.0*x28 < 0
    )
    : (
        x38 < 0
    )) ? -1.0 : 1.0;
    const double x40 = x25 ? 1.0 : 0.0;
    const double x41 = 1 - x40;

    double out;
    out = atan2(x39*(x29*x40 + x38*x41), x39*(x29*x41 + x40*(-x30 - x31 - x32 + x33 + x34 + x35 + x37)));
    return out;

}
