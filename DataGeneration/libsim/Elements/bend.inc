double symx(double a, double c, double d, double k, double s, double threshold, const VectorXd& x) {
    const double x0 = cos(d);
    const double x1 = s + 1;
    const double x2 = cos(a);
    const double x3 = pow(x2, 2);
    const double x4 = sin(a);
    const double x5 = pow(x4, 2);
    const double x6 = 1 - c;
    const double x7 = sin(d);
    const double x8 = x1*x2*x4 - x2*x4*x6;
    const double x9 = (x0*x8 + x7*(x1*x5 + x3*x6))*x[1] + (x0*(x1*x3 + x5*x6) + x7*x8)*x[0];

    double out;
    if (threshold < fabs(k)) {
        out = -(x[2] - 1/k)*sin(k*x9);
    }
    else {
        out = -x9*(k*x[2] - 1);
    }
    return out;

}

void dsymxdx(double a, double c, double d, double k, double s, double threshold, const VectorXd& x, VectorXd& out) {
    const double x0 = cos(d);
    const double x1 = s + 1;
    const double x2 = cos(a);
    const double x3 = pow(x2, 2);
    const double x4 = sin(a);
    const double x5 = pow(x4, 2);
    const double x6 = 1 - c;
    const double x7 = sin(d);
    const double x8 = x1*x2*x4 - x2*x4*x6;
    const double x9 = x0*(x1*x3 + x5*x6) + x7*x8;
    const double x10 = x0*x8 + x7*(x1*x5 + x3*x6);
    const double x11 = x10*x[1] + x9*x[0];
    const double x12 = k*x11;
    const double x13 = k*(-x[2] + 1.0/k)*cos(x12);
    const double x14 = threshold < fabs(k);
    const double x15 = k*x[2] - 1;
    out.resize(3);
    if (x14) {
    
        out[0] = x13*x9;
    }
    else {
        out[0] = -x15*x9;
    }
    if (x14) {
        out[1] = x10*x13;
    }
    else {
        out[1] = -x10*x15;
    }
    if (x14) {
        out[2] = -sin(x12);
    }
    else {
        out[2] = -k*x11;
    }

}

void d2symxdx2(double a, double c, double d, double k, double s, double threshold, const VectorXd& x, MatrixXd& out) {
    const double x0 = cos(d);
    const double x1 = s + 1;
    const double x2 = cos(a);
    const double x3 = pow(x2, 2);
    const double x4 = sin(a);
    const double x5 = pow(x4, 2);
    const double x6 = 1 - c;
    const double x7 = sin(d);
    const double x8 = x1*x2*x4 - x2*x4*x6;
    const double x9 = x0*(x1*x3 + x5*x6) + x7*x8;
    const double x10 = x0*x8 + x7*(x1*x5 + x3*x6);
    const double x11 = k*(x10*x[1] + x9*x[0]);
    const double x12 = pow(k, 2)*(-x[2] + 1.0/k)*sin(x11);
    const double x13 = threshold < fabs(k);
    const double x14 = x13 ? -x10*x12*x9 : 0;
    const double x15 = k*cos(x11);
    const double x16 = x13 ? -x15*x9 : -k*x9;
    const double x17 = x13 ? -x10*x15 : -k*x10;
    out.resize(3, 3);
    if (x13) {
    
        out(0, 0) = -x12*pow(x9, 2);
    }
    else {
        out(0, 0) = 0;
    }
    out(0, 1) = x14;
    out(0, 2) = x16;
    out(1, 0) = x14;
    if (x13) {
        out(1, 1) = -pow(x10, 2)*x12;
    }
    else {
        out(1, 1) = 0;
    }
    out(1, 2) = x17;
    out(2, 0) = x16;
    out(2, 1) = x17;
    out(2, 2) = 0;

}

double symy(double a, double c, double d, double s, const VectorXd& x) {
    const double x0 = sin(d);
    const double x1 = s + 1;
    const double x2 = cos(a);
    const double x3 = pow(x2, 2);
    const double x4 = sin(a);
    const double x5 = pow(x4, 2);
    const double x6 = 1 - c;
    const double x7 = cos(d);
    const double x8 = x1*x2*x4 - x2*x4*x6;

    double out;
    out = (-x0*x8 + x7*(x1*x5 + x3*x6))*x[1] + (-x0*(x1*x3 + x5*x6) + x7*x8)*x[0];
    return out;

}

void dsymydx(double a, double c, double d, double s, VectorXd& out) {
    const double x0 = sin(d);
    const double x1 = s + 1;
    const double x2 = cos(a);
    const double x3 = pow(x2, 2);
    const double x4 = sin(a);
    const double x5 = pow(x4, 2);
    const double x6 = 1 - c;
    const double x7 = cos(d);
    const double x8 = x1*x2*x4 - x2*x4*x6;

    out.resize(3);
    out[0] = -x0*(x1*x3 + x5*x6) + x7*x8;
    out[1] = -x0*x8 + x7*(x1*x5 + x3*x6);
    out[2] = 0;

}

void d2symydx2(MatrixXd& out) {

    out.resize(3, 3);
    out(0, 0) = 0;
    out(0, 1) = 0;
    out(0, 2) = 0;
    out(1, 0) = 0;
    out(1, 1) = 0;
    out(1, 2) = 0;
    out(2, 0) = 0;
    out(2, 1) = 0;
    out(2, 2) = 0;

}

double symz(double a, double c, double d, double k, double s, double threshold, const VectorXd& x) {
    const double x0 = 1.0/k;
    const double x1 = cos(d);
    const double x2 = s + 1;
    const double x3 = cos(a);
    const double x4 = pow(x3, 2);
    const double x5 = sin(a);
    const double x6 = pow(x5, 2);
    const double x7 = 1 - c;
    const double x8 = sin(d);
    const double x9 = x2*x3*x5 - x3*x5*x7;
    const double x10 = (x1*x9 + x8*(x2*x6 + x4*x7))*x[1] + (x1*(x2*x4 + x6*x7) + x8*x9)*x[0];
    const double x11 = (1.0/2.0)*pow(x10, 2);

    double out;
    if (threshold < fabs(k)) {
        out = x0 + (-x0 + x[2])*cos(k*x10);
    }
    else {
        out = -pow(k, 2)*x11*x[2] + k*x11 + x[2];
    }
    return out;

}

void dsymzdx(double a, double c, double d, double k, double s, double threshold, const VectorXd& x, VectorXd& out) {
    const double x0 = cos(d);
    const double x1 = s + 1;
    const double x2 = cos(a);
    const double x3 = pow(x2, 2);
    const double x4 = sin(a);
    const double x5 = pow(x4, 2);
    const double x6 = 1 - c;
    const double x7 = x0*(x1*x3 + x5*x6);
    const double x8 = sin(d);
    const double x9 = x1*x2*x4 - x2*x4*x6;
    const double x10 = x8*x9;
    const double x11 = x10 + x7;
    const double x12 = x8*(x1*x5 + x3*x6);
    const double x13 = x0*x9;
    const double x14 = x12 + x13;
    const double x15 = x11*x[0] + x14*x[1];
    const double x16 = k*x15;
    const double x17 = k*(x[2] - 1/k)*sin(x16);
    const double x18 = threshold < fabs(k);
    const double x19 = 2*x10 + 2*x7;
    const double x20 = (1.0/2.0)*pow(k, 2);
    const double x21 = x15*x20*x[2];
    const double x22 = 2*x12 + 2*x13;
    out.resize(3);
    if (x18) {
    
        out[0] = -x11*x17;
    }
    else {
        out[0] = (1.0/2.0)*k*x15*x19 - x19*x21;
    }
    if (x18) {
        out[1] = -x14*x17;
    }
    else {
        out[1] = (1.0/2.0)*k*x15*x22 - x21*x22;
    }
    if (x18) {
        out[2] = cos(x16);
    }
    else {
        out[2] = -pow(x15, 2)*x20 + 1;
    }

}

void d2symzdx2(double a, double c, double d, double k, double s, double threshold, const VectorXd& x, MatrixXd& out) {
    const double x0 = cos(d);
    const double x1 = s + 1;
    const double x2 = cos(a);
    const double x3 = pow(x2, 2);
    const double x4 = sin(a);
    const double x5 = pow(x4, 2);
    const double x6 = 1 - c;
    const double x7 = x0*(x1*x3 + x5*x6);
    const double x8 = sin(d);
    const double x9 = x1*x2*x4 - x2*x4*x6;
    const double x10 = x8*x9;
    const double x11 = x10 + x7;
    const double x12 = pow(k, 2);
    const double x13 = x8*(x1*x5 + x3*x6);
    const double x14 = x0*x9;
    const double x15 = x13 + x14;
    const double x16 = x11*x[0] + x15*x[1];
    const double x17 = k*x16;
    const double x18 = x12*(x[2] - 1/k)*cos(x17);
    const double x19 = threshold < fabs(k);
    const double x20 = 2*x10 + 2*x7;
    const double x21 = (1.0/2.0)*x20;
    const double x22 = x12*x[2];
    const double x23 = x21*x22;
    const double x24 = x19 ? -x11*x15*x18 : (1.0/2.0)*k*x15*x20 - x15*x23;
    const double x25 = k*sin(x17);
    const double x26 = x12*x16;
    const double x27 = x19 ? -x11*x25 : -x21*x26;
    const double x28 = 2*x13 + 2*x14;
    const double x29 = (1.0/2.0)*x28;
    const double x30 = x19 ? -x15*x25 : -x26*x29;
    out.resize(3, 3);
    if (x19) {
    
        out(0, 0) = -pow(x11, 2)*x18;
    }
    else {
        out(0, 0) = (1.0/2.0)*k*x11*x20 - x11*x23;
    }
    out(0, 1) = x24;
    out(0, 2) = x27;
    out(1, 0) = x24;
    if (x19) {
        out(1, 1) = -pow(x15, 2)*x18;
    }
    else {
        out(1, 1) = (1.0/2.0)*k*x15*x28 - x15*x22*x29;
    }
    out(1, 2) = x30;
    out(2, 0) = x27;
    out(2, 1) = x30;
    out(2, 2) = 0;

}
