#include "CrossSection.h"
#include "NRKinematics.h"

void test(){

double m = 200.0;
double M = 1;
double v = 1.0;
double Lambda = 10.0;
double gs = 1.0;
double s = 1;

DM::CrossSectionModel model(Lambda,gs,s);
model.set_M(M);
model.set_m(m);
model.set_v(v);

model.xsec_costheta(0.5);

}

