#ifndef _CrossSection_h_
#define _CrossSection_h_

#include "NRKinematics.h"

namespace DM {

const double pi = 3.14159;

// Speed of light in km/s
const double c = 3e5; 

// Reduced Planck's constant in GeVs
const double hbar = 6.582e-19;


// Legendre Polynomial calculator
class LegendrePolynomial {

  public:
    LegendrePolynomial(int n);

    double (LegendrePolynomial::*eval)(double) const;
    double Eval(double x) const;       

  private:
    int _n;

    double Eval0(double x) const;
    double Eval1(double x) const;
    double Eval2(double x) const;
    
};

LegendrePolynomial::LegendrePolynomial(int n) : _n(n)
{
  if(_n == 0) eval = &LegendrePolynomial::Eval0;
  else if(_n == 1) eval = &LegendrePolynomial::Eval1;
  else if(_n == 2) eval = &LegendrePolynomial::Eval2;
  else std::cout << "Invalid spin, spin of 0,1,2 supported" << std::endl;
}

double LegendrePolynomial::Eval(double x) const
{
  return (this->*eval)(x);
}

double LegendrePolynomial::Eval0(double x) const
{
  return 1;
}

double LegendrePolynomial::Eval1(double x) const
{
  return x;
}

double LegendrePolynomial::Eval2(double x) const
{
  return 0.5*(3*x*x-1);
}

class CrossSectionModel
{

  public:
    CrossSectionModel(double Lambda, double gs, int s);
    
    double xsec_Er(const double& costheta) const;
    double xsec_costheta(const double& costheta) const;
    void set_m(double m){ _m = m; };
    void set_M(double M){ _M = M; };
    void set_v(double v){ _v = v; };
  
  private:
    const double _Lambda;
    const double _gs;
    const int _s;
    const LegendrePolynomial _lp;
 
    double _m;
    double _M;
    double _v;

};

CrossSectionModel::CrossSectionModel(double Lambda, double gs, int s) : _Lambda(Lambda), _gs(gs), _s(s), _lp(s)
{

}

// Cross section with respect to costheta (costheta in the lab frame)
double CrossSectionModel::xsec_costheta(const double& costheta) const
{

  if(costheta < -1 || costheta > 1.0){
    std::cout << "Invalid costheta of " << costheta << " returning xsec of 0" << std::endl;
    return 0;
  }

  const NRKinematics kin(_m,_M,_v);
  double jacobian = 2*_M*kin.v_N(costheta)*kin.v_N(costheta)/costheta; 
  //std::cout << "jacobian = " << jacobian << std::endl;

  return jacobian*xsec_Er(costheta);

}

// Cross section with respect to nuclear recoil energy in the lab frame 
double xsec_Er(const double& costheta) const
{

  const NRKinematics kin(_m,_M,_v);
  TVector3 p = kin.p(costheta);
  TVector3 k = kin.k(costheta);
  TVector3 pprime = kin.pprime(costheta);
  TVector3 kprime = kin.kprime(costheta);

}

}

#endif
