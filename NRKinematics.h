#ifndef _NRKinematics_h_
#define _NRKinematics_h_

#include "TVector3.h"

namespace DM {

class NRKinematics 
{

  public:

    NRKinematics(double m, double M, double v);

    double v_N(const double& costheta) const;
    double vprime(const double& costheta) const;
    double cosalpha(const double& costheta) const;
    double Er(const double& costheta) const;
    TVector3 p(const double& costheta) const;
    TVector3 k(const double& costheta) const;
    TVector3 pprime(const double& costheta) const;
    TVector3 kprime(const double& costheta) const;
    TVector3 phat(const double& costheta) const;
    TVector3 khat(const double& costheta) const;
    TVector3 pprimehat(const double& costheta) const;
    TVector3 kprimehat(const double& costheta) const;

    void gallilean_boost(TVector3& p,TVector3 boost, double m) const;
    TVector3 cms_velocity(TVector3& p1,TVector3& p2, double m1, double m2) const;

  private:
    
    const double _m;
    const double _M;
    const double _v;    

};	

NRKinematics::NRKinematics(double m, double M, double v) : _m(m), _M(M), _v(v)
{

}

double NRKinematics::v_N(const double& costheta) const 
{
  return 2*(_M*_v*costheta)/(_M + _M*_M/_m);
}


double NRKinematics::vprime(const double& costheta) const
{
  double vn = v_N(costheta);
  return sqrt(_v*_v - _M/_m*vn*vn);
}

double NRKinematics::cosalpha(const double& costheta) const
{
  return (_m*_v - _M*v_N(costheta)*costheta)/(_m*vprime(costheta));
}

double NRKinematics::Er(const double& costheta) const
{
  return 0.5*_M*v_N(costheta)*v_N(costheta);
}

TVector3 NRKinematics::p(const double& costheta) const
{
  return TVector3(0,0,0);
}

TVector3 NRKinematics::k(const double& costheta) const
{
  return TVector3(0,0,_m*_v);
}

TVector3 NRKinematics::pprime(const double& costheta) const
{
  return TVector3(0,_M*v_N(costheta)*sqrt(1-costheta*costheta),_M*v_N(costheta)*costheta);
}

TVector3 NRKinematics::kprime(const double& costheta) const
{
  double ca = cosalpha(costheta);
  double sa = _M*v_N(costheta)*sqrt(1-costheta*costheta)/_m/vprime(costheta);
  return TVector3(0,-_m*vprime(costheta)*sa,_m*vprime(costheta)*ca);
}

TVector3 NRKinematics::phat(const double& costheta) const
{
  return TVector3(0,0,-1);
}

TVector3 NRKinematics::khat(const double& costheta) const 
{
  return TVector3(0,0,1);
}

TVector3 NRKinematics::pprimehat(const double& costheta) const
{
  return TVector3(0,sqrt(1-costheta*costheta),costheta);
}

TVector3 NRKinematics::kprimehat(const double& costheta) const
{
  double calpha = cosalpha(costheta);
  return TVector3(0,sqrt(1-calpha*calpha),calpha);
}

void NRKinematics::gallilean_boost(TVector3& p,TVector3 boost, double m) const
{
  p[0] += boost.X()*m;
  p[1] += boost.Y()*m;
  p[2] += boost.Z()*m;
}

TVector3 NRKinematics::cms_velocity(TVector3& p1,TVector3& p2, double m1, double m2) const
{
  return (1.0/(m1 + m2))*(p1 + p2);
}

}

#endif
