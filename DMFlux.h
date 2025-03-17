#ifndef _DMFlux_h_
#define _DMFlux_h_

// NDimensional guassian

namespace DM {

class NDGaussian {

  public:
  
    NDGaussian(unsigned int ndim,std::vector<double> mean,std::vector<double> sigma);
    double Prob(std::vector<double> pos) const;

  private:

    const unsigned int _ndim;    
    const std::vector<double> _mean;
    const std::vector<double> _sigma;
    double _norm; 

};

// ND Guassian - note sigma is the uncertainty (not the covariance)
NDGaussian::NDGaussian(unsigned int ndim,std::vector<double> mean,std::vector<double> sigma) : _ndim(ndim), _mean(mean), _sigma(sigma)
{
  assert(_ndim == _mean.size() && _ndim == _sigma.size());
  _norm = 1.0;
  for(int i=0;i<_ndim;i++) _norm *= _sigma.at(i);
  _norm = 1.0/pow(2*3.1415,_ndim/2)/_norm;
}

double NDGaussian::Prob(std::vector<double> pos) const
{
  double arg = 0.0;
  for(int i=0;i<_ndim;i++) arg += (pos[i] - _mean[i])*(pos[i] - _mean[i])/_sigma[i]/_sigma[i];
  return exp(-arg/2.0)*_norm;
}


// DM velocity model - all units are km/s
class DMFlux {

  public:

    DMFlux(std::vector<double> vE);
    double GetFlux(std::vector<double> v) const;

  private:

    const double _vc = 220.0;
    const std::vector<double> _vE; 
    NDGaussian* _p_prob; 

};

DMFlux::DMFlux(std::vector<double> vE) : _vE(vE)
{

  std::vector<double> neg_ve,sigma;

  for(size_t i=0;i<_vE.size();i++){
    neg_ve.push_back(-_vE.at(i));
    sigma.push_back(_vc/sqrt(2));
  }

  _p_prob = new NDGaussian(3,neg_ve,sigma);

}

double DMFlux::GetFlux(std::vector<double> v) const 
{
  return _p_prob->Prob(v);
}

}

#endif
