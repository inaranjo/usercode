#include <math.h> 
#include "TMath.h" 
#include <limits>


class ratioEfficiencyMu {
public:
  ratioEfficiencyMu(){} ;
  double efficiency(double m, double m0, double sigma, double alpha, double n, double norm) const ;
  double mcEfficiency(double pt, bool EB=true) const ;
} ;

double ratioEfficiencyMu::efficiency(double m, double m0, double sigma, double alpha, double n, double norm) const 
 { 
   const double sqrtPiOver2 = 1.2533141373;
   const double sqrt2 = 1.4142135624;

   double sig = fabs((double) sigma);
   
   double t = (m - m0)/sig ;
   
   if (alpha < 0)
     t = -t;

   double absAlpha = fabs(alpha / sig);
   double a = TMath::Power(n/absAlpha,n)*exp(-0.5*absAlpha*absAlpha);
   double b = absAlpha - n/absAlpha;

   if (a>=std::numeric_limits<double>::max()) return -1. ;

   double ApproxErf ;
   double arg = absAlpha / sqrt2 ;
   if (arg > 5.) ApproxErf = 1 ;
   else if (arg < -5.) ApproxErf = -1 ;
   else ApproxErf = erf(arg) ;

   double leftArea = (1 + ApproxErf) * sqrtPiOver2 ;
   double rightArea = ( a * 1/TMath::Power(absAlpha - b,n-1)) / (n - 1);
   double area = leftArea + rightArea;


   if ( t <= absAlpha ){
     arg = t / sqrt2 ;
     if (arg > 5.) ApproxErf = 1 ;
     else if (arg < -5.) ApproxErf = -1 ;
     else ApproxErf = erf(arg) ;
     return norm * (1 + ApproxErf) * sqrtPiOver2 / area ;
   }
   else{
     return norm * (leftArea +  a * (1/TMath::Power(t-b,n-1) - 1/TMath::Power(absAlpha - b,n-1)) / (1 - n)) / area ;
   }
  
 } 


double ratioEfficiencyMu::mcEfficiency(double pt, bool EB) const
{
  double ratio_ = 0 ;
  if (pt<15) return ratio_ ;

  double meanMC,sigmaMC,alphaMC,nMC,normMC ;

  if (EB) {
    meanMC=16.99389526 ; sigmaMC=-0.04080190 ; alphaMC=0.00794730 ; nMC=1.60377906 ; normMC=0.99626161;
  } else {
    meanMC=16.99065795 ; sigmaMC=-0.11993730 ; alphaMC=0.01384991 ; nMC=2.38867304 ; normMC=0.86552275;
  }
  double effMC = efficiency(pt,meanMC,sigmaMC,alphaMC,nMC,normMC) ;

  ratio_ = effMC ;
  return ratio_ ;
}
