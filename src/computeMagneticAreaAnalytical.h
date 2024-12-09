#include<iostream>
#include<fstream>
#include<cmath>
#include <vector>
#include <algorithm>
#include <string>
#include <stdexcept>
#include <iostream>


struct FieldParameters{
  double amplitude;
  double frequency;
};

double computeMaxArea(double badim) {
  double weight = exp(-0.5 / pow(badim, 1.8));
  double logarg = log(1 + 0.1 * pow(badim, 0.9) + 0.056 * pow(badim, 1.65) + 1.54e-05 * pow(badim, 3.21));
  double maxA = M_PI * badim / 6 + weight * (4 - M_PI * badim / 6 - 4 / (1 + logarg));
  return maxA;
}

double computeAreaAdim(double wadim, double badim){
  double p0        = - 0.45 * exp(-8.5 / pow(badim, 1.75)) + 1;
  double amax      = computeMaxArea(badim);
  double wmax      = 1 + 0.287 * pow(badim, 0.8915) / (1 + 5.52 * pow(badim, -1.16));
  double slope_inf = M_PI * badim / (3 * wmax);
  double slope_0   = M_PI * badim / (3 * pow(1 + 0.16 * pow(badim, 2.05), 1.0 / 2.05));

  double p1        = log(slope_inf / amax) / log((p0 + 1) / p0);
  double p2        = (p0 + 1) / p1;
  double p3        = 1.8;
  double p4        = (1 - p0) / p3;
  
  double c1        = p0 * pow(amax, -1 / p1) / (p1 * p2);
  double c2        = c1 / p0;
  double c3        = pow(1 / (wmax * slope_0 * pow(c2, p1)), 1 / p4);
  double c4        = 1660 * exp(-6.33 * pow(wadim, 0.18)) * exp(-11 / badim) * pow(badim, -0.5);
  double num       = (wadim / wmax) / pow(c3 + (1 + c4 ) * pow(wadim / wmax, p3), p4);
  
  return num / pow((c2 + c1 * pow(wadim / wmax, p2)), p1);
}

double computeArea(FieldParameters fb,
                   double coreRadius,
                   double msat,
                   double coatingWidth,
                   double viscosity,
                   double kBT){

  double hydroRadius  = coreRadius + coatingWidth;
  double hydroRadius3 = hydroRadius * hydroRadius * hydroRadius;
  double coreRadius3  = coreRadius * coreRadius * coreRadius;
  double m0           = 4./3.*M_PI*coreRadius3 * msat;
  double wadim        = 2*M_PI*fb.frequency*4*M_PI*viscosity*hydroRadius3/kBT;
  double badim        = m0*fb.amplitude/kBT;

  return computeAreaAdim(wadim, badim);
}

std::vector<double> computeArea(std::vector<FieldParameters>& fb,
                                double coreRadius,
                                double msat,
                                double coatingWidth,
                                double viscosity,
                                double kBT) {
    std::vector<double> areas(fb.size());
    
    std::transform(fb.begin(), fb.end(), areas.begin(),
                   [&](const FieldParameters& fbp) {
                       return computeArea(fbp, coreRadius, msat, coatingWidth, viscosity, kBT);
                   });
    
    return areas;
}
  
//   /* if (std_coreRadius >0.0){ */
//   /*   std_coreRadius = std::min(std_coreRadius, 0.9); */
//   /*   std_coreRadius = std::max(std_coreRadius, 0.01); */
//   /* } */
  
//   /* std_coreRadius*=coreRadius; */
  
//   double f = fb.frequency;
//   double b = fb.amplitude;
//   double area;
  
//   static double last_coreRadius     = std::numeric_limits<double>::quiet_NaN();
//   static double last_std_coreRadius = std::numeric_limits<double>::quiet_NaN();
//   static double coreRadius3_mean    = 0.0;
  
//   if (std_coatingWidth == 0.0 && std_coreRadius == 0.0) {
//     double coreRadius3  = coreRadius*coreRadius*coreRadius;
//     double hydroRadius3 = (coreRadius + coatingWidth)*(coreRadius + coatingWidth)*(coreRadius + coatingWidth);
//     double wadim        = 8 * M_PI * f * M_PI * viscosity * hydroRadius3 / temperature;
//     double badim        = 4 * M_PI * b * msat * coreRadius3 / (3 * temperature);
//     area  = computeAreaAdim(wadim, badim);
//   } else if (std_coatingWidth == 0 && std_coreRadius != 0.0) {
//     std::cout<<"TODO\n";
//     return area*normalization;
// }
