#pragma once
#include<iostream>
#include<fstream>
#include<cmath>
#include <vector>
#include <algorithm>
#include <string>
#include <stdexcept>
#include <iostream>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/distributions/lognormal.hpp>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

struct FieldParameters{
  double amplitude;
  double frequency;
};

struct IntegrandParams {
  double sigma;
  double mu;
  double coatingWidth;
  double viscosity;
  double kBT;
  double msat;
  FieldParameters fb;
};

double computeMaxArea(double badim) {
  double weight = exp(-0.5 / pow(badim, 1.8));
  double logarg = log(1 + 0.1 * pow(badim, 0.9) + 0.056 * pow(badim, 1.65) + 1.54e-05 * pow(badim, 3.21));
  double maxA   = M_PI * badim / 6 + weight * (4 - M_PI * badim / 6 - 4 / (1 + logarg));
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
  double coreRadius3  = coreRadius  *  coreRadius * coreRadius;
  double m0           = 4./3. * M_PI * coreRadius3 * msat;
  double wadim        = 8*M_PI*fb.frequency*M_PI*viscosity*hydroRadius3/kBT;
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

double lognormal_pdf(double x, double sigma, double mu) {
  if (x <= 0.0) return 0.0;
  return (1.0 / (x * sigma * std::sqrt(2.0 * M_PI))) *
    std::exp(-std::pow(std::log(x) - mu, 2) / (2.0 * sigma * sigma));
}

double area_integrand(double coreRadius, void* paramsIn) {
  const IntegrandParams* p = static_cast<const IntegrandParams*>(paramsIn);
  double area = computeArea(p->fb, coreRadius, p->msat, p->coatingWidth,
                            p->viscosity, p->kBT);
  double rc3 = coreRadius * coreRadius * coreRadius;
  return  area * lognormal_pdf(coreRadius, p->sigma, p->mu) * rc3;
}

double computeArea(FieldParameters fb, double meanCoreRadius, double msat,
                   double coatingWidth, double viscosity, double kBT, double stdCoreRadius) {

  
    double variance        = stdCoreRadius * stdCoreRadius;
    double sigma           = std::sqrt(std::log(1 + variance / (meanCoreRadius * meanCoreRadius)));
    double mu              = std::log(meanCoreRadius) - 0.5 * sigma * sigma;

    double meanCoreRadius3 = std::exp(3 * mu + 4.5 * sigma * sigma);
    IntegrandParams params = {sigma, mu, coatingWidth, viscosity, kBT, msat, fb};
    gsl_function F;

    F.function = [](double x, void* p) { return area_integrand(x, p); };
    F.params   = &params;

    gsl_integration_workspace* workspace = gsl_integration_workspace_alloc(1000);
    
    double result, error;
    
    if (stdCoreRadius/meanCoreRadius<0.01){
      result = computeArea(fb, meanCoreRadius, msat, coatingWidth, viscosity, kBT)*meanCoreRadius3;
    } else if (stdCoreRadius/meanCoreRadius>0.75){
      double lowerLimit = 0;
      double upperLimit = meanCoreRadius + 5 * stdCoreRadius;
      gsl_integration_qag(&F, lowerLimit, upperLimit, 1e-8, 1e-8, 1000, GSL_INTEG_GAUSS15, workspace, &result, &error);
    } else {
      gsl_integration_qagiu(&F, 0.0, 1e-8, 1e-8, 1000, workspace, &result, &error);
    }
    gsl_integration_workspace_free(workspace);
    return result / meanCoreRadius3;
}
