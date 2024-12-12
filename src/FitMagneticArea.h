#include <iostream>
#include <map>
#include "computeMagneticAreaAnalytical.h"
#include "FittingAlgorithms.h"

using namespace FittingAlgorithms;

// Helper function to retrieve a parameter with extra and validation
double getParameter(const std::map<std::string, double>& fittingParams,
                    const std::map<std::string, double>& extraParams,
                    const std::string& key,
                    std::optional<double> defaultValue = std::nullopt) {
    bool inMainParams = fittingParams.count(key) > 0;
    bool inExtraParams = extraParams.count(key) > 0;

    if (inMainParams && inExtraParams) {
        throw std::runtime_error("Parameter \"" + key + "\" is present in both fittingParams and extraParams.");
    }

    if (!inMainParams && !inExtraParams) {
        // If key is not found in either map, use the default value if provided
        if (defaultValue.has_value()) {
            return defaultValue.value();
        }
        throw std::runtime_error("Parameter \"" + key + "\" is missing in both fittingParams and extraParams, and no default value is provided.");
    }

    return inMainParams ? fittingParams.at(key) : extraParams.at(key);
}

double computeAreaFit(FieldParameters fb,
                      std::map<std::string, double> fittingParameters,
                      std::map<std::string, double> extraParameters){ 
  
  // Retrieve parameters from fittingParameters or extra to extraParameters if not found
  double coatingWidth  = getParameter(fittingParameters, extraParameters, "coatingWidth");
  double coreRadius    = getParameter(fittingParameters, extraParameters, "coreRadius");
  double msat          = getParameter(fittingParameters, extraParameters, "msat");
  double numParticles  = getParameter(fittingParameters, extraParameters, "numParticles");
  double viscosity     = getParameter(fittingParameters, extraParameters, "viscosity");
  double kBT           = getParameter(fittingParameters, extraParameters, "kBT");
  double stdCoreRadius = getParameter(fittingParameters, extraParameters, "stdCoreRadius", 0.0);

  if (stdCoreRadius <= 0)
    return numParticles * computeArea(fb, coreRadius, msat, coatingWidth, viscosity, kBT);
  else
    return numParticles * computeArea(fb, coreRadius, msat, coatingWidth, viscosity, kBT, stdCoreRadius);
}


GaussNewton::FitResult fitAreas(std::vector<FieldParameters> &fb,
                                std::vector<double> &areas,
                                ParallelTempering::Parameters ptParams,
                                GaussNewton::Parameters gnParams,
                                std::vector<StringDoubleMap> &initialGuesses,
                                StringDoubleMap extraParameters){
  
  StringDoubleMap fit1 = ParallelTempering::fit<FieldParameters>(fb, areas, computeAreaFit,
                                                                 initialGuesses, ptParams,
                                                                 squaredRelativeError,
                                                                 extraParameters);
  
  GaussNewton::FitResult fitDef = GaussNewton::fit<FieldParameters>(fb, areas, computeAreaFit,
                                                                    fit1, gnParams,
                                                                    squaredRelativeError,
                                                                    extraParameters);
  
  return fitDef;
}

GaussNewton::FitResult fitAreas(std::vector<FieldParameters> &fb,
                                std::vector<double> &areas,
                                ParallelTempering::Parameters ptParams,
                                GaussNewton::Parameters gnParams,
                                StringDoubleMap &initialGuess,
                                StringDoubleMap &extraParameters){

  std::vector<StringDoubleMap> initialGuesses(ptParams.temperatures.size(), initialGuess);
  return fitAreas(fb, areas, ptParams, gnParams, initialGuesses, extraParameters); 
}
