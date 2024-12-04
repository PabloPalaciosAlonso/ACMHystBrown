#include<iostream>
#include<cmath>
#include <algorithm>
#include <numeric>
#include"../../ParallelTemperingCode/src/parallelTempering.h"
#include"../../Gauss-Newton-Code/src/GaussNewton.h"
#include"../src/laguerreQuadrature.h"
#include"../src/usefulFunctions.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <nlohmann/json.hpp> // Incluye la biblioteca JSON

double computeMaximumArea(double badim){

  double b = 0.0324349;
  double c = 1.88238;
  double d = 1.17e-5;
  double e = 3.683;

  // Definición de las constantes
  
  // Cálculo de los términos intermedios
  double term1 = badim * M_PI / 24.0;
  double term2 = b * pow(badim, c);
  double term3 = d * pow(badim, e);
  
  // Cálculo del denominador dentro del logaritmo
  double log_argument = 1.0 + term1 + term2 + term3;
  
  // Verificación para evitar el logaritmo de un número no positivo
  if (log_argument <= 0.0) {
    fprintf(stderr, "Error: El argumento del logaritmo debe ser positivo.\n");
    return -1.0; // Valor de error
  }
    
  // Cálculo de fmax según la fórmula proporcionada
  return 4.0 * (1.0 - 1.0 / (1.0 + log(log_argument)));
}

//This function returns the parameters that ensure the position of the maximum and its value is correct
std::tuple<double, double, double, double> calculate_parameters(double badim, double p) {

  double amax     = computeMaximumArea(badim);
  double slopeInf = M_PI*badim/3;
  double wmax     = 1 + 0.289*pow(badim, 0.89)/(1+11.83*pow(badim, -1.633));
  
  // Verificar que los valores ingresados sean válidos
  if (p <= 0 || amax <= 0 || wmax <= 0 || slopeInf <= 0) {
    throw std::invalid_argument("Todos los parámetros deben ser positivos y mayores que cero.");
  }

  // Paso 1: Calcular beta
  double numerator_beta = (p + 1) * (std::log(p) - std::log(1.0 + p));
  double denominator_beta = -std::log(slopeInf / (wmax * amax));

  // Manejo de división por cero o logaritmos de números no positivos
  if (denominator_beta == 0.0) {
    throw std::runtime_error("La expresión en el denominador para beta es cero. Ajusta los parámetros para evitar división por cero.");
  }

  if (denominator_beta < 0.0) {
    // Si el denominador es negativo, invertimos el signo del numerador y del denominador
    numerator_beta = -numerator_beta;
    denominator_beta = -denominator_beta;
  }

  double beta = numerator_beta / denominator_beta;

  // Paso 2: Calcular q
  double q = - (p + 1) / beta;

  // Paso 3: Calcular u
  double u = (std::log(amax) - q * std::log(1.0 + p)) / p;

  // Paso 4: Calcular x0
  double x0 = wmax * std::exp(-u);

  // Paso 5: Calcular alpha
  double ln_alpha = std::log(p) - beta * u;
  double alpha = std::exp(ln_alpha);

  // Retornar los parámetros calculados
  return std::make_tuple(beta, q, alpha, x0);
}

// compute_max_area function implementation
double compute_max_area(double badim, double c1, double c2, double c3, double c4, double c5, double c6, double c7, double c8) {
    double weight = exp(-c1 / pow(badim, c2));
    double logarg = log(1 + c3 * pow(badim, c4) + c5 * pow(badim, c6) + c7 * pow(badim, c8));
    double maxA = M_PI * badim / 6 + weight * (4 - M_PI * badim / 6 - 4 / (1 + logarg));
    return maxA;
}

double computeArea(double wadim, double badim){
  double c1 = 8.5;
  double c2 = 1.75;
  double c3 = 0.55;
  double c4 = 0.16;
  double c5 = 2.05;
  
  double c7 = 0.5;
  double c10 = 1.8;
  double c11 = 1660;
  double c13 = 11;
  double c18 = 0.18;
  double c19 = 1;
  double c21 = 0.8915;
  double c31 = 0.5;
  double c32 = 1.8;
  double c33 = 0.1;
  double c34 = 0.9;
  double c35 = 0.056;
  double c38 = 3.21;
  
  double c12= 6.33;
  double c20= 0.287;
  double c22= 5.52;
  double c23= 1.16;
  double c36= 1.65;
  double c37= 1.54e-05;
  
  double p = (c3 - 1) * exp(-c1 / pow(badim, c2)) + 1;
  
  double amax = compute_max_area(badim, c31, c32, c33, c34, c35, c36, c37, c38);
  
  double wmax = 1 + c20 * pow(badim, c21) / (1 + c22 * pow(badim, -c23));
  double slope_inf = M_PI * badim / (3 * wmax);
  double slope_0 = M_PI * badim / (3 * pow(1 + c4 * pow(badim, c5), 1.0 / c5));
  
  double q = log(slope_inf / amax) / log((p + 1) / p);
  double r = (p + 1) / q;
  double b = p * pow(amax, -1 / q) / (q * r);
  double a = b / p;
  
  double l1 = c10;
  double l2 = (1 - p) / l1;
  
  double kk = pow(1 / (wmax * slope_0 * pow(a, q)), 1 / l2);
  
  double a1 = (wadim / wmax) / pow(kk + (1 + c11 * exp(-c12 * pow(wadim, c18)) * exp(-c13 / pow(badim, c19)) * pow(badim, -c7)) * pow(wadim / wmax, l1), l2);
  
  return a1 / pow((a + b * pow(wadim / wmax, r)), q);
}

bool compareMaps(const Gauss_Newton::FitResult& fits) {

  std::map<std::string, double> fittingErrors = fits.errors;
  std::map<std::string, double> fittingParams = fits.parameters;
  
  
  for (const auto& errorEntry : fittingErrors) {
    const std::string& key = errorEntry.first;
    double errorValue = errorEntry.second;
    
    // Buscar la clave en el segundo mapa
    auto paramIt = fittingParams.find(key);
    if (paramIt != fittingParams.end()) {
      double paramValue = paramIt->second;
      
      // Comparar valores
      if (errorValue > paramValue) {
        return false;
      }
    }
  }
  return true;
}

void GSL_Error_Handler(const char *reason, const char *file, int line, int gsl_errno) {
  std::cerr << "GSL Error: " << reason << " in " << file << ":" << line << ", error code: " << gsl_errno << std::endl;
}

struct FieldParameters{
  double field;
  double freq;
};

// Wrapper for the integrand function to use with GSL
double integrandr3(double r, void *params) {
    double *p = (double *)params;
    double rmean = p[0];
    double rstd = p[1];
    return pow(r, 3) * lognormal(r, rmean, rstd);
}

// Function to compute the average cubic radius using GSL for integration
double computeAverageCubicRadius(double rmean, double rstd) {

  gsl_integration_workspace *w = gsl_integration_workspace_alloc(5000);
  
  double result, error;
  double params[2] = {rmean, rstd};
  gsl_function F;
  F.function = &integrandr3;
  F.params = &params;
  
  // Integrate from 0 to infinity
  int status = gsl_integration_qagiu(&F, 0.0, 0, 1e-9, 1000, w, &result, &error);
  if (status != 0) {
      throw std::runtime_error("Error en la integración con GSL.");
  }
  
  gsl_integration_workspace_free(w);
  
  if (result == 0.0) result = rmean*rmean*rmean;
  return result;
}

// Wrapper for the integrand function to use with GSL
double integrandr6(double r, void *params) {
    double *p = (double *)params;
    double rmean = p[0];
    double rstd = p[1];
    return pow(r, 6) * lognormal(r, rmean, rstd);
}

// Function to compute the average cubic radius using GSL for integration
double computeAverageRadius6(double rmean, double rstd) {
    // GSL integration setup
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(5000);
    double result, error;
    double params[2] = {rmean, rstd};
    gsl_function F;
    F.function = &integrandr6;
    F.params = &params;

    // Integrate from 0 to infinity
    int status = gsl_integration_qagiu(&F, 0.0, 0, 1e-9, 1000, w, &result, &error);
    if (status != 0) {
      throw std::runtime_error("Error en la integración con GSL.");
    }
    gsl_integration_workspace_free(w);
    return result;
}


double computeAreaOld(double wadim, double badim) {

  
  double c1= 0.73;
  double c10= 1.92418;
  double c11= 0.127154;
  double c12= 2.148;
  double c13= 0.453308;
  double c14= 0.923318;
  double c15= 0.163892;
  double c16= 2.05725;
  double c17= 0.88;
  double c18= 0.93;
  double c2= 0.625683;
  double c3= 0.916712;
  double c4= 0.00406689;
  double c5= 24.2396;
  double c6= 0.272884;
  double c7= 1.69641;
  double c8= 1.59651;
  double c9= 3.19173;


  
  double t1 = wadim*wadim;
  double t2 = pow(1+c15*pow(badim,c16), 1./c16);
  double t3 = c11 * badim * badim * wadim / (c12 + pow(badim*wadim, -c13 + 1))*exp(-c14/(pow(badim, c17)*pow(wadim,c18)));
  
  double locCorrMax = c1 + c2*pow(badim, c3)/(1+pow(badim, c4));
  double hmax = c5;
  double amax = c6;
  double b = (c7-c8)/c8*pow(hmax,c7);
  double a = amax*c7/c8*pow(hmax,c7-c8);
  double ampCorrMax = a * pow(badim,c8)/ (b + pow(badim,c7));

  double b2 = (c9-c10)/c10*pow(locCorrMax, c9);
  double a2 = ampCorrMax*c9/c10*pow(locCorrMax, c9-c10);

  double corr = a2 * pow(wadim,c10)/ (b2 + pow(wadim,c9));

  double area1 = badim*M_PI*wadim/(3*(t1 + t2 + t3));
  
  double area = area1/(1-corr);

  // if (area>4)
  //   std::cout<<wadim<<" "<<badim<<" "<<area<<std::endl;
  
  return area;
}

// Wrapper for the integrand function to use with GSL
double integrandArea(double r, void *params) {
  double *p = (double *)params;
  double coreRadius         = p[0];
  double std_coreRadius     = p[1];
  double wadim_nocoreRadius = p[2];
  double badim_nocoreRadius = p[3];
  double coatingWidth       = p[4];
  //std::cout<<wadim_nocoreRadius * std::pow(coatingWidth + r, 3)<<std::endl;
  //std::cout<<"    "<<r<<std::endl;
  return pow(r, 3) *
    lognormal(r, coreRadius, std_coreRadius) *
    computeArea(wadim_nocoreRadius * std::pow(coatingWidth + r, 3), badim_nocoreRadius * std::pow(r, 3));
}

// Function to compute the area
double computeArea(double coreRadius, double std_coreRadius,
                   double wadim_nocoreRadius, double badim_nocoreRadius,
                   double coatingWidth, double coreRadius3_mean,
                   gsl_integration_workspace *w) {
   
    double result, error;
    double params[5] = {coreRadius, std_coreRadius, wadim_nocoreRadius, badim_nocoreRadius, coatingWidth};
    gsl_function F;
    F.function = &integrandArea;
    F.params = &params;

    // Integrate from 0 to infinity
    //int status = gsl_integration_qags(&F, 0.0, 500.0, 0, 1e-9, 1000, w, &result, &error);
    int status = gsl_integration_qagiu(&F, 0.0, 0, 1e-6, 1000, w, &result, &error);
    if (status != 0) {
      throw std::runtime_error("Error en la integración con GSL.");
      gsl_integration_workspace_free(w);
    }
    return result / coreRadius3_mean;
}

double computeAreasDim(FieldParameters fb,
                       std::map<std::string, double> &fittingParameters,
                       std::map<std::string, double> &extraParameters,
                       gsl_integration_workspace *w){
  // Extract required parameters from extraParameters
  double viscosity   = extraParameters["viscosity"];
  double temperature = extraParameters["temperature"];

  // Retrieve parameters from fittingParameters or fallback to extraParameters if not found
  double coatingWidth     = fittingParameters.count("coatingWidth") ? fittingParameters["coatingWidth"] : extraParameters["coatingWidth"];
  double std_coatingWidth = fittingParameters.count("std_coatingWidth") ? fittingParameters["std_coatingWidth"] : extraParameters["std_coatingWidth"];
  double coreRadius       = fittingParameters.count("coreRadius") ? fittingParameters["coreRadius"] : extraParameters["coreRadius"];
  double std_coreRadius   = fittingParameters.count("std_coreRadius") ? fittingParameters["std_coreRadius"] : extraParameters["std_coreRadius"];
  double msat             = fittingParameters.count("msat") ? fittingParameters["msat"] : extraParameters["msat"];
  double normalization    = fittingParameters.count("normalization") ? fittingParameters["normalization"] : extraParameters["normalization"];
  
  // Check if any parameter is not found in both fittingParameters and extraParameters
  if (std::isnan(coatingWidth) || std::isnan(std_coatingWidth) || std::isnan(coreRadius) || std::isnan(std_coreRadius) || std::isnan(msat)) {
    std::cerr << "Parameter not found in both fittingParameters and extraParameters" << std::endl;
    throw std::runtime_error("Parameter not found in both fittingParameters and extraParameters");
  }

  if (std_coreRadius >0.0){
    std_coreRadius = std::min(std_coreRadius, 0.9);
    std_coreRadius = std::max(std_coreRadius, 0.01);
  }

  std_coreRadius*=coreRadius;

  double f = fb.freq;
  double b = fb.field;
  double area;

  static double last_coreRadius     = std::numeric_limits<double>::quiet_NaN();
  static double last_std_coreRadius = std::numeric_limits<double>::quiet_NaN();
  static double coreRadius3_mean    = 0.0;
    
  if (std_coatingWidth == 0.0 && std_coreRadius == 0.0) {
    double coreRadius3 = coreRadius*coreRadius*coreRadius;
    double hydroRadius3 = (coreRadius + coatingWidth)*(coreRadius + coatingWidth)*(coreRadius + coatingWidth);
    double wadim = 8 * M_PI * f * M_PI * viscosity * hydroRadius3 / temperature;
    double badim = 4 * M_PI * b * msat * coreRadius3 / (3 * temperature);
    area  = computeArea(wadim, badim);
  } else if (std_coatingWidth == 0 && std_coreRadius != 0.0) {
    std_coreRadius = std::max(std_coreRadius, 1.0);
    double wadim_nocoreRadius = 8 * std::pow(M_PI, 2) * f * viscosity / temperature;
    double badim_nocoreRadius = 4 * M_PI * b * msat / (3 * temperature);
    if (std::isnan(last_coreRadius) || coreRadius != last_coreRadius || std_coreRadius != last_std_coreRadius) {
      coreRadius3_mean    = computeAverageCubicRadius(coreRadius, std_coreRadius);
      last_coreRadius     = coreRadius;
      last_std_coreRadius = std_coreRadius;
    }

    area = computeArea(coreRadius, std_coreRadius, wadim_nocoreRadius,
                       badim_nocoreRadius, coatingWidth, coreRadius3_mean, w);
  }
  return area*normalization;
}

std::vector<double> computeAreasDimVec(std::vector<FieldParameters> &fb,
                                       std::map<std::string, double> &fittingParameters,
                                       std::map<std::string, double> &extraParameters){
  std::vector<double> areas(fb.size());
  gsl_integration_workspace *w = gsl_integration_workspace_alloc(100000);

  //#pragma omp parallel for
  fori(0, fb.size()){
    areas[i] = computeAreasDim(fb[i], fittingParameters, extraParameters, w);
  }
  
  gsl_integration_workspace_free(w);
  return areas;
}

std::vector<double> computeAreasDimVecLog(std::vector<FieldParameters> &fb,
                                          std::map<std::string, double> &fittingParameters,
                                          std::map<std::string, double> &extraParameters){
  std::vector<double> areasLog(fb.size());
  gsl_integration_workspace *w = gsl_integration_workspace_alloc(100000);

  //#pragma omp parallel for
  fori(0, fb.size()){
    areasLog[i] = log(computeAreasDim(fb[i], fittingParameters, extraParameters, w));
  }
  
  gsl_integration_workspace_free(w);
  return areasLog;
}

std::vector<double> computeNormalizedAreasVec(std::vector<FieldParameters> &fb,
                                              std::vector<double> &ydata_exp,
                                              std::map<std::string, double> &fittingParameters,
                                              std::map<std::string, double> &extraParameters){
  std::vector<double> normalizedAreas(fb.size());
  gsl_integration_workspace *w = gsl_integration_workspace_alloc(5000);

  //#pragma omp parallel for
  fori(0, fb.size()){
    double ypred = computeAreasDim(fb[i], fittingParameters, extraParameters, w);
    normalizedAreas[i] = sqrt(fabs(ypred-ydata_exp[i])/ydata_exp[i]);
  }
  
  gsl_integration_workspace_free(w);
  return normalizedAreas;
}

double costFunctionMax(std::vector<FieldParameters> &xdata,
                       std::vector<double> &ydata,
                       std::map<std::string, double> &fittingParameters,
                       std::map<std::string, double> &extraParameters) {

  std::vector<double> ypred     = computeAreasDimVec(xdata, fittingParameters, extraParameters);
  if (ypred.size() != ydata.size()) {
    throw std::invalid_argument("Los vectores deben tener el mismo tamaño");
  }
  
  double max_diff = 0.0;
  
  // Calcular (ypred - ydata) / ydata y encontrar el máximo valor absoluto en el mismo bucle
  for (size_t i = 0; i < ypred.size(); ++i) {
    double diff = fabs(ypred[i] - ydata[i]);
    if (diff > max_diff) {
      max_diff = diff;
    }
  }
  
  return max_diff;  // Devolver el máximo
}

double costFunctionNormalizedInverse(std::vector<FieldParameters> &xdata,
                                     std::vector<double> &ydata,
                                     std::map<std::string, double> &fittingParameters,
                                     std::map<std::string, double> &extraParameters) {

  // Get the predictions using the provided function
  std::vector<double> ypred     = computeAreasDimVec(xdata, fittingParameters, extraParameters);
  std::vector<double> ydataCopy = ydata;
  
  // Check if sizes of ypred and ydata are equal
    if (ypred.size() != ydata.size()) {
        throw std::runtime_error("Size mismatch between predictions and actual data");
    }

    // // Find the max values in ypred and ydata
    // double max_ypred = *std::max_element(ypred.begin(), ypred.end());
    // double max_ydata = *std::max_element(ydata.begin(), ydata.end());

    // // Normalize ypred and ydata
    // for (auto& val : ypred) {
    //     val /= max_ypred;
    // }
    // for (auto& val : ydataCopy) {
    //     val /= max_ydata;
    // }

    // Compute the mean squared error
    double mse = std::inner_product(ypred.begin(), ypred.end(), ydataCopy.begin(), 0.0,
                                    std::plus<double>(),
                                    [](double a, double b) { return (1/(b))*(a - b)*(a - b); }) / ydata.size();

    double coreRadius = fittingParameters.count("coreRadius") ? fittingParameters["coreRadius"] : extraParameters["coreRadius"];
    double std_coreRadius = fittingParameters.count("std_coreRadius") ? fittingParameters["std_coreRadius"] : extraParameters["std_coreRadius"];
    double msat = fittingParameters.count("msat") ? fittingParameters["msat"] : extraParameters["msat"];
    // mse += 50*(coreRadius-13.63)*(coreRadius-13.63)+100000000*(msat-0.000166)*(msat-0.000166)
    //   +5*(std_coreRadius-2.36)*(std_coreRadius-2.36);
    return mse;
}


double costFunctionNormalizedLog(std::vector<FieldParameters> &xdata,
                                 std::vector<double> &ydata,
                                 std::map<std::string, double> &fittingParameters,
                                 std::map<std::string, double> &extraParameters) {
  
  // Get the predictions using the provided function


  std::vector<double> ypred(xdata.size());//     = computeAreasDimVec(xdata, fittingParameters, extraParameters);
  try {
    ypred = computeAreasDimVec(xdata, fittingParameters, extraParameters);
  } catch (const std::exception& e) {
    ypred = std::vector<double>(xdata.size(), 4.0);
  }

  std::vector<double> ydataCopy = ydata;
  
  // Check if sizes of ypred and ydata are equal
  if (ypred.size() != ydata.size()) {
    throw std::runtime_error("Size mismatch between predictions and actual data");
  }
  
  
  // Compute the mean squared error
  double mse = std::inner_product(ypred.begin(), ypred.end(), ydataCopy.begin(), 0.0,
                                  std::plus<double>(),
                                  [](double a, double b) {return log(a/b)*log(a/b); }) / ydata.size();

  return mse;
}

double costFunctionNormalized(std::vector<FieldParameters> &xdata,
                              std::vector<double> &ydata,
                              std::map<std::string, double> &fittingParameters,
                              std::map<std::string, double> &extraParameters) {
  
  // Get the predictions using the provided function


  std::vector<double> ypred(xdata.size());//     = computeAreasDimVec(xdata, fittingParameters, extraParameters);
  try {
    ypred = computeAreasDimVec(xdata, fittingParameters, extraParameters);
  } catch (const std::exception& e) {
    ypred = std::vector<double>(xdata.size(), 4.0);
  }

  std::vector<double> ydataCopy = ydata;
  
  // Check if sizes of ypred and ydata are equal
  if (ypred.size() != ydata.size()) {
    throw std::runtime_error("Size mismatch between predictions and actual data");
  }
  
  
  // Compute the mean squared error
  double mse = std::inner_product(ypred.begin(), ypred.end(), ydataCopy.begin(), 0.0,
                                  std::plus<double>(),
                                  [](double a, double b) {return fabs(a-b)/b; }) / ydata.size();

  return mse;
}

std::vector<double> costFunctionNormalizedVec(std::vector<FieldParameters> &xdata,
                                              std::vector<double> &ydata,
                                              std::map<std::string, double> &fittingParameters,
                                              std::map<std::string, double> &extraParameters) {
  
  // Get the predictions using the provided function

  int nPoints = xdata.size();
  std::vector<double> ypred(nPoints);
  ypred = computeAreasDimVec(xdata, fittingParameters, extraParameters);
  
  // Check if sizes of ypred and ydata are equal
  if (ypred.size() != ydata.size()) {
    throw std::runtime_error("Size mismatch between predictions and actual data");
  }

  std::vector<double> mse(nPoints);

  for (int i = 0; i<nPoints; i++)
    mse[i] = sqrt(fabs(ypred[i]-ydata[i])/ydata[i]);
  return mse;
}


void readAreas(const std::string& filename, std::vector<FieldParameters>& fb, std::vector<double>& areas) {

  double mu0       = 4*M_PI*1e-7*1e9;
  std::ifstream infile(filename);
  if (!infile) {
    throw std::runtime_error("Could not open file " + filename);
  }
  
  std::string line;
  while (std::getline(infile, line)) {
    // Ignore comment lines starting with #
    if (line.empty() || line[0] == '#') {
      continue;
    }
    
    std::istringstream iss(line);
    double freq, field, area;
    if (!(iss >> freq >> field >> area)) {
      throw std::runtime_error("Error reading line: " + line);
    }
    
    FieldParameters fp; fp.freq = freq*1e-6; fp.field = field*1e-6*mu0;
    fb.push_back(fp);
    areas.push_back(area/(field));
  }
  
  infile.close();
}

std::string replaceExtension(const std::string& filename) {
    size_t lastDot = filename.find_last_of(".");
    if (lastDot == std::string::npos) {
        // No extension found, just append ".json"
        return filename + "Fit.json";
    } else {
        // Replace the extension with ".json"
        return filename.substr(0, lastDot) + "Fit.json";
    }
}

void writeResults(const std::string& in_filename,
                  const std::map<std::string, double>& knownParameters,
                  const Gauss_Newton::FitResult& fits) {
    // Crear un objeto JSON
    nlohmann::json jsonOutput;

    // Agregar parámetros conocidos
    for (const auto& [key, value] : knownParameters) {
        jsonOutput["knownParameters"][key] = value;
    }

    // Agregar parámetros ajustados
    for (const auto& [key, value] : fits.parameters) {
        jsonOutput["fittingParameters"][key] = value;
    }

    // Agregar errores de ajuste
    for (const auto& [key, value] : fits.errors) {
      jsonOutput["fittingErrors"][key] = value;
    }

    // Determinar el nombre del archivo de salida
    std::string filename = replaceExtension(in_filename);

    // Escribir el archivo JSON
    std::ofstream file(filename);
    if (file.is_open()) {
        file << jsonOutput.dump(4); // Formato legible con una sangría de 4 espacios
        file.close();
    } else {
        throw std::ios_base::failure("No se pudo abrir el archivo: " + filename);
    }
}


void writeResults(const std::string& in_filename,
                  const std::map<std::string, double>& knownParameters,
                  const std::map<std::string, double>& fittingParameters) {
    // Crear un objeto JSON
    nlohmann::json jsonOutput;

    // Agregar parámetros conocidos
    for (const auto& [key, value] : knownParameters) {
        jsonOutput["knownParameters"][key] = value;
    }

    // Agregar parámetros ajustados
    for (const auto& [key, value] : fittingParameters) {
        jsonOutput["fittingParameters"][key] = value;
    }

    // Determinar el nombre del archivo de salida
    std::string filename = replaceExtension(in_filename);

    // Escribir el archivo JSON
    std::ofstream file(filename);
    if (file.is_open()) {
        file << jsonOutput.dump(4); // Formato legible con una sangría de 4 espacios
        file.close();
    } else {
        throw std::ios_base::failure("No se pudo abrir el archivo: " + filename);
    }
}



int main(int argc, char* argv[]){
  gsl_set_error_handler(&GSL_Error_Handler);

  double minRc = 10;
  double maxRc = 30;

  double minMsat = 10000*1e-9;
  double maxMsat = 1000000*1e-9;

  double minWidth = 1;
  double maxWidth = 50;

  double minRcStd = 0.05;
  double maxRcStd = 1;

  double viscosity = 0.000918;
  double kb        = 1.38e-23*1e18;
  double temp      = 300;
  
  int nsteps      = 3000;
  int nstepsSwap  = 200;
  int nstepsBreak = 2000;
  
  std::vector<double> temperatures = {0.00001,0.00005, 0.0001, 0.0005, 0.0025, 0.0125, 0.0625};
  std::vector<double> jumpSize = {0.01,0.01, 0.01, 0.01, 0.01, 0.01, 0.01};

  auto initialGuesses = generateInitialGuesses(minRc, maxRc, minMsat, maxMsat, minWidth, maxWidth,
                                               minRcStd, maxRcStd, temperatures.size());



  // std::map<std::string, double> initialGuess = {
  //   {"coreRadius", 3.},
  //   {"coatingWidth", 3.1925},
  //   {"msat", 2.1},
  //   {"std_coreRadius", 0.1625}
  // };

  
  std::map<std::string, double> knownParameters = {
    {"temperature", 0.0041124},
    {"viscosity", viscosity},
    {"std_coatingWidth", 0.0},
    //{"coatingWidth", 0.0},
    {"coreRadius", 15},
    {"msat", 0.00015},
    {"normalization", 1},
    {"std_coreRadius", 0.1}
  };
  
  std::string filename = std::string(argv[1]);
  
  std::vector<FieldParameters> fb;
  std::vector<double> areas;
  readAreas(filename, fb, areas);

  std::cout<<knownParameters["temperature"]<<" "<<knownParameters["viscosity"]<<std::endl;

  parallelTempering::MCParameters par;
  par.temperatures = temperatures;
  par.numberStepsFinish  = nstepsBreak;
  par.numberStepsSwap    = nstepsSwap;
  par.maximumNumberSteps = nsteps;
  par.jumpSize           = jumpSize;
  par.errorBreak         = 0.0075;

  
  auto fittingParameters = parallelTempering::fitData(fb, areas,
                                                      costFunctionNormalized,
                                                      par,
                                                      initialGuesses, knownParameters);
  
  // Gauss-Newton algorithm parameters
  Gauss_Newton::GNParameters gnParams;
  gnParams.maxIterations     = 500;
  gnParams.tolerance         = 1e-4;
  gnParams.regularization    = 1e-6;
  gnParams.printSteps        = 5;
  
  std::vector<double> areasLog(areas.size());
  
  fori(0, areasLog.size()) areasLog[i] = log(areas[i]);

  if (areas.size()<= fittingParameters.size()){
    writeResults(filename, knownParameters, fittingParameters);
    
  } else {


    
    try {
      auto results = Gauss_Newton::fitParams(fb, areas, gnParams,
                                             costFunctionNormalizedVec,
                                             fittingParameters,
                                             knownParameters);
      bool smallError = compareMaps(results);

      int count = 0;
      while (not smallError and count <5){
        std::cout<<"Not valid solution, starting again:\n";
        auto initialGuesses = generateInitialGuesses(minRc, maxRc, minMsat, maxMsat, minWidth, maxWidth,
                                               minRcStd, maxRcStd, temperatures.size());
        fittingParameters = parallelTempering::fitData(fb, areas,
                                                      costFunctionNormalized,
                                                      par,
                                                      initialGuesses, knownParameters);

        results = Gauss_Newton::fitParams(fb, areas, gnParams,
                                             costFunctionNormalizedVec,
                                             fittingParameters,
                                             knownParameters);
        smallError = compareMaps(results);
        count++;
        
      }
      writeResults(filename, knownParameters, results);
    } catch (...){
      writeResults(filename, knownParameters, fittingParameters);
    }
  }  
  return 0;
}
