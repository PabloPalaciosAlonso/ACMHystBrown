#include "ACMHystBrownWrapper.h"

// Registrar FieldParameters
void registerFieldParameters(py::module_& module) {
    py::class_<FieldParameters>(module, "FieldParameters")
        .def(py::init<>())
        .def(py::init<double, double>(), py::arg("amplitude"), py::arg("frequency"))
        .def_readwrite("amplitude", &FieldParameters::amplitude)
        .def_readwrite("frequency", &FieldParameters::frequency);
}

// Registrar parámetros de ajuste
void registerFittingParameters(py::module_& module) {
    using PTParams = FittingAlgorithms::ParallelTempering::Parameters;
    using GNParams = FittingAlgorithms::GaussNewton::Parameters;

    py::class_<PTParams>(module, "PTParameters")
        .def(py::init<>())
        .def_readwrite("maxIterations",  &PTParams::maxIterations)
        .def_readwrite("numStepsSwap",   &PTParams::numStepsSwap)
        .def_readwrite("printSteps",     &PTParams::printSteps)
        .def_readwrite("numStepsFinish", &PTParams::numStepsFinish)
        .def_readwrite("tolerance",      &PTParams::tolerance)
        .def_readwrite("temperatures",   &PTParams::temperatures)
        .def_readwrite("jumpSize",       &PTParams::jumpSize);

    py::class_<GNParams>(module, "GNParameters")
        .def(py::init<>())
        .def_readwrite("maxIterations",  &GNParams::maxIterations)
        .def_readwrite("tolerance",      &GNParams::tolerance)
        .def_readwrite("regularization", &GNParams::regularization)
        .def_readwrite("printSteps",     &GNParams::printSteps);
}

void registerFitAreas(py::module_& module) {
    module.def("fitAreas", [](std::vector<FieldParameters>& fb,
                              std::vector<double>& areas,
                              FittingAlgorithms::ParallelTempering::Parameters& ptParams,
                              FittingAlgorithms::GaussNewton::Parameters& gnParams,
                              std::vector<StringDoubleMap>& initialGuesses,
                              StringDoubleMap& extraParameters){
      auto fit = fitAreas(fb, areas, ptParams, gnParams, initialGuesses, extraParameters);
      return py::make_tuple(fit.parameters, fit.errors);
    },
               py::arg("fb"),
               py::arg("areas"),
               py::arg("ptParams"),
               py::arg("gnParams"),
               py::arg("initialGuesses"),
               py::arg("extraParameters"));
    
    module.def("fitAreas", [](std::vector<FieldParameters>& fb,
                              std::vector<double>& areas,
                              FittingAlgorithms::ParallelTempering::Parameters& ptParams,
                              FittingAlgorithms::GaussNewton::Parameters& gnParams,
                              StringDoubleMap& initialGuesses,
                              StringDoubleMap& extraParameters){
      auto fit = fitAreas(fb, areas, ptParams, gnParams, initialGuesses, extraParameters);
      return py::make_tuple(fit.parameters, fit.errors);
    },
               py::arg("fb"),
               py::arg("areas"),
               py::arg("ptParams"),
               py::arg("gnParams"),
               py::arg("initialGuesses"),
               py::arg("extraParameters"));

}

void registerComputeArea(py::module_& module) {
    module.def("computeArea", [](FieldParameters fb,
                                 double meanCoreRadius,
                                 double msat,
                                 double coatingWidth,
                                 double viscosity,
                                 double kBT,
                                 double stdCoreRadius) {
        return computeArea(fb, meanCoreRadius, msat, coatingWidth, viscosity, kBT, stdCoreRadius);
    },
               py::arg("fb"),
               py::arg("meanCoreRadius"),
               py::arg("msat"),
               py::arg("coatingWidth"),
               py::arg("viscosity"),
               py::arg("kBT"),
               py::arg("stdCoreRadius"),
               "Computa el área bajo la curva de magnetización.");
}


// Módulo principal
PYBIND11_MODULE(ACMHystBrown, m) {
  registerFieldParameters(m);
  registerFittingParameters(m);
  registerFitAreas(m);
  registerComputeArea(m);
 }
