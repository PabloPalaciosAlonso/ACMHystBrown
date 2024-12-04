#ifndef ACMHYSTAREABROWN_H
#define ACMHYSTAREABROWN_H

#include <pybind11/pybind11.h>
#include <pybind11/stl.h> // Para conversiones de std::vector y std::map
#include <pybind11/stl_bind.h> // Para soporte adicional de contenedores STL (opcional)
#include <vector>
#include <tuple>
#include "FitMagneticArea.h"
#include "GaussNewton.h"
#include "ParallelTempering.h"
#include "utils/defines.h"

namespace py = pybind11;
void registerFieldParameters(py::module_& module);
void registerFittingParameters(py::module_& module);
void registerFitAreas(py::module_& module);

#endif // ACMHYSTAREABROWN_H
