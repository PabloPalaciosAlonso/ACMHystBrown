from setuptools import setup, Extension, find_packages
import os
import pybind11

# Directorio base
base_dir = os.path.abspath(os.path.dirname(__file__))

# Directorios de inclusión
include_dirs = [
    pybind11.get_include(),
    os.path.join(base_dir, "src"),
    os.path.join(base_dir, "wrappers"),
    os.path.join(base_dir, "external", "FittingAlgorithms", "src"),
    "/usr/include",  # Asegura que se incluya la cabecera de GSL
    "/usr/include/gsl",  # Algunos sistemas la tienen aquí
]

# Definición de la extensión
extension = Extension(
    name="ACMHystBrown.ACMHystBrown",
    sources=["wrappers/ACMHystBrownWrapper.cpp"],
    include_dirs=include_dirs,
    language="c++",
    extra_compile_args=["-std=c++17", "-fopenmp", "-O3"],
    extra_link_args=["-fopenmp", "-lgsl", "-lgslcblas", "-lm"],
)

# Configuración de setup
setup(
    name="ACMHystBrown",
    version="1.0.0",
    packages=find_packages(include=["ACMHystBrown", "ACMHystBrown.*"]),
    package_dir={"ACMHystBrown": "ACMHystBrown"},
    ext_modules=[extension],
    zip_safe=False,
    python_requires=">=3.6",
    install_requires=["pybind11>=2.6.0", "numpy"],
)
