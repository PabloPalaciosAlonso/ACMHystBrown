from setuptools import setup, Extension, find_packages
import os
import pybind11

# Directorio base
base_dir = os.path.abspath(os.path.dirname(__file__))

# Directorios de inclusión
include_dirs = [
    pybind11.get_include(),
    os.path.join(base_dir, "ACMHystBrown", "src"),  # Ruta actualizada
    os.path.join(base_dir, "wrappers"),
    os.path.join(base_dir, "external", "FittingAlgorithms", "src"),
]

# Definición de la extensión
extension = Extension(
    name="ACMHystBrown.ACMHystBrown",  # Nombre del módulo dentro del paquete
    sources=["wrappers/ACMHystBrownWrapper.cpp"],  # Código fuente del módulo C++
    include_dirs=include_dirs,
    language="c++",
    extra_compile_args=["-std=c++17", "-fopenmp", "-O3"],  # Opciones de compilación
    extra_link_args=["-fopenmp"],  # Opciones de enlace
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
    install_requires=["pybind11>=2.6.0", "numpy"],  # Dependencias necesarias
)
