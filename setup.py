from setuptools import setup, Extension
import pybind11
import os

base_dir = os.path.abspath(os.path.dirname(__file__))

include_dirs = [
    pybind11.get_include(),
    os.path.join(base_dir, "src"),
    os.path.join(base_dir, "wrappers"),
    os.path.join(f"{base_dir}/external/FittingAlgorithms", "src"),
    
]

extensions = [
    Extension(
        name="ACMHystAreaBrown",
        sources=[
            "wrappers/ACMHystAreaBrownWrapper.cpp",
        ],
        include_dirs=include_dirs,
        language="c++",
        extra_compile_args=["-std=c++17", "-fopenmp", "-O3"],
        extra_link_args=["-fopenmp"],
    )
]

setup(
    name="ACMHystAreaBrown",
    version="1.0.0",
    author="Pablo Palacios-Alonso",
    author_email="pablo.palaciosa@uam.es",
    description="Python module with some fitting algorithms.",
    ext_modules=extensions,
    zip_safe=False,
    python_requires=">=3.6",
    install_requires=[
        "pybind11>=2.6.0",
    ],
)
