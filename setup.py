# setup.py
from setuptools import setup, find_packages

setup(
    name="PV_Circuit_Model",      # Name of your package
    version="0.1.0",               # Version number
    packages=find_packages(),      # Automatically find all packages in your source
    install_requires=[],           # List any dependencies here, like ["numpy", "requests"]
    author="Johnson Wong",
    description="PV Circuit model for cells, modules and components with additional tools for simulations and model fits to measurement data",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
    ],
)