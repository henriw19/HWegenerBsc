# Import the setuptools module
import setuptools

# Configure the package with the setup() function from setuptools
setuptools.setup(

    # Set the name of the package
    name="hwbsc",

    # Set the version number of the package
    version="0.0.0",

    # Specify the package to be included in the distribution
    packages=["hwbsc"],

    # Specify the package dependencies required for the package to work properly
    install_requires=[
        "numpy",
        "matplotlib",
        "jupyter",
    ],

    # Specify the minimum Python version required to run the package
    python_requires=">=3.10",

    # Specify the directory where the source files are located
    package_dir={"": "src"},
)