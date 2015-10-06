# freeDegreeEst

This is a python wrapper, based on Boost Python, for code contained in 
Willet et al. 2007, "Multiscale Poisson Intensity and Density Estimation". 
The original Matlab code can be downloaded from:

http://willett.ece.wisc.edu/code/FreeDegree.zip 

## Installation

You need Boost Python installed. If you have a custom installation of Boost
Python, you can specify its root by using the BOOSTROOT environment variable.

Then you can install the code as:

> python setup.py install

## Usage

See the example.py script under Examples. Note that you need numpy, matplotlib
and scipy installed to run the example.
