# freeDegreeEst

This is a python wrapper, based on Boost Python, for code contained in 
Willet et al. 2007, "Multiscale Poisson Intensity and Density Estimation". 
The original Matlab code can be downloaded from:

http://willett.ece.wisc.edu/code/FreeDegree.zip 

*** The original code is property of the respective authors. I do not own
it and I didn't have any part in the research which produced that code. I simply
wrote this wrapper because I found that code useful for my research.

## Installation

You need Boost Python installed. If you have a custom installation of Boost
Python, you can specify its root by using the BOOSTROOT environment variable.

Then you can install the code as:

> python setup.py install

## Usage

See the example.py script under Examples. Note that you need numpy, matplotlib
and scipy installed to run the example.
