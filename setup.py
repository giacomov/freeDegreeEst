#!/usr/bin/env python
 
from distutils.core import setup
from distutils.extension import Extension
from distutils.command.install_headers import install_headers
import os

boost_root = os.environ.get("BOOSTROOT")

if(boost_root):
  #The user want to override pre-defined location of boost
  
  print("\n\n **** Using boost.python from the env. variable $BOOSTROOT (%s)" %(boost_root))
  
  include_dirs = [ os.path.join(boost_root,'include')]
  library_dirs = [ os.path.join(boost_root,'lib') ]
  
  print("     Include dir: %s" %(include_dirs))
  print("     Library dir: %s" %(library_dirs))

else:

  include_dirs = []
  library_dirs = []

pass
 
setup(
    
    name="freeDegreeEst",
    
    packages = [],
    
    version = 'v1.0.0',
    
    description = "Multiscale Poisson Intensity and Density Estimation",
    
    author = 'Giacomo Vianello, a wrapper around code from Willet et al. 2007',
    
    author_email = 'giacomo.vianello@gmail.com',
        
    ext_modules=[
        
        Extension("freeDegreeEst", 
                  
                  ["freeDegreeEst/freeDegreeEst.cpp"],
        
        
        
        libraries = ["boost_python"],
        
        include_dirs=include_dirs,
        
        library_dirs=library_dirs),
    ],
    
    
    headers=[],
    
    install_requires=[])

