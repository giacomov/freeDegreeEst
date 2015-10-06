import numpy
import matplotlib.pyplot as plt
from math import log, floor
import scipy.stats

import freeDegreeEst

#Generate simulated data

Nevents = 1000 #a thousand events

pdf = scipy.stats.beta(2,9)

events = sorted( pdf.rvs( size = Nevents ) )

#Now reconstruct the PDF by using the events

#Default value for the estimation length (see Willed et al. 2007)
estLength = int(pow(2.0,floor(log(Nevents)/log(2.0))))

#Default value for the penalization for the log-likelihood
penalization = log( Nevents ) / 2.0

#Maximum order of polynomial. Default value is equal to the number
#of points
maxDegree = Nevents

res = freeDegreeEst.getDensity( events, len(events), 
                                estLength, penalization,
                                maxDegree )

#Now plot the result, as well as the true distribution
x = numpy.linspace( 0, 1, estLength)
truePdf = pdf.pdf( x )

plt.plot( x, res, label='free-degree est.')
plt.plot( x, truePdf, '--', label='true distrib.')
plt.hist( events, 30, label='data', normed=True)

plt.xlabel("Independent variable")
plt.ylabel("Density")

plt.legend()

plt.show()

raw_input("Press ENTER when you are done")
