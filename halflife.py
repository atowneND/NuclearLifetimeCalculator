#!/usr/bin/python
# halflife.py

import scipy.signal
import numpy
import csv
import matplotlib.pyplot as plt

class calculateHalfLife:
    # __init__

    def __init__(self, **kwargs):
        self.__dict__.update(**kwargs)

    def read_CSV(self, conversionFactorFilename):
        """
            read in the data
            normalizes by dividing each data point by the sum of all data
            points times the bin width
        """
        channel = []
        counts = []

        legit_channels = []
        legit_counts = []
        with open(conversionFactorFilename) as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                channel.append(float(row['Channel']))
                counts.append(float(row['Counts']))

        return channel, counts

    def find_time_per_channel(self,counts):
        """
            find normalization factor
        """
        x=numpy.arange(3,5)
        foo=scipy.signal.find_peaks_cwt(counts,x)
        y=numpy.diff(foo)
        avg_channel_diff = numpy.mean(y)
        calibration_constant = 0.16E-6 # 0.16 usec
        time_per_channel = calibration_constant/avg_channel_diff

        return time_per_channel

    def convert_channels_to_time(self,channels,time_per_channel):
        """
            normalize the data
        """
        time_axis = [i * time_per_channel for i in channels]
        return time_axis
        
    def polynomial_fit(self, x, y):
        """
            call the other methods and compute half-life
        """
        # linearize the data
        ylog = []
        xlog = []
        for i in range(len(y)):
            if (y[i]!=0):
                ylog.append(numpy.log(y[i]))
                xlog.append(x[i])

        xlin = xlog[40:1200]
        ylin = ylog[40:1200]
            
        # polynomial fit, order 1
        pcoeffs = numpy.polyfit(xlin, ylin, 1)
        tau = -1./pcoeffs[0]
        yfit = [i*pcoeffs[0]+pcoeffs[1] for i in xlog]

        # plot polynomial fit
        # plot the log data and the whole fit
        plt.plot(
            xlog,ylog,'r.',
            xlog,yfit,'k--',linewidth=2.0
        )
        plt.axis([0,.15E-5,0,9])
        plt.show()

        # plot only the linear part of the log data with the fit
        yfit = [i*pcoeffs[0]+pcoeffs[1] for i in xlin]
        plt.plot(
            xlin,ylin,'r.',
            xlin,yfit,'k--',linewidth=2.0
        )
        plt.show()
        return tau

    def get_lifetime_from_tau(self,tau):
        """
            given input tau, compute lifetime
        """
        print "tau =",tau
        print "T =",tau*numpy.log(2)

if __name__ == '__main__':
    conversionFactorFilename = "conversion_factor.csv"
    # initialize
    foo = calculateHalfLife()

    # read file for normalization
    channels, counts = foo.read_CSV( 
        conversionFactorFilename,
    )

    # find normalization
    time_per_channel = foo.find_time_per_channel(
        counts,
    )

    # read data file
    dataFilename = "Take2_readable.csv"
    channels2, counts2 = foo.read_CSV(
        dataFilename,
    )
    dataFilename = "Take3_readable.csv"
    channels3, counts3 = foo.read_CSV(
        dataFilename,
    )
    dataFilename = "Take4_readable.csv"
    channels4, counts4 = foo.read_CSV(
        dataFilename,
    )

    # join data files on channel
    if (channels4==channels2) and (channels3==channels2):
        data_counts = [counts2[i]+counts3[i]+counts4[i] for i in range(len(counts2))]
        # normalize data file
        time_axis = foo.convert_channels_to_time(
            channels2,
            time_per_channel,
        )
        # polynomial fit
        tau = foo.polynomial_fit(
            time_axis,
            data_counts,
        )
        foo.get_lifetime_from_tau(tau)
    else:
        print "Couldn't join data files. Try again."
