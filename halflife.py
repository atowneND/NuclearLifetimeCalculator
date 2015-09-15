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

        xlin = xlog[30:1200]
        ylin = ylog[30:1200]
            
        # polynomial fit, order 1
        pcoeffs = numpy.polyfit(xlin, ylin, 1)
        tau = -1./pcoeffs[0]
        print pcoeffs
        print "tau =",tau
        print "T =",tau*numpy.log(2)
        yfit = [i*pcoeffs[0]+pcoeffs[1] for i in xlog]
        plt.plot(
            xlog,ylog,'r.',
            xlog,yfit,'k--',linewidth=2.0
        )
        plt.axis([0,.15E-5,0,8])
        plt.show()

    def plot_data(self, xvals, yvals):
        """
            plot any dataset given by x and y
        """
        plt.plot(
            xvals,
            yvals,
            'r-',
        )
        plt.show()

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
    # normalize data file
    time_axis = foo.convert_channels_to_time(
        channels2,
        time_per_channel,
    )
    # plot the data and analyze
#    foo.plot_data(
#        channels2,
#        counts2,
#    )

    foo.polynomial_fit(
        time_axis,
        counts2,
    )

