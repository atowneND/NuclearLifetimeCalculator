#!/usr/bin/python
# halflife.py

import scipy.signal
import numpy
import csv
import matplotlib.pyplot as plt

class calculateHalfLife:
    # __init__
    # normalize_data
    # plot_data
    # calculate_quantum_bins_histogram
    # process_semiclassical_data

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
        
    def calculate_halflife(self, x, y):
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
#        self.plot_data(
#            x,
#            y,
#        )
#        self.plot_data(
#            xlog,
#            ylog,
#        )
        self.plot_data(
            xlog[30:1200],
            ylog[30:1200],
        )
            
        # least squares fit
        def lin_fit(x, a, b):
            return a * x + b
        def log_fit(x, a, b, t):
            return a*numpy.exp(-x/t)+b
        popt, pcov = scipy.optimize.curve_fit(
            lin_fit,
            xlog[30:1200],
            ylog[30:1200],
        )
        print popt, pcov
#        popt, pcov = scipy.optimize.curve_fit(
#            log_fit,
#            x,
#            y,
#        )
#        self.plot_data(
#            xlog,
#            ylog,
#        )
#        print popt, pcov

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

    foo.calculate_halflife(
        channels2,
        counts2,
    )

