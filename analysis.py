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

    def find_time_per_channel(self,channels,counts):
        """
            find normalization factor
        """
        # plotting lists
        peak_channels = []

        # sums and counts
        ctr = 0
        peak_sum = 0
        count_sum = 0
        last_zero = 0

        # loop
        for i in numpy.arange(len(channels)):
            if (counts[i]!=0):
                if (channels[i]-channels[last_zero]==1):
                    peak_sum = peak_sum + channels[i]*counts[i]
                    count_sum = count_sum + counts[i]
                    ctr = ctr + 1
                else:
                    if (ctr!=0):
                        peak_channels.append(peak_sum/count_sum)
                    count_sum = 0
                    peak_sum = 0
                    ctr = 0
                last_zero = i

        # make this an argument passed to this method!!!!
        calibration_constant = 0.16E-6 # 0.16 usec
        calibration_constant = 0.01E-6 # 0.01 usec

        # known time between peaks
        time_list = numpy.arange(len(peak_channels))*calibration_constant

        # linear regression fit
        slope, intercept, r_val, p_val, std_err = scipy.stats.linregress(peak_channels,time_list)

        return slope, std_err

    def convert_channels_to_time(self,channels,time_per_channel):
        """
            normalize the data
        """
        time_axis = [i * time_per_channel for i in channels]
        return time_axis

    def log_fit(self, x, y):
        """ 
            from data, compute half-life
        """
        def func(x, a, b):
            return a*numpy.exp(-b*x)

        popt, pcov = scipy.optimize.curve_fit(func, numpy.asarray(x), numpy.asarray(y))
        tau = 1./popt[1]
        self.get_lifetime_from_tau(tau)

        #print numpy.std(numpy.diag(pcov))
        return popt, pcov

    def plot_log_fit(self, x, y, popt, pcov, noise):
        """
            plot data and fit
        """
        yfit = popt[0]*numpy.exp(-(popt[1]*x)) + noise
        stdev = numpy.sqrt(sum([i**2 for i in numpy.diag(pcov)]))
        stdev = numpy.std(numpy.diag(pcov))
        print stdev
        plt.plot(
                x,y,'b.',
                x,yfit,'r-',linewidth=2.5,
                )
        plt.show()
        chisq,p =scipy.stats.chisquare(y,yfit)
        print chisq,p

    def set_window(self, x, y):
        """
            return window for fit
        """
        channels = x[150:1200]
        counts = y[150:1200]
        return channels, counts

    def get_noise_level(self, y):
        """
            return window for noise data
        """
        counts = y[4500:7500]
        avg_noise = numpy.mean(counts)
        return avg_noise 

    def get_lifetime_from_tau(self,tau):
        """
            given input tau, compute lifetime
        """
        print "tau =",tau
        print "T =",tau*numpy.log(2)

    def get_data(self):
        """
            read data and join
        """
        # read data file
        dataFilename = "nuclife_data/Take2_readable.csv"
        channels2, counts2 = foo.read_CSV(
                dataFilename,
                )
        dataFilename = "nuclife_data/Take3_readable.csv"
        channels3, counts3 = foo.read_CSV(
                dataFilename,
                )
        dataFilename = "nuclife_data/Take4_readable.csv"
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
        else:
            print "Couldn't join data files. Try again."
            exit

        return time_axis, data_counts

class speedOfLight:
    def __init__(self, **kwargs):
        self.__dict__.update(**kwargs)

    def run_light_speed(self,data):
        print "here"
        # read calibration file
        conversionFactorFilename = "lightspeed_data/conversion_factor.csv"
        channel, counts = data.read_CSV(conversionFactorFilename)

        # find conversion
        time_per_channel, calibration_error = data.find_time_per_channel(channel,counts)

        # read data
        datafilename = "lightspeed_data/Five_Peaks_readable.csv"
        data_channels, data_counts = data.read_CSV(datafilename)
        data_time = data.convert_channels_to_time(data_channels, time_per_channel)

        # isolate peaks
        # error

if __name__ == '__main__':
    conversionFactorFilename = "nuclife_data/conversion_factor.csv"
    # initialize
    # can only run light speed, not nuclear lifetime
    foo = calculateHalfLife()
    bar = speedOfLight()
    bar.run_light_speed(foo)

###############################################################
# MOVE TO FUNCTION
#    # read file for normalization
#    channels, counts = foo.read_CSV( 
#            conversionFactorFilename,
#            )
#
#    # find normalization
#    time_per_channel, calibration_error = foo.find_time_per_channel(
#            channels,
#            counts,
#            )
#
#    time_axis, data_counts = foo.get_data()
#    noise = foo.get_noise_level(data_counts)
#    t_norm, counts_norm = foo.set_window(time_axis,data_counts)
#    popt , pcov = foo.log_fit(t_norm,counts_norm - noise)
##    foo.plot_log_fit(time_axis,data_counts,popt,pcov,noise)
