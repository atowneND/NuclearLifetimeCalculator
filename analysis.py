#!/usr/bin/python
# halflife.py

import scipy.signal
import numpy
import csv
import matplotlib.pyplot as plt

class modernPhysicsLabData:
    # __init__
    # read_CSV - read channels,counts
    # peaks2line_weighted_avg - using a weighted average, find center of peaks and slope of line through those peaks
    # convert_channels_to_time - scale channels to time

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

    def peaks2line_weighted_avg(self,channels,counts,calibration_constant):
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
        #calibration_constant = 0.16E-6 # 0.16 usec
        #calibration_constant = 0.01E-6 # 0.01 usec

        # known time between peaks
        time_list = numpy.arange(len(peak_channels))*calibration_constant

        # linear regression fit
        slope, intercept, r_val, p_val, std_err = scipy.stats.linregress(peak_channels,time_list)

        return slope, intercept, std_err, peak_channels

    def convert_channels_to_time(self,channels,time_per_channel):
        """
            normalize the data
        """
        time_axis = [i * time_per_channel for i in channels]
        return time_axis

class nuclearLifetime:
    # __init__
    # log_fit
    # plot_log_fit
    # set_window
    # get_noise_level
    # get_lifetime_from_tau
    # get_data

    def __init__(self, **kwargs):
        self.__dict__.update(**kwargs)

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

    def get_data(self,lab):
        """
            read data and join
        """
        # read data file
        dataFilename = "nuclife_data/Take2_readable.csv"
        channels2, counts2 = lab.read_CSV(
                dataFilename,
                )
        dataFilename = "nuclife_data/Take3_readable.csv"
        channels3, counts3 = lab.read_CSV(
                dataFilename,
                )
        dataFilename = "nuclife_data/Take4_readable.csv"
        channels4, counts4 = lab.read_CSV(
                dataFilename,
                )

        # join data files on channel
        if (channels4==channels2) and (channels3==channels2):
            data_counts = [counts2[i]+counts3[i]+counts4[i] for i in range(len(counts2))]
            # normalize data file
            time_axis = lab.convert_channels_to_time(
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

    def run_light_speed(self,lab):
        # read calibration file
        conversionFactorFilename = "lightspeed_data/conversion_factor.csv"
        channel, counts = lab.read_CSV(conversionFactorFilename)

        # find conversion
        mult_const = 0.01E-6 # 0.01E-6 usec
        time_per_channel, calibration_intercept, calibration_error, calibration_peaks = lab.peaks2line_weighted_avg(channel,counts,mult_const)

        # read data
        datafilename = "lightspeed_data/Five_Peaks_readable.csv"
        data_channels, data_counts = lab.read_CSV(datafilename)
        data_time = lab.convert_channels_to_time(data_channels, time_per_channel)
        #plt.plot(data_time,data_counts)
        #plt.show()

        # weighted average of five peaks
        t_center, data_intercept, data_error, data_channel_peaks = lab.peaks2line_weighted_avg(data_channels,data_counts,1)
        data_time_peaks = lab.convert_channels_to_time(data_channel_peaks,time_per_channel)
        plt.plot(data_time_peaks)
        plt.show()
        x = numpy.arange(100.0)
        y = t_center*x + data_intercept
        plt.plot(x,y)
        plt.show()

        # channel v time
#        plt.plot(data_channel_peaks, 
        # channel v position
        # error

if __name__ == '__main__':
    conversionFactorFilename = "nuclife_data/conversion_factor.csv"
    # initialize
    # can only run light speed, not nuclear lifetime
#    nuclife = nuclearLifetime()
    lab = modernPhysicsLabData()
    light = speedOfLight()
    light.run_light_speed(lab)

################################################################
## MOVE TO FUNCTION
#    # read file for normalization
#    channels, counts = lab.read_CSV( 
#            conversionFactorFilename,
#            )
#
#    # find normalization
#    mult_const = 0.16E-6 # 0.16 usec
#    time_per_channel, calibration_intercept, calibration_error, calibration_peaks = lab.peaks2line_weighted_avg(channels,counts,mult_const)
#    time_axis, data_counts = nuclife.get_data(lab)
#    noise = nuclife.get_noise_level(data_counts)
#    t_norm, counts_norm = nuclife.set_window(time_axis,data_counts)
#    popt , pcov = nuclife.log_fit(t_norm,counts_norm - noise)
#    nuclife.plot_log_fit(time_axis,data_counts,popt,pcov,noise)
