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
        true_peak_channels = []
        true_peak_counts = []
        all_peak_channels = []
        all_peak_counts = []

        # sums and counts
        ctr = 0
        peak_sum = 0
        count_sum = 0
        last_nonzero = 0

        # loop
        for i in numpy.arange(len(channels)):
            # if anything was detected
            if (counts[i]!=0):
                if ctr==0: # first bin for the peak
                    peak_sum = channels[i]*counts[i]
                    count_sum = counts[i]
                elif (channels[i]-channels[last_nonzero]==1): # for adjacent bins
                    peak_sum = peak_sum + channels[i]*counts[i]
                    count_sum = count_sum + counts[i]
                # increment number of adjacent nonzero bins
                ctr = ctr + 1

                last_nonzero = i

            # for zero counts
            else:
                # if the counter isn't 0 => a peak was detected - ignore unfinished peaks
                if (ctr!=0):
                    # add last peak
                    all_peak_channels.append(peak_sum/count_sum)
                    all_peak_counts.append(count_sum)

                    # reset counters
                    ctr = 0
                    count_sum = 0
                    peak_sum = 0
            
        peak_max = max(all_peak_counts)
        for i in numpy.arange(len(all_peak_channels)):
            if (all_peak_counts[i] > 0.4 * peak_max):
                true_peak_channels.append(all_peak_channels[i])
                true_peak_counts.append(all_peak_counts[i])

        npoints = len(true_peak_counts)
        time_list = numpy.arange(len(true_peak_channels))*calibration_constant

#        print npoints
        if (npoints < 3):
            if (npoints == 2):
                slope = (time_list[1]-time_list[0])/(true_peak_channels[1]-true_peak_channels[0])
                intercept = time_list[1]-true_peak_channels[1]*slope
                r_val = 1
                std_err = 0
            else:
                print "error"
        else:
            # linear regression fit
            slope, intercept, r_val, p_val, std_err = scipy.stats.linregress(true_peak_channels,time_list)

        return slope, intercept, std_err, true_peak_channels

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
#        conversionFactorFilename = "foo.csv"
        channel, counts = lab.read_CSV(conversionFactorFilename)

        # read distance data
        distanceFilename = "lightspeed_data/distance.csv"
        peak, distance = lab.read_CSV(distanceFilename)

        # find conversion
        mult_const = 0.01E-6 # 0.01E-6 usec
        time_per_channel, calibration_intercept, calibration_error, calibration_peaks = lab.peaks2line_weighted_avg(channel,counts,mult_const)
        m_D = 1/time_per_channel

        # plot channel vs time
        x = numpy.arange(2000)
        y = time_per_channel*x + calibration_intercept
        t = numpy.arange(len(calibration_peaks)) * mult_const
        chvt_fit = plt.plot(y,x,'k-',label='fit')
        chvt_data = plt.plot(t,calibration_peaks,'r.',label='data')
        plt.xlabel('Time (s)')
        plt.ylabel('Channel')
        plt.title('Channel vs. Time')
        plt.legend(loc=2)
        plt.xlim(-0.1*1.2E-8,1.2E-8)
        plt.ylim(1000,2000)
#        plt.show()

#        # read data
        datafilename = "lightspeed_data/Five_Peaks_readable.csv"
        data_channels, data_counts = lab.read_CSV(datafilename)
        data_time = lab.convert_channels_to_time(data_channels, time_per_channel)

        # weighted average of five peaks
        t_center, data_intercept, data_error, data_channel_peaks = lab.peaks2line_weighted_avg(data_channels,data_counts,1)
        data_time_peaks = lab.convert_channels_to_time(data_channel_peaks,time_per_channel)
        x = numpy.arange(100.0)
        y = t_center*x + data_intercept

#        plt.plot(x,y,'r-')
#        plt.show()

        # channel v position
        m_ch_v_pos, a,b,c,d = scipy.stats.linregress(data_channel_peaks,distance)
        m_S = m_ch_v_pos/100
        print m_D, m_S, 2 * m_D / m_S
        plt.plot(distance,data_channel_peaks,'r.')
        plt.ylabel('Channel')
        plt.xlabel('Distance (cm)')
#        plt.show()
        
#        #plt.plot(data_channel_peaks,distance)
#        # divide slope of ^ by time_per_channel
#        # error - fit gaussian to five peaks

if __name__ == '__main__':
    conversionFactorFilename = "nuclife_data/conversion_factor.csv"
    conversionFactorFilename = "foo"
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
