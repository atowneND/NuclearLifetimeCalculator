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

#        print "chan:",true_peak_channels
#        print "diff:",numpy.diff(true_peak_channels)
        npoints = len(true_peak_counts)
        time_list = numpy.arange(len(true_peak_channels))*calibration_constant

#        print npoints
        if (npoints < 3):
            if (npoints == 2):
                slope = (time_list[1]-time_list[0])/(true_peak_channels[1]-true_peak_channels[0])
                errfac = 3.1937389/(true_peak_channels[1]-true_peak_channels[0])
                print "calibration slope error =",errfac*slope
                print "     slope =",slope
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
        channel, counts = lab.read_CSV(conversionFactorFilename)
        plt.plot(channel,counts)
        plt.xlabel('Channel')
        plt.ylabel('Counts')
        plt.title('Calibration Data')
        plt.show()

        # read distance data
        distanceFilename = "lightspeed_data/distance.csv"
        peak, distance_cm = lab.read_CSV(distanceFilename)
        distance = map(lambda distance_cm: distance_cm/100., distance_cm)

        # find conversion
        mult_const = 0.01E-6 # 0.01E-6 usec
        time_per_channel, calibration_intercept, calibration_error, calibration_peaks = lab.peaks2line_weighted_avg(channel,counts,mult_const)
        m_D = 1/time_per_channel
        print m_D

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
        plt.text(5E-9,1100,"y = 8.124E10 * x + 1.345E-11")
        plt.show()

        # read data
        datafilename = "lightspeed_data/Five_Peaks_readable.csv"
        data_channels, data_counts = lab.read_CSV(datafilename)
        data_time = lab.convert_channels_to_time(data_channels, time_per_channel)
        plt.plot(data_channels,data_counts)
        plt.xlabel('Channel')
        plt.ylabel('Counts')
        plt.title('Raw Data')
        plt.show()

        # weighted average of five peaks
        t_center, data_intercept, data_error, data_channel_peaks = lab.peaks2line_weighted_avg(data_channels,data_counts,1)
        data_time_peaks = lab.convert_channels_to_time(data_channel_peaks,time_per_channel)
        x = numpy.arange(100.0)
        y = t_center*x + data_intercept

        # channel v position
        m_ch_v_pos, a,b,c,data_err = scipy.stats.linregress(distance,data_channel_peaks)
        m_S = m_ch_v_pos
        print m_S
        print "error =",data_err
        c = 2 * m_D / m_S
        print "speed of light =",c
        c_ref = 2.99792458E8
        print "diff =",c-c_ref
        print "%diff =",(c-c_ref)/c_ref
        
        x = numpy.linspace(0,max(distance)+0.5,100)
        y = m_ch_v_pos * x + a
        plt.plot(distance,data_channel_peaks,'r.',label='Data')
        plt.plot(x,y,'k-',label='Fit')
        plt.legend(loc=2)
        plt.text(0.05,1400,"y = 532.22 * x + 134.77")
        plt.title('Channel vs. Distance')
        plt.ylabel('Channel')
        plt.xlabel('Distance (m)')
        plt.show()

class gammaSpectroscopy:
    def __init__(self, **kwargs):
        self.__dict__.update(**kwargs)

    def run_gamma_spec(self,lab):
        filename_list = ("gammaspec_data/Eu.csv", "gammaspec_data/Co.csv", "gammaspec_data/Ho.csv", "gammaspec_data/KCl.csv")
        for fd in xrange(len(filename_list)):
            channel, counts = lab.read_CSV(filename_list[fd])
            print ""
            print filename_list[fd]
            if (fd == 0): #Eu
                x=0
                energy_per_channel = self.calibrateEu(channel,counts) # keV
                print energy_per_channel
                plt.title("$^{152}$Eu")
            elif (fd == 1): #Co
#                newcounts = []
#                for x in xrange(len(counts)):
#                    if (counts[x] < 0.5*max(counts)):
#                        newcounts.append(0)
#                    else:
#                        newcounts.append(counts[x])
#
#                peakind = scipy.signal.find_peaks_cwt(newcounts,numpy.arange(75,80),noise_perc = 90)
#                peak_centers = []
#                for i in peakind:
#                    print i
#                    max_counts, foo = (max(newcounts[i-40:i+40]))
#                    #peak_centers.append(channel[max_counts])

                peakinds = (4736,5382)
                resolution,res_err = self.compute_resolution(channel,counts,peakinds)
                peak_energies = [i*energy_per_channel for i in peakinds]
                print "peaks = ",peak_energies
                print "resolution =",resolution," +-",res_err
                peak2compton_ratio = self.peak_to_compton(energy_per_channel,counts,peak_energies)
                print "peak to compton =",peak2compton_ratio
                plt.title("$^{60}$Co")
            elif (fd == 2): #Ho
                peakind = scipy.signal.find_peaks_cwt(counts,numpy.arange(30,50))
                E_peaks = [i*energy_per_channel for i in peakind]

                peakinds = (733,1121)
                
                print "peaks = ",[i*energy_per_channel for i in peakinds]
                peak_energies = [i*energy_per_channel for i in peakinds]
                all_energies = [i*energy_per_channel for i in channel]

                peak_energies = [87.62,180.9,276.7,707.4,805.6,946.1,1236.3]
                red_number = [1,2,3,4,6,8,9]
                peak_energies = [87.62, 180.9, 276.7]
                red_number = [1,2,3]
                self.moment_of_inertia(peak_energies,red_number)
                plt.title("$^{166}$Ho")
#                plt.plot(all_energies,counts)
#                plt.xlabel(E_peaks)
#                plt.show()
            elif (fd == 3): #K
                #peakinds = scipy.signal.find_peaks_cwt(counts,numpy.arange(80,130))
                peakinds = 5890
                print "peaks =",peakinds*energy_per_channel
                resolution, res_err = self.compute_single_resolution(channel,counts,peakinds)
                print "resolution =",resolution," +-",res_err
                snr = self.peak_to_background(channel,counts,peakinds)
                print "snr =", snr
                plt.title("$^{40}$K")

#            font = {'family' : 'normal',
#                    'weight' : 'bold',
#                    'size'   : 22}
#            plt.rc('font',**font)
#            energies = 0
#            energies = [i*energy_per_channel for i in channel]
#            plt.plot(energies,counts,linewidth=2.0)
#            plt.xlabel("Energy (keV)")
#            plt.ylabel("Counts")
#            plt.show()
            
    def calibrateEu(self,channel,counts):
        known_peaks = [121.7817, 244.6974, 344.2785] # keV
        uncertainties = [.0003, .0008, .0012]
        channel_uncertainties = [1.0, 1.0, 1.0]

        peakind = scipy.signal.find_peaks_cwt(counts,numpy.arange(10,50))

        dchannel = numpy.diff(peakind[0:3])
        denergy = numpy.diff(known_peaks)

        calibration_constant = denergy[0]/dchannel[0]

        energy_error = uncertainties[0] + uncertainties[1]
        channel_error = channel_uncertainties[0] + channel_uncertainties[1]
        percent_energy_error = energy_error / denergy[0]
        percent_channel_error = channel_error / dchannel[0]
        total_error = (percent_energy_error + percent_channel_error) * calibration_constant
        print "calibration uncertainty =",total_error, percent_channel_error+percent_energy_error

        return calibration_constant

    def compute_resolution(self,channel,counts,peakind):
        print "resolution:"
        resolution = []
        errors = []
        for x in peakind:
            print "x =",x
            half_max = 0.5 * counts[x]
            i=1
            done = 0
            left = 0
            right = 0
            while (done!=1):
                left_counts = counts[x-i]
                right_counts = counts[x+i]
                if (left_counts<half_max) & (left != 1):
                    left_channel = channel[x-i+1]
                    left = 1
                if (right_counts<half_max) & (right != 1):
                    right_channel = channel[x+i-1]
                    right = 1
                done = left & right
                i = i + 1
                if (i>10):
                    break
            fwhm = right_channel-left_channel+1
            resolution.append(100. * fwhm / x)
            print "fwhm =",fwhm
            channel_uncertainties = [1.0, 1.0, 1.0]
            fwhm_channel_error = (channel_uncertainties[0] + channel_uncertainties[1]) / (fwhm)
            peak_channel_error = channel_uncertainties[2] / x
            errors.append((fwhm_channel_error+peak_channel_error)*100.*fwhm/x)

        return resolution,errors

    def compute_single_resolution(self,channel,counts,peakind):
#        plt.plot(channel,counts)
#        plt.xlabel(peakind)
        x = peakind
        print "x =",x
        half_max = 0.5 * counts[x]
        i=1
        done = 0
        left = 0
        right = 0
        while (done!=1):
            left_counts = counts[x-i]
            right_counts = counts[x+i]
            if (left_counts<half_max) & (left != 1):
                left_channel = channel[x-i+1]
                left = 1
            if (right_counts<half_max) & (right != 1):
                right_channel = channel[x+i-1]
                right = 1
            done = left & right
            i = i + 1
            if (i>50):
                break
        fwhm = right_channel-left_channel+1
        resolution = (100. * fwhm / x)
        print "fwhm =",fwhm,left_channel,right_channel
        channel_uncertainties = [1.0, 1.0, 1.0]
        fwhm_channel_error = (channel_uncertainties[0] + channel_uncertainties[1]) / (fwhm)
        peak_channel_error = channel_uncertainties[2] / x
        errors = (fwhm_channel_error+peak_channel_error)*100.*fwhm/x
#        print "error (percent, abs)=",fwhm_channel_error+peak_channel_error,(fwhm_channel_error+peak_channel_error)*100.*fwhm/x

        return resolution,errors

    def peak_to_background(self,channel,counts,peakind):
        noiserange = [2000,4000]
        noise_avg = numpy.mean(counts[noiserange[0]:noiserange[1]])
        return counts[peakind]/noise_avg

    def peak_to_compton(self,energy_per_channel,counts,peak_energies):
        compton_energy_range = [1040,1096]
        channel_range = [int(i / energy_per_channel) for i in compton_energy_range]
        compton_level = numpy.mean(counts[channel_range[0]:channel_range[1]])
        peak2_to_compton = peak_energies[1] / compton_level
        return peak2_to_compton

    def moment_of_inertia(self,peak_energies,red_number):
        print "moment of inertia"
        peak = peak_energies[0]
        print "peak =",peak
        hbar = 1.055E-34 # J*s
        I_exp = 6. * (hbar**2) / ( 2. * (peak * 1000.)*1.602E-19 ) # J
        print "ANS:",I_exp,"J;",I_exp * 1.602E-19,"eV"

        ll1 = []
        E_diffs = []
        hbar = 6.582E-16
        E_peaks = [i*1000 for i in peak_energies]
        print red_number,E_peaks
        for i in xrange(len(red_number)):
            ll1.append(red_number[i] * (red_number[i] + 1) * (hbar ** 2) / 2)

        E_diffs = E_peaks
        slope, intercept, r_val, p_val, std_err  = scipy.stats.linregress(ll1,E_diffs)
        print "slope =",slope
        I_exp = (1./slope) * 1.602E-19
        electron_mass = 0.00054858 # amu
        mass_atomic = 164.93031 # amu
        kg_per_amu = 1.6605E-27 # kg/amu
        r_o = 176E-12 # m

        mass = (mass_atomic - 67 * electron_mass) * kg_per_amu # kg
        r_avg = r_o * (166 ** (1./3.)) # m
        beta = 0
        print "mass, ravg =",mass,r_avg
        I_rigid = (2./5.) * mass * (r_avg ** 2) * (1 + 0.31 * beta) # J

        beta = 0.29
        I_fluid = (9./(8.*numpy.pi)) * mass * (r_avg ** 2) * (beta ** 2 ) # J

        print I_rigid," >",I_exp," >",I_fluid

#,        plt.plot(ll1,E_diffs)
#        plt.show()

class Compton:
    def __init__(self, **kwargs):
        self.__dict__.update(**kwargs)

    def runCompton(self,lab):
        date_list = [1110]
        theta_list = [20,30,45,60,75,90,105,120,135]
        calfile = "compton/calibration_" + str(date_list[0]) + ".csv"
        cal_files = ["compton/calibration_" + str(i) + ".csv" for i in date_list]

        for d in date_list:
            calfile = "compton/calibration_" + str(d) + ".csv"
            cal_channel, cal_counts = lab.read_CSV(calfile)
            if (d==1110):
                print "Calibrating:"
#                peak1 = scipy.signal.find_peaks_cwt(cal_counts[0:2000],numpy.arange(200,400))
#                peak2 = scipy.signal.find_peaks_cwt(cal_counts[5000:8000],numpy.arange(1000,1500))
                peak1=[328,1836]
                peak2 = [5117, 6448, 7821] 
                peaks = [328, 6448]
                energy_per_channel = self.calibrateCs(peaks) #keV
            for t in theta_list:
                channel, counts = self.filterNoise(t, d)
                energies = [ch * energy_per_channel for ch in channel]
                peaks = scipy.signal.find_peaks_cwt(counts,numpy.arange(600,1500))
#                print all_peaks
#                peaks = [i for i in peaks]
#                peaks = [6006, 6024, 6031, 6034]
#                plt.plot(channel,counts)
#                plt.xlabel(peaks)
#                plt.show()

                # find energy for each peak
#                peaks = [6006, 6024, 6031, 6034]
                peak_energies = [ch * float(energy_per_channel) for ch in peaks] #keV
                # plug into compton equation
                output_filename = "Cu_" + str(t) + "_peaks.data"
                f = open(output_filename,'w')
                f.write("E'(keV)\t E'/E\tcross section(m^-2)\n")
                for e in peak_energies:
                    f.write(str(e) + ' ')

                    e_gamma = self.comptonEquation(t, e)
                    ratio = e / e_gamma
                    f.write(str(ratio) + ' ')

                    cross_section = self.klein_nishina(t,e)
                    f.write(str(cross_section))
                    f.write('\n')

    def klein_nishina(self, theta, energy):
        m_e = 9.11E-31 #kg
        c = 3E8 #m/s
        alpha = energy / (m_e * c**2 * 1000/1.602E-19)
        r0 = 2.818E-15 #m

        dsigma = (r0**2 / (1+ alpha*(1-numpy.cos(theta)))**3) * (1+ (alpha*(1-numpy.cos(theta)))**2/((1+numpy.cos(theta)**2)*(1+alpha*(1-numpy.cos(theta)))))
        return dsigma

    def comptonEquation(self, theta, energy):
        m_e = 9.11E-31 #kg
        c = 3E8 #m/s
        mc2_kev = (m_e * c**2) / (1.602E-19 / 1000) #keV
        E_out = energy / (1- (energy / mc2_kev)*(1-numpy.cos(theta)))
        return E_out

    def filterNoise(self,angle,date):
        data_file = "compton/Cu_" + str(angle) + "_data_" + str(date) + ".csv"
        data_channel, data_counts = lab.read_CSV(data_file)
        noise_file = "compton/Cu_" + str(angle) + "_noise_" + str(date) + ".csv"
        noise_channel, noise_counts = lab.read_CSV(noise_file)
        channel = data_channel
        counts = [data_counts[index] - noise_counts[index] for index in xrange(len(data_counts))]
        return channel, counts

    def calibrateCs(self,peaks):
        known_energies = [32., 662.] #keV
        print peaks
        print known_energies
        energy_per_channel = numpy.diff(known_energies) / numpy.diff(peaks)
        return energy_per_channel

if __name__ == '__main__':
    conversionFactorFilename = "nuclife_data/conversion_factor.csv"
    conversionFactorFilename = "foo"
    # initialize
    # can only run light speed, not nuclear lifetime
#    nuclife = nuclearLifetime()
    lab = modernPhysicsLabData()
#    light = speedOfLight()
#    light.run_light_speed(lab)
#    gamma = gammaSpectroscopy()
#    gamma.run_gamma_spec(lab)
    scatter = Compton()
    scatter.runCompton(lab)

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
