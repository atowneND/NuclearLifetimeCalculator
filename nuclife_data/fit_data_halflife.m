time_per_channel = 2.411E-10;

filename2 = '/home/ashley/Documents/SeniorYear/ModPhysILab/NucLifetime/halflife/Take2_readable.csv';
filename3 = '/home/ashley/Documents/SeniorYear/ModPhysILab/NucLifetime/halflife/Take3_readable.csv';
filename4 = '/home/ashley/Documents/SeniorYear/ModPhysILab/NucLifetime/halflife/Take4_readable.csv';

x2 = csvread(filename2,1,0);
x3 = csvread(filename3,1,0);
x4 = csvread(filename4,1,0);

counts = x2(:,2) + x3(:,2) + x4(:,2) - 365.592;
channels = x2(:,1);
plot(channels,counts);
t_tic = channels * time_per_channel;
a = 5472;
b = 7.29E6;
x = linspace(0,2E-6,1000);
y = a*exp(-b*x)+365.592;

plot(t_tic,counts+365.592,'k.',x,y,'r-')
xlabel('Time \mus');
ylabel('Counts');
title('Counts vs. Time');
legend('Measured Data','Exponential Fit');
TeXString = texlabel('Fit Equation:');
text(.5E-6,4500,TeXString)
TeXString = texlabel('f(t)=5472e^{-7.29E6t}+365.592');
text(.5E-6,4200,TeXString)

TeXString = texlabel('Chi-Squared:');
text(.5E-6,3500,TeXString)
TeXString = texlabel('3.25E6');
text(.5E-6,3200,TeXString)

TeXString = texlabel('Standard Deviation:');
text(.5E-6,2500,TeXString)
TeXString = texlabel('55.23');
text(.5E-6,2200,TeXString)
