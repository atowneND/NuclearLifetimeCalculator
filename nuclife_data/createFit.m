function [fitresult, gof] = createFit(x, y)
%CREATEFIT(X,Y)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : x
%      Y Output: y
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 21-Sep-2015 17:40:41


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( x, y );

% Set up fittype and options.
ft = fittype( 'exp1' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [5486.21168067942 -7313223.68616421];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name', 'Nuclear Lifetime' );
h = plot( fitresult, xData, yData );
legend( h, 'Counts vs. Time', 'Computed counts', 'Counts', 'Location', 'NorthEast' );
% Label axes
xlabel('time micro seconds')
ylabel('counts')
grid on


