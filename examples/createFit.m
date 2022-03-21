function [fitresult, gof] = createFit(tc, eind, x0)
%CREATEFIT1(TC,EIND)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : tc
%      Y Output: eind
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 12-Mar-2022 14:55:17


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( tc, eind );

% Set up fittype and options.
ft = fittype( 'm * x + i0 + (c * sin(b * x - x1) + d * cos(b * x - x2)) * exp(-(x - x0)^2/sigma^2)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
% opts.StartPoint = [0.758272677307835 0.896825321017695 0.278498218867048 0.957166948242946 0.970592781760616 0.546881519204984 0.957506835434298 0.964888535199277 0.157613081677548];
opts.StartPoint = x0;
% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, xData, yData );
% legend( h, 'eind vs. tc', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
% % Label axes
% xlabel( 'tc', 'Interpreter', 'none' );
% ylabel( 'eind', 'Interpreter', 'none' );
grid on


