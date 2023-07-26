
% function myLeastSquaresRoutine

%**************************************************************************
%**************************************************************************
%**************************************************************************

%** Description of the program **
% This program is an implementation for estimating a set of parameters
% of a linear model "M(x) = mx+b" via the least squares approach

%** Input variables **
% None

%** Output variables **
% None

%**************************************************************************
%**************************************************************************
%**************************************************************************


%**************************************************************************
%************ Declaration and initialization of parameters ****************
%**************************************************************************

clc       % Clear the display in the command window.
close all % Clear all stored variables from the Workspace.
clear all % Close all figures.

%** Deterministic Model Parameters **
m = 3;        % slope of the line
b = 5;       % vertical or y-intercept
q = [m b]'; % Vector of model parameters
p = length(q);     % dimension of vector q

%** Independent variables **
% X variable in _ _ _ Units
xMin = 0.0;            % lower bound for the x variable
dx = 1;              % step size in x
xMax = 20;       % upper bound for the x variable
x = [xMin:dx:xMax]';   % x vector
n = length(x); % Total number of points in x

%**************************************************************************
%**************************************************************************
%**************************************************************************


%**************************************************************************
%************************ generating synthetic data ***********************
%**************************************************************************

nl = 1.5;%0.05; % noise level or standard error of epsilon
% nl = 0.05;%1.0509;%5.0325;%0.05; % noise level or standard error of epsilon

% % noise structure
xi = 0; % Absolute noise.

%** Output **

detModel = deterministicModel(x,q);

epsilon = (nl*randn(n,1)); % epsilon ~ N(0,nl^2)

E = detModel.^(xi).*epsilon; % Error
observedData = detModel + E;  % General Statistical Model

% data plots
figure(1)
a = plot(x,observedData,'ro',x,detModel,'bo--');
legend(a,'Observated Data','Deterministic Model')
xlabel('X (Units)')
ylabel('Y (Units)')
xlim([xMin xMax])
ylim([0 ceil(max([max(observedData),max(detModel)]))+10])
title('Title: X vs. Y')

%**************************************************************************
%**************************************************************************
%**************************************************************************


%**************************************************************************
%******************** Least Squares Procedure *****************************
%**************************************************************************

statModelParameters = [xi,p]';

q_opt = 5*ones(p,1); % initial guess
q = q_opt;

w = ones(n,1);
    
[q_opt, J_q_opt] = fminsearch(@(q) objectiveFunctional(x,observedData,q,statModelParameters,w),q);
    
    %** Output **
    detModel = deterministicModel(x,q_opt);

    % weight formulation
    w = detModel.^(-2*xi);
    
    figure(2)
    subplot(2,2,[1 3])
    plot(5,5,'xr') % initial guess
    hold on
    plot(m,b,'ob'); % "true values"
    hold on
    plot(q_opt(1),q_opt(2),'*r');
    hold on
    xlabel('slope m')
    ylabel('y-intercept b')

    figure(3)
    a = plot(x,observedData,'ro',x,detModel,'bo-');
    legend(a,'Observed Data','Deterministic model')
    xlabel('X (Units)')
    ylabel('Y (Units)')
    xlim([xMin xMax])
    ylim([0 ceil(max([max(observedData),max(detModel)]))+10])
    title('Title: X vs. Y')
    
    %** Residual plots **
    
    % Residual Plots (CV & NCV)
    
    % 1) Constant variance xi = 0
    
    figure(4)
    subplot(1,2,1)
	plot(x,observedData-detModel,'o')
    title('Residuals over X values plot assuming constant variance')
%    title('Residual over X values plot assuming non-constant variance')
    xlabel('X (Units)')
    ylabel('Residuals')
    xlim([xMin xMax])
    subplot(1,2,2)
    plot(detModel,observedData-detModel,'o')
    title('Deterministic Model vs. Residual assuming constant variance')
%    title('Deterministic Model vs. Residual assuming constant variance')    
    xlabel('Deterministic Model')
    xlim([min(detModel) max(detModel)]);
    ylabel('Residuals')

%**************************************************************************
%*************** end function MyParameterEstimationSIR ********************
%**************************************************************************


%**************************************************************************
%************************** List of functions *****************************
%**************************************************************************

%**************************************************************************
%**************** Deterministic Model ***************************
%**************************************************************************
function output = deterministicModel(x,q)
% m = q(1); slope of the line
% b = q(2); vertical or y-intercept 
output = q(1).*x + q(2);
end
%**************************************************************************
%**************************************************************************
%**************************************************************************

%**************************************************************************
%*********************** Objective Functional (or Cost Function)*****************************
%**************************************************************************
function output = objectiveFunctional(x,observedData,q,statModelParameters,w)

%** Output **
detModel = deterministicModel(x,q);

% weights formulation
xi = statModelParameters(1);
w = detModel.^(-2*xi);

% Euclidean norm
output = norm( sqrt(w).*(observedData - detModel) )^2;
end
%**************************************************************************
%**************************************************************************
%**************************************************************************
