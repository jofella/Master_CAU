%===================| Tutorial III: Numerical Methods |====================
%======================| Computer-Based Exercise 1 |=======================
% This script applies importance sampling to simulate the moments of a
% truncated standard normal distribution.  

% Clear workspace memory and command window:
clear
clc

% Set truncation bounds of the target pdf: -- given in the exercise
a = -1;  
b = 1.5;

% Set number of replications: -- tool for researcher
S = 100000;

% Draw S standard normal random numbers and compute weights:
theta = randn(S, 1);             %column vector, most efficent way
w     = (theta>a).*(theta<b);   %if-and-check, output 1 - must lie within the bound(truncated normal)

% Compute mean of the target pdf:
theta_bar = mean(w.*theta)./mean(w);

% Compute variance of the target pdf:
S2_theta = mean(w.*(theta-theta_bar).^2)./mean(w);

% Analytical results: (posterior moments) - truncated standard normal
ma = (normpdf(a,0,1)-normpdf(b,0,1))/(normcdf(b,0,1)-normcdf(a,0,1));
va = 1 + (a*normpdf(a,0,1)-b*normpdf(b,0,1))/(normcdf(b,0,1)-normcdf(a,0,1))...
    - ((normpdf(a,0,1)-normpdf(b,0,1))/(normcdf(b,0,1)-normcdf(a,0,1)))^2;

% Display results:
fprintf(1,'\n\n Results after %8.0f replications: \n\n',S)
fprintf(1,' Theoretical mean     = %7.4f, Simulated mean     = %7.4f \n',[ma, theta_bar])
fprintf(1,' Theoretical variance = %7.4f, Simulated variance = %7.4f \n',[va, S2_theta])
