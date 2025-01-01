%=====| Tutorial I: Introduction to Bayesian Computation |=====
%================| Computer-Based Exercise 1 |=================
% This script estimates pi by means of Monte Carlo integration.

% Clear workspace memory and command window:
clear
clc

% Set number of replications S:
S = 1000000;                  

% Draw uniform random points (x^(s),y^(s)) inside the square A:
x = 2*rand(S,1)-1;  % x coordinate draw
y = 2*rand(S,1)-1;  % y coordinate draw

% Check whether z^(s)=(x^(s),y^(s)) lies inside the circle:
P_z = x.^2 + y.^2;     % distance of each point z_s from the center (0,0)
isinside = P_z <= 1;   % 1 if inside, 0 else

% Estimate pi by means of Monte Carlo integration:
pi_hat = 4*mean(isinside);  % estimate of pi

% Report results:
disp(['Replications: ' int2str(S) ', pi_hat : ' num2str(pi_hat)...
    ', Absolute Error: ' num2str(abs(pi-pi_hat))])

% Part (b):
% Is the theoretical numerical standard error correct?
% If yes, then in 95% of repeated MC simulations (as the one we performed
% above) should yield an estimate of pi_hat such that the approximation
% error |pi-pi_hat|<0.01 if we choose S=103,599.

S   = 103599;                % number of replications for each Monte Carlo integration
rep = 10000;                 % number of Monte Carlo integrations
pi_hat_mean = zeros(rep,1);  % pre-allocate column-vector to store the estimates for pi

tic;
parfor i = 1:rep
    x = 2*rand(S,1)-1;       % x coordinate drawn from uniform distribution
    y = 2*rand(S,1)-1;       % y coordinate drawn from uniform distribution
    P_z = x.^2 + y.^2;       % distance from center
    isinside = P_z <= 1;     % 1 if inside, 0 else
    pi_hat_mean(i) = 4*mean(isinside);  % estimate of pi at iteration i   
end
toc;

pi_hat_ok = mean(abs(pi_hat_mean-pi)<0.01);
disp([int2str(rep) ' repeated Monte Carlo integrations. '...
    num2str(100*pi_hat_ok) ' percent of the estimates less than 0.01 away from pi.'])
