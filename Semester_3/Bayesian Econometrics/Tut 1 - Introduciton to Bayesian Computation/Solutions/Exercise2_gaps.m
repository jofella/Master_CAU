%=====| Tutorial I: Introduction to Bayesian Computation |=====
%================| Computer-Based Exercise 2 |=================
% This script estimates posterior moments of theta ~ N(0,1) by Monte Carlo
% Integration.

% Clear workspace memory and command window:
clear
clc

%% Part (a)

% Set number of replications S:
S = 1000000; 

% Generate S random draws and estimate posterior moments:
theta  = randn(S,1);     % draw from N(0,1)
mu     = mean(theta);    % calculate posterior mean
sigma2 = var(theta);     % calculate posterior variance

% Report results:
disp(['Replications: ' int2str(S) ', estimate of mean: ' num2str(mu)...
    ', estimate of variance: ' num2str(sigma2)])

%% Part (b)

% Theoretical number of necessary replications S:
% (i)  P(|theta_bar-E(theta)|<0.001) = 0.95      =>  S_mu  = 3,841,600
% (ii) P(|sigma2_hat-Var(theta)|<0.001) = 0.95   =>  S_var = 7,683,199

S_mu   = 3841600;        % number of replications for each Monte Carlo integration
S_var  = 7683199;        % number of replications for each Monte Carlo integration
rep    = 1000;           % number of Monte Carlo integrations
mu     = zeros(rep,1);   % store posterior mean estimates of each Monte Carlo integration
sigma2 = zeros(rep,1);   % store posterior variance estimates of each Monte Carlo integration

tic;
parfor i = 1:rep
    theta_mu  = randn(S_mu,1);   % draw from N(0,1)
    theta_var = randn(S_var,1);  % draw from N(0,1)
    mu(i)     = mean(theta_mu);  % calculate posterior mean
    sigma2(i) = var(theta_var);  % calculate posterior variance
end
toc;

mu_ok     = mean(abs(mu-0)<0.001);      % fraction of times |theta_bar-E(theta)|<0.001
sigma2_ok = mean(abs(sigma2-1)<0.001);  % fraction of times |sigma2_hat-Var(theta)|<0.001

disp([num2str(mu_ok*100) '% of the times |theta_bar-E(theta)|<0.001'])
disp([num2str(sigma2_ok*100) '% of the times |sigma2_hat-Var(theta)|<0.001'])
