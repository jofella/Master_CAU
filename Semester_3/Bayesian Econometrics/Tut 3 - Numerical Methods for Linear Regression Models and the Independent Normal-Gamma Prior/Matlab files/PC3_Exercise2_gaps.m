%===================| Tutorial III: Numerical Methods |====================
%======================| Computer-Based Exercise 2 |=======================
% This script applies numerical methods to simulate posterior moments of a 
% logistic distribution using uniform and standard normal random numbers.

% Clear workspace memory and command window:
clear
clc

% Set parameters of the posterior logistic distribution: (alpha upper and beta upper; given)
alpha = 4;
beta  = 2;

%% ===| (a) Probability integral transform

S         = 100000;                     % number of replications
u         = rand(S,1);                  % draw uniform random numbers
theta     = alpha - beta*log(1./u-1);   % compute theta^(S) by the inverse cdf transformation
theta_bar = mean(theta);                % simulated posterior mean
NSE_mean  = std(theta)/sqrt(S);         % NSE for posterior mean
S2_theta  = var(theta);                 % simulated posterior variance
NSE_var   = std((theta-theta_bar).^2)/sqrt(S);   % NSE for posterior variance
g_theta   = mean(exp(sqrt(abs(theta))-1));   % simulated posterior value g(theta)
NSE_g     = std(exp(sqrt(abs(theta))-1))/sqrt(S);    % NSE for g(theta)

% Display results:
fprintf('\n Probability integral transform with %d replications \n\n',S)
fprintf(' Mean:     theoretical = %8.4f, simulated = %7.4f, NSE = %7.4f \n',[alpha, theta_bar, NSE_mean])
fprintf(' Var:      theoretical = %8.4f, simulated = %7.4f, NSE = %7.4f \n',[beta^2*pi^2/3, S2_theta, NSE_var])
fprintf(' g(theta): theoretical = %8.4f, simulated = %7.4f, NSE = %7.4f \n\n',[NaN, g_theta, NSE_g])

%% ===| (b) Acceptance-rejection method

S       = 100000;                   % number of replications
k       = 7;                        % degrees of freedom for the proposal t-distribution (get fatter tail)
mu      = alpha;                    % mean for the proposal t-distribution (=alpha)
sigma   = sqrt((k-2)/k)*beta*pi/sqrt(3);  % std. dev. for the proposal t-distribution (=std. dev. of the logistic)
y       = mu + sigma.*trnd(k, [S, 1]); % draw t-distributed random numbers (draws from standard t)

p       = (1/beta)*exp((y-alpha)/beta)./(1+exp((y-alpha)/beta)).^2; % compute target pdf (logistic) for each y
c       = 1/(sqrt(pi)*k^(-k/2)*gamma(k/2)/gamma((k+1)/2));  % compute integrating constant for the proposal pdf
q       = c*(1/sigma)*(k+(y-mu).^2/sigma^2).^(-(k+1)/2);    % compute proposal pdf (t-distribution) for each y
M       = max(p./q);    % compute factor M

u       = rand(S, 1);       % draw uniform random numbers
theta   = y(u<(p./(M*q)));  % acceptance-rejection decision
S_tilde = length(theta);    % number of accepted draws

theta_bar = mean(theta);                % simulated posterior mean
NSE_mean  = std(theta)/sqrt(S_tilde);         % NSE for posterior mean
S2_theta  = var(theta);                 % simulated posterior variance
NSE_var   = std((theta-theta_bar).^2)/sqrt(S_tilde);   % NSE for posterior variance
g_theta   = mean(exp(sqrt(abs(theta))-1));   % simulated posterior value g(theta)
NSE_g     = std(exp(sqrt(abs(theta))-1))/sqrt(S_tilde);    % NSE for g(theta

% Display results:
fprintf('\n Acceptance-rejection method with %d replications \n\n',S)
fprintf(' Acceptance Rate: theoretical = %8.4f, simulated = %7.4f \n',[1/M, S_tilde/S])
fprintf(' Mean:            theoretical = %8.4f, simulated = %7.4f, NSE = %7.4f \n',[alpha, theta_bar, NSE_mean])
fprintf(' Var:             theoretical = %8.4f, simulated = %7.4f, NSE = %7.4f \n',[beta^2*pi^2/3, S2_theta, NSE_var])
fprintf(' g(theta):        theoretical = %8.4f, simulated = %7.4f, NSE = %7.4f \n\n',[NaN, g_theta, NSE_g])

% Figure: target and proposal pdfs
figure
area(y,p,'Facecolor','#0072BD','edgecolor','none','FaceAlpha',0.75)
xlim([-20 30])
hold on
area(y,M*q,'Facecolor','#A2142F','edgecolor','none','FaceAlpha',0.25)
grid on
set(gca,'FontSize',16)
title('Logistic target pdf $p(\theta)$ (blue) and proposal pdf $q(\theta)$ (red)','interpreter', 'latex')

%% ===| (c) Importance sampling - using the t random numbers y generated in (b)  

w              = p./q;      % importance sampling weights
g              = exp(sqrt(abs(theta))-1);   % MC function of interest g(theta) for the mean
theta_bar      = sum(w.*y)/sum(w);          % simulated posterior mean
theta_bar_Avar = "...";          % asymptotic variance of the mean estimate
NSE_mean       = "...";          % NSE for posterior mean

g              = "...";          % MC function of interest g(theta) for the variance
S2_theta       = "...";          % simulated posterior variance
S2_theta_Avar  = "...";          % asymptotic variance of the variance estimate
NSE_S2         = "...";          % NSE for posterior variance

g              = "...";          % MC function of interest g(theta) in (iii)
g_theta        = "...";          % simulated posterior mean of g(theta)
g_theta_Avar   = "...";          % asymptotic variance of the g(theta) estimate
NSE_g          = "...";          % NSE for posterior function of interest g(theta)

% Display results:
fprintf('\n Importance sampling with %d replications \n\n',S)
fprintf(' Average Weight        = %8.4f, Std. Dev. = %7.4f \n',[mean(w), std(w)])
fprintf(' Mean:     theoretical = %8.4f, simulated = %7.4f, NSE = %7.4f \n',[alpha, theta_bar, NSE_mean])
fprintf(' Var:      theoretical = %8.4f, simulated = %7.4f, NSE = %7.4f \n',[beta^2*pi^2/3, S2_theta, NSE_S2])
fprintf(' g(theta): theoretical = %8.4f, simulated = %7.4f, NSE = %7.4f \n\n',[NaN, g_theta, NSE_g])