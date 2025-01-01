%=====================| Tutorial III: Gibbs Sampling |=====================
%======================| Computer-Based Exercise 3 |=======================
% This script replicates the house price example in Chapter 4 (p.73-77)
% using the independent Normal-Gamma prior.

% -------------------------------------------------------------------------
% Application of the Gibbs sampler to simulate 
% beta|h ~ N(beta1,V1)
% h|beta ~ Gamma(1/ssq1,nu1)

% Remark (*) - drawing Gamma random numbers in MATLAB:
% Note that MATLAB's definition of the Gamma pdf differs from the one 
% presented in Koop's textbook!
% Textbook's Gamma(mu,nu) is Matlab's Gamma(nu/2,2*mu/nu).
% -------------------------------------------------------------------------

%% The Gibbs sampler

% Clear workspace memory and command window:
clear
clc

% Load data and set regression parameters/variables: 
load HPRICE.txt;
N = size(HPRICE,1);  % Sample size
y = HPRICE(:,1);     % Dependent variable
X = HPRICE(:,2:5);   % Explanatory variables
X = [ones(N,1) X];   % Add intercept
k = size(X,2);       % Number of regressors

% Set informative prior:
beta0    = [0, 10, 5000, 10000, 10000]';             % Prior mean for beta
ssq0     = 5000^2;                                   % Inverse of the prior mean for h
nu0      = 5;                                        % Prior degrees of freedom for h
V0       = diag([10000^2 5^2 2500^2 5000^2 5000^2]); % Prior variance for beta 
kappa0   = inv(V0);                                  % Prior precision for beta

% OLS estimates:
XX       = "...";    % Compute X'X
Xy       = "...";    % Compute X'y
beta_OLS = "...";    % OLS estimates
ssq_OLS  = "...";    % Variance of error term
V_OLS    = "...";    % Variance-covariance matrix of OLS estimates

% Set simulation parameters:
S0 = "...";          % Burn-in replications
S1 = "...";          % Monte Carlo posterior replications 

% Initialize matrices to store posterior draws:
beta = "...";        % Matrix of beta posterior draws
h    = "...";        % Vector of h posterior draws

% Initial draw beta(0) using OLS estimates:
beta_initial = "...";

% First draw h(1)|beta_initial:
nu1  = "...";        % Posterior degrees of freedom for h|beta
ssq1 = "...";        % Inverse of the posterior mean for h|beta
h(1) = "...";        % Draw from the conditional Gamma for h|beta

% First draw beta(1)|h(1):
V1        = "...";   % Posterior variance for beta|h
V1        = "...";   % Enforce symmetry of V1
beta1     = "...";   % Posterior mean for beta|h
beta(1,:) = "...";   % Draw from the conditional Normal for beta|h

% Start the Gibbs sampler:
for i=2:S0+S1    
    % Draw h(i)|beta(i-1):
    ssq1 = "...";
    h(i) = "...";   
    % Draw beta(i)|h(i)
    V1        = "...";
    V1        = "...";
    beta1     = "..."; 
    beta(i,:) = "...";   
end

% Save burn-in draws S0 separately and then discard them in beta and h:
beta_S0 = "...";
h_S0    = "...";
beta    = "...";
h       = "...";

% Compute posterior point estimates using the MC posterior draws S1:
beta_mean = "...";   % Posterior mean of beta
beta_std  = "...";   % Posterior std. deviation of beta
beta_prc  = "...";   % Posterior 2.5/97.5 percentiles of beta
h_mean    = "...";   % Posterior mean of h
h_std     = "...";   % Posterior std. deviation of h
h_prc     = "...";   % Posterior 2.5/97.5 percentiles of h

% Trace plots:
names_params = {'$\beta_1$','$\beta_2$','$\beta_3$','$\beta_4$','$\beta_5$'};
figure
sgtitle('Markov Chain Trace Plots','Interpreter','Latex','FontSize',18)
for i = 1:k
    subplot(2,3,i)
    plot([beta_S0(:,i);beta(:,i)])
    title(names_params{i},'Interpreter','Latex','FontSize',16)
end
subplot(2,3,6)
plot([h_S0(:);h(:)])
title('h','Interpreter','Latex','FontSize',16)

%% Numerical standard errors (NSE) 

% NSE based on Newey-West long-run variances:
beta_NSE = "...";             % Initialize row-vector to store NSEs for all beta's 
for j = 1:k
    beta_NSE(j) = "...";      % Compute NSEs for beta
end
h_NSE = "...";                % Compute NSEs for h

%% CD test for convergence of the Markov chain

% Divide the Markov chain into different subsamples: 
A = "..."; B = "..."; C = "...";

% Compute NSE of posterior mean estimates in subsamples A and C: 
betaA_NSE = "...";     % Initialize row-vector to store NSEs for subsample A
betaC_NSE = "...";     % Initialize row-vector to store NSEs for subsample C
for j = 1:k
    betaA_NSE(j) = "...";  % Compute NSEs for beta in subsample A
    betaC_NSE(j) = "...";  % Compute NSEs for beta in subsample C
end
hA_NSE = "...";            % Compute NSEs for h in subsample A
hC_NSE = "...";            % Compute NSEs for h in subsample C

% Compute CD convergence statistic:
CD_beta = "..."; 
CD_h    = "..."; 


%%%%%%%%%% ---------- DISPLAY RESULTS ---------- %%%%%%%%%%
% display Table 4.1 (p. 75)
fprintf(1,'\n\n')
fprintf(1,'--------------------------------------------------------------------------\n')
fprintf(1,'        Prior and posterior results for beta (std. dev. in brackets) \n\n')
disp('        Prior                      Posterior             NSE         CD')
disp(' --------------------        ---------------------     -------    -------')
fprintf(1,'%9.2f (%9.2f)       %9.2f  (%9.2f)  %9.2f  %9.2f  \n',[beta0 sqrt(diag(V0)) beta_mean' beta_std' beta_NSE' CD_beta']');
fprintf(1,'---------------------------------------------------------------------------\n')

%% (e) (*) Prediction

X_star = "...";        % Information on X^* to predict y^*

% Monte Carlo integration to predict y^*:
y_star = "...";        % Initialize vector to store y^* predictions 
for i=1:S1
    y_star(i) = "..."; % Draw from the predictive density
end

% Point estimates:
y_star_mean = "...";   % Mean of predictive density
y_star_sd   = "...";   % Standard deviation of predictive density
y_star_prc  = "...";   % Posterior 2.5/97.5 percentiles of beta

% Histogram of the predictive density:
figure
histfit(y_star,100)
hold on
line([y_star_mean, y_star_mean], ylim, 'LineWidth', 2, 'Color', 'r');        % Add predictive mean
line([y_star_prc(1), y_star_prc(1)], ylim, 'LineWidth', 0.5, 'Color', 'r');  % Add predictive 2.5 percentile
line([y_star_prc(2), y_star_prc(2)], ylim, 'LineWidth', 0.5, 'Color', 'r');  % Add predictive 97.5 percentile
title('Predictive density for $y^*$ and fitted normal pdf','Interpreter','Latex','FontSize',16)