%======================| Tutorial IV: MH algorithm |=======================
%======================| Computer-Based Exercise 1 |=======================
% -------------------------------------------------------------------------
% This script replicates the nonlinear CES regression in Chapter 5 (p.107)
% using the independence chain MH algorithm with candidate draws coming 
% from a multivariate normal distribution with mean = OLS estimates and 
% variance = asymptotic variance times a scale factor (to account for 
% potential problems of underdispersion).
% -------------------------------------------------------------------------

clear
clc

% Load data and set regression parameters/variables: 
load ch5data.out;
y   = "...";         % Dependent variable
lab = "...";         % Labor input
cap = "...";         % Capital input
N   = "...";         % Sample size
X   = "...";         % Regressors
k   = "...";         % Number of regressors

% Set simulation parameters:
S0 = "...";          % Burn-in replications
S1 = "...";          % MC posterior replications
c  = "...";          % Scale factor for the proposal variance

% Nonlinear conditional mean function (CES):
f = "...";

% Initial estimate for mean and variance of gamma: 
gammahat = "..."; gammahat(4) = "...";  % OLS estimates
e    = "...";                           % OLS residuals     
sig2 = "...";                           % Residual variance

% Scores of the NLS estimator evaluated at the OLS estimator:
score(:,1) = ones(N,1);
score(:,2) = lab;
score(:,3) = cap;
score(:,4) = gammahat(2)*lab.*log(lab) + gammahat(3)*cap.*log(cap)...
    - (gammahat(2)*lab + gammahat(3)*cap).*log(gammahat(2)*lab + gammahat(3)*cap);

% Variance matrix based on asymptotic normal distribution of OLS estimator:
V         = sig2*eye(k)/(score'*score/N); % Asymptotic variance matrix
Sigma     = c*V/N;                        % Avar(gammahat) scaled by c
cholSigma = chol(Sigma);                  % Choleski factor 
iSigma    = inv(Sigma);                   % Inverse of variance matrix

% Initialize storage matrices to speed up computation:
gamma = "...";  % Pre-allocate matrix of posterior estimates
alpha = "...";  % Pre-allocate vector of acceptance probabilities 
acc   = "...";  % Pre-allocate vector of acceptance decisions

% Starting value (use OLS point estimates as first draw):
gamma(:,1) = "...";

% Fix seed for replication purposes:
rng("...")

% Start iteration:
for s = 2:(S0+S1)
    
    % Proposal draw from the multivariate normal distribution with 
    % mean gammahat and variance Sigma:
    z_star = "...";         % Draw from N(0,I)
    g_star = "...";         % Transform to a draw from N(gammahat,Sigma)
    
    % Ratio of posterior kernels:
    u_star = "...";
    u_last = "...";
    p_ratio = "...";
        
    % Ratio of proposal kernels:
    % To speed up execution, use the fact that: (g_star-gammahat)'*inv(Sigma)*(g_star-gammahat) = z_star'*z_star
    q_ratio = "...";

    % Acceptance probability:
    alpha(s) = "...";
    
    % Acceptance-rejection decision:
    u = "...";
    acc(1,s) = "...";
    gamma(:,s) = "...";
end

% Save burn-in draws S0 separately and then discard them:
gamma_S0 = "...";          % Save burn-in draws of gamma
gamma    = "...";          % Save S1 posterior draws 

% Convergence diagnostics using MCMC trace plots:
names_params = {'$\gamma_1$','$\gamma_2$','$\gamma_3$','$\gamma_4$'};
figure
tl = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
for i = 1:k
    nexttile
    plot([gamma_S0(i,:) gamma(i,:)])
    title(names_params{i},'Interpreter','Latex','FontSize',16)
end
title(tl,'Markov Chain Trace Plots','Interpreter','Latex','FontSize',18)

% Display results:
simcov = cov(gamma');
disp('---------------------------------------------------------------')
disp('Results of the independence chain Metropolis-Hastings algorithm')
disp('---------------------------------------------------------------')
disp('       gamma1      gamma2      gamma3      gamma4              ')
disp('---------------------------------------------------------------')
disp(['Means: ' num2str(mean(gamma'))]) 
disp(['s.e.:  ' num2str(std(gamma'))])
disp('---------------------------------------------------------------')
disp(['Acceptance rate = ' num2str(mean(acc))])
disp('---------------------------------------------------------------')
disp(['Elasticity of substitution = ' num2str(1/(1-mean(gamma(end,:))))])
disp('---------------------------------------------------------------')