%======================| Tutorial IV: MH algorithm |=======================
%======================| Computer-Based Exercise 1 |=======================
% -------------------------------------------------------------------------
% This script replicates the nonlinear CES regression in Chapter 5 (p.107)
% using the independence chain MH algorithm with candidate draws coming 
% from a multivariate normal distribution with mean = OLS estimates and 
% variance = asymptotic variance times a scale factor (to account for 
% potential problems of underdispersion). In addition, it also estimates
% the error precision h by brute force MH and assuming mutual independence
% between gamma and h. 
% (ii) Candidate draws for h are taken from a conditional Gamma density
% with parameters nu=N and mu = N/[y-f(X,\gamma^*]'[y-f(X,\gamma^*)], hence
% this proposal draw is given gamma^*.
% -------------------------------------------------------------------------

clear
clc

% Load data and set regression parameters/variables: 
load ch5data.out;
y   = ch5data(:,1);         % Dependent variable
lab = ch5data(:,2);         % Labor input
cap = ch5data(:,3);         % Capital input
N   = size(y,1);            % Sample size
X   = [ones(N,1) lab cap];  % Regressors
k   = size(X,2)+1;          % Number of regressors

% Set simulation parameters:
S0 = 10000;                 % Burn-in replications
S1 = 100000;                % MC posterior replications
c  = 2;                     % Scale factor for the proposal variance

% Nonlinear conditional mean function (CES):
f = @(x,g) g(1)*x(:,1)+(g(2)*x(:,2).^g(4)+g(3)*x(:,3).^g(4)).^(1/g(4));

% Initial estimate for mean and variance of gamma: 
gammahat = X\y; gammahat(4) = 1;    % OLS estimates
e    = y - f(X,gammahat);           % OLS residuals     
sig2 = (e'*e)/(N-k);                % Residual variance
uu_hat = e'*e;                      % Residual squares

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
gamma = zeros(k,S0+S1);  % Pre-allocate matrix of posterior estimates for gamma
h     = zeros(1,S0+S1);  % Pre-allocate vector of posterior estimates for h
alpha = zeros(1,S0+S1);  % Pre-allocate vector of acceptance probabilities 
acc   = zeros(1,S0+S1);  % Pre-allocate vector of acceptance decisions

% Starting value (use OLS point estimates as first draw):
gamma(:,1) = gammahat;
h(1) = N/uu_hat;

% Fix seed for replication purposes:
rng(123)

% Start iteration:
for s = 2:(S0+S1)
    
    % Proposal draw from the multivariate normal distribution with 
    % mean gammahat and variance Sigma:
    z_star = randn(4,1);                    % Draw from N(0,I)
    g_star = cholSigma'*z_star + gammahat;  % Transform to a draw from N(gammahat,Sigma)
    
    % Update residuals conditional on g_star and gamma(s-1):
    u_star = y - f(X,g_star);
    u_last = y - f(X,gamma(:,s-1));  

    % Proposal draw from the Gamma distribution for h:
    h_star = gamrnd(N/2,2/(u_star'*u_star)); % Proposal draw from conditional Gamma

    % Ratios of p/q for (g_star,h_star) and (gamma(s-1),h(s-1)):
    pqstar_ratio = exp(0.5*z_star'*z_star) / (u_star'*u_star)^(N/2);
    pqlast_ratio = exp(0.5*(gamma(:,s-1)-gammahat)'*iSigma*...
        (gamma(:,s-1)-gammahat)) / (u_last'*u_last)^(N/2);

    % Acceptance probability:
    alpha(s) = min(pqstar_ratio/pqlast_ratio,1);
    
    % Acceptance-rejection decision:
    u = rand(1);
    acc(1,s) = u<=alpha(s);
    gamma(:,s) = acc(1,s)*g_star + (1-acc(1,s))*gamma(:,s-1);
    h(:,s) = acc(1,s)*h_star + (1-acc(1,s))*h(:,s-1);
end

% Discard burn-in replications S0:
gamma = gamma(:,S0+1:end);
h = h(:,S0+1:end);

% Display results:
disp('---------------------------------------------------------------')
disp('Results of the independence chain Metropolis-Hastings algorithm')
disp('---------------------------------------------------------------')
disp('       gamma1      gamma2      gamma3      gamma4      h       ')
disp('---------------------------------------------------------------')
disp(['Means: ' num2str([mean(gamma') mean(h')])]) 
disp(['s.e.:  ' num2str([std(gamma') std(h')])])
disp('---------------------------------------------------------------')
disp(['Acceptance rate = ' num2str(mean(acc))])
disp('---------------------------------------------------------------')
disp(['Elasticity of substitution = ' num2str(1/(1-mean(gamma(end,:))))])
disp('---------------------------------------------------------------')