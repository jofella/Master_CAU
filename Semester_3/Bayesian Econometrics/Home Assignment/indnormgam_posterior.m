function [params, moments] = indnormgam_posterior(y, X, beta0, kappa0, ssq1_0, nu1_0, ssq2_0, nu2_0, S0, S1)
% -------------------------------------------------------------------------
% [Ref Textbook - p.64 (Gibbs sampling)]
% 
% This function calculates the posterior of the multiple
% regression model with matrix form given by:
%
% y(unemployment) = X(features) * beta(priors) + e,   e ~ N(0,1/h)
%
%
% beta ~ N(beta0,1/(h*kappa0))  - Coefficients (independent of H)
% h1      ~ G(1/s0^2,nu0)       - Error precision (calm period)
% h2      ~ G(1/s0^2,nu0)       - Error precision (recession)
%
%
% The posterior conditionals:
%
% p(beta|y, h1, h2) ~ N(beta1, 1/kappa1)
% p(h1|y, beta, h2) ~ N(1/s1^2, nu1)
% p(h2|y, beta, h1) ~ N(1/s2^2, nu2)
%
% -------------------------------------------------------------------------
% Input: --> 0 --> prior
% y      (N x 1) vector of the dependent variable
% X      (N x k) matrix of exogenous variables

% -------------------------------------------------------------------------
% Output (structures): --> 1 is posterior
% params.beta1    (k x 1) vector: posterior beta

% -------------------------------------------------------------------------


%%%%-- FOR Reference: _0 = prior, _1 = posterior --%%%%


% ===| Set regression parameters:
% Recall: ...

k   =    size(X,2);      % Number of regressors
N   =    size(y, 1);     % Sample size
T   =    N;
T2  =    sum(diag(d));   % Number of recession periods, diag(d) is re-transorming Matrix into vector
T1  =    T-T2;           % Number of calm periods


% ===| OLS estimates:
% Recall: ...

XX       = X'*X;    % Compute X'X
Xy       = X'*y;    % Compute X'y
beta_OLS = X\y;     % OLS estimates
ssq_OLS  = ((y-X*beta_OLS)'*(y-X*beta_OLS))/(N-k) ;    % Variance of error term
V_OLS    = ssq_OLS*inv(XX);    % Variance-covariance matrix of OLS estimates


% ===| Gibbs Sampler:
% Recall: ...

% Intitialize storage matrices for posterior draws:
beta    = zeros(S0+S1, k);
h1      = zeros(S0+S1, 1);
h2      = zeros(S0+S1, 1);

% Initial draw beta(0) using OLS estimates:
beta_initial = mvnrnd(beta_OLS, V_OLS, 1); % random draw from OLS --> since conditional is multnormal

% First draw h(1)|beta_intial:
nu1_1 = T1 + nu1_0;
ssq1_1 = ((y-X*beta_initial')'*d*(y-X*beta_initial')+nu1_0*ssq1_0)/nu1_1;
h1(1) = gamrnd(nu1_1/2, 2/(nu1_1*ssq1_1));

% First draw h(2)|beta_intial:
nu2_1 = T2 + nu2_0;
ssq2_1 = ((y-X*beta_initial')'*d*(y-X*beta_initial')+nu2_0*ssq2_0)/nu2_1;
h2(1) = gamrnd(nu2_1/2, 2/(nu2_1*ssq2_1));

% First draw beta(1)|h1(1), h1(2):
H         = diag(h1(1) * (1 - diag(d)) + h2(1) * diag(d));
V1        = inv(V_OLS + X'*H*X);        % Posterior variance for beta|h1,h2, inverse cause V^-1
beta1     = V1*(X'*H*y);                % Posterior mean for beta|h1,h2 - "*"V1 cause its a matrix so mult by inverse
beta(1,:) = mvnrnd(beta1',V1,1);        % Draw from the conditional Normal for beta|h


% Gibbs sampler loop:
for i=2:S0+S1
    % Draw h1(i)|beta(i-1):
    ssq1_1 = ((y-X*beta(i-1,:)')'*(y-X*beta(i-1,:)')+nu0*ssq0)/nu1;
    h1(i) = gamrnd(nu1/2, 2/(nu1*ssq1_1));
    
    % Draw h2(i)|beta(i-1):
    ssq2_1 = ((y-X*beta(i-1,:)')'*(y-X*beta(i-1,:)')+nu0*ssq0)/nu1;
    h2(i) = gamrnd(nu1/2, 2/(nu1*ssq2_1));
    
    % Draw beta(i)|h1(i),h2(i)
    H         = diag(h1(i) * (1 - diag(d)) + h2(i) * diag(d));
    V1        = inv(V_OLS + X'*H*X);
    V1        = (V1 + V1')/2;
    beta1     = V1*(X'*H*y); 
    beta(i,:) = mvnrnd(beta1, V1, 1);   
end

% ===| Store posterior parameters: 
% Recall: ...

% Save burn-in draws S0 separately and then discard them in beta and h:

% "Burnt"-values - Analyze convergence MC-Chain
params.beta_S0  = beta(1:S0,:); 
params.h1_S0    = h1(1:S0,:);
params.h2_S0    = h2(1:S0,:);

% Actual sample
params.beta     = beta(S0+1:end,:);
params.h1       = h1(S0+1:end,:);
params.h2       = h2(S0+1:end,:);


% ===| Get moments: 
% Recall: ...

% Posterior mean and variance of beta 
moments.beta_mean = mean(beta);
moments.beta_std = std(beta); 

% Posterior mean and variance of h1
moments.h1_mean = mean(h1);
moments.h1_std = std(h1);

% Posterior mean and variance of h2
moments.h2_mean = mean(h2);
moments.h2_std = std(h2);