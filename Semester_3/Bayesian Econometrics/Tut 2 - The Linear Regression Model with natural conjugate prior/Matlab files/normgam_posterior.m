function [params,moments] = normgam_posterior(y, X, beta0, kappa0, ssq0, nu0)
% -------------------------------------------------------------------------
% [Ref textbook 3.9.]
% 
% This function calculates the normal-gamma posterior of the multiple
% regression model with matrix form given by:
%
% y (house prices) = X(features)*beta(priors) + e,   e ~ N(0,1/h)
%
% where the natural conjugate normal-gamma prior follows
%
% beta|h ~ N(beta0,1/(h*kappa0))
% h      ~ G(1/s0^2,nu0)
%
% The posterior distribution is normal-gamma:
%
% p(beta,h|y) ~ NG(beta1, 1/kappa1, 1/s1^2, nu1) --> describe the shape
% -------------------------------------------------------------------------
% Input: --> 0 --> prior
% y      (N x 1) vector of the dependent variable
% X      (N x k) matrix of exogenous variables
% beta0  (k x 1) vector of prior means for beta
% kappa0 (k x k) matrix of prior precision for beta
% ssq0   inverse of the scalar prior mean for h: s0^2
% nu0    scalar prior degrees of freedom for h
% -------------------------------------------------------------------------
% Output (structures): --> 1 is posterior
% params.beta1    (k x 1) vector: posterior beta
% params.kappa1   (k x k) matrix: posterior kappa
% params.ssq1     scalar: posterior s^2 
% params.nu1      scalar: posterior nu
% moments.mbeta   (k x 1) vector of posterior means for beta
% moments.vbeta   (k x k) matrix of posterior variance matrix for beta
% moments.mh      scalar posterior mean for h
% moments.vh      scalar posterior variance for h
% -------------------------------------------------------------------------

% Get data-based statistics (OLS quantities): // Ref.p.37
k       = size(X, 2);                   % number of regressors
N       = size(y, 1);                   % sample size
nu      = N-k;                          % OLS degrees of freedom 
betahat = (X'*X)\(X'*y);                % OLS estimator - slower alternative: inv(X'*X)*X'*y
ssq     = (1/nu)*(y-X*betahat)'*(y-X*betahat);   % OLS estimator of the error variance
kappa   = X'*X;                        % OLS precision

% Posterior parameters of the normal-gamma distribution:
params.beta1  = (kappa0+kappa)\(kappa0*beta0 + kappa*betahat); % slower alternative: inv(kappa0+kappa)*(kappa0*beta0 + kappa*betahat), since inverse backslash else we use normal "/"
params.kappa1 = kappa0+kappa;
H = inv(eye(k) + kappa0*inv(kappa))*kappa0; %This also works for kappa0 -> 0, eye(k) = identity matrix with dimension k
% This is more efficient than compution: inv(inv(kappa0)+inv(kappa))
params.ssq1   = (nu0*ssq0 + nu*ssq + (betahat-beta0)'*H*(betahat-beta0)) / (nu0+N); %posterior variance p.41
params.nu1    = nu0+N;


% Getting moments p.42

% Posterior mean and variance of beta 
moments.mbeta = params.beta1;
moments.vbeta = (params.nu1/(params.nu1-2)) * params.ssq1*inv(params.kappa1); % or slow alternative: (params.nu1/(params.nu1-2)) * params.ssq1*inv(params.kappa1)

% Posterior mean and variance of h
moments.mh = 1/params.ssq1;
moments.vh = 2/(params.nu1*params.ssq1^2); % from results of posterior dist. / variance of error precision
