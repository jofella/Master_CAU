%======| Tutorial II: Linear Regression with Natural Conjugate Prior |=====
%======================| Computer-Based Exercise 1 |=======================
% This script replicates the house price example in Chapter 3 (p.47-52),
% also discussed in the lecture slides.

% Clear workspace memory and command window:
clear
clc

% Load data and set regression parameters/variables: 
load HPRICE.txt;
N = size(HPRICE,1);
y = HPRICE(:,1);
X = HPRICE(:,2:5); %pre-defined by exercise (data set is bigger)
X = [ones(N,1) X];
k = size(X,2);

% Set informative prior and compute posterior parameters/moments: --
% Ref.textbook --> e.g. educated guesses (asking (local) real estate
% agents) - p.48 textbook
beta0    = [0, 10, 5000, 10000, 10000]';             % prior mean for beta --> add 1 unit increases price by ...-units
ssq0     = 5000^2;                                   % inverse of the prior mean for h
nu0      = 5;                                        % prior degrees of freedom for h
Varbeta0 = [10000^2, 5^2, 2500^2, 5000^2, 5000^2]';  % prior variance for beta (non-normalized)
V0       = diag(Varbeta0)*(nu0-2)/(nu0*ssq0);        % prior variance for beta
kappa0   = inv(V0);                                  % prior precision for beta
[p_inf, m_inf] = normgam_posterior(y,X,beta0,kappa0,ssq0,nu0); %copy paste form function

% Set noninformative prior and compute posterior parameters/moments:
beta0_non  = zeros(5,1);
ssq0_non   = 1; %get fully non-inf prior to get infinitly high variance --> set parameters 
nu0_non    = 0;
kappa0_non = zeros(k,k);
[p_non,m_non] = normgam_posterior(y,X,beta0_non,kappa0_non,ssq0_non,nu0);

% Display results: Table 3.1 (p. 51)
fprintf(1,'\n\n')
fprintf(1,'--------------------------------------------------------------------------------------\n')
fprintf(1,'        Prior and posterior means for beta (std.dev. in brackets) \n\n')
disp('        Prior                               Posterior based on')
disp('     Informative             noninformative prior           informative prior')
disp(' --------------------        ---------------------        ---------------------')
fprintf(1,'%9.2f (%9.2f)  --   %9.2f  (%9.2f)  --   %9.2f  (%9.2f) \n',...
    [beta0 sqrt(Varbeta0) m_non.mbeta sqrt(diag(m_non.vbeta)) m_inf.mbeta sqrt(diag(m_inf.vbeta))]');
fprintf(1,'--------------------------------------------------------------------------------------\n')

% Display results: Table 3.2 (p. 51)
fprintf(1,'\n\n')
fprintf(1,'--------------------------------------------------------------------------------------\n')
fprintf(1,'       Prior and posterior means for h (std.dev. in brackets \n\n')
disp('        Prior                               Posterior basend on')
disp('     Informative             noninformative prior           informative prior')
disp(' --------------------        ---------------------        ---------------------')
fprintf(1,'%9.2e (%9.2e)  --   %9.2e  (%9.2e)  --   %9.2e  (%9.2e) \n',...
    [1/ssq0 sqrt(2/(ssq0^2*nu0)) m_non.mh sqrt(m_non.vh) m_inf.mh sqrt(m_inf.vh)]');
fprintf(1,'--------------------------------------------------------------------------------------\n')

%% Credible Intervals (HPDI Intervals) and Posterior Odds Ratios - Infrence/Hypothesis testing in frequentist


% Compute 95% credible intervals:
intv_inf95 = t_interval(p_inf.beta1, p_inf.ssq1*inv(p_inf.kappa1), p_inf.nu1, 0.05);
intv_non95 = t_interval(p_non.beta1, p_non.ssq1*inv(p_non.kappa1), p_non.nu1, 0.05);

% Compute 99% credible intervals:
intv_inf99 = t_interval(p_inf.beta1, p_inf.ssq1*inv(p_inf.kappa1), p_inf.nu1, 0.01);
intv_non99 = t_interval(p_non.beta1, p_non.ssq1*inv(p_non.kappa1), p_non.nu1, 0.01);

% Compute probabilities p(beta_i>0|y):
probs_inf = t_nonzeroprob(p_inf.beta1, p_inf.ssq1*inv(p_inf.kappa1), p_inf.nu1);
probs_non = t_nonzeroprob(p_non.beta1, p_non.ssq1*inv(p_non.kappa1), p_non.nu1);

% Compute log of marginal likelihood (to circumvent infty's):
ic = gammaln(p_inf.nu1/2) + (nu0/2)*log(nu0*ssq0) + (-N/2)*log(pi) - gammaln(nu0/2);
marglik_inf = ic - 0.5*log(det(p_inf.kappa1)) - 0.5*log(det(V0)) - (p_inf.nu1/2)*log(p_inf.nu1*p_inf.ssq1);

% [Reference p. 3.6.2] - Equality restrictions

% Compute posterior odds for restricting one beta at a time to zero
PO = zeros(5,1);
for i = 1:5
    subset = setdiff(1:5,i);
    Xr = X(:,subset);
    beta0r = beta0(subset);
    ssq0r = ssq0;
    nu0r = nu0;
    Varbeta0r = Varbeta0(subset);
    V0r = diag(Varbeta0r)*(nu0r-2)/(nu0r*ssq0r);
    kappa0r = inv(V0r);
    [p_r,m_r] = normgam_posterior(y,Xr,beta0r,kappa0r,ssq0r,nu0r);
    icr = gammaln(p_r.nu1/2) + (nu0r/2)*log(nu0r*ssq0r) + (-N/2)*log(pi) - gammaln(nu0r/2);
    marglik_r = icr - 0.5*log(det(p_r.kappa1)) - 0.5*log(det(V0r)) - (p_r.nu1/2)*log(p_r.nu1*p_r.ssq1);
    PO(i,1) = exp(marglik_r-marglik_inf);
end

% Display results: Table 3.3 (p. 52)
fprintf(1,'\n\n')
fprintf(1,'--------------------------------------------------------------------------------------\n')
fprintf(1,'                                 Informative prior \n\n')
disp(' p(beta_j>0|y)           95% HPDI                  99% HPDI          PO for beta_j=0   ')
disp(' -------------    ----------------------    ----------------------   ---------------')
fprintf(1,'%9.2f         [%9.2f, %9.2f]    [%9.2f, %9.2f]      %9.2e \n',...
    [probs_inf(:,2) intv_inf95 intv_inf99 PO]');
fprintf(1,'--------------------------------------------------------------------------------------\n\n')
fprintf(1,'                               Noninformative prior \n\n')
disp(' p(beta_j>0|y)           95% HPDI                  99% HPDI          PO for beta_j=0   ')
disp(' -------------    ----------------------    ----------------------   ---------------')
fprintf(1,'%9.2f         [%9.2f, %9.2f]    [%9.2f, %9.2f]      %6.2e \n',...
    [probs_non(:,2) intv_non95 intv_non99 ones(5,1)*NaN]');
fprintf(1,'--------------------------------------------------------------------------------------\n')