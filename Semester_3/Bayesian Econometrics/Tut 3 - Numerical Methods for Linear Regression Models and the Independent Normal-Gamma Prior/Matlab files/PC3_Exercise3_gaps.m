%=====================| Tutorial III: Gibbs Sampling |=====================
%======================| Computer-Based Exercise 3 |=======================
% This script replicates the house price example in Chapter 4 (p.73-77)
% using the independent Normal-Gamma prior.

% -------------------------------------------------------------------------
% Application of the Gibbs sampler to simulate 
% beta|h ~ N(beta1,V1) --> posterior is marked with 1
% h|beta ~ Gamma(1/ssq1,nu1)

% Remark (*) - drawing Gamma random numbers in MATLAB:
% Note that MATLAB's definition of the Gamma pdf differs from the one 
% presented in Koop's textbook!
% Textbook's Gamma(mu,nu) is Matlab's Gamma(nu/2, 2*mu/nu). --> IMPORTANT
% TO KEEP IN MIND!
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
X = [ones(N,1) X];   % Add intercept - first column full of ones
k = size(X,2);       % Number of regressors

% Set informative prior:
beta0    = [0, 10, 5000, 10000, 10000]';             % Prior mean for beta
ssq0     = 5000^2;                                   % Inverse of the prior mean for h
nu0      = 5;                                        % Prior degrees of freedom for h
V0       = diag([10000^2 5^2 2500^2 5000^2 5000^2]); % Prior variance for beta 
kappa0   = inv(V0);                                  % Prior precision for beta

% OLS estimates:
XX       = X'*X;    % Compute X'X
Xy       = X'*y;    % Compute X'y
beta_OLS = X\y;    % OLS estimates
ssq_OLS  = ((y-X*beta_OLS)'*(y-X*beta_OLS))/(N-k) ;    % Variance of error term
V_OLS    = ssq_OLS*inv(XX);    % Variance-covariance matrix of OLS estimates

% Set simulation parameters:
S0 = 1000;          % Burn-in replications
S1 = 10000;          % Monte Carlo posterior replications 

% Initialize matrices to store posterior draws:
beta = zeros(S0+S1, k);        % Matrix of beta posterior draws
h    = zeros(S0+S1, 1);        % Vector of h posterior draws - due to homoscadicity assumption

% Initial draw beta(0) using OLS estimates:
beta_initial = mvnrnd(beta_OLS, V_OLS, 1); %takes random draw from OLS --> since conditional is multnormal

% First draw h(1)|beta_initial: (solutions from lecture slides)
nu1  = N + nu0;        % Posterior degrees of freedom for h|beta
ssq1 = ((y-X*beta_initial')'*(y-X*beta_initial')+nu0*ssq0)/nu1;        % Inverse of the posterior mean for h|beta
h(1) = gamrnd(nu1/2, 2/(nu1*ssq1));        % Draw from the conditional Gamma for h|beta --> DO Calculations!!

% First draw beta(1)|h(1):
V1        = inv(kappa0 + h(1)*XX);   % Posterior variance for beta|h
V1        = (V1 + V1')/2;   % Enforce symmetry of V1 - either close to 0 or too large
beta1     = V1*(kappa0*beta0 + h(1)*Xy);   % Posterior mean for beta|h
beta(1,:) = mvnrnd(beta1,V1,1);   % Draw from the conditional Normal for beta|h

% Start the Gibbs sampler:
for i=2:S0+S1    
    % Draw h(i)|beta(i-1):
    ssq1 = ((y-X*beta(i-1,:)')'*(y-X*beta(i-1,:)')+nu0*ssq0)/nu1;
    h(i) = gamrnd(nu1/2, 2/(nu1*ssq1));
    % Draw beta(i)|h(i)
    V1        = inv(kappa0 + h(i)*XX);
    V1        = (V1 + V1')/2;
    beta1     = V1*(kappa0*beta0 + h(i)*Xy); 
    beta(i,:) = mvnrnd(beta1, V1, 1);   
end

% Save burn-in draws S0 separately and then discard them in beta and h:
beta_S0 = beta(1:S0,:); %keep to analyze info of convergence of MC-Chain
h_S0    = h(1:S0,:);
beta    = beta(S0+1:end,:);
h       = h(S0+1:end,:);

% Compute posterior point estimates using the MC posterior draws S1:
beta_mean = mean(beta);   % Posterior mean of beta
beta_std  = std(beta);   % Posterior std. deviation of beta
beta_prc  = prctile(beta, [2.5, 97.5]);   % Posterior 2.5/97.5 percentiles of beta
h_mean    = mean(h);   % Posterior mean of h
h_std     = std(h);   % Posterior std. deviation of h
h_prc     = prctile(h, [2.5, 97.5]);   % Posterior 2.5/97.5 percentiles of h

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
beta_NSE = zeros(1,k);             % Initialize row-vector to store NSEs for all beta's 
for j = 1:k
    beta_NSE(j) = sqrt(NeweyWest(beta(:,j),10))/sqrt(S1);      % Compute NSEs for beta
end
h_NSE = sqrt(NeweyWest(h,10))/sqrt(S1);                % Compute NSEs for h

%% CD test for convergence of the Markov chain

% Divide the Markov chain into different subsamples: 
A = round(0.1*S1); B = round(0.5*S1); C = round(0.4*S1);

% Compute NSE of posterior mean estimates in subsamples A and C: 
betaA_NSE = zeros(1,k);     % Initialize row-vector to store NSEs for subsample A
betaC_NSE = zeros(1,k);     % Initialize row-vector to store NSEs for subsample C
for j = 1:k
    betaA_NSE(j) = sqrt(NeweyWest(beta(1:A, j),10))/sqrt(A);  % Compute NSEs for beta in subsample A
    betaC_NSE(j) = sqrt(NeweyWest(beta(1+A+B:A+B+C, j),10))/sqrt(C);  % Compute NSEs for beta in subsample C
end
hA_NSE = sqrt(NeweyWest(h(1:A),10))/sqrt(A);            % Compute NSEs for h in subsample A
hC_NSE = sqrt(NeweyWest(h(1+A+B:A+B+C),10))/sqrt(C);            % Compute NSEs for h in subsample C

% Compute CD convergence statistic:
CD_beta = (mean(beta(1:A,:))-mean(beta(1+A+B:A+B+C,:)))./(betaA_NSE + betaC_NSE); 
CD_h    = (mean(h(1:A))-mean(h(1+A+B:A+B+C)))./(hA_NSE + hC_NSE); 


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

X_star = [1 5000 2 2 1];        % Information on X^* to predict y^*

% Monte Carlo integration to predict y^*:
y_star = zeros(S1,1);        % Initialize vector to store y^* predictions 
for i=1:S1
    y_star(i) = normrnd(X_star*beta(i,:)',1/sqrt(h(i))); % Draw from the predictive density
end

% Point estimates:
y_star_mean = mean(y_star);   % Mean of predictive density
y_star_sd   = std(y_star);   % Standard deviation of predictive density
y_star_prc  = prctile(y_star, [2.5 97.5]);   % Posterior 2.5/97.5 percentiles of beta

% Histogram of the predictive density:
figure
histfit(y_star,100)
hold on
line([y_star_mean, y_star_mean], ylim, 'LineWidth', 2, 'Color', 'r');        % Add predictive mean
line([y_star_prc(1), y_star_prc(1)], ylim, 'LineWidth', 0.5, 'Color', 'r');  % Add predictive 2.5 percentile
line([y_star_prc(2), y_star_prc(2)], ylim, 'LineWidth', 0.5, 'Color', 'r');  % Add predictive 97.5 percentile
title('Predictive density for $y^*$ and fitted normal pdf','Interpreter','Latex','FontSize',16)