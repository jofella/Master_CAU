%===========| Home Assignment Econometrics III (WiSe 2024/25) |============
% ===| Participants:
% Josef Fella 1186179
% 
% -------------------------------------------------------------------------
% ===| Bayesian estimation of the dynamic linear regression model (1):
% UNEMP_t = mu + alpha_1 * UNEMP_t-1 
%          + beta_1 * INPRO_t-1 + ... + beta_q * INPRO_t-q
%          + gamma_1 * CPI_t-1 + ... + gamma_q * CPI_t-q
%          + phi_1 * BCONF_t-1 + ... + phi_q * BCONF_t-q + lambda * COVID_t + e_t, e_t ~ (0,1/h_i)
% 
% where h_i = h_1 if t is a calm period and h_i = h_2 if t is a recession.
%
% In matrix form:
%
% y = X*beta + e,   e ~ N(0,H^{-1}),
%
% with independent Normal-Gamma prior as follows:
%
% beta ~ N(beta0,1/kappa0)
% h1   ~ G(1/s0_1^2,nu0_1)
% h2   ~ G(1/s0_2^2,nu0_2)
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% ===| Data description:
% dates      (T-by-1) vector of monthly dates
% data       (T-by-4) matrix of macroeconomic variables for the US economy (stationary series after transformations)
% recessions (J-by-2) matrix of US recessions with start and finish points, where J is the number of dated recessions in the sample 
% -------------------------------------------------------------------------

% ===| Clear data memory:
clear, clc

% ===| Load the data set:
load US_macro.mat

% ===| Plot the transformed (stationary) data:
figure
tiledlayout('flow','TileSpacing','compact','Padding','compact');
for k = 1:size(data,2)
    nexttile, hold on
    plot(dates,data(:,k),'Color','#0072BD','LineWidth',1.5)
    grid on
    yline(0,'LineWidth',1)
    recessionplot
    axis tight
    set(gca,'FontSize',14)
    set(gca,'TickLabelInterpreter','latex','FontSize',18)
    title(VARnames{k},'Interpreter','latex')
end
set(gcf,'Position',[100 100 1400 1000])

% ===| Set baseline specs and adjust the dataset:
q = 4;                                                                      % Set distributed lags
X = [lagmatrix(data(:,1),1) lagmatrix(data(:,2),1:q)...                     % Construct matrix of covariates
      lagmatrix(data(:,3),1:q) lagmatrix(data(:,4),1:q)];
X = X(q+1:end,:);                                                           % Adjust matrix of covariates
y = data(q+1:end,1);                                                        % Set dependent variable and match sample
T = size(y,1);  % sample size N                                                            % Effective sample size
dates = dates(q+1:end);                                                     % Match sample size of dates
COVID = dates >= datetime(2020,04,01) & dates <= datetime(2020,06,01);      % Define COVID-19 dummy variable
X = [ones(T,1) X COVID];                                                    % Matrix of covariates
k = size(X,2);                                                              % Number of covariates
J = size(recessions,1);                                                     % Number of US recessions in the sample
tRecession = zeros(T,1);                                                    % Pre-allocate dummy-vector for recession periods
for j = 1:J
    tRecession = tRecession | (dates >= recessions(j, 1) & dates <= recessions(j, 2)); % Identify dates within the start and end period of each recession
end


%% ===| Excercise 3 - Write a function for indpnormgam :

% Set up input parameters
d = diag(tRecession); % diagonal Matrix indicating calm/recession period (H-Matrix-Mechanism)


%% Testing Zone
k =     size(X,2);      % Number of regressors
N =     size(y, 1);     % Sample size
T =     N;
T2 =    sum(diag(d));   % Number of recession periods, diag(d) is re-transorming Matrix into vector
T1 =    T-T2;           % Number of calm periods
kappa0 = zeros(k,k);
nu1_0 = 5;
ssq1_0 = 1;
nu2_0 = 5;
ssq2_0 = 1;
I_t = eye(T);
S0 = 10;
S1 = 10;


% Intitialize storage matrices for posterior draws:
beta    = zeros(S0+S1, k);
h1      = zeros(S0+S1, 1);
h2      = zeros(S0+S1, 1);


%OLS estimates
XX       = X'*X;    % Compute X'X
Xy       = X'*y;    % Compute X'y
beta_OLS = X\y;     % OLS estimates
ssq_OLS  = ((y-X*beta_OLS)'*(y-X*beta_OLS))/(N-k) ;    % Variance of error term
V_OLS    = ssq_OLS*inv(XX);    % Variance-covariance matrix of OLS estimates


% Initial draw beta(0) using OLS estimates:
beta_initial = mvnrnd(beta_OLS, V_OLS, 1);


% === Posterior parameters:

% First draw h(1)|beta_intial:
nu1_1 = T1 + nu1_0;
ssq1_1 = ((y-X*beta_initial')'*d*(y-X*beta_initial')+nu1_0*ssq1_0)/nu1_1;
h1(1) = gamrnd(nu1_1/2, 2/(nu1_1*ssq1_1));

% First draw h(2)|beta_intial:
nu2_1 = T2 + nu2_0;
ssq2_1 = ((y-X*beta_initial')'*d*(y-X*beta_initial')+nu2_0*ssq2_0)/nu2_1;
h2(1) = gamrnd(nu2_1/2, 2/(nu2_1*ssq2_1));



%% ===| Excercise 4 - Apply Function + Non-informative OLS prior :

% First draw beta(1)|h1(1), h1(2):
H         = diag(h1(1) * (1 - diag(d)) + h2(1) * diag(d));
V1        = inv(V_OLS + X'*H*X);     % Posterior variance for beta|h1,h2, inverse cause V^-1
beta1     = V1*(X'*H*y);   % Posterior mean for beta|h1,h2 - "*"V1 cause its a matrix so mult by inverse
beta(1,:) = mvnrnd(beta1',V1,1);     % Draw from the conditional Normal for beta|h



%% ===| Test area:

H_diag = h1(1) * (1 - diag(d)) + h2(1) * diag(d);
H = diag(H_diag);

H

%% ===| Excercise ... :
% ===| Set regression parameters:
% Recall: ...
ssq0 = 1; % Prior variance
nu0 = 2;  % Degrees of freedom
nu1 = 2; 
beta0 = zeros(k, 1); % Prior mean
kappa0 = 0.1; % Prior precision
d = diag(tRecession);
S0 = 1000;
S1 = 10000;

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
ssq_OLS  = ((y-X*beta_OLS)'*(y-X*beta_OLS))/(N-k);    % Variance of error term
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

