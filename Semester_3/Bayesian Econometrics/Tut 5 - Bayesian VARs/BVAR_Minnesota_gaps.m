%--------------------------------------------------------------------------
% Bayesian VAR estimation using the Minnesota Prior
% The VAR model is
%
%     y(t) = A0 + A1 x y(t-1) + ... + Ap x y(t-p) + e(t)
%
% The matrix representation of the VAR is
%
%                   Y = X x A + E
%
% where E ~ N(0,SIGMA), E[e(t)e(s)]=0 for any t not equal to s. 
% Note that A summarizes all VAR coefficients and its vectorized form will 
% be given by a=vec(A). 
%
% Written by Richard Schnorrenberger (2022): richard.schn@stat-econ.uni-kiel.de
%--------------------------------------------------------------------------

%----------------| LOAD DATA AND GET STATIONARY VARIABLES |----------------
% ===| Clear data memory and Matlab environment:
clear;
clc;

% ===| Load monthly data set on the US macroeconomy from 1982M1 to 2020M2:
load US_macrodata.mat

% ===| Plot raw data:
figure
for i = 1:8
    subplot(3,3,i)
    plot(dates,data(:,i),'Color','#0072BD','LineWidth',1.5);
    title(VARnames{i});
    recessionplot;
end

% ===| Transform raw data into stationary time series (or almost stationary):
CPI    = diff(log(data(:,1)))*100;    % CPI inflation (m.o.m. percent growth rate)
PROD   = diff(log(data(:,2)))*100;    % Industrial production (m.o.m. percent growth rate)
UNEMP  = data(2:end,3);               % Unemployment rate (percent levels)
EXP    = diff(log(data(:,4)))*100;    % Real expenditure (m.o.m. percent growth rate)
FED3M  = diff(data(:,5));             % 3-Month rate (m.o.m. difference in levels)
CRED   = diff(log(data(:,6)))*100;    % Bank credit (m.o.m. percent growth rate)
SPREAD = data(2:end,7);               % Yield spread rate (percent levels)
SP500  = diff(log(data(:,8)))*100;    % S&P 500 (m.o.m. percent growth rate)
Yraw   = [CPI PROD UNEMP EXP FED3M CRED SPREAD SP500]; % Stationary data Y
dates  = dates(2:end);

% If you want to use this code with your data set, name the data set as
% 'Yraw'. Note that 'Yraw' is a matrix with T rows by M columns, where T 
% is the time series sample size and M is the number of endogenous variables
% or number of equations in the VAR(p) model.

%------------------------| PRELIMINARIES |---------------------------------
% ===| Define VAR lags and intercept specification:
constant = "...";           % 1: if you desire intercepts, 0: otherwise 
p = "...";                  % Number of lags on dependent variables

% Note that the number of VAR coefficients equals (pM+1)*M, including the
% constant. Hence it easily gets overparameterized.

%------------------------| DATA HANDLING |---------------------------------
% ===| Get initial dimensions of dependent variables:
[Traw, M] = "...";
        
% ===| Generate lagged Y matrix, which will be part of the X matrix:
Ylag = "...";     % Y is T-by-M while ylag is T-by-(pM)

% ===| Define matrix of VAR regressors X which has all the R.H.S. variables 
% including the constant and lags of the dependent variable:
if constant == 1
    X = "...";
else
    X = "...";  
end

% ===| Obtain number of regressors in each equation of the VAR:  
K = "...";        % Also equal to the column dimension of matrix X 

% ===| Create final Y matrix by removing those initial rows of Y that refer to 
% the number of choosen lags:  
Y = "...";        % This is the final Y matrix used for the VAR

% ===| Traw was the dimension of the initial data. T is the number of actual 
% time series observations of Y and X:
T = "...";

%-----------------------| OLS QUANTITIES |---------------------------------
% ===| Get OLS estimators which are part of the posterior formulae:
A_OLS = "...";                   % Matrix of OLS estimates of VAR coefficients
a_OLS = "...";                   % Vectorized form of A_OLS such that a_OLS=vec(A_OLS)
SSE = "...";                     % OLS sum of squared residuals
SIGMA_OLS = "...";               % OLS residual variance

%---------------------| PRIOR SPECIFICATION |------------------------------
% ===| Prior mean:
    A_prior = "...";             % Matrix of prior means
    a_prior = "...";             % Vector of prior means
    
% ===| Get residual variances of univariate p-lag autoregressions, which are
% used in the Minnesota prior covariance matrix Sigma:   
    sigma_sq = "...";            % Vector to store residual variances
    for i = 1:M
        % Create lags of dependent variable in i-th equation:
        Ylag_i = "...";
        Ylag_i = "...";
        % Dependent variable in i-th equation:
        Y_i = "...";
        % OLS estimates of i-th equation:
        alpha_i = "...";
        sigma_sq(i,1) = "...";
    end 
    
% ===| Hyperparameters on the Minnesota variance of alpha:
    a_bar_1 = "...";      % for coefficients on own lags (if i=j): lambda^2
    a_bar_2 = "...";      % for coefficients on lags of other variables: theta_1^2 x lambda^2
    a_bar_3 = "...";      % for coefficients on intercepts/exogenous variables: theta_2^2 x lambda^2
    
% ===| Create a matrix of dimensions K x M, which will contain the K diagonal
% elements of the prior covariance matrix, in each of the M equations.
    V_i = zeros(K,M);
    % Index in each equation which are the own lags
    ind = zeros(M,p);
    for i=1:M
        ind(i,:) = constant+i:M:K;
    end
    for i = 1:M  % for each i-th equation
        for j = 1:K   % for each j-th RHS variable
            if constant==1
                if j==1 % if there is constant, use this code
                    V_i(j,i) = a_bar_3*sigma_sq(i,1); % variance on constant                
                elseif find(j==ind(i,:))>0
                    V_i(j,i) = a_bar_1./(ceil((j-1)/M)^2); % variance on own lags           
                    % Note: the "ceil((j-1)/M)" command finds the associated lag 
                    % number for each parameter
                else
                    for kj=1:M
                        if find(j==ind(kj,:))>0
                            ll = kj;                   
                        end
                    end                 % variance on other lags  
                    V_i(j,i) = (a_bar_2*sigma_sq(i,1))./((ceil((j-1)/M)^2)*sigma_sq(ll,1));           
                end
            else   % if no constant is defined, then use this code
                if find(j==ind(i,:))>0
                    V_i(j,i) = a_bar_1./(ceil(j/M)^2); % variance on own lags
                else
                    for kj=1:M
                        if find(j==ind(kj,:))>0
                            ll = kj;
                        end                        
                    end                 % variance on other lags  
                    V_i(j,i) = (a_bar_2*sigma_sq(i,1))./((ceil(j/M)^2)*sigma_sq(ll,1));            
                end
            end
        end
    end
 
    
    % The prior covariance matrix is diagonal with elements V_i:
    V_prior = diag(V_i(:));     % Prior covariance matrix of vector a_prior  
    
%----------------------| POSTERIOR COMPUTATION |---------------------------
% ===| Posterior variance of VAR coefficients:
V_post = "...";

% ===| Posterior mean of VAR coefficients:
a_post = "...";

% ===| Reshape posterior results:
A_post     = "...";
A_post_std = "...";
% Note that each column here stores posterior estimates for each equation
% of the VAR(p) model.
