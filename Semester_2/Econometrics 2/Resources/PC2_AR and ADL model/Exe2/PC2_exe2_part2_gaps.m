% Load the exercise 2, part 1
PC2_exe2_gaps;
close all;
clc;

%=========================================================================
% ------------------------  ADL model -----------------------------------
%=========================================================================
%% Part 2a) Autocorrelation test for the residuals
%  -------------------- Estimation -------------------------
% Model 2:
%      GDPgr=a+\sum^M_{i=0} b_i Tax_{t-i} + \sum^N_{j=1} c_j GDPgr_{t-j} + u_t
% In matrix form
%                         Y = X B + u
% X=[1 Xlags Ylags]; 
% where: 
% Xlags = [Tax_t Tax_{t-1} Tax_{t-2} ... Tax_{t-M}], with M=12;
% Ylags = [GDPgr_{t-1} GDPgr_{t-2} ... GDPgr_{t-N}], with N=11;

% ===== Construct the Xlags and Ylags: use the function lagmatrix
M = 12;                        % lags number of Tax shock
N = 11;                        % lags number of GDP growth

Xlmat = ...;                   % lag matrix of tax shock
Ylmat = ...;                   % lag matrix of GDP growth

% ===== construct X matrix in model 2.
Xlags = Xlmat(M+1:end,:);      
Ylags = Ylmat(M+1:end,:);
XM2   = [Xlags, Ylags];
% ===== use function fitlm to estimate model 2 
ADLmdl = ...;
disp(ADLmdl)
% store estimated values
BhatM2     = ...;             % estimators
resM2      = ...;             % residuals
covBhatM2  = ...;             % estimated variance of estimators

% -------------------- Residual diagnostic -----------------
% ===== plot autocorrelogram: autocorr
...

% ===== Ljung-Box Q-test for residual autocorrelation: lbqtest
...

%% Part 2b) Calculate the (cumulated) impulse responses
% --------------------- (cumulated) impulse responses ----------------
% B = [a,b0,b1,b2,..b12,c1,c2,..., c11]';
...

% ------------------- BOOTSTRAP -----------------------------------
% Using boostrap to calculate the estimated variance of the responses.
% ===== boostrap
ndraws=100000;
P=chol(covBhatM2,'lower');
gamM2_bstr=zeros(ndraws,M+1);
gamM2_2bstr=zeros(ndraws,M+1);
for j=1:ndraws
    Bhat_bstr=BhatM2+P*mvnrnd(0,1,M+N+2); % draw beta from the normal distribution
    bhat=Bhat_bstr(2:M+2);
    ahat=[Bhat_bstr(M+3:end); 0];
    Del=bhat(1);
    del_i=zeros(M+1,1);
    del_i(1)=Del(1);
    for i=1:12
        del_i(i+1,1)=bhat(i+1)+ahat(1:i)'*Del;
        Del=[del_i(i+1); Del];
    end
    bsum=cumsum(del_i);
    gamM2_bstr(j,:)=bsum'; 
end
IRFmean_M2= mean(gamM2_bstr,1)';
IRFse_M2 = std(gamM2_bstr,1)';

% 68\% confidence interval (1 standard deviation)
respM2_lb=respM2+tinv(0.16,T).*IRFse_M2;
respM2_ub=respM2+tinv(0.84,T)*IRFse_M2;

%------------------- END BOOSTRAP ------------------------------
% Replicate figure 5
...


%% Part 2c) Information criteria (IC)
% ADL model:
% Y_t=a+b0 X_t+b1 X_{t-1}+...+bp X_{t-p}+c1 Y_{t-1}+...+cq Y_{t-q}+u_t
lagmax    = 10; % maximum lag number in the model (p=q)
Ynew      = GDPgr(2:end); % because of missing values of the first obs
Xnew      = Tax(2:end);
AIC       = zeros(lagmax,1);
AICc      = zeros(lagmax,1);
BIC       = zeros(lagmax,1);
HQ        = zeros(lagmax,1);
for p=1:lagmax
    % re-estimate model for each case of p
    Ylags        = ...;               % lagmatrix for GDP growth
    Xlags        = ...;               % lagmatrix for exogenous tax changes
    Ytemp        = ...;               % dependent variable Y
    Xtemp        = ...;               % X=[Xlags Ylags]
    [nobs,npara] = size(Xtemp);       % number of observations and number of parameters
    mdl_p        = fitlm(Xtemp,Ytemp);% function: fitlm to estimate the model 
    logL         = ...;               % store value of loglikelihood
    [aic,bic,ic] = ...;               % function: aicbic for calculate information criteria
    % save the results
    AIC(p)  = aic;
    AICc(p) = ic.aicc;
    BIC(p)  = bic;
    HQ(p)   = ic.hqc;   
end
[~,AIClags]    = min(AIC);
[~,AICclags]   = min(AICc);
[~,BIClags]    = min(BIC);
[~,HQlags]     = min(HQ);

% display
ICtable = array2table([(1:lagmax)',AIC,AICc,BIC,HQ],"VariableNames",["lags","AIC","AICc","BIC","HQ"]);
disp(ICtable)
disp(['AIC selects the lag number = ', num2str(AIClags)])
disp(['AICc selects the lag number = ', num2str(AICclags)])
disp(['BIC selects the lag number = ', num2str(BIClags)])
disp(['HQ selects the lag number = ', num2str(HQlags)])





