clear; clc;
rng(123456);
% ========================== PC2: Exe2 ===============================
% ====================================================================
% Import data
data      = readtable('RRdata.xlsx');  % function: readtable
GDPgr     = data.PCGDP1;               % GDP growth (PCGDP1)
Tax       = data.EXOGENRRATIO;         % Tax shock (EXOGENRRATIO) 
timeraw      = (datetime(1947,01,01):calquarters(1):datetime(2007,12,01))'; 
%% ------------------------ part 1.a) -------------------------------
% Estimate the model: GDPgr= a + sum^{M=12}_{i=0} b_i Tax_{t-i} + u_t;
% In the matrix form
%                     Y_t = X_t B + u_t
% where: 
%          Y_t = GDPgr
%          X_t = [1 X_t X_{t-1} .. X_{t-M}]
%          B = [a b_i]'
% ------------------------------ Set up -----------------------------
M        = 12; % lag length
Xlags    = ...;                                    % create lag matrix: lagmatrix
X        = Xlags(M+1:end,:);                       % sufficient observations of X
Y        = GDPgr(M+1:end);                         % sufficient observations of Y
time     = timeraw(M+1:end);  
hor      = 0:1:M;                                  % number of response horizon
K        = size(Xlags,2)+1;                        % number of estimators
T        = size(Xlags,1);                          % number of sufficient observations
% --------------------------- Estimation --------------------------------
mdlest = ...;                                      % estimate model using fitlm
disp(mdlest)                                       % display regression output
% ===== store the estimated values
Bhat    = mdlest.Coefficients.Estimate;            % store the estimators
res     = mdlest.Residuals.Raw;                    % store the resiudals
VarBhat = mdlest.CoefficientCovariance;
% ===== cumulated effects
beta    = Bhat(2:end);                             % beta=[b0 b1 ...b12]' 
resp    = ...;                                     % cumsum
% use bootstrap to calculate the standard error estimation for 
% impulse responses (cumulative effects)
% ======================== boostrap =====================
ndraws=10000;
P = chol(VarBhat,'lower');  % Cholesky decomposition.
resp_bstr=zeros(ndraws,M+1);
% start boostrappingg
for j=1:ndraws
    coeff_bstr=Bhat+P*mvnrnd(0,1,K); % draw beta from the normal distribution
    b_bstr=coeff_bstr(2:end);         % estimated coefficients of tax shocks
    bsum=cumsum(b_bstr);              % the cumulated effects 
    resp_bstr(j,:)=bsum';            
end
%=========================== end boostrap ================
% mean and estimated standard error 
resp_mean = mean(resp_bstr,1)';            % mean 
resp_std  = std(resp_bstr,1)';             % std
% t-value
gam_tval  = resp./resp_std;
% lower and upper bound of 68% CI
resp_lb   = resp+tinv(0.16,T-M)*resp_std;  % lower bound 
resp_ub   = resp+tinv(0.84,T-M)*resp_std;  % upper bound 
% ===== plot 
...

%% ------------------------ Part 1b:  Autocorrelation -------------------
% ===== plot autocorrelogram: autocorr
...

% ===== Ljung-Box Q-test for residual autocorrelation: lbqtest
...

%% --------------- Part 1c: Newey West standard error --------------------
% ===== Newey-West covariance matrix estimation: hac
covNW = ...; 

% =========================== boostrap ==================================
ndraws=10000; 
P=chol(covNW,'lower');  % Cholesky decomposition.
resp_bstr=zeros(ndraws,M+1);
% start boostrapping
for j=1:ndraws
    ABhat_bstr=Bhat+P*mvnrnd(0,1,K); % draw beta from the normal distribution
    b_bstr=ABhat_bstr(2:end);         % estimated coefficients of tax shocks
    bsum=cumsum(b_bstr);              % the cumulated effects 
    resp_bstr(j,:)=bsum';            
end
% ===== Mean and Estimated standard error for the responses
respmean     = mean(resp_bstr,1)';
respseNW     = std(resp_bstr,1)';

% ======================== End boostrap =============================
% ===== Newey-West t-value
gam_NWtval=resp./respseNW;
% construct table 
results3=array2table([hor', resp, gam_tval, gam_NWtval],"VariableNames",...
    ["horizon","responses","conventional t-value","Newey-West t-value"]);
disp(results3)
