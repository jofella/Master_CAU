% =========================================================================
% ======================   PC4: Cointegration =============================
% ========================================================================
%% ------------------------- Exercise 1 ----------------------------------
clear; clc;
% ===== import data from BPdata.xlsx
data  = ...;
tfp   = ...;      % total factor productivity
sp    = ...;      % stock prices
dates = ...;      % sequence of datetime

%% a) plot time series of tfp and sp
...

%% b) ADF test for tfp and sp: function adftest
lags        = ...;
vars        = ...;
[~,pvalAR]  = ...;  % test WITHOUT intercept and trend
[~,pvalARD] = ...;  % test WITH intercept
[~,pvalTS]  = ...;  % test WITH intercept and trend
% display and save the results
...

%% c) Engel-Granger cointegration test: function egcitest
lags          = ...;
Y             = ...;
[~,pegci_nc]  = ...;    % WITHOUT constant and trend
[~,pegci_c]   = ...;    % WITH constant and NO trend
[~,pegci_ct]  = ...;    % WITH constant and trend
[~,pegci_ctt] = ...;    % WITH constant, trend, and quadratic trend

% display and save the results
...
