% =========================================================================
% ======================   PC4: Cointegration =============================
% ========================================================================
%% ------------------------- Exercise 2 ----------------------------------
clear; clc;
data = ...;
%% part a) Plot time series of LM3real, LGDPreal, and LT
...
  
%% part b) ADF test
lags        = ...;
vars        = ...;
[~,pvalAR]  = ...; % ADF test WITHOUT intercept and trend
[~,pvalARD] = ...; % ADF test WITH intercept
[~,pvalTS]  = ...; % ADF test WITH intercept and trend

% display and save results
...
%% part c) Engle-Granger test
lags          = ...;
Y             = ...;
[~,pegci_nc]  = ...;    % WITHOUT constant and trend
[~,pegci_c]   = ...;    % WITH constant and NO trend
[~,pegci_ct]  = ...;    % WITH constant and trend
[~,pegci_ctt] = ...;    % WITH constant, trend, and quadratic trend

% display and save the results
...
