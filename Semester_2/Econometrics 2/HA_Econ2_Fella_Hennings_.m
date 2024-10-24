clear; clc;
%% Group
%{
Josef Fella / QF / stu245231
Robert Hennings / QF / stu...
...
...
...
%}

%% Exercise 1 - Plot ts

% Importing files and structure data
GDP = xlsread('HAdata');
[T, N] = size(GDP);
varnames = {'GDP';'Term Spread'};
unit     = {'index';'rate'};

% ------------------------

% Build timeseries (quarterly)
time = (datetime(1985,01,01):calmonths(3):datetime(2024,02,01))';

% ------------------------

% Plotting results in one figure object
figure
for i = 1:N
    subplot(N, 1, i)
    plot(time, GDP(:, i), 'color', 'r', 'LineWidth', 2);
    recessionplot; % nice to see US business cycle
    
    xlabel('quarterly');
    ylabel(['' unit{i}])
    title(['' varnames{i}])
end

%% Exercise 2 - Perform ADF test
clear, clc;

% Prepare dataset
data = readtable("HAdata.xlsx");
dates = data.Quarterly;

% Define variables
GDP = data.GDP;
Spread = data.SPREAD;


% ------------------------

% Compute ADF test

% Prepare data
var_GDP = log(GDP);
var_Spread = log(abs(Spread));
lags = 8; % Max. lag order --> empirically shown for qrtly


% 1. Test GDP ------------------------
[GDP_AR, GDP_pVal_AR] = adftest(var_GDP, "Model", "AR", "Lags", 1:lags); % 1.ADF
[GDP_ARD, GDP_pVal_ARD] = adftest(var_GDP, "Model", "ARD", "Lags", 1:lags); % 2. ADF with intercept
[GDP_TS, GDP_pVal_TS] = adftest(var_GDP, "Model", "TS", "Lags", 1:lags); % 3. ADF with intercept + trend

% Display results
 results_GDP = table((1:lags)', GDP_pVal_AR', GDP_pVal_ARD', GDP_pVal_TS','VariableNames',["Lags",...
 "GDP_pval-AR","GDP_pval-ARD","GDP_pval-TS"]);
 disp('Augmented Dickey-Fuller Tests (p-values)')
 disp(results_GDP)


% 2. Test Spread ------------------------
[Spread_AR, Spread_pVal_AR] = adftest(var_Spread, "Model", "AR", "Lags", 1:lags); % 1.ADF
[Spread_ARD, Spread_pVal_ARD] = adftest(var_Spread, "Model", "ARD", "Lags", 1:lags); % 2. ADF with intercept
[Spread_TS, Spread_pVal_TS] = adftest(var_Spread, "Model", "TS", "Lags", 1:lags); % 3. ADF with intercept + trend


% Display results
 results_Spread = table((1:lags)', Spread_AR', Spread_ARD', Spread_pVal_TS','VariableNames',["Lags",...
 "Spread_pval-AR","Spread_pval-ARD","Spread_pval-TS"]);
 disp('Augmented Dickey-Fuller Tests (p-values)')
 disp(results_Spread)


% ------------------------

% Interpretation of ADF hypothesis

%{
1. GDP:
has not unit root



2. Spread:
has a unit root


%}


%% Exercise 3
clear, clc;

% Load the data from HAdata.xlsx
data = readtable('HAdata.xlsx');
GDP = data.GDP;
SPREAD = data.SPREAD;
GDPgr = diff(log(GDP)); % stationary GDP


% Construct ADL model ---------------------------

% ADL model parameters
lagmax = 10; % Maximum lag number in the model (p = q)

% Containers to store results
AIC = zeros(lagmax, 1);
AICc = zeros(lagmax, 1);
BIC = zeros(lagmax, 1);
HQ = zeros(lagmax, 1);

% Loop over possible values of p/q
for p = 1:lagmax
    % Re-estimate model for each case of p
    Ylags = lagmatrix(GDPgr, 1:p); % Lagmatrix for GDP growth
    Xlags = lagmatrix(SPREAD(2:end), 0:p); % Lagmatrix for exogenous spread changes
    Ytemp = GDPgr(p+1:end); % Dependent variable Y
    Xtemp = [Xlags(p+1:end,:), Ylags(p+1:end,:)]; % X = [Xlags Ylags]
    
    % Remove NaNs (resulting from lagging)
    validRows = all(~isnan(Xtemp), 2);
    Ytemp = Ytemp(validRows);
    Xtemp = Xtemp(validRows, :);
    
    % Estimate the model
    mdl_p = fitlm(Xtemp, Ytemp); % Function: fitlm to estimate the model
    logL = mdl_p.LogLikelihood; % Store value of log-likelihood
    nobs = mdl_p.NumObservations; % Number of observations
    npara = mdl_p.NumEstimatedCoefficients; % Number of parameters
    
    % Calculate information criteria
    [aic, bic, ic] = aicbic(logL, npara, nobs); % calculates information criterion
    
    % Save the results
    AIC(p) = aic;
    AICc(p) = ic.aicc;
    BIC(p) = bic;
    HQ(p) = ic.hqc;   
end

% Find the lag length with the minimum AIC, AICc, BIC, and HQ
[~, AIClags] = min(AIC);
[~, AICclags] = min(AICc);
[~, BIClags] = min(BIC);
[~, HQlags] = min(HQ);


% Plot results ---------------------------

% Display the information criteria for all lags
ICtable = array2table([(1:lagmax)', AIC, AICc, BIC, HQ], ...
    "VariableNames", ["lags", "AIC", "AICc", "BIC", "HQ"]);
disp(ICtable);

% Display the selected lags
disp(['AIC selects the lag number = ', num2str(AIClags)]);
disp(['AICc selects the lag number = ', num2str(AICclags)]);
disp(['BIC selects the lag number = ', num2str(BIClags)]);
disp(['HQ selects the lag number = ', num2str(HQlags)]);


%% Exercise 4 -- Tut 2 Ex.2
% a)

clear, clc;

% Prepare data
data = readtable('HAdata.xlsx');
GDP = data.GDP;
SPREAD = data.SPREAD;
GDPgr = diff(log(GDP));  % This is ∆Y_t

% Lagged variables for the ADL(1,1) model
Ylag1 = lagmatrix(GDPgr, 1);           % ∆Y_{t-1}
SPREADlag1 = lagmatrix(SPREAD(2:end), 1); % SPREAD_{t-1}

% Trim data to account for lags
validRows = ~isnan(Ylag1) & ~isnan(SPREADlag1);
Y = GDPgr(validRows);
X = [ones(sum(validRows), 1), Ylag1(validRows), SPREAD(2:end-1), SPREADlag1(validRows)];

% ------------------------------------------------

% Estimate the ADL using OLS
ADLmdl = fitlm(X, Y);

% Display the model summary
disp(ADLmdl);

% Get residuals
residuals = ADLmdl.Residuals.Raw;


% ------------------------------------------------

% Plot the autocorrelogram (ACF) of the residuals
figure;
autocorr(residuals);
title('ACF of Residuals for ADL(1,1) Model');

% Perform Ljung-Box Q-test for residuals
[h, pValue] = lbqtest(residuals);
disp(['Ljung-Box Q-test p-value = ', num2str(pValue)]);
if h == 0
    disp('Residuals are not autocorrelated.');
else
    disp('Residuals are autocorrelated.');
end



% b)Wald test

if h ~= 0
    % Re-estimate the model with Newey-West standard errors
    [b, ~, ~, covb] = hac(X, Y, 'Type', 'NeweyWest', 'Lags', 1);
    
    % Display the coefficients and Newey-West standard errors
    se = sqrt(diag(covb));
    tStats = b ./ se;
    pValues = 2 * (1 - tcdf(abs(tStats), length(Y) - length(b)));
    resultsTable = table(b, se, tStats, pValues, ...
        'VariableNames', {'Coefficients', 'NW_SE', 'tStats', 'pValues'});
    disp('Newey-West Results:');
    disp(resultsTable);
    
    % Perform Wald test on the coefficients (excluding intercept)
    R = [0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1]; % Restriction matrix for α1, β0, β1
    q = zeros(size(R, 1), 1); % Null hypothesis: coefficients are zero
    W = (R * b - q)' * inv(R * covb * R') * (R * b - q); % Wald statistic
    pWald = 1 - chi2cdf(W, size(R, 1)); % p-value of the Wald test
    disp(['Wald test p-value = ', num2str(pWald)]);
end