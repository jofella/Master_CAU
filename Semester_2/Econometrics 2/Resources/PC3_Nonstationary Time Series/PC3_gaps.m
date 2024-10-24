% =========================== PC Tutorial 03 ==============================
clear, clc
 %% Part (a)
% ===== Import the dataset:
data = readtable("Exchange_Rates.xlsx"); % readtable

% ===== Construct a sequence of dates
dates = data.Dates;

% ===== Define variables
E_UK = data.E_UK; % Define the UK/US nominal exchange rate
E_CH = data.E_CH; % Define the CH/US nominal exchange rate

% ===== Plot nominal exchange rate
figure
subplot(2,1,1)

plot(dates,E_UK,'LineWidth',1.5);
title("nominal exchange rate between UK Pound and US dollar")
subplot(2,1,2)

plot(dates,E_CH,'LineWidth',1.5);
title("nominal exchange rate between Chinese Yuan and US dollar")

% ------------------------------------------------------------------
%% Part (b)
% ===== Re-define the nominal exchange rate of CH/US from 2005Q3
CH_float = find(dates == datetime(2005,07,01));                  
% Find date row that initiates the managed float Chinese system

E_CH = data.E_CH(CH_float:end,1); %  Re-define the CH/US nominal exchange rate

% ===== Compute ADF tests: adftest
var  = log(E_UK); 

% Set variable to compute unit root tests
lags = 8;                        

% Maximum lag order
[h_AR,pVal_AR]   = adftest(var,"Model","AR","Lags",1:lags); % ADF without intercept and trend
[h_ARD,pVal_ARD] = adftest(var,"Model","ARD","Lags",1:lags); % ADF with intercept
[h_TS,pVal_TS]   = adftest(var,"Model","TS","Lags",1:lags); % ADF with intercept and trend

% ===| Display results:
 results = table((1:lags)',pVal_AR',pVal_ARD',pVal_TS','VariableNames',["Lags",...
 "pval-AR","pval-ARD","pval-TS"]);
 disp('Augmented Dickey-Fuller Tests (p-values)')
 disp(results)
% ----------------------------------------------------------------------
%% Part (c)

%% Part (c)

RER_UK   = data.E_UK.*(data.CPI_US./data.CPI_UK);                            
%Define the UK/US real exchange rate

RER_CH   = E_CH.*(data.CPI_US(CH_float:end,1)./ data.CPI_CH(CH_float:end,1));                            

 % ===== Plot real exchange rates:

figure

subplot(2,1,1)
plot(dates,RER_UK,'LineWidth',1.5);
title("real exchange rate between UK Pound and US dollar")
subplot(2,1,2)
plot(dates(CH_float:end,1),RER_CH,'LineWidth',1.5);
title("real exchange rate between Chinese Yuan and US dollar")
