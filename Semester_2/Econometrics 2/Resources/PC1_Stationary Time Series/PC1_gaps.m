clear; clc;
%% -------------------------- EXE1 ------------------------
% Part a
 % ===== import data: xlsread
 Y        = xlsread('USdata.xlsx');     % load dataset from excel file
 [T, N]   = size(Y); % dimensions of Y
 varnames = {'Industrial production';'Unemployment rate';'Money supply';'Stock prices (% change)'};
 unit     = {'index';'rate';'index';'rate'};
 % to construct a sequence of time: datetime
 time = (datetime(1985,01,01):calmonths(1):datetime(2024,02,01))';
 % ===== plot the ts: plot
 figure
 plot(time,Y(:,1),'color','r','LineWidth',1.5);
 title('Industrial production')

 % multiple graphs in one figure: subplot
...

%% Part b
% ===== Data transformation
tcode = [3;2;3;1];             % 3: first diff of log, 2: first diff, 1: no transform
Yn    = zeros(T-1,N);          % create zeros matrix to store data
for i=1:N
    if tcode(i)==3      
        Yn(:,i) = ...;         % log(Y_{t}) - log(Y_{t-1})
    elseif tcode(i)==2
        Yn(:,i) = ...;         % Y_{t} - Y_{t-1}   
    else
        Yn(:,i) = Y(2:end,i);        
    end
end
time2 = time(2:end);
% plot time series
...

%% ------------------------------- EXE2---------------------------------
% part a
clear; clc;
% Import the dataset
Y            = ...; % xlsread
PCGDP        = ...; % real GDP
PCGDPl       = ...; % first diff of log GDP
EXOGENRRATIO = ...; % exogenous tax changes
time         = ...; % generate a sequence of quarterly datetime
% plot
figure
colororder({'k','k'});
yyaxis left
plot(time,PCGDPl,'-r','LineWidth',0.7); hold on;
yyaxis right
plot(time,log(PCGDP),'-b','LineWidth',0.7);
legend('GDP growth rate','GDP');
%% part b
% sample autocorrelation function: autocorr
[acf_s,lags,CI] = ...;
disp([lags,acf_s])
disp(CI)
% plot the autocorrelogram for the fiscal shock and GDP growth rate 
...
















