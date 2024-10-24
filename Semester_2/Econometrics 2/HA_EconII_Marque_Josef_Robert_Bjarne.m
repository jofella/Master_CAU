% Josef Fella stu245231 QF
% Robert Hennings stu236320 QF
% Marque Mollenhauer stu227420 QE
% Bjarne Ehmcke-Kasch stu245454 Economics

%% Clear all
clear % to remove all variables from the workspace
clc % to clear the command window
%% Excercise 1
data = readtable('HAdata.xlsx');
time = (datetime(1985,01,01):calquarters(1):datetime(2024,01,01))';
GDP = data.GDP;
SPREAD = data.SPREAD;
figure
subplot(2,1,1)
plot(time,GDP,'LineWidth',1.5);
recessionplot;
title("GDP")
xlabel('quarterly');
ylabel('GDP');

subplot(2,1,2)
plot(time,SPREAD,'LineWidth',1.5);
recessionplot;
title("Spread")
xlabel('quarterly');
ylabel('Spread');

%% Exercise 2
%GDP
lags = 8;
vars = GDP;
[~,pvalAR] = adftest(vars,"Model","AR","Lags",1:lags); % test WITHOUT interceptand trend
[~,pvalARD] = adftest(vars,"Model","ARD","Lags",1:lags); % test WITH intercept
[~,pvalTS] = adftest(vars,"Model","TS","Lags",1:lags); % test WITH intercept and trend
% display and save the results
resultsGDP = table((1:lags)',pvalAR',pvalARD',pvalTS','VariableNames',["lag","AR","ARD","TS"])
disp('Augmented Dickey-Fuller Test GDP (p-values)');
disp(resultsGDP)

% GDP contains the unit root as the p-values are close to one so we do not
% reject H0 of the ADF-test --> non-stationary


%SPREAD
lags = 8;
vars = SPREAD;
[~,pvalAR] = adftest(vars,"Model","AR","Lags",1:lags); % test WITHOUT intercept and trend
[~,pvalARD] = adftest(vars,"Model","ARD","Lags",1:lags); % test WITH intercept
[~,pvalTS] = adftest(vars,"Model","TS","Lags",1:lags); % test WITH intercept and trend
% display and save the results
resultsSPREAD = table((1:lags)',pvalAR',pvalARD',pvalTS','VariableNames',["lag","AR","ARD","TS"]);
disp('Augmented Dickey-Fuller Test SPREAD (p-values)');
disp(resultsSPREAD)

% For the variable Spread we reject H0 = no unit root --> likely stationary
%% Exercise 3
dlogGDP = diff(log(GDP));
SPREADshift = SPREAD +1.55; % assure positive values
min(SPREADshift);
lagmax = 10;
Ynew = dlogGDP(1:end);
Xnew = SPREAD(2:end);
AIC = zeros(lagmax,1);
BIC = zeros(lagmax,1);
for p=1:lagmax
    Ylags = lagmatrix(Ynew,1:p); % lagmatrix for GDP 
    Xlags = lagmatrix(Xnew,0:p); % lagmatrix for Spread
    Ytemp = Ynew(p+1:end,1); % dependent variable Y
    Xtemp = [Xlags(p+1:end,:),Ylags(p+1:end,:)];
    [nobs,npara] = size(Xtemp); % number of observations and number of parameters
    mdl_p = fitlm(Xtemp,Ytemp); % function: fitlm to estimate the model
    logL = mdl_p.LogLikelihood; % store value of loglikelihood
    [aic,bic,ic] = aicbic(logL,npara,nobs); % function: aicbic for calculate information criteria
 % save the results
    AIC(p) = aic;
    
    BIC(p) = bic;
   
end
[~,AIClags] = min(AIC);
[~,BIClags] = min(BIC);

% display
ICtable = array2table([(1:lagmax)',AIC,BIC],"VariableNames",["lags","AIC","BIC",]);
disp(ICtable)

% Using the information criteria AIC (-933.09 ) and BIC (-923.96) the minimizing values for p
% and q is one

%% Exercise 4
Ylags=lagmatrix(dlogGDP,1);
Xlags=lagmatrix(SPREAD,0:1);
Xtemp4=[Ylags(2:end,:),Xlags(3:end,:)];
Ydep4=dlogGDP(2:end);

mdlest=fitlm(Xtemp4, Ydep4);

disp(mdlest)

resmd4 = mdlest.Residuals.Raw;
figure;
autocorr(resmd4, 20);


%a) Based on the autocorrelogram of the residuals the residuals seem to be
%not (too heavily) autocorrelated, because the lie within the confidence
%bounds. In the magnitude they are rather small and there is no overall
%structure to be seen and therefor a high degree of dispersion.


% ==== L jung-Box Q-test for residual autocorrelation: llbqtest
[~,pval_lbq] = lbqtest(resmd4,"Lags",1:12);
disp(pval_lbq)

%Based on the results for the p-values for all lag lengths we cannot reject
%the H0 and therefore we do not assume autocorrelation. The p-value
%increases with the lag order.


%b) ===== Newey-West covariance matrix estimation: hac
T= size(Xlags,1);
q = round(T^(1/4));
[EstCoeffCov,se,coeff]= hac(Xtemp4,Ydep4,"Bandwidth",q+1);
se

r = [0;0];
R = [0 0 1 0; 0 0 0 1 ];

[h,p,Wstat,crit] = waldtest (r,R, EstCoeffCov,0.05)
p

% We reject the H0