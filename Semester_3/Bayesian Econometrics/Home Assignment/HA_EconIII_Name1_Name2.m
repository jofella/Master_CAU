%===========| Home Assignment Econometrics III (WiSe 2024/25) |============
% ===| Participants:
% Name 1 (matriculation number)
% Name 2 (matriculation number)
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
T = size(y,1);                                                              % Effective sample size
dates = dates(q+1:end);                                                     % Match sample size of dates
COVID = dates >= datetime(2020,04,01) & dates <= datetime(2020,06,01);      % Define COVID-19 dummy variable
X = [ones(T,1) X COVID];                                                    % Matrix of covariates
k = size(X,2);                                                              % Number of covariates
J = size(recessions,1);                                                     % Number of US recessions in the sample
tRecession = zeros(T,1);                                                    % Pre-allocate dummy-vector for recession periods
for j = 1:J
    tRecession = tRecession | (dates >= recessions(j, 1) & dates <= recessions(j, 2)); % Identify dates within the start and end period of each recession
end