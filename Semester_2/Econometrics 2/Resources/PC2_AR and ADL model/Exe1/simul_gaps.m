clear; clc;
% ====== a)
% Open the function ar1simul.m and fill in the gaps

% ====== b)
% set parameter values
sig  = 1;                    % variance of disturbance
rep  = 1000;                 % number of replications
avec = [0;0.5;0.9;0.99];     % AR(1) coefficient
Tvec = [20; 40; 80];         % number of observations

% set matrices to store values in simulation
ahatvar_theo = zeros(length(avec),length(Tvec));     % =zeros(4,3)
Ahat         = zeros(rep,length(avec),length(Tvec)); % zeros(1000,4,3)

% start simulation and estimation
for i=1:length(avec)
    a=avec(i);
    for j=1:length(Tvec)
        T=Tvec(j);
        for s=1:rep
            Y           = ar1simul(a,T);                      % simulate data set
            ylag        = ...;                                % compute y(t-1)
            yt          = ...;                                % compute y(t)
            ahat        = ...;                                % OLS estimator
            Ahat(s,i,j) = ahat;
        end     

        % compute theoretical asymptotic variance
        ahatvar_theo(i,j) = ...;                  % theoretical variance 

    end
end

% plot the simulated and theoretical distribution
for j = 1:length(Tvec)
    figure(j)
    for i = 1:length(avec)
        [pdf_sim, as] = ksdensity(Ahat(:,i,j));                      % estimate empirical density "pdf_sim" at points "as"  
        pdf_theo = normpdf(as,avec(i),sqrt(ahatvar_theo(i,j)));      % compute asymptotic density "pdf_theo" at points "as"
        subplot(2,2,i); plot(as,pdf_sim,'k-');                       % plot empirical density
        subplot(2,2,i); hold on; plot(as,pdf_theo,'r-');             % plot theoretical density
        subplot(2,2,i); hold on; xline(mean(Ahat(:,i,j)),'-.k')
        legend('Simulated','Asymptotic','Location','NorthWest')      % annotate figure
        xlabel('\alpha')
        ylabel('density')
        title(['\alpha = ', num2str(avec(i)), ', T = ', num2str(Tvec(j))])
    end
    
end




