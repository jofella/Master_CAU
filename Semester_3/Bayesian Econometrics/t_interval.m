function interval = t_interval(mu,Sigma,nu,alpha)
% -------------------------------------------------------------------------
% This function constructs symmetric 100*(1-alpha)% intervals for the 
% multivariate t distribution with parameters mu, Sigma, and nu.
% -------------------------------------------------------------------------
% Input:
% mu          (k x 1) mean vector
% Sigma       (k x k) scale variance matrix
% nu          scalar degrees of freedom
% -------------------------------------------------------------------------
% Output:
% interval    (k x 2) matrix of lower and upper interval bounds
% -------------------------------------------------------------------------

% Quantile from the t distribution:
cv = tinv(1-alpha/2,nu); % inverse cdf of student-t and getting CI, two sided interval

% Vector of standard deviations provided by Sigma:
sd = sqrt(diag(Sigma));

% Credible interval:
interval = [mu - cv*sd, mu + cv*sd]; %computes interval


