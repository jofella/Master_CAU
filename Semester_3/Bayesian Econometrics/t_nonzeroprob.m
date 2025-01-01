function probs = t_nonzeroprob(mu,Sigma,nu)
% -------------------------------------------------------------------------
% This function computes probabilities P(Y<0) and P(Y>0) for the random 
% vector Y that has multivariate t distribution with parameters mu, Sigma, 
% and nu.
% -------------------------------------------------------------------------
% Input:
% mu          (k x 1) mean vector
% Sigma       (k x k) scale variance matrix
% nu          scalar degrees of freedom
% -------------------------------------------------------------------------
% Output:
% probs       (k x 2) matrix of P(Y<0) and P(Y>0) probabilities
% -------------------------------------------------------------------------

% Vector of standard deviations provided by Sigma:
sd = sqrt(diag(Sigma));

% Probabilities from the t distribution:
probs = tcdf(-mu./sd,nu);
probs(:,2) = 1-probs(:,1);



