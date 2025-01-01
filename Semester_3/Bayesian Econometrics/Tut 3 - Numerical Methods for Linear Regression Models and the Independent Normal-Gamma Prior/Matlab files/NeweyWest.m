function [lrvar, gamma] = NeweyWest(y,k)
% -------------------------------------------------------------------------
% This function computes Newey-West long run variance using a prespecified
% lag window 
% -------------------------------------------------------------------------
% Input:    
% y        vector of data
% k        scalar number of lags
% -------------------------------------------------------------------------
% Output:   
% lrvar    long-run variance
% gamma    vector of autocovariances 0 to k
% -------------------------------------------------------------------------

T = length(y);          % dimension of y
y = reshape(y,T,1);     % ensures that y is a column vector

% Compute autocovariances:
gamma = zeros(k+1,1);
for i = 0:k
    gamma(i+1) = (y(1+i:end)-mean(y(1+i:end)))'*(y(1:end-i)-mean(y(1:end-i)))/(T-k-1);
end

% Compute long-run variance:
w = [1, 2*(1-(1:k)/(k+1))]';
lrvar = sum(gamma.*w);

    