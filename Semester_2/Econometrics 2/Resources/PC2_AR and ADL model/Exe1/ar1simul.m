function [y] = ar1simul(a,T)
%%   This function simulates an AR(1) process.
% --------- AR(1) model: y(t) = a*y(t-1) + e(t)
% --------- Input:
%                       a = coefficient
%                       T = sample size
% ---------- Output
%                       Y = simulated data with T observations

e    = randn(T,1);  % error term
y    = zeros(T,1);  % store
y(1) = e(1);

%  use 'for' loop to generate observations y(t) 
...


