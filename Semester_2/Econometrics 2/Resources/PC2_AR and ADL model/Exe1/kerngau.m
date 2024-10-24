function [pdf,grd] = kerngau(x,srt)

% [pdf,grd,cdf,kd] = kerngau(x,h,srt)
%
% Kernel density estimation for data in vector x with Gaussian kernel.
%
% smoothing constant h = 0.79*IQR*n^(-1/5)
%
%    s: standard deviation of data
%  IQR: interquartile range of data
%    n: sample size
%  srt: Schrittweite, d.h. für srt=2 wird bei der Berechnung nur jede 2. Beobachtung
%       des geordneten samples verwendet, um die Rechenzeit zu verringern.


nobs=length(x);
h=0.9*min([std(x) iqr(x)/1.34])*nobs^(-.2); % s. Silverman (1986), p. 48
x=sort(x);
spdf=zeros(length(x),1);
kk = fix(nobs/srt)*srt;
for k=1:srt:kk
    xx=(x-x(k))/h;
    pdf=(1/sqrt(2*pi))*exp(-(xx.^2)/2);
    spdf=spdf+pdf;
end
pdf=spdf/(kk*h/srt);
grd=x;

