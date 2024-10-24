clear, clc

%% ============ Exe 1 ==============
% a)
m = [1 2 3 4 5];
m = 1:5;
n = m';
% b)
O = ones(2);
Z = zeros(2,3);
I= eye(3);
% c )
Q = [1 2 3 4;
    5 6 7 8;
    9 10 11 12];

%% ============= Exe 2 ===============
% a)
P = Q';
% b)
c = P(:,end);
r = P(3,:);
m = P(3:end,2:end)

%% ============== Exe 3 ================
% a) 
T = 100;
u = randn(T,1);
% b)
y = zeros(T,1); % common way to initalize the vector/matrix
y0 = 0; a = 0.5; b = 0.2;
y(1) = a + b*y0 + u(1);
for t = 2:T
    y(t) = a+b*y(t-1)+u(t);
end

figure
plot(y)

%% =========== Exe 4 ===============
ylag1 = y(1:end-1);
X = [ones(T-1,1) ylag1];
Y = y(2:end);
Bhat = (X'*X)\(X'*Y)
disp(Bhat)
disp([a; b])



