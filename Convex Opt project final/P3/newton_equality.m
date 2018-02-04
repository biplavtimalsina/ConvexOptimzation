clear all
clc
n=100;
p=30;
% A=randn(p,n);
% x0=rand(n,1);
load A
load x0
b=A*x0;
MAXITERS = 100;
ALPHA = 0.01;
BETA = 0.5;
NTTOL = 1e-7;
x = x0;
%find p_star
% cvx_begin
% variables x(n)
% minimize sum(x'*log(x));
% cvx_end
% p_star = cvx_optval;
p_star=-34.347345215327020;
for iter=1:MAXITERS
% calculate value function
val = x'*log(x);
% calculate gradient of objective fcn
grad = 1+log(x);
% calculate hessian of objective fcn
hess = diag(1./x);
% solve the equation Ax=b based pn KKT matrix
dxw = -[hess A'; A zeros(p,p)] \ [grad; zeros(p,1)];
% calculate direction
dx = dxw(1:n);
% calculate lambda for stopping criterion
lambda = grad'*dx;
% check stopping criterion
if (abs(lambda) < NTTOL), break; end;
% implement backtracking line search
t=1;
while (min(x+t*dx) <= 0), t = BETA*t; end;
while ((x+t*dx)'*log(x+t*dx) >= val + t*ALPHA*lambda), t=BETA*t; end;
% update x
x = x + t*dx;
% calculate error
f=x'*log(x);
error(iter)=f-p_star;
end;
figure(1)
loglog(error,'ko-')