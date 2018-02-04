clear all
clc
n=100;
p=30;
%load A
A=rand(100,30)
x0=ones(n,1);
%%
MAXITERS = 100;
ALPHA = 0.01;
BETA = 0.5;
RESTOL = 1e-7;
x=x0; 
nu=zeros(p,1);
p_star=-34.347345215642460;
for i=1:MAXITERS
% build the residual vector
r = [1+log(x)+A'*nu; A*x-b]; 
% solve the equation Ax=b based on KKT matrix
sol = -[diag(1./x) A'; A zeros(p,p)] \ r;
Dx = sol(1:n); Dnu = sol(n+[1:p]);
% check stopping criterion
if (norm(r) < RESTOL), break; end;
% implement backtracking line search
t=1;
while (min(x+t*Dx) <= 0), t = BETA*t; end;
while norm([1+log(x+t*Dx)+A'*(nu+Dnu); A*(x+Dx)-b]) > ...
(1-ALPHA*t)*norm(r), t=BETA*t; end;
% update x and v
x = x + t*Dx; nu = nu + t*Dnu;
% calculate error
f=x'*log(x);
error(i)=f-p_star;
res_dual(i)=norm(r(1:n));
res_pri(i)=norm(r(n+[1:p]));
end;
figure(1)
loglog(error,'ko-')

i=[1:6]
loglog(i,res_dual,i,res_pri);
