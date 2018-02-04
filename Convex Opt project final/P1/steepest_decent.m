% Steepest Descend Algorithm 
clear all
clc
% input dimension of Q matrix
% n=input('Enter the dimension of the Q:');
n=100;
%% generate PSD matrix with definite condition numebr
con_num=20; % con_num=2,5,10,25,50 ==> condition number is con_num
QP=randi(con_num,[n,n]);
% QP=QP'*QP;
P=diag(diag(QP));
 ma=randi(n);
 mi=randi(n);
 while mi==ma
     mi=randi(n);
 end
 P(ma,ma)=con_num;
 P(mi,mi)=1;
%% randomly generated matrix Q, vector q and initial point x0
Q=P;
q=(randn(n,1));
x0 = (randn(n,1));
P=diag(diag(Q));

 % make sure Q is PSD
 eigen_Q=eig(Q);
 if (any(eigen_Q)<0)
     disp('Matrix Q is not a PSD matrix');
 end 

%% run cvx to find p_star
cvx_begin
variable x(n);
dual variable y;
minimize ((1/2)*quad_form(x,Q)+q'*x);
cvx_end
p_star=cvx_optval;

%%
% termination tolerance
tol = 1e-6;

% maximum number of allowed iterations
maxiter = 10000;

% minimum allowed perturbation
dxmin = 1e-6;

% step size r)
alpha = 0.2;
beta=0.5;
% initialize gradient norm, optimization vector, iteration counter, perturbation
gnorm = inf;
x = x0;
niter = 0;
Dx = inf;

% define the objective function:
f=@(x) (1/2)*x'*Q*x+q'*x;
fp=feval(f,x);


%% Gradient descent algorithm:

% and(and(gnorm>=tol,fp-p_star>=tol), and(niter <= maxiter, Dx >= dxmin))
% and(gnorm>=tol, and(niter <= maxiter, Dx >= dxmin))
while and(and(gnorm>=tol,fp-p_star>=tol), and(niter <= maxiter, Dx >= dxmin))
    % calculate gradient:
    g = (1/2)*(Q'+Q)*x+q;
    gnorm = norm(g);
    dx=-P\g;%diag(diag(Q))
%     take step: two approaches
%         t=exact_line(x,dx,Q,q);
    t=backtrack_linesearch(f,dx,x,beta,alpha);

    % update x
   xnew = x - t*g;
    % check step
    if ~isfinite(xnew)
        display(['Number of iterations: ' num2str(niter)])
        error('x is inf or NaN')
    end
    
    % plot current point for 2D
    %     plot([x(1) xnew(1)],[x(2) xnew(2)],'ko-')
    %     refresh
    
    % update termination metrics
    niter = niter + 1;
    Dx = norm(xnew-x);
    x = xnew;
    fp=feval(f,x);
    error(niter)=fp-p_star;
end
xopt = x;
fopt = f(xopt);
niter = niter - 1;
%% calculating condition number for quadratic optimization problem
e=eig(Q);
ymax=max(e);
ymin=min(e);
cond=ymax/ymin;