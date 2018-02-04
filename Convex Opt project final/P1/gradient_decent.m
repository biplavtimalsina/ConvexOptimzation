clear all
clc
% input dimension of Q matrix
% n=input('Enter the dimension of the Q:');
n=5;%n=25,50,100,200,500,1000 for size analysis

%% randomly generated matrix Q, vector q and initial point x0
% example 1 with n=2;
%  load ex1_n2 % generate Q,q,x0 with saved data
%  load ex2_n2_cond
%  P=randn(n);
%% generate PSD matrix with definite condition numebr
con_num=10; % con_num=2,5,10,25,50 ==> condition number is con_num
P=diag(diag(randi(con_num,[n,n])));
ma=randi(n);
mi=randi(n);
while mi==ma
    mi=randi(n);
end
P(ma,ma)=con_num;
P(mi,mi)=1;
%%

% make sure Q is PSD
Q=P;
q=randn(n,1);
x0=round(randn(n,1));
eigen_Q=eig(Q);
if (any(eigen_Q)<0)
    disp('Matrix Q is not a PSD matrix');
end

% example 2 with n=8; generatinig sparse symmetric matrix
% v = sparse([1,-.25,0,0,0.25,0,0,-.25]);
%  v= sparse(rand(1,n));
%  H = gallery('circul',v);
%  Q=H'*H;
% f=-4:3;

%% run cvx to find p_star
cvx_begin
variable x(n);
dual variable y;
minimize ((1/2)*quad_form(x,Q)+q'*x);
cvx_end
p_star=cvx_optval;

%% Parameters initialization

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

%% Gradient decent Algorithm

% and(and(gnorm>=tol,fp-p_star>=tol), and(niter <= maxiter, Dx >= dxmin))
% and(gnorm>=tol, and(niter <= maxiter, Dx >= dxmin))
while and(and(gnorm>=tol,fp-p_star>=tol), and(niter <= maxiter, Dx >= dxmin))
    % calculate gradient:
    g = (1/2)*(Q'+Q)*x+q;
    gnorm = norm(g);
    dx=-g;
    %     take step: two approaches
     %      t=exact_line(x,dx,Q,q);
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
plot(error)