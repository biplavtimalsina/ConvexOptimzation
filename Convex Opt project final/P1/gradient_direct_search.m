 %function[xopt,error]= gradient_Descent(Q,q,n)
% input dimension of Q matrix
% n=input('Enter the dimension of the Q:');

 n=5;
 %n=25,50,100,200,500,1000 for size analysis

%% randomly generated matrix Q, vector q and initial point x0
% example 1 with n=2;
%  load ex1_n2 % generate Q,q,x0 with saved data
%  load ex2_n2_cond
 P=randn(n);
%% generate PSD matrix with definite condition numebr
con_num=2; % con_num=2,5,10,25,50 ==> condition number is con_num
P=diag(diag(randi(con_num,[n,n])));
ma=randi(n);
mi=randi(n);
while mi==ma
    mi=randi(n);
end
P(ma,ma)=con_num;
P(mi,mi)=1;
%%
%  Q=[3,0,0,0,0;0,10,0,0,0;0,0,1,0,0;0,0,0,7,0;0,0,0,0,1]
%  q=[0.877628740871287;1.03360564858667;0.419791320369650;0.601069227881255;-0.674018313717088]
% % make sure Q is PSD
 Q=P;
 q=randn(n,1);
% x0=[1;1;1;1;1]
%x0=[5;5;5;5;-15]
x0=round(randn(n,1));
eigen_Q=eig(Q);
if (any(eigen_Q)<0)
    disp('Matrix Q is not a PSD matrix');
end



%% run cvx to find p_opt
cvx_begin
variable x(n);
dual variable y;
minimize ((1/2)*quad_form(x,Q)+q'*x);

cvx_end
p_opt=cvx_optval;

%% Parameters initialization

% termination tolerance
tol = 1e-6;

% maximum number of allowed iterations
num_iter = 10000;

% minimum allowed perturbation
del_x_min = 1e-6;

% step size r)
alpha = 0.2;
beta=.5;
% initialize gradient norm, optimization vector, iteration counter, perturbation
norm_g = inf;
x = x0;
num_iter = 0;
Dx = inf;

% define the objective function:
f=@(x) (1/2)*x'*Q*x+q'*x;
fp=feval(f,x);

%% Gradient decent Algorithm

while( norm_g>=tol && fp-p_opt>=tol  && num_iter <= num_iter && Dx >= del_x_min)
    % calculate gradient:
    gradient = (1/2)*(Q'+Q)*x+q;
    norm_g = norm(gradient);
    dx=-gradient;
    
  %t=exact_line(x,dx,Q,q);
  t=backtrack_linesearch(f,dx,x,beta,alpha);
    
    % update x
    xnew = x - t*gradient;
    % check step
       
    % update termination metrics
    num_iter = num_iter + 1;
    Dx = norm(xnew-x);
    x = xnew;
    fp=feval(f,x);
    error(num_iter)=abs(fp-p_opt);
end
xopt = x;
fopt = f(xopt);
num_iter = num_iter - 1;
%% calculating condition number for quadratic optimization problem
e=eig(Q);
ymax=max(e);
ymin=min(e);
cond=ymax/ymin;
%%

plot(error)

% end