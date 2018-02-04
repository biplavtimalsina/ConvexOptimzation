clear all
clc

% define starting point
% x0 = round(3*rand(2,1));
x0=[-1;2];
% termination tolerance
tol = 1e-6;

% maximum number of allowed iterations
maxiter = 1000;

% minimum allowed perturbation
dxmin = 1e-6;

% step size ( beta=0.8 causes instability, beta=0.6 accurate with lowest iteration number)
alpha = 0.01;
beta=0.5;
% initialize gradient norm, optimization vector, iteration counter, perturbation
gnorm = inf;
x = x0;
niter = 0;
dx = inf;

% define the objective function:
f=@(x1,x2) exp(x1+3*x2-0.1)+exp(x1-3*x2-0.1)+exp(-x1-0.1);
% plot objective function contours for visualization for 2D:
figure(1); clf; fcontour(f); %axis equal; 
hold on

% redefine objective function syntax for use with optimization:
f2 = @(x) f(x(1),x(2));

% gradient descent algorithm:
while and(gnorm>=tol, and(niter <= maxiter, dx >= dxmin))
    % calculate gradient:
    g = grad(x(1),x(2));
    gnorm = norm(g);
    d=-g;
    % take step:
    t=backtrack_linesearch(f2,d,x,beta,alpha);
    xnew = x + t*d;
    % check step
    if ~isfinite(xnew)
        display(['Number of iterations: ' num2str(niter)])
        error('x is inf or NaN')
    end
    
    % plot current point for 2D
    plot([x(1) xnew(1)],[x(2) xnew(2)],'ko-')
    refresh

    % update termination metrics
    niter = niter + 1;
    dx = norm(xnew-x);
    x = xnew;
    
end
xopt = x;
fopt = f2(xopt);
niter = niter - 1;