function t = backtrack_linesearch(f,dx,x,beta,alpha)
t = 1;
fk = feval(f,x);
xx = x;
x = x + t*dx;
fk1 = feval(f,x);
while fk1 > fk - alpha*t*(dx'*dx),
%   ff=fk - alpha*t*(dx'*dx);
  t = t*beta;
  x = xx + t*dx;
  fk1 = feval(f,x);
end