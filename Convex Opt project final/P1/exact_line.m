function t=exact_line(x,dx,Q,q)
t0=1;
ff=@(t) (1/2)*(x+t*dx)'*Q*(x+t*dx)+q'*(x+t*dx);
t=fminsearch(ff,t0);
if t>1
    t=1;
end
end