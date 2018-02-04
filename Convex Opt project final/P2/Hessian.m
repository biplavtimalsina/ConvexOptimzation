% define the Hesssian of the objective
function h = Hessian(x1,x2)
h = [exp(x1+3.*x2-0.1)+exp(x1-3.*x2-0.1)+exp(-x1-0.1) 3.*exp(x1+3.*x2-0.1)-3.*exp(x1-3.*x2-0.1)
    3.*exp(x1+3.*x2-0.1)-3.*exp(x1-3.*x2-0.1) 9.*exp(x1+3.*x2-0.1)+9.*exp(x1-3.*x2-0.1)];
end