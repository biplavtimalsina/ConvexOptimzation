% define the gradient of the objective
function g = grad(x1,x2)
g = [exp(x1+3.*x2-0.1)+exp(x1-3.*x2-0.1)-exp(-x1-0.1);3.*exp(x1+3.*x2-0.1)-3.*exp(x1-3.*x2-0.1)];
end