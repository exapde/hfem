function p=LagrangePolynomial(x,xn,f)
% The Lagrange interpolating polynomial is the polynomial P(x) of degree <=(n-1) 
% that passes through the n points (x_1,f(x_1)), (x_2,f(x_2)), ...,
% (x_n,f(x_n)), and is given by

% This function computes Lagrange polynomials at x when the pair (xn,f) are given 

n = length(x);
m = length(xn);
p = zeros(n,1);
for k=1:n
    for i=1:m
        y=xn;y(i)=[];
        sy=prod(xn(i)-y);
        sx=prod(x(k)-y);
        p(k)=p(k)+sx*f(i)/sy;
    end
end
