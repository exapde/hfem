syms F11 F21 F31 F12 F22 F32 F13 F23 F33
syms mu lambda 

if nd==2
    F = [F11 F12; F21 F22];    
else
    F = [F11 F12 F13; F21 F22 F23; F31 F32 F33];    
end

I = eye(nd);       % identity matrix
J = det(F);        % Jacobian of the deformation gradient
C = F.'*F;         % Cauchy-Green tensor
E = 0.5*(C - I);   % Green Tensor
e = 0.5*(F + F.'); % small strain tensor
D = inv(F).';

wcase=2;
switch (wcase)
    case 0 % linear elasticity model        
        % W = mu*trace(e.'*e) + 0.5*lambda*(trace(e))^2
        P = 2*mu*e + lambda*trace(e)*I;  
    case 1 % Saint Venant?Kirchhoff model      
        % W  = mu*trace(E.'*E) + 0.5*lambda*(trace(E))^2 
        P = F*(2*mu*E + lambda*trace(E)*I);
    case 2 % Neo-Hookean moedel
        % W  = (mu/2)*(trace(CF)+Cd-3) - mu*log(JF) + (lambda/2)*(log(JF))^2;    
        P = mu*F + (lambda*log(J)-mu)*D;
    otherwise
        error('Not a valid model');
end

f = -simplify(P(:));

%  ((F22*F33 - F23*F32)*(mu - lambda*log(F11*F22*F33 - F11*F23*F32 - F12*F21*F33 + F12*F23*F31 + F13*F21*F32 - F13*F22*F31)))/(F11*F22*F33 - F11*F23*F32 - F12*F21*F33 + F12*F23*F31 + F13*F21*F32 - F13*F22*F31) - F11*mu
%  - F21*mu - ((F12*F33 - F13*F32)*(mu - lambda*log(F11*F22*F33 - F11*F23*F32 - F12*F21*F33 + F12*F23*F31 + F13*F21*F32 - F13*F22*F31)))/(F11*F22*F33 - F11*F23*F32 - F12*F21*F33 + F12*F23*F31 + F13*F21*F32 - F13*F22*F31)
%    ((F12*F23 - F13*F22)*(mu - lambda*log(F11*F22*F33 - F11*F23*F32 - F12*F21*F33 + F12*F23*F31 + F13*F21*F32 - F13*F22*F31)))/(F11*F22*F33 - F11*F23*F32 - F12*F21*F33 + F12*F23*F31 + F13*F21*F32 - F13*F22*F31) - F31*mu
%  - F12*mu - ((F21*F33 - F23*F31)*(mu - lambda*log(F11*F22*F33 - F11*F23*F32 - F12*F21*F33 + F12*F23*F31 + F13*F21*F32 - F13*F22*F31)))/(F11*F22*F33 - F11*F23*F32 - F12*F21*F33 + F12*F23*F31 + F13*F21*F32 - F13*F22*F31)
%    ((F11*F33 - F13*F31)*(mu - lambda*log(F11*F22*F33 - F11*F23*F32 - F12*F21*F33 + F12*F23*F31 + F13*F21*F32 - F13*F22*F31)))/(F11*F22*F33 - F11*F23*F32 - F12*F21*F33 + F12*F23*F31 + F13*F21*F32 - F13*F22*F31) - F22*mu
%  - F32*mu - ((F11*F23 - F13*F21)*(mu - lambda*log(F11*F22*F33 - F11*F23*F32 - F12*F21*F33 + F12*F23*F31 + F13*F21*F32 - F13*F22*F31)))/(F11*F22*F33 - F11*F23*F32 - F12*F21*F33 + F12*F23*F31 + F13*F21*F32 - F13*F22*F31)
%    ((F21*F32 - F22*F31)*(mu - lambda*log(F11*F22*F33 - F11*F23*F32 - F12*F21*F33 + F12*F23*F31 + F13*F21*F32 - F13*F22*F31)))/(F11*F22*F33 - F11*F23*F32 - F12*F21*F33 + F12*F23*F31 + F13*F21*F32 - F13*F22*F31) - F13*mu
%  - F23*mu - ((F11*F32 - F12*F31)*(mu - lambda*log(F11*F22*F33 - F11*F23*F32 - F12*F21*F33 + F12*F23*F31 + F13*F21*F32 - F13*F22*F31)))/(F11*F22*F33 - F11*F23*F32 - F12*F21*F33 + F12*F23*F31 + F13*F21*F32 - F13*F22*F31)
%    ((F11*F22 - F12*F21)*(mu - lambda*log(F11*F22*F33 - F11*F23*F32 - F12*F21*F33 + F12*F23*F31 + F13*F21*F32 - F13*F22*F31)))/(F11*F22*F33 - F11*F23*F32 - F12*F21*F33 + F12*F23*F31 + F13*F21*F32 - F13*F22*F31) - F33*mu
 
%  mu*q11 - ((mu - lambda*log(q11*q23*q32 - q11*q22*q33 + q12*q21*q33 - q12*q23*q31 - q13*q21*q32 + q13*q22*q31))*(q22*q33 - q23*q32))/(q11*q22*q33 - q11*q23*q32 - q12*q21*q33 + q12*q23*q31 + q13*q21*q32 - q13*q22*q31)
%  mu*q21 + ((mu - lambda*log(q11*q23*q32 - q11*q22*q33 + q12*q21*q33 - q12*q23*q31 - q13*q21*q32 + q13*q22*q31))*(q12*q33 - q13*q32))/(q11*q22*q33 - q11*q23*q32 - q12*q21*q33 + q12*q23*q31 + q13*q21*q32 - q13*q22*q31)
%  mu*q31 - ((mu - lambda*log(q11*q23*q32 - q11*q22*q33 + q12*q21*q33 - q12*q23*q31 - q13*q21*q32 + q13*q22*q31))*(q12*q23 - q13*q22))/(q11*q22*q33 - q11*q23*q32 - q12*q21*q33 + q12*q23*q31 + q13*q21*q32 - q13*q22*q31)
%  mu*q12 + ((mu - lambda*log(q11*q23*q32 - q11*q22*q33 + q12*q21*q33 - q12*q23*q31 - q13*q21*q32 + q13*q22*q31))*(q21*q33 - q23*q31))/(q11*q22*q33 - q11*q23*q32 - q12*q21*q33 + q12*q23*q31 + q13*q21*q32 - q13*q22*q31)
%  mu*q22 - ((mu - lambda*log(q11*q23*q32 - q11*q22*q33 + q12*q21*q33 - q12*q23*q31 - q13*q21*q32 + q13*q22*q31))*(q11*q33 - q13*q31))/(q11*q22*q33 - q11*q23*q32 - q12*q21*q33 + q12*q23*q31 + q13*q21*q32 - q13*q22*q31)
%  mu*q32 + ((mu - lambda*log(q11*q23*q32 - q11*q22*q33 + q12*q21*q33 - q12*q23*q31 - q13*q21*q32 + q13*q22*q31))*(q11*q23 - q13*q21))/(q11*q22*q33 - q11*q23*q32 - q12*q21*q33 + q12*q23*q31 + q13*q21*q32 - q13*q22*q31)
%  mu*q13 - ((mu - lambda*log(q11*q23*q32 - q11*q22*q33 + q12*q21*q33 - q12*q23*q31 - q13*q21*q32 + q13*q22*q31))*(q21*q32 - q22*q31))/(q11*q22*q33 - q11*q23*q32 - q12*q21*q33 + q12*q23*q31 + q13*q21*q32 - q13*q22*q31)
%  mu*q23 + ((mu - lambda*log(q11*q23*q32 - q11*q22*q33 + q12*q21*q33 - q12*q23*q31 - q13*q21*q32 + q13*q22*q31))*(q11*q32 - q12*q31))/(q11*q22*q33 - q11*q23*q32 - q12*q21*q33 + q12*q23*q31 + q13*q21*q32 - q13*q22*q31)
%  mu*q33 - ((mu - lambda*log(q11*q23*q32 - q11*q22*q33 + q12*q21*q33 - q12*q23*q31 - q13*q21*q32 + q13*q22*q31))*(q11*q22 - q12*q21))/(q11*q22*q33 - q11*q23*q32 - q12*q21*q33 + q12*q23*q31 + q13*q21*q32 - q13*q22*q31)
 