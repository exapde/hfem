d0=fileparts([pwd,filesep]);

syms u1 u2 u3 u4 u5 u6 u7 u8 u9 u10 u11 u12
syms uh1 uh2 uh3 uinf1 uinf2 uinf3
syms x1 x2 x3 nl1 nl2 nl3
syms time
syms param1 param2 param3 
ui = [uinf1 uinf2 uinf3];
param = [param1 param2 param3];
    
for wcase=0:0
    
switch (wcase)
    case 0 % linear elasticity model        
        mkdir([d0,'/LE_uq']);
        cd([d0,'/LE_uq']);        
    case 1 % Saint Venant?Kirchhoff model      
        mkdir([d0,'/SVK_uq']);
        cd([d0,'/SVK_uq']);        
    case 2 % Neo-Hookean moedel
        % W  = (mu/2)*(trace(CF)+Cd-3) - mu*log(JF) + (lambda/2)*(log(JF))^2;    
        mkdir([d0,'/NeoHookean_uq']);
        cd([d0,'/NeoHookean_uq']);        
    otherwise
        error('Not a valid model');
end    

for nd=2:3    
    
clear f fh fb;

if nd==2
    F = -[u3 u5; u4 u6];    
    udg = [u1 u2 u3 u4 u5 u6];
    uh = [uh1 uh2];
    uinf = [uinf1 uinf2];
    pg = [x1 x2];  
    nl = [nl1 nl2];
    Cd = 1;    
else
    F = -[u4 u7 u10; u5 u8 u11; u6 u9 u12];    
    udg = [u1 u2 u3 u4 u5 u6 u7 u8 u9 u10 u11 u12];
    uh = [uh1 uh2 uh3];
    uinf = [uinf1 uinf2 uinf3];
    pg  = [x1 x2 x3];
    nl  = [nl1 nl2 nl3];    
    Cd = 0;
end

nch = length(uh);

I = eye(nd);       % identity matrix
J = det(F);        % Jacobian of the deformation gradient
C = F.'*F;         % Cauchy-Green tensor
E = 0.5*(C - I);   % Green Tensor
e = 0.5*(F + F.' - 2*I); % small strain tensor
D = inv(F).';      % transpose of deformation inverse 

switch (wcase)
    case 0 % linear elasticity model        
        % W = mu*trace(e.'*e) + 0.5*lambda*(trace(e))^2
        P = 2*param1*e + param2*trace(F-I)*I;  
    case 1 % Saint Venant?Kirchhoff model      
        % W  = mu*trace(E.'*E) + 0.5*lambda*(trace(E))^2 
        P = F*(2*param1*E + param2*trace(E)*I);
    case 2 % Neo-Hookean moedel
        % W  = (mu/2)*(trace(CF)+Cd-3) - mu*log(JF) + (lambda/2)*(log(JF))^2;    
        P = param1*F + (param2*log(J)-param1)*D;
    otherwise
        error('Not a valid model');
end

f = -simplify(P(:));

if nd==2    
    fh = [f(1)*nl1 + f(3)*nl2 + param3*(u1 - uh1);... 
          f(2)*nl1 + f(4)*nl2 + param3*(u2 - uh2)];
      
    fb{1} = [x1+uinf1 - uh1; x2+uinf2 - uh2];
    fb{2} = [f(1)*nl1 + f(3)*nl2 + param3*(u1 - uh1) - uinf1;...
             f(2)*nl1 + f(4)*nl2 + param3*(u2 - uh2) - uinf2];
    fb{3} = [x1+uinf1 - uh1; ...
             f(2)*nl1 + f(4)*nl2 + param3*(u2 - uh2) - uinf2];
    fb{4} = [f(1)*nl1 + f(3)*nl2 + param3*(u1 - uh1) - uinf1;...
             x2+uinf2 - uh2];     
else    
    fh = [f(1)*nl1 + f(4)*nl2 + f(7)*nl3 + param3*(u1 - uh1);... 
          f(2)*nl1 + f(5)*nl2 + f(8)*nl3 + param3*(u2 - uh2);...
          f(3)*nl1 + f(6)*nl2 + f(9)*nl3 + param3*(u3 - uh3)];
          
    fb{1} = [x1+uinf1 - uh1; x2+uinf2 - uh2; x3+uinf3 - uh3];
    fb{2} = [f(1)*nl1 + f(4)*nl2 + f(7)*nl3 + param3*(u1 - uh1) - uinf1;... 
             f(2)*nl1 + f(5)*nl2 + f(8)*nl3 + param3*(u2 - uh2) - uinf2;...
             f(3)*nl1 + f(6)*nl2 + f(9)*nl3 + param3*(u3 - uh3) - uinf3];
    fb{3} = [x1+uinf1 - uh1; ...
             f(2)*nl1 + f(5)*nl2 + f(8)*nl3 + param3*(u2 - uh2) - uinf2;...
             f(3)*nl1 + f(6)*nl2 + f(9)*nl3 + param3*(u3 - uh3) - uinf3];
    fb{4} = [f(1)*nl1 + f(4)*nl2 + f(7)*nl3 + param3*(u1 - uh1) - uinf1;... 
             x2+uinf2 - uh2;...
             f(3)*nl1 + f(6)*nl2 + f(9)*nl3 + param3*(u3 - uh3) - uinf3];
    fb{5} = [f(1)*nl1 + f(4)*nl2 + f(7)*nl3 + param3*(u1 - uh1) - uinf1;... 
             f(2)*nl1 + f(5)*nl2 + f(8)*nl3 + param3*(u2 - uh2) - uinf2;...
             x3+uinf3 - uh3];     
    fb{6} = [x1+uinf1 - uh1;...
             x2+uinf2 - uh2;...
             f(3)*nl1 + f(6)*nl2 + f(9)*nl3 + param3*(u3 - uh3) - uinf3];
    fb{7} = [f(1)*nl1 + f(4)*nl2 + f(7)*nl3 + param3*(u1 - uh1) - uinf1;... 
             x2+uinf2 - uh2;...
             x3+uinf3 - uh3];
    fb{8} = [x1+uinf1 - uh1;...
             f(2)*nl1 + f(5)*nl2 + f(8)*nl3 + param3*(u2 - uh2) - uinf2;...
             x3+uinf3 - uh3];          
end

% this function will generate c code 
nc = nd*nd+nd;
ncu = nd;
appname = 'leuq';
filename1 = ['flux_leuq' num2str(nd) 'd' '.c'];
filename2 = ['source_leuq' num2str(nd) 'd'  '.c'];
filename3 = ['fhat_leuq' num2str(nd) 'd' '.c'];
filename4 = ['fbou_leuq' num2str(nd) 'd' '.c'];
fhat = fh;
genccode;
end
cd('..');
end



