syms u1 u2 u3 pres q11 q21 q31 q12 q22 q32 q13 q23 q33
syms uh1 uh2 uh3 uinf1 uinf2 uinf3
syms x1 x2 x3 nl1 nl2 nl3
syms time
syms tau mu kappa 

param = [mu kappa tau];
if nd==2
    q = [q11 q12; q21 q22];    
    udg = [u1 u2 pres q11 q21 q12 q22];
    uh = [uh1 uh2];
    uinf = [uinf1 uinf2];
    pg = [x1 x2];  
    nl = [nl1 nl2];
    C = 1;
else
    q = -[q11 q12 q13; q21 q22 q23; q31 q32 q33];    
    udg = [u1 u2 u3 pres q11 q21 q31 q12 q22 q32 q13 q23 q33];
    uh = [uh1 uh2 uh3];
    uinf = [uinf1 uinf2 uinf3];
    pg  = [x1 x2 x3];
    nl  = [nl1 nl2 nl3];    
    C = 0;
end

nch = length(uh);

JF=det(q);
IF=inv(q);
TF=q.';
VF=IF.';
CF=q.'*q;

wcase=1;
switch (wcase)
    case 1
        W  = (mu/2)*(trace(CF)+C-3) - mu*log(JF) + (kappa/2)*(JF-1)^2;
        Wi = (mu/2)*(trace(CF)+C-3) - mu*log(JF);
        Wv = (kappa/2)*(JF-1)^2;
        Pv = pres*JF*VF;
        g  = JF-1; 
    case 2
        W  = (mu/2)*(trace(CF)+C-3) - mu*log(JF) + (kappa/2)*(log(JF))^2;
        Wi = (mu/2)*(trace(CF)+C-3) - mu*log(JF);
        Wv = (kappa/2)*(log(JF))^2;
        Pv = pres*VF;
        g  = log(JF);
    case 3
        W  = (mu/2)*(((JF^2)^(-1/3))*(trace(CF)+C)-3) + (kappa/2)*(JF-1)^2;
        Wi = (mu/2)*(((JF^2)^(-1/3))*(trace(CF)+C)-3);
        Wv = (kappa/2)*(JF-1)^2;
        Pv = pres*JF*VF;
        g  = JF-1; 
    case 4
        W  = (mu/2)*(JF^(-2/3)*(trace(CF)+C)-3) + (kappa/2)*(log(JF))^2;
        Wi = (mu/2)*(JF^(-2/3)*(trace(CF)+C)-3);
        Wv = (kappa/2)*(log(JF))^2;
        Pv = pres*VF;
        g  = log(JF);
    otherwise
        error('Not a valid model');
end

if nd==2    
    Pi = jacobian(Wi,udg);
    f = Pi(4:end).' + Pv(:);
        
    fh = [f(1)*nl1 + f(3)*nl2 + tau*(u1 - uh1);... 
          f(2)*nl1 + f(4)*nl2 + tau*(u2 - uh2)];
      
    %s = u1*u1 + u2*u2;
        
    fb{1} = [x1+uinf1 - uh1; x2+uinf2 - uh2];
    fb{2} = 'neumann';
    fb{3} = [x1+uinf1 - uh1; ...
             f(2)*nl1 + f(4)*nl2 + tau*(u2 - uh2) + uinf2];
    fb{4} = [f(1)*nl1 + f(3)*nl2 + tau*(u1 - uh1) + uinf1;
             x2+uinf2 - uh2];     
else    
    Pi = jacobian(Wi,udg);
    f = Pi(5:end).' - Pv(:);
    
    fh = [f(1)*nl1 + f(4)*nl2 + f(7)*nl3 + tau*(u1 - uh1);... 
          f(2)*nl1 + f(5)*nl2 + f(8)*nl3 + tau*(u2 - uh2);...
          f(3)*nl1 + f(6)*nl2 + f(9)*nl3 + tau*(u3 - uh3)];
      
    %s = u1*u1 + u2*u2 + u3*u3;  
    
    %fb{1} = 'dirichlet';
    fb{1} = [x1+uinf1 - uh1; x2+uinf2 - uh2; x3+uinf3 - uh3];
    fb{2} = 'neumann';
    fb{3} = [x1+uinf1 - uh1; ...
             f(2)*nl1 + f(5)*nl2 + f(8)*nl3 + tau*(u2 - uh2) - uinf2;...
             f(3)*nl1 + f(6)*nl2 + f(9)*nl3 + tau*(u3 - uh3) - uinf3];
    fb{4} = [f(1)*nl1 + f(4)*nl2 + f(7)*nl3 + tau*(u1 - uh1) - uinf1;... 
             x2+uinf2 - uh2;...
             f(3)*nl1 + f(6)*nl2 + f(9)*nl3 + tau*(u3 - uh3) - uinf3];
    fb{5} = [f(1)*nl1 + f(4)*nl2 + f(7)*nl3 + tau*(u1 - uh1) - uinf1;... 
             f(2)*nl1 + f(5)*nl2 + f(8)*nl3 + tau*(u2 - uh2) - uinf2;...
             x3+uinf3 - uh3];     
    fb{6} = [x1+uinf1 - uh1;...
             x2+uinf2 - uh2;...
             f(3)*nl1 + f(6)*nl2 + f(9)*nl3 + tau*(u3 - uh3) - uinf3];
    fb{7} = [f(1)*nl1 + f(4)*nl2 + f(7)*nl3 + tau*(u1 - uh1) - uinf1;... 
             x2+uinf2 - uh2;...
             x3+uinf3 - uh3];
    fb{8} = [x1+uinf1 - uh1;...
             f(2)*nl1 + f(5)*nl2 + f(8)*nl3 + tau*(u2 - uh2) - uinf2;...
             x3+uinf3 - uh3];          
end

filename1 = ['flux' num2str(nd) 'd' '.m'];
filename2 = ['source' num2str(nd) 'd'  '.m'];
filename3 = ['fhat' num2str(nd) 'd' '.m'];
filename4 = ['fbou' num2str(nd) 'd' '.m'];

% this function will generate matlab code 
genmatlabcode;




