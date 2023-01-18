syms u1 u2 u3 u4 u5
syms uh1 uh2 uh3 uh4 uh5 
syms uinf1 uinf2 uinf3 uinf4 uinf5
syms x1 x2 x3 
syms nl1 nl2 nl3
syms time
syms param1 param2

param = [param1,param2];
tau = param(2);
if nd==2    
    udg = [u1 u2 u3 u4];
    uh = [uh1 uh2 uh3 uh4];
    uinf = [uinf1 uinf2 uinf3 uinf4];
    pg = [x1 x2];  
    nl = [nl1 nl2];    
else
    udg = [u1 u2 u3 u4 u5];
    uh = [uh1 uh2 uh3 uh4 uh5];
    uinf = [uinf1 uinf2 uinf3 uinf4 uinf5];
    pg = [x1 x2 x3];  
    nl = [nl1 nl2 nl3];    
end

[An,Anm] = symAn(nl,uh,param(1),1);


