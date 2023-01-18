d0=fileparts([pwd,filesep]);

syms u1 u2 u3 u4 u5 u6 u7 u8 u9 u10 u11 u12
syms uh1 uh2 uh3 uinf1 uinf2 uinf3
syms x1 x2 x3 nl1 nl2 nl3
syms time
syms param1 param2 param3 
ui = [uinf1 uinf2 uinf3];
param = [param1 param2 param3];
    
  
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

% linear elasticity model (displacement/displacement gradient)               
P = param1*(F+F.') + param2*trace(F)*I;  
f = -simplify(P(:));
if nd==2    
    fh = [f(1)*nl1 + f(3)*nl2 + param3*(u1 - uh1);... 
          f(2)*nl1 + f(4)*nl2 + param3*(u2 - uh2)];
      
    fb{1} = [uinf1 - uh1; uinf2 - uh2];
    fb{2} = [f(1)*nl1 + f(3)*nl2 + param3*(u1 - uh1) - uinf1;...
             f(2)*nl1 + f(4)*nl2 + param3*(u2 - uh2) - uinf2];
    fb{3} = [uinf1 - uh1; ...
             f(2)*nl1 + f(4)*nl2 + param3*(u2 - uh2) - uinf2];
    fb{4} = [f(1)*nl1 + f(3)*nl2 + param3*(u1 - uh1) - uinf1;...
             uinf2 - uh2];     
else    
    fh = [f(1)*nl1 + f(4)*nl2 + f(7)*nl3 + param3*(u1 - uh1);... 
          f(2)*nl1 + f(5)*nl2 + f(8)*nl3 + param3*(u2 - uh2);...
          f(3)*nl1 + f(6)*nl2 + f(9)*nl3 + param3*(u3 - uh3)];
          
    fb{1} = [uinf1 - uh1; uinf2 - uh2; uinf3 - uh3];
    fb{2} = [f(1)*nl1 + f(4)*nl2 + f(7)*nl3 + param3*(u1 - uh1) - uinf1;... 
             f(2)*nl1 + f(5)*nl2 + f(8)*nl3 + param3*(u2 - uh2) - uinf2;...
             f(3)*nl1 + f(6)*nl2 + f(9)*nl3 + param3*(u3 - uh3) - uinf3];
    fb{3} = [uinf1 - uh1; ...
             f(2)*nl1 + f(5)*nl2 + f(8)*nl3 + param3*(u2 - uh2) - uinf2;...
             f(3)*nl1 + f(6)*nl2 + f(9)*nl3 + param3*(u3 - uh3) - uinf3];
    fb{4} = [f(1)*nl1 + f(4)*nl2 + f(7)*nl3 + param3*(u1 - uh1) - uinf1;... 
             uinf2 - uh2;...
             f(3)*nl1 + f(6)*nl2 + f(9)*nl3 + param3*(u3 - uh3) - uinf3];
    fb{5} = [f(1)*nl1 + f(4)*nl2 + f(7)*nl3 + param3*(u1 - uh1) - uinf1;... 
             f(2)*nl1 + f(5)*nl2 + f(8)*nl3 + param3*(u2 - uh2) - uinf2;...
             uinf3 - uh3];     
    fb{6} = [uinf1 - uh1;...
             uinf2 - uh2;...
             f(3)*nl1 + f(6)*nl2 + f(9)*nl3 + param3*(u3 - uh3) - uinf3];
    fb{7} = [f(1)*nl1 + f(4)*nl2 + f(7)*nl3 + param3*(u1 - uh1) - uinf1;... 
             uinf2 - uh2;...
             uinf3 - uh3];
    fb{8} = [uinf1 - uh1;...
             f(2)*nl1 + f(5)*nl2 + f(8)*nl3 + param3*(u2 - uh2) - uinf2;...
             uinf3 - uh3];          
end

filename1 = ['flux' num2str(nd) 'd' '.m'];
filename2 = ['source' num2str(nd) 'd'  '.m'];
filename3 = ['fhat' num2str(nd) 'd' '.m'];
filename4 = ['fbou' num2str(nd) 'd' '.m'];

% this function will generate matlab code 
genmatlabcode;

% % this function will generate c code 
% nc = nd*nd+nd;
% ncu = nd;
% appname = 'ledisp';
% filename1 = ['flux_ledisp' num2str(nd) 'd' '.c'];
% filename2 = ['source_ledisp' num2str(nd) 'd'  '.c'];
% filename3 = ['fhat_ledisp' num2str(nd) 'd' '.c'];
% filename4 = ['fbou_ledisp' num2str(nd) 'd' '.c'];
% fhat = fh;
% genccode;

end




