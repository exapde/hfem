syms u1 u2 u3 u4 u5
syms uh1 uh2 uh3 uh4 uh5 
syms uinf1 uinf2 uinf3 uinf4 uinf5
syms x1 x2 x3 
syms nl1 nl2 nl3
syms time
syms param1 param2 param3 param4

param = [param1,param2,param3,param4];
gam = param(1);
tau = param(3);
rmin = param(4);
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

ncu = length(uh);
nc = length(udg);
nch = ncu;

if nd==2        
    gam1 = gam - 1.0;                                             
    r    = udg(1);
    ru   = udg(2);
    rv   = udg(3);
    rE   = udg(4);
    
    rt   = (1/20)*log(1+exp(20*(r/rmin)))*rmin;
    r1   = 1/rt;
    uv   = ru*r1;
    vv   = rv*r1;
    E    = rE*r1;
    af   = 0.5*(uv*uv+vv*vv);
    p    = gam1*(rE-rt*af);
    h    = E+p*r1;
    
    f = [ru, ru*uv+p, rv*uv,   ru*h, ...
         rv, ru*vv,   rv*vv+p, rv*h];
     
    %fhat = f(uhat) * n + An(uhat)*(u-uhat) 
    
    rh    = uh(1);
    rth   = (1/20)*log(1+exp(20*(rh/rmin)))*rmin;
    ruh   = uh(2);
    rvh   = uh(3);
    rEh   = uh(4);
    r1h   = 1/rth;
    uvh   = ruh*r1h;
    vvh   = rvh*r1h;
    Eh    = rEh*r1h;
    afh   = 0.5*(uvh*uvh+vvh*vvh);
    ph    = gam1*(rEh-rth*afh);
    hh    = Eh+ph*r1h;
    
    fh = [ruh, ruh*uvh+ph, rvh*uvh,    ruh*hh, ...
          rvh, ruh*vvh,    rvh*vvh+ph, rvh*hh]; 
      
    An = tau;%symAn(uh,nl,param(1),2);
    sh = simplify(An*(udg.'-uh.'));
    fh = simplify(fh(1:ncu)*nl(1) + fh(ncu+1:end)*nl(2) + sh.');

%     % boundary fluxes
%     Ap = symAn(uh,nl,param(1),3);
%     Am = symAn(uh,nl,param(1),4);
%     fb{1} = simplify(Ap*(udg.'-uh.')) + simplify(Am*(uinf.'-uh.'));
%     
%     un = udg(2)*nl(1) + udg(3)*nl(2);        
%     uimf = udg;
%     uimf(2) = uimf(2) - nl(1)*un;
%     uimf(3) = uimf(3) - nl(2)*un;
%     fb{2} = (uimf - uh).';
%     
%     p  = uinf(1);
%     uimf = [udg(1), udg(2), udg(3), p/(param(1)-1) + 0.5*udg(2)*udg(2)/udg(1) + 0.5*udg(3)*udg(3)/udg(1)];
%     fb{3} = simplify(Ap*(udg.'-uh.')) + simplify(Am*(uimf.'-uh.'));
%         
%     fb{4} = simplify(An*(udg.'-uh.'));
else    
    gam1 = param(1) - 1.0;                                             
    r    = udg(1);
    ru   = udg(2);
    rv   = udg(3);
    rw   = udg(4);
    rE   = udg(5);

    r1   = 1/r;
    uv   = ru*r1;
    vv   = rv*r1;
    wv   = rw*r1;
    E    = rE*r1;
    af   = 0.5*(uv*uv+vv*vv+wv*wv);
    p    = gam1*(rE-r*af);
    h    = E+p*r1;
    
    f = [ru, ru*uv+p, rv*uv,   rw*uv,   ru*h,  ...
         rv, ru*vv,   rv*vv+p, rw*vv,   rv*h, ...
         rw, ru*wv,   rv*wv,   rw*wv+p, rw*h];    
     
    rh    = uh(1);
    ruh   = uh(2);
    rvh   = uh(3);
    rwh   = uh(4);
    rEh   = uh(5);
    r1h   = 1/rh;
    uvh   = ruh*r1h;
    vvh   = rvh*r1h;
    wvh   = rwh*r1h;
    Eh    = rEh*r1h;
    afh   = 0.5*(uvh*uvh+vvh*vvh+wvh*wvh);
    ph    = gam1*(rEh-rh*afh);
    hh    = Eh+ph*r1h;
          
    fh = [ruh, ruh*uvh+ph, rvh*uvh,    rwh*uvh,    ruh*hh,  ...
          rvh, ruh*vvh,    rvh*vvh+ph, rwh*vvh,    rvh*hh, ...
          rwh, ruh*wvh,    rvh*wvh,    rwh*wvh+ph, rwh*hh];    
      
    An = tau;% symAn(uh,nl,param(1),2);
    sh = simplify(An*(udg.'-uh.'));
    fh = simplify(fh(1:ncu)*nl(1) + fh(ncu+1:2*ncu)*nl(2) + fh(2*ncu+1:end)*nl(3) + sh.');
   
    % boundary fluxes
%     Ap = symAn(uh,nl,param(1),3);
%     Am = symAn(uh,nl,param(1),4);
%     fb{1} = simplify(Ap*(udg.'-uh.')) + simplify(Am*(uinf.'-uh.'));
%     
%     un = udg(2)*nl(1) + udg(3)*nl(2) + udg(4)*nl(3);        
%     uimf = udg;
%     uimf(2) = uimf(2) - nl(1)*un;
%     uimf(3) = uimf(3) - nl(2)*un;
%     uimf(4) = uimf(4) - nl(3)*un;
%     fb{2} = (uimf - uh).';
%     
%     p  = uinf(1);
%     uimf = [udg(1), udg(2), udg(3), udg(4), p/(gamma-1) + 0.5*udg(2)*udg(2)/udg(1) + 0.5*udg(3)*udg(3)/udg(1) + 0.5*udg(4)*udg(4)/udg(1)];
%     fb{3} = simplify(Ap*(udg.'-uh.')) + simplify(Am*(uimf.'-uh.'));
%         
%     fb{4} = simplify(An*(udg.'-uh.'));            
end


filename1 = ['flux' num2str(nd) 'd' '.m'];
filename2 = ['source' num2str(nd) 'd'  '.m'];
filename3 = ['fhat' num2str(nd) 'd' '.m'];
filename4 = ['fbou' num2str(nd) 'd' '.m'];

genmatlabcode;




