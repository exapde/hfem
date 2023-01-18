syms u1 u2 u3 u4 u5 
syms uh1 uh2 uh3 uh4 uh5
syms uinf1 uinf2 uinf3 uinf4 uinf5
syms x1 x2 x3 
syms nl1 nl2 nl3
syms time
syms param1 param2 param3 param4 param5 param6

includegetan = 1;
appname = 'euler';

param = [param1 param2 param3 param4 param5 param6];
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

gam  = param(1);
gam1 = param(1) - 1.0;          
Re   = param(3);
Pr   = param(4);
Minf = param(5);
tau  = param(6);
Re1  = 1/Re;
M2   = Minf^2;
c23  = 2.0/3.0;
fc   = 1/(gam1*M2*Re*Pr);

ncu = length(uh);
nc = length(udg);
nch = ncu;

if nd==2                                               
    r    = udg(1);
    ru   = udg(2);
    rv   = udg(3);
    rE   = udg(4);

    r1   = 1/r;
    uv   = ru*r1;
    vv   = rv*r1;
    E    = rE*r1;
    q   = 0.5*(uv*uv+vv*vv);
    p    = gam1*(rE-r*q);
    h    = E+p*r1;
    
    f    = [ru, ru*uv+p, rv*uv,   ru*h, ...
            rv, ru*vv,   rv*vv+p, rv*h];
     
    rh    = uh(1);
    ruh   = uh(2);
    rvh   = uh(3);
    rEh   = uh(4);
    r1h   = 1/rh;
    uvh   = ruh*r1h;
    vvh   = rvh*r1h;
    Eh    = rEh*r1h;
    qh   = 0.5*(uvh*uvh+vvh*vvh);
    ph    = gam1*(rEh-rh*qh);
    hh    = Eh+ph*r1h;
    
    fh  = [ruh, ruh*uvh+ph, rvh*uvh,    ruh*hh, ...
           rvh, ruh*vvh,    rvh*vvh+ph, rvh*hh]; 
      
    An = tau;%symAn(uh,nl,param(1),2);
    sh = simplify(An*(udg(1:ncu).'-uh.'));
    fhat = simplify(fh(1:ncu)*nl(1) + fh(ncu+1:end)*nl(2) + sh.');
else                                         
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
    q   = 0.5*(uv*uv+vv*vv+wv*wv);
    p    = gam1*(rE-r*q);
    h    = E+p*r1;
    
    f  = [ru, ru*uv+p, rv*uv,   rw*uv,   ru*h,  ...
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
    qh   = 0.5*(uvh*uvh+vvh*vvh+wvh*wvh);
    ph    = gam1*(rEh-rh*qh);
    hh    = Eh+ph*r1h;
          
    fh = [ruh, ruh*uvh+ph, rvh*uvh,    rwh*uvh,    ruh*hh,  ...
           rvh, ruh*vvh,    rvh*vvh+ph, rwh*vvh,    rvh*hh, ...
           rwh, ruh*wvh,    rvh*wvh,    rwh*wvh+ph, rwh*hh];    
      
    An = tau;% symAn(uh,nl,param(1),2);
    sh = simplify(An*(udg(1:ncu).'-uh.'));
    fhat = simplify(fh(1:ncu)*nl(1) + fh(ncu+1:2*ncu)*nl(2) + fh(2*ncu+1:end)*nl(3) + sh.');   
end

filename1 = ['flux_' appname num2str(nd) 'd' '.c'];
filename2 = ['source_' appname num2str(nd) 'd'  '.c'];
filename3 = ['fhat_' appname num2str(nd) 'd' '.c'];

genccode; % generate source codes

% generate an application file
gid = fopen('fluxes.c','w');
if includegetan==1
    str = '#include "../getan.c"';
    fprintf(gid, '%s\n', str);    
end
fid = fopen('../../fluxes_template.c','r');    
tline = fgetl(fid); 
while ischar(tline)        
    str = tline;
    str = strrep(str, 'appname', appname);    
    fprintf(gid, '%s\n', str);    
    tline = fgetl(fid);            
end
fclose(fid);
fclose(gid);




