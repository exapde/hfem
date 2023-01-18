clear 
nd = 2;

syms u1 u2 u3 u4 u5 u6 u7 u8 u9
syms uh1 uh2 uh3 uh4 uh5 uh6 uh7 uh8 uh9
syms uinf1 uinf2 uinf3 uinf4 uinf5 uinf6 uinf7 uinf8 uinf9
syms x1 x2 x3 
syms nl1 nl2 nl3
syms time
syms param1 param2 param3 param4 param5 param6

includegetan = 1;
appname = 'mhd';

if nd == 23    
    param = [param1 param2];
    udg  = [u1 u2 u3 u4 u5 u6 u7];
    uh   = [uh1 uh2 uh3 uh4 uh5 uh6 uh7];
    uinf = [uinf1 uinf2 uinf3 uinf4 uinf5 uinf6 uinf7];
    pg   = [x1 x2];  
    nl   = [nl1 nl2];   
elseif nd == 2
    param = [param1 param2];
    udg  = [u1 u2 u3 u4 u5 u6 u7 u8];
    uh   = [uh1 uh2 uh3 uh4 uh5 uh6 uh7 uh8];
    uinf = [uinf1 uinf2 uinf3 uinf4 uinf5 uinf6 uinf7 uinf8];
    pg   = [x1 x2];  
    nl   = [nl1 nl2]; 
else
    param = [param1 param2];
    udg  = [u1 u2 u3 u4 u5 u6 u7 u8 u9];
    uh   = [uh1 uh2 uh3 uh4 uh5 uh6 uh7 uh8 uh9];
    uinf = [uinf1 uinf2 uinf3 uinf4 uinf5 uinf6 uinf7 uinf8 uinf9];
    pg   = [x1 x2 x3];  
    nl   = [nl1 nl2 nl3];    
end

gam  = param(1);
gam1 = param(1) - 1.0;        
tau  = param(2);

ncu = length(uh);
nc  = length(udg);
nch = ncu;

if nd == 23     

    
elseif nd == 2
    
    r    = udg(1);
    ru   = udg(2);
    rv   = udg(3);
    rw   = udg(4);
    rE   = udg(5);
    bx   = udg(6);
    by   = udg(7);
    bz   = udg(8);
    
    r1   = 1/r;
    uv   = ru*r1;
    vv   = rv*r1;
    wv   = rw*r1;
    E    = rE*r1;
    q    = 0.5*(uv*uv+vv*vv+wv*wv);
    b    = (bx*bx + by*by + bz*bz)*0.5;
    p    = gam1*(rE-r*q-b);
    h    = E + p*r1 + b*r1;
    uvb  = uv*bx + vv*by + wv*bz;

    f    = [ru, ru*uv+p+b-bx*bx, rv*uv-by*bx    , rw*uv-bz*bx, ru*h-uvb*bx, 0            , -(vv*bx-by*uv) , -(wv*bx-bz*uv) , ...
            rv, ru*vv-bx*by    , rv*vv+p+b-by*by, rw*vv-bz*by, rv*h-uvb*by, -(uv*by-bx*vv) , 0            , -(wv*by-bz*vv)];    
      
    rh    = uh(1);
    ruh   = uh(2);
    rvh   = uh(3);
    rwh   = uh(4);
    rEh   = uh(5);
    bxh   = uh(6);
    byh   = uh(7);
    bzh   = uh(8);
    
    r1h   = 1/rh;
    uvh   = ruh*r1h;
    vvh   = rvh*r1h;
    wvh   = rwh*r1h;
    Eh    = rEh*r1h;
    bh    = (bxh*bxh + byh*byh + bzh*bzh)*0.5;
    qh    = 0.5*(uvh*uvh+vvh*vvh+wvh*wvh);
    ph    = gam1*(rEh-rh*qh-bh);
    hh    = Eh+ph*r1h+bh*r1h;
    uvbh  = uvh*bxh + vvh*byh + wvh*bzh;
             
    fh    = [ruh, ruh*uvh+ph+bh-bxh*bxh, rvh*uvh-byh*bxh      , rwh*uvh-bzh*bxh, ruh*hh-uvbh*bxh, 0            , -(vvh*bxh-byh*uvh), -(wvh*bxh-bzh*uvh) , ...
             rvh, ruh*vvh-bxh*byh      , rvh*vvh+ph+bh-byh*byh, rwh*vvh-bzh*byh, rvh*hh-uvbh*byh, -(uvh*byh-bxh*vvh), 0            , -(wvh*byh-bzh*vvh) , ];  
      
    An = tau; 
    sh   = simplify(An*(udg(1:ncu).'-uh.'));    
    fhat = simplify(fh(1:ncu)*nl(1) + fh(ncu+1:end)*nl(2) + sh.');
    fh = fhat;
    
    filename1 = ['flux_mhd23di.c'];
    filename2 = ['source_mhd23di.c'];
    filename3 = ['fhat_mhd23di.c'];
    
else                                         


end


genccode; % generate source codes


