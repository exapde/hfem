clear 
nd = 2;

syms u1 u2 u3 u4 u5 u6 u7 u8 u9 u10 u11 u12 u13 u14 u15 u16 u17 u18 u19 u20 u21 u22 u23 u24 u25 u26 u27
syms uh1 uh2 uh3 uh4 uh5 uh6 uh7 uh8 uh9
syms uinf1 uinf2 uinf3 uinf4 uinf5 uinf6 uinf7 uinf8 uinf9
syms x1 x2 x3 
syms nl1 nl2 nl3
syms time
syms param1 param2 param3

includegetan = 1;
appname = 'mhd';

if nd == 2    
    param = [param1 param2 param3];
    udg  = [u1 u2 u3 u4 u5 u6 u7 u8 u9 u10 u11 u12 u13 u14 u15 u16 u17 u18 u19 u20 u21 u22 u23 u24 u25 u26 u27];
    uh   = [uh1 uh2 uh3 uh4 uh5 uh6 uh7 uh8 uh9];
    uinf = [uinf1 uinf2 uinf3 uinf4 uinf5 uinf6 uinf7 uinf8 uinf9];
    pg   = [x1 x2];  
    nl   = [nl1 nl2];   
else
    disp('3D MHD q-u formulation not yet implemented');
    return;
%    param = [param1 param2 param3];
%    udg  = [u1 u2 u3 u4 u5 u6 u7 u8 u9];
%    uh   = [uh1 uh2 uh3 uh4 uh5 uh6 uh7 uh8 uh9];
%    uinf = [uinf1 uinf2 uinf3 uinf4 uinf5 uinf6 uinf7 uinf8 uinf9];
%    pg   = [x1 x2 x3];  
%    nl   = [nl1 nl2 nl3];    
end

gam  = param(1);
gam1 = param(1) - 1.0;        
tau  = param(2);
alpha1 = param(3);
alpha = 2.0;

ncu = length(uh);
nc  = length(udg);
nch = ncu;

if nd == 2                                               
    r    = udg(1);
    ru   = udg(2);
    rv   = udg(3);
    rw   = udg(4);
    rE   = udg(5);
    bx   = udg(6);
    by   = udg(7);
    bz   = udg(8);
    phi  = udg(9);
    
    rx   = udg(10);
    rux  = udg(11);
    rvx  = udg(12);
    rwx  = udg(13);
    rEx  = udg(14);
    bxx  = udg(15);
    byx  = udg(16);
    bzx  = udg(17);
    phix = udg(18);
    
    ry   = udg(19);
    ruy  = udg(20);
    rvy  = udg(21);
    rwy  = udg(22);
    rEy  = udg(23);
    bxy  = udg(24);
    byy  = udg(25);
    bzy  = udg(26);
    phiy = udg(27);
    
    r1   = 1/r;
    uv   = ru*r1;
    vv   = rv*r1;
    wv   = rw*r1;
    E    = rE*r1;
    q    = 0.5*(uv*uv+vv*vv+wv*wv);
    b    = (bx*bx + by*by + bz*bz)*0.5;
    p    = gam1*(rE-r*q-b-0.5*phi^2);
    h    = E + p*r1 + b*r1;
    uvb  = uv*bx + vv*by + wv*bz;
    divB = bxx+byy;

    f    = [ru, ru*uv+p+b-bx*bx, rv*uv-by*bx    , rw*uv-bz*bx, ru*h-uvb*bx+alpha1*phi*bx, alpha1*phi            , -(vv*bx-by*uv) , -(wv*bx-bz*uv) , alpha1*bx, ...
            rv, ru*vv-bx*by    , rv*vv+p+b-by*by, rw*vv-bz*by, rv*h-uvb*by+alpha1*phi*by, -(uv*by-bx*vv) , alpha1*phi            , -(wv*by-bz*vv) , alpha1*by];    
      
    rh    = uh(1);
    ruh   = uh(2);
    rvh   = uh(3);
    rwh   = uh(4);
    rEh   = uh(5);
    bxh   = uh(6);
    byh   = uh(7);
    bzh   = uh(8);
    phih  = uh(9);
    
    r1h   = 1/rh;
    uvh   = ruh*r1h;
    vvh   = rvh*r1h;
    wvh   = rwh*r1h;
    Eh    = rEh*r1h;
    bh    = (bxh*bxh + byh*byh + bzh*bzh)*0.5;
    qh    = 0.5*(uvh*uvh+vvh*vvh+wvh*wvh);
    ph    = gam1*(rEh-rh*qh-bh-0.5*phih^2);
    hh    = Eh+ph*r1h+bh*r1h;
    uvbh  = uvh*bxh + vvh*byh + wvh*bzh;
             
    fh    = [ruh, ruh*uvh+ph+bh-bxh*bxh, rvh*uvh-byh*bxh      , rwh*uvh-bzh*bxh, ruh*hh-uvbh*bxh+alpha1*phih*bxh, alpha1*phih            , -(vvh*bxh-byh*uvh), -(wvh*bxh-bzh*uvh) , alpha1*bxh, ...
             rvh, ruh*vvh-bxh*byh      , rvh*vvh+ph+bh-byh*byh, rwh*vvh-bzh*byh, rvh*hh-uvbh*byh+alpha1*phih*byh, -(uvh*byh-bxh*vvh), alpha1*phih            , -(wvh*byh-bzh*vvh) , alpha1*byh];  
      
    An = tau; 
    sh   = simplify(An*(udg(1:ncu).'-uh.'));    
    fhat = simplify(fh(1:ncu)*nl(1) + fh(ncu+1:end)*nl(2) + sh.');
    fh = fhat;
    
%    s     = -[0.0, divB*bx, divB*by, divB*bz, divB*uvb, divB*uv, divB*vv, divB*wv, alpha*phi];
    s     = -[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, alpha*phi];
    
    filename1 = ['test_flux_mhd23dq.c'];
    filename2 = ['test_source_mhd23dq.c'];
    filename3 = ['test_fhat_mhd23dq.c'];
end

genccode;
