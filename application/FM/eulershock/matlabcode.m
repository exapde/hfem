syms u1 u2 u3 u4 u5 u6 u7 u8 u9 u10 u11 u12 u13 u14 u15 u16 u17 u18 u19 u20
syms uh1 uh2 uh3 uh4 uh5 
syms uinf1 uinf2 uinf3 uinf4 uinf5
syms x1 x2 x3 
syms nl1 nl2 nl3
syms time
syms param1 param2 param3 param4 param5 param6

param = [param1 param2 param3 param4 param5 param6];
if nd==2    
    udg = [u1 u2 u3 u4 u5 u6 u7 u8 u9 u10 u11 u12];
    uh = [uh1 uh2 uh3 uh4];
    uinf = [uinf1 uinf2 uinf3 uinf4];
    pg = [x1 x2];  
    nl = [nl1 nl2];    
else
    udg = [u1 u2 u3 u4 u5 u6 u7 u8 u9 u10 u11 u12 u13 u14 u15 u16 u17 u18 u19 u20];
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

    rx   = udg(5);
    rux  = udg(6);
    rvx  = udg(7);
    rEx  = udg(8);

    ry   = udg(9);
    ruy  = udg(10);
    rvy  = udg(11);
    rEy  = udg(12);
    
    r1   = 1/r;
    uv   = ru*r1;
    vv   = rv*r1;
    E    = rE*r1;
    q   = 0.5*(uv*uv+vv*vv);
    p    = gam1*(rE-r*q);
    h    = E+p*r1;
    
    fi   = [ru, ru*uv+p, rv*uv,   ru*h, ...
            rv, ru*vv,   rv*vv+p, rv*h];
     
    uv_r  = -uv*r1;
    uv_ru =  r1;
    vv_r  = -vv*r1;
    vv_rv =  r1;

    ux  = (rux - rx*uv)*r1;
    vx  = (rvx - rx*vv)*r1;    
    qx  = uv*ux + vv*vx;
    px  = gam1*(rEx - rx*q - r*qx);
    Tx  = gam*M2*(px*r - p*rx)*r1^2;

    uy  = (ruy - ry*uv)*r1;
    vy  = (rvy - ry*vv)*r1;    
    qy  = uv*uy + vv*vy;
    py  = gam1*(rEy - ry*q - r*qy);
    Ty  = gam*M2*(py*r - p*ry)*r1^2;

    txx = Re1*c23*(2*ux - vy);
    txy = Re1*(uy + vx);
    tyy = Re1*c23*(2*vy - ux);
    
    fv = [0, txx, txy, uv*txx + vv*txy + fc*Tx, ...
          0, txy, tyy, uv*txy + vv*tyy + fc*Ty];

    f = fi+fv;
    
    %fhat = f(uhat) * n + An(uhat)*(u-uhat) 
    
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
    
    fhi = [ruh, ruh*uvh+ph, rvh*uvh,    ruh*hh, ...
           rvh, ruh*vvh,    rvh*vvh+ph, rvh*hh]; 
      
    uvh_r  = -uvh*r1h;
    uvh_ru =  r1h;
    vvh_r  = -vvh*r1h;
    vvh_rv =  r1h;

    uxh  = (rux - rx*uvh)*r1h;
    vxh  = (rvx - rx*vvh)*r1h;    
    qxh  = uvh*uxh + vvh*vxh;
    pxh  = gam1*(rEx - rx*qh - rh*qxh);
    Txh  = gam*M2*(pxh*rh - ph*rx)*r1h^2;

    uyh  = (ruy - ry*uvh)*r1h;
    vyh  = (rvy - ry*vvh)*r1h;    
    qyh  = uvh*uyh + vvh*vyh;
    pyh  = gam1*(rEy - ry*qh - rh*qyh);
    Tyh  = gam*M2*(pyh*rh - ph*ry)*r1h^2;

    txxh = Re1*c23*(2*uxh - vyh);
    txyh = Re1*(uyh + vxh);
    tyyh = Re1*c23*(2*vyh - uxh);
    
    fhv = [0, txxh, txyh, uvh*txxh + vvh*txyh + fc*Txh, ...
           0, txyh, tyyh, uvh*txyh + vvh*tyyh + fc*Tyh];

    fh = fhi + fhv;
    
    An = tau;%symAn(uh,nl,param(1),2);
    sh = simplify(An*(udg(1:ncu).'-uh.'));
    fhat = simplify(fh(1:ncu)*nl(1) + fh(ncu+1:end)*nl(2) + sh.');
else                                         
    r    = udg(1);
    ru   = udg(2);
    rv   = udg(3);
    rw   = udg(4);
    rE   = udg(5);

    rx   = udg(6);
    rux  = udg(7);
    rvx  = udg(8);
    rwx  = udg(9);
    rEx  = udg(10);

    ry   = udg(11);
    ruy  = udg(12);
    rvy  = udg(13);
    rwy  = udg(14);
    rEy  = udg(15);
    
    rz   = udg(16);
    ruz  = udg(17);
    rvz  = udg(18);
    rwz  = udg(19);
    rEz  = udg(20);
    
    r1   = 1/r;
    uv   = ru*r1;
    vv   = rv*r1;
    wv   = rw*r1;
    E    = rE*r1;
    q   = 0.5*(uv*uv+vv*vv+wv*wv);
    p    = gam1*(rE-r*q);
    h    = E+p*r1;
    
    fi = [ru, ru*uv+p, rv*uv,   rw*uv,   ru*h,  ...
          rv, ru*vv,   rv*vv+p, rw*vv,   rv*h, ...
          rw, ru*wv,   rv*wv,   rw*wv+p, rw*h];    
      
    uv_r  = -uv*r1;
    uv_ru =  r1;
    vv_r  = -vv*r1;
    vv_rv =  r1;
    wv_r  = -wv*r1;
    wv_rv =  r1;

    ux  = (rux - rx*uv)*r1;
    vx  = (rvx - rx*vv)*r1;
    wx  = (rwx - rx*wv)*r1;
    qx  = uv*ux + vv*vx + wv*wx;
    px  = gam1*(rEx - rx*q - r*qx);
    Tx  = gam*M2*(px*r - p*rx)*r1^2;

    uy  = (ruy - ry*uv)*r1;
    vy  = (rvy - ry*vv)*r1;
    wy  = (rwy - ry*wv)*r1;
    qy  = uv*uy + vv*vy + wv*wy;
    py  = gam1*(rEy - ry*q - r*qy);
    Ty  = gam*M2*(py*r - p*ry)*r1^2;

    uz  = (ruz - rz*uv)*r1;
    vz  = (rvz - rz*vv)*r1;
    wz  = (rwz - rz*wv)*r1;
    qz  = uv*uz + vv*vz + wv*wz;
    pz  = gam1*(rEz - rz*q - r*qz);
    Tz  = gam*M2*(pz*r - p*rz)*r1^2;
    
    txx = Re1*c23*(2*ux - vy - wz);
    txy = Re1*(uy + vx);
    txz = Re1*(uz + wx);    
    tyy = Re1*c23*(2*vy - ux - wz);
    tyz = Re1*(vz + wy);
    tzz = Re1*c23*(2*wz - ux - vy);
    
    fv = [0, txx, txy, txz, uv*txx + vv*txy + wv*txz + fc*Tx, ...
          0, txy, tyy, tyz, uv*txy + vv*tyy + wv*tyz + fc*Ty,...
          0, txz, tyz, tzz, uv*txz + vv*tyz + wv*tzz + fc*Tz];

    f = fi+fv;         
     
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
          
    fhi = [ruh, ruh*uvh+ph, rvh*uvh,    rwh*uvh,    ruh*hh,  ...
           rvh, ruh*vvh,    rvh*vvh+ph, rwh*vvh,    rvh*hh, ...
           rwh, ruh*wvh,    rvh*wvh,    rwh*wvh+ph, rwh*hh];    
      
    uvh_r  = -uvh*r1h;
    uvh_ru =  r1h;
    vvh_r  = -vvh*r1h;
    vvh_rv =  r1h;
    wvh_r  = -wvh*r1h;
    wvh_rv =  r1h;

    uxh  = (rux - rx*uvh)*r1h;
    vxh  = (rvx - rx*vvh)*r1h;
    wxh  = (rwx - rx*wvh)*r1h;
    qxh  = uvh*uxh + vvh*vxh + wvh*wxh;
    pxh  = gam1*(rEx - rx*qh - rh*qxh);
    Txh  = gam*M2*(pxh*rh - ph*rx)*r1h^2;

    uyh  = (ruy - ry*uvh)*r1h;
    vyh  = (rvy - ry*vvh)*r1h;
    wyh  = (rwy - ry*wvh)*r1h;
    qyh  = uvh*uyh + vvh*vyh + wvh*wyh;
    pyh  = gam1*(rEy - ry*qh - rh*qyh);
    Tyh  = gam*M2*(pyh*rh - ph*ry)*r1h^2;

    uzh  = (ruz - rz*uvh)*r1h;
    vzh  = (rvz - rz*vvh)*r1h;
    wzh  = (rwz - rz*wvh)*r1h;
    qzh  = uvh*uzh + vvh*vzh + wvh*wzh;
    pzh  = gam1*(rEz - rz*qh - rh*qzh);
    Tzh  = gam*M2*(pzh*rh - ph*rz)*r1h^2;
    
    txxh = Re1*c23*(2*uxh - vyh - wzh);
    txyh = Re1*(uyh + vxh);
    txzh = Re1*(uzh + wxh);    
    tyyh = Re1*c23*(2*vyh - uxh - wzh);
    tyzh = Re1*(vzh + wyh);
    tzzh = Re1*c23*(2*wzh - uxh - vyh);
    
    fhv = [0, txxh, txyh, txzh, uvh*txxh + vvh*txyh + wvh*txzh + fc*Txh, ...
           0, txyh, tyyh, tyzh, uvh*txyh + vvh*tyyh + wvh*tyzh + fc*Tyh,...
           0, txzh, tyzh, tzzh, uvh*txzh + vvh*tyzh + wvh*tzzh + fc*Tzh];
    
    fh = fhi+fhv;
    
    An = tau;% symAn(uh,nl,param(1),2);
    sh = simplify(An*(udg(1:ncu).'-uh.'));
    fhat = simplify(fh(1:ncu)*nl(1) + fh(ncu+1:2*ncu)*nl(2) + fh(2*ncu+1:end)*nl(3) + sh.');
   
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

fh = fhat;
filename1 = ['flux' num2str(nd) 'd' '.m'];
filename2 = ['source' num2str(nd) 'd'  '.m'];
filename3 = ['fhat' num2str(nd) 'd' '.m'];
filename4 = ['fbou' num2str(nd) 'd' '.m'];
% 
% % this function will generate matlab code 
genmatlabcode;


