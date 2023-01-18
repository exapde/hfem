syms u1 u2 u3 u4 u5 u6 u7 u8 u9 u10 u11 u12 u13 u14 u15 u16 u17 u18 u19 u20
syms uh1 uh2 uh3 uh4 uh5 
syms uinf1 uinf2 uinf3 uinf4 uinf5
syms x1 x2 x3 x4
syms nl1 nl2 nl3
syms time
syms param1 param2 param3 param4 param5 param6 param7

param = [param1 param2 param3 param4 param5 param6 param7];
if nd==2    
    udg = [u1 u2 u3 u4 u5 u6 u7 u8 u9 u10 u11 u12];
    uh = [uh1 uh2 uh3 uh4];
    uinf = [uinf1 uinf2 uinf3 uinf4];
    pg = [x1 x2 x3];  
    nl = [nl1 nl2];    
else
    udg = [u1 u2 u3 u4 u5 u6 u7 u8 u9 u10 u11 u12 u13 u14 u15 u16 u17 u18 u19 u20];
    uh = [uh1 uh2 uh3 uh4 uh5];
    uinf = [uinf1 uinf2 uinf3 uinf4 uinf5];
    pg = [x1 x2 x3 x4];  
    nl = [nl1 nl2 nl3];    
end

apar = 20;
bpar = 0.0;    

gam  = param(1);
gam1 = param(1) - 1.0;          
Re   = param(3);
Pr   = param(4);
Minf = param(5);
tau  = param(6);
epar = param(7);
%epar = param(8);
Re1  = 1/Re;
M2   = Minf^2;
c23  = 2.0/3.0;
fc   = 1/(gam1*M2*Re*Pr);

ncu = length(uh);
nc = length(udg);
nch = ncu;

if nd==2           
    hpar = pg(3);
    
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
    
    rHx = rEx+px;
    rHy = rEy+py;          
    Div = (ux+vy)*hpar;
    
    %av  = (epar/apar)*log(1+exp(apar*(Div-bpar)));       
    av  = epar*(Div.*(atan(100*Div)/pi + 1/2) + 1/2 - atan(100)/pi);
    f1  = [av*rx, av*rux, av*rvx, av*rHx,...
           av*ry, av*ruy, av*rvy, av*rHy];    
       
else       
    hpar = pg(4);
    
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
    
    rHx = rEx+px;
    rHy = rEy+py;      
    rHz = rEz+pz;
    Div = (ux+vy+wz)*hpar;
    
    %av = (epar/apar)*log(1+exp(apar*(Div-bpar)));    
    av  = epar*(Div.*(atan(100*Div)/pi + 1/2) + 1/2 - atan(100)/pi);
    f1 = [av*rx, av*rux, av*rvx, av*rwx, av*rHx,...
          av*ry, av*ruy, av*rvy, av*rwy, av*rHy,...
          av*rz, av*ruz, av*rvz, av*rwz, av*rHz];      
end

f = f1;
filename1 = ['avflux' num2str(nd) 'd' '.py'];
genpythoncode;

f = f1(:);
filename1 = ['avflux' num2str(nd) 'd' '.m'];
genmatlabcode;



