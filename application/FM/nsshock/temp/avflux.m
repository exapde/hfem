function [f,f_udg] = avflux(pg,udg,param,time)

if nargout<2
    f = sensor(pg,udg,param);
    return;
end

gam  = param{1};
gam1 = gam-1.0;

[ng,nc] = size(udg);
nch = nc/3;

f_u = zeros(ng,nch);
f_q = zeros(ng,nch,2);
   
r    = udg(:,1);
ru   = udg(:,2);
rv   = udg(:,3);
rE   = udg(:,4);

rx   = udg(:,5);
rux  = udg(:,6);

ry   = udg(:,9);    
rvy  = udg(:,11);    

r1   = 1./r;
u    = ru.*r1;
v    = rv.*r1;

ux  = (rux - rx.*u).*r1;    

ux_r   = 2*rx.*ru.*r1.^3 - rux.*r1.^2;
ux_ru  = -rx.*r1.^2;
ux_rv  = 0*r;
ux_rE  = 0*r;

ux_rx  = -u.*r1;
ux_rux = r1;
ux_rvx = 0*r;
ux_rEx = 0*r;

ux_ry  = 0*r;
ux_ruy = 0*r;
ux_rvy = 0*r;
ux_rEy = 0*r;    

vy  = (rvy - ry.*v).*r1;            

vy_r  = 2*ry.*rv.*r1.^3 - rvy.*r1.^2;
vy_ru = 0*r;
vy_rv = -ry.*r1.^2;
vy_rE = 0*r;

vy_rx  = 0*r;
vy_rux = 0*r;
vy_rvx = 0*r;
vy_rEx = 0*r;

vy_ry  = -v.*r1;
vy_ruy = 0*r;
vy_rvy = r1;
vy_rEy = 0*r;        

Duv    = ux+vy;    
Duv_r  = ux_r+vy_r;    
Duv_ru = ux_ru+vy_ru;    
Duv_rv = ux_rv+vy_rv;    
Duv_rE = ux_rE+vy_rE;
 
Duv_rx = ux_rx+vy_rx;
Duv_rux= ux_rux+vy_rux;
Duv_rvx= ux_rvx+vy_rvx;
Duv_rEx= ux_rEx+vy_rEx;

Duv_ry = ux_ry+vy_ry;
Duv_ruy= ux_ruy+vy_ruy;
Duv_rvy= ux_rvy+vy_rvy;
Duv_rEy= ux_rEy+vy_rEy;


r1   = 1./r;
r1m1 = -1./(r.^2);
r1m2 = 0*r;
r1m3 = 0*r;
r1m4 = 0*r;

uv   = ru.*r1;
uvm1 =          ru.*r1m1;
uvm2 =     r1 + ru.*r1m2;
uvm3 =          ru.*r1m3;
uvm4 =          ru.*r1m4;

vv   = rv.*r1;
vvm1 =          rv.*r1m1;
vvm2 =          rv.*r1m2;
vvm3 =     r1 + rv.*r1m3;
vvm4 =          rv.*r1m4;

af   = 0.5*(uv.*uv+vv.*vv);
afm1 = uv.*uvm1 + vv.*vvm1;
afm2 = uv.*uvm2 + vv.*vvm2;
afm3 = uv.*uvm3 + vv.*vvm3;
afm4 = uv.*uvm4 + vv.*vvm4;

p    = gam1*(rE -r.*af);
pm1  = gam1*(   -   af - r.*afm1);
pm2  = gam1*(          - r.*afm2);
pm3  = gam1*(          - r.*afm3);
pm4  = gam1*(1       - r.*afm4);

c2   = gam* p.*r1;
c2m1 = gam*(pm1.*r1 + p.*r1m1);
c2m2 = gam*(pm2.*r1 + p.*r1m2);
c2m3 = gam*(pm3.*r1 + p.*r1m3);
c2m4 = gam*(pm4.*r1 + p.*r1m4);

if max(c2)<0.1
    cm = 0.1;  % minimum sound speed
    a1 = 2/(cm^2)
    a2 = cm^2;
    c3 = log(1 + exp(a1*(c2-a2)))/a1 + a2;
    g3 = exp(a1*(c3-a2))./(exp(a1*(c3-a2)) + 1);
    c3m1 = g3.*c2m1;
    c3m2 = g3.*c2m2;
    c3m3 = g3.*c2m3;
    c3m4 = g3.*c2m4;
else
    c3 = c2;
    c3m1 = c2m1;
    c3m2 = c2m2;
    c3m3 = c2m3;
    c3m4 = c2m4;
end

c    = sqrt(c3);
c_r  = 0.5*c3m1./c;
c_ru  = 0.5*c3m2./c;
c_rv  = 0.5*c3m3./c;
c_rE  = 0.5*c3m4./c;

hs = pg(:,3);
cmax0    = 1./c;
cmax1    = hs./(c);
cmax1_r  = -hs.*cmax0.^2.*(c_r);
cmax1_ru = -hs.*cmax0.^2.*(c_ru);
cmax1_rv = -hs.*cmax0.^2.*(c_rv);
cmax1_rE = -hs.*cmax0.^2.*(c_rE);

Dc1     = Duv.*cmax1;
Dc1_r   = Duv_r.*cmax1+Duv.*cmax1_r; 
Dc1_ru  = Duv_ru.*cmax1+Duv.*cmax1_ru; 
Dc1_rv  = Duv_rv.*cmax1+Duv.*cmax1_rv; 
Dc1_rE  = Duv_rE.*cmax1+Duv.*cmax1_rE;  
Dc1_rx  = Duv_rx.*cmax1;                    
Dc1_rux = Duv_rux.*cmax1;
Dc1_rvx = Duv_rvx.*cmax1;
Dc1_rEx = Duv_rEx.*cmax1;
Dc1_ry  = Duv_ry.*cmax1;                    
Dc1_ruy = Duv_ruy.*cmax1;
Dc1_rvy = Duv_rvy.*cmax1;
Dc1_rEy = Duv_rEy.*cmax1;    

% artificial viscosity parameters
a1=12;
a2=0.5;

f = log(1 + exp(a1*(Dc1-a2)))/a1;
g = exp(a1*(Dc1-a2))./(exp(a1*(Dc1-a2)) + 1);
in = a1*(Dc1-a2)>20;
f(in) = Dc1(in)-a2;
g(in) = 1;

f_u(:,1) = g.*Dc1_r; 
f_u(:,2) = g.*Dc1_ru; 
f_u(:,3) = g.*Dc1_rv; 
f_u(:,4) = g.*Dc1_rE; 

f_q(:,1,1) = g.*Dc1_rx;                    
f_q(:,2,1) = g.*Dc1_rux;
f_q(:,3,1) = g.*Dc1_rvx;
f_q(:,4,1) = g.*Dc1_rEx;

f_q(:,1,2) = g.*Dc1_ry;                    
f_q(:,2,2) = g.*Dc1_ruy;
f_q(:,3,2) = g.*Dc1_rvy;
f_q(:,4,2) = g.*Dc1_rEy;    

f_udg = cat(2,f_u,reshape(f_q,ng,2*nch));

f     = param{6}*f;
f_udg = param{6}*f_udg;


function f = sensor(pg, UDG, param)

gam  = 1.4;
gam1 = gam-1;

r    = UDG(:,1,:);
ru   = UDG(:,2,:);
rv   = UDG(:,3,:);
rE   = UDG(:,4,:);

rx   = UDG(:,5,:);
rux  = UDG(:,6,:);

ry   = UDG(:,9,:);    
rvy  = UDG(:,11,:);    

r1   = 1./r;
u    = ru.*r1;
v    = rv.*r1;

ux  = (rux - rx.*u).*r1;
vy  = (rvy - ry.*v).*r1;

af    = 0.5*(u.*u+v.*v);
Duv   = ux+vy; 
p     = gam1*(rE -r.*af);
c2    = gam* p.*r1;

cm = 0.1;  % minimum sound speed
a1 = 2/(cm^2);
a2 = cm^2;
c3 = log(1 + exp(a1*(c2-a2)))/a1 + a2;
in = a1*(c2-a2)>20;
c3(in) = c2(in);
c     = sqrt(c3);

hs    = pg(:,3,:);
cmax1 = hs./(c);
Dc1   = Duv.*cmax1;

a1 = 12;
a2 = 0.5;
f = log(1 + exp(a1*(Dc1-a2)))/a1;
in = a1*(Dc1-a2)>20;
f(in) = Dc1(in)-a2;
f = param{6}*f;








