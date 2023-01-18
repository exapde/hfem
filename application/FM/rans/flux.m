function [f,f_udg] = flux(pg,udg,param,time)
% FLUX Volume flux function
%   [f,fu,fq] = flux(p,u,q,param)
% 
%      P(N,ND)              Coordinates for N points
%      U(N,NC)              Unknown vector for N points with NC components
%      Q(N,NC,ND)           Flux vector for N points with NC components in the
%                           coordinate directions
%      PARAM                Parameter list
%      F(N,NC,ND):          Volume flux at N points
%      FU(N,NC,ND,NC):      Jacobian of the flux flux vector w.r.t. U
%      FQ(N,NC,ND,NC,ND):   Jacobian of the flux flux vector w.r.t. Q

[ng,nc] = size(udg);
nch = nc/3;
zero = zeros(ng,1);
one  = ones(ng,1);

gam  = param{1};
gam1 = gam-1.0;
Re   = param{3};
Re1  = 1/Re;
Pr   = param{4};
Minf = param{5};
M2   = Minf^2;

c23  = 2/3;
cv1  = 7.1;
sigm = 2/3;
%fc   = 1/(gam1*M2*Re*Pr);
                                             
r    = udg(:,1);
ru   = udg(:,2);
rv   = udg(:,3);
rE   = udg(:,4);
rN   = udg(:,5);

rx   = udg(:,6);
rux  = udg(:,7);
rvx  = udg(:,8);
rEx  = udg(:,9);
rNx  = udg(:,10);

ry   = udg(:,11);
ruy  = udg(:,12);
rvy  = udg(:,13);
rEy  = udg(:,14);
rNy  = udg(:,15);

r1   = 1./r;
u    = ru.*r1;
v    = rv.*r1;
E    = rE.*r1;
N    = rN.*r1;
q    = 0.5*(u.*u+v.*v);
p    = gam1*(rE-r.*q);
h    = E+p.*r1;
                                        
f = zeros(ng,nch,2);
f(:,:,1) = [ru, ru.*u+p, rv.*u,   ru.*h, rN.*u];
f(:,:,2) = [rv, ru.*v,   rv.*v+p, rv.*h, rN.*v];


f_u = zeros(ng,nch,2,nch);
f_u(:,:,1,1) = [zero, 0.5*((gam-3)*u.*u+gam1*v.*v), -u.*v,         2*gam1*u.*q-gam*E.*u, -N.*u];
f_u(:,:,1,2) = [ one,                    (3-gam)*u,     v, gam*E-0.5*gam1*(3*u.*u+v.*v),     N];
f_u(:,:,1,3) = [zero,                      -gam1*v,     u,                   -gam1*u.*v,  zero];
f_u(:,:,1,4) = [zero,                     gam1*one,  zero,                        gam*u,  zero];
f_u(:,:,1,5) = [zero,                         zero,  zero,                         zero,     u];

f_u(:,:,2,1) = [zero, -u.*v, 0.5*((gam-3)*v.*v+gam1*u.*u),         2*gam1*v.*q-gam*E.*v, -N.*v];
f_u(:,:,2,2) = [zero,     v,                      -gam1*u,                   -gam1*u.*v,  zero];
f_u(:,:,2,3) = [ one,     u,                    (3-gam)*v, gam*E-0.5*gam1*(3*v.*v+u.*u),     N];
f_u(:,:,2,4) = [zero,  zero,                     gam1*one,                        gam*v,  zero];
f_u(:,:,2,5) = [zero,  zero,                         zero,                         zero,     v];

u_r  = -u.*r1;
u_ru =  r1;
v_r  = -v.*r1;
v_rv =  r1;

ux  = (rux - rx.*u).*r1;
vx  = (rvx - rx.*v).*r1;
Ex  = (rEx - rx.*E).*r1;
Nx  = (rNx - rx.*N).*r1;
qx  = u.*ux + v.*vx;
px  = gam1*(rEx - rx.*q - r.*qx);
Tx  = gam*M2*(px.*r - p.*rx).*r1.^2;

uy  = (ruy - ry.*u).*r1;
vy  = (rvy - ry.*v).*r1;
Ey  = (rEy - ry.*E).*r1;
Ny  = (rNy - ry.*N).*r1;
qy  = u.*uy + v.*vy;
py  = gam1*(rEy - ry.*q - r.*qy);
Ty  = gam*M2*(py.*r - p.*ry).*r1.^2;

% chi = rN*Re;
% fv1 = chi.^3./(chi.^3+cv1^3);
% muN = rN.*fv1;
% in  = muN<0;
% muT = muN + Re1;
% muT(in) = Re1;
% fc  = muT/(gam1*M2*Pr);
% in  = rN<0;
% muS = rN + Re1;
% muS(in) = Re1;
% fs  = muS/sigm;

a = 0.0;
b = 20;
chi = rN*Re;
shi = log(1 + exp(b*(chi-a)))/b;
in  = b*(chi-a)>b;
shi(in) = chi(in)-a;
fv1 = shi.^3./(shi.^3+cv1^3);
muN = Re1*shi.*fv1;
muT = muN + Re1;
fc  = muT/(gam1*M2*Pr);
muS = Re1*shi + Re1;
fs  = muS/sigm;

txx = c23*muT.*(2*ux - vy);
txy = muT.*(uy + vx);
tyy = c23*muT.*(2*vy - ux);

fv = zeros(ng,nch,2);
fv(:,:,1) = [zero, txx, txy, u.*txx + v.*txy + fc.*Tx, fs.*Nx];
fv(:,:,2) = [zero, txy, tyy, u.*txy + v.*tyy + fc.*Ty, fs.*Ny];

% in = muN<0;
% muT_rN = (4*rN.^3*Re^3)./(rN.^3*Re^3 + cv1^3) - (3*rN.^6*Re^6)./(rN.^3*Re^3 + cv1^3).^2;
% muT_rN(in) = 0.0; %-((4*rN(in).^3*Re^3)./(rN(in).^3*Re^3 + cv1^3) - (3*rN(in).^6*Re^6)./(rN(in).^3*Re^3 + cv1^3).^2);
% in = rN<0;
% muS_rN = one;
% muS_rN(in) = 0.0; %-1.0;

chi_rN = Re;
shi_rN = chi_rN./(exp(b*(a - chi)).*(1./exp(b*(a - chi)) + 1));
shi_rN(in) = Re;
fv1_rN = shi_rN.*((3*shi.^2)./(cv1^3 + shi.^3) - (3*shi.^5)./(cv1^3 + shi.^3).^2);
muT_rN = Re1*(shi.*fv1_rN + shi_rN.*fv1);
muS_rN = Re1*shi_rN;

fc_rN  = muT_rN/(gam1*M2*Pr);
fs_rN  = muS_rN/sigm;

txx_r  =  c23*muT.*((4*ru.*rx-2*rv.*ry)-r.*(2*rux-rvy)).*r1.^3;
txx_ru = -c23*2*muT.*rx.*r1.^2;
txx_rv =  c23*muT.*ry.*r1.^2;
txx_rE =  zero;
txx_rN =  c23*muT_rN.*(2*ux - vy);

txx_rx  = -c23*2*muT.*ru.*r1.^2;
txx_rux =  c23*2*muT.*r1;
txx_rvx =  zero;
txx_rEx =  zero;
%txx_rNx =  zero;

txx_ry  =  c23*muT.*rv.*r1.^2;
txx_ruy =  zero;
txx_rvy = -c23*muT.*r1;
txx_rEy =  zero;
%txx_rNy =  zero;

txy_r  =  muT.*(2*(ru.*ry+rv.*rx)-r.*(ruy+rvx)).*r1.^3;
txy_ru = -muT.*ry.*r1.^2;
txy_rv = -muT.*rx.*r1.^2;
txy_rE =  zero;
txy_rN =  muT_rN.*(uy + vx);

txy_rx  = -muT.*rv.*r1.^2;
txy_rux =  zero;
txy_rvx =  muT.*r1;
txy_rEx =  zero;
%txy_rNx =  zero;

txy_ry  = -muT.*ru.*r1.^2;
txy_ruy =  muT.*r1;
txy_rvy =  zero;
txy_rEy =  zero;
%txy_rNy =  zero;

tyy_r  =  c23*muT.*((4*rv.*ry-2*ru.*rx)-r.*(2*rvy-rux)).*r1.^3;
tyy_ru =  c23*muT.*rx.*r1.^2;
tyy_rv = -c23*2*muT.*ry.*r1.^2;
tyy_rE =  zero;
tyy_rN =  c23*muT_rN.*(2*vy - ux);

tyy_rx  =  c23*muT.*ru.*r1.^2;
tyy_rux = -c23*muT.*r1;
tyy_rvx =  zero;
tyy_rEx =  zero;
%tyy_rNx =  zero;

tyy_ry  = -c23*2*muT.*rv.*r1.^2;
tyy_ruy =  zero;
tyy_rvy =  c23*2*muT.*r1;
tyy_rEy =  zero;
%tyy_rNy =  zero;

Tx_r  = -M2*gam*gam1*(rEx.*r.^2-2*rux.*r.*ru-2*rvx.*r.*rv-2*rE.*rx.*r+3*rx.*(ru.^2+rv.^2)).*r1.^4;
Tx_ru = -M2*gam*gam1*(r.*rux-2*ru.*rx).*r1.^3;
Tx_rv = -M2*gam*gam1*(r.*rvx-2*rv.*rx).*r1.^3;
Tx_rE = -M2*gam*gam1*rx.*r1.^2;
%Tx_rN =  zero;

Tx_rx  =  M2*gam*gam1*(ru.^2+rv.^2-r.*rE).*r1.^3;
Tx_rux = -M2*gam*gam1*ru.*r1.^2;
Tx_rvx = -M2*gam*gam1*rv.*r1.^2;
Tx_rEx =  M2*gam*gam1*r1;
%Tx_rNx =  zero;

Ty_r  = -M2*gam*gam1*(rEy.*r.^2-2*ruy.*r.*ru-2*rvy.*r.*rv-2*rE.*ry.*r+3*ry.*(ru.^2+rv.^2)).*r1.^4;
Ty_ru = -M2*gam*gam1*(r.*ruy-2*ru.*ry).*r1.^3;
Ty_rv = -M2*gam*gam1*(r.*rvy-2*rv.*ry).*r1.^3;
Ty_rE = -M2*gam*gam1*ry.*r1.^2;
%Ty_rN =  zero;

Ty_ry  = Tx_rx;
Ty_ruy = Tx_rux;
Ty_rvy = Tx_rvx;
Ty_rEy = Tx_rEx;
%Ty_rNy = zero;

%Nx  = (rNx - rx.*N).*r1 = (rNx - rx.*rN/r).*r1;
Nx_r   = -rNx.*r1.^2 + 2*rN.*rx.*r1.^3;
Nx_ru  = zero;
Nx_rv  = zero;
Nx_rE  = zero;
Nx_rN  = -rx.*r1.^2;

Nx_rx  = -rN.*r1.^2;
Nx_rux = zero;
Nx_rvx = zero;
Nx_rEx = zero;
Nx_rNx = r1;

%Ny  = (rNy - ry.*N).*r1 = (rNy - ry.*rN/r).*r1;
Ny_r   = -rNy.*r1.^2 + 2*rN.*ry.*r1.^3;
Ny_ru  = zero;
Ny_rv  = zero;
Ny_rE  = zero;
Ny_rN  = -ry.*r1.^2;

Ny_ry  = -rN.*r1.^2;
Ny_ruy = zero;
Ny_rvy = zero;
Ny_rEy = zero;
Ny_rNy = r1;


% fv(:,:,1) = [zero, txx, txy, u.*txx + v.*txy + fc.*Tx, fs.*Nx];
% fv(:,:,2) = [zero, txy, tyy, u.*txy + v.*tyy + fc.*Ty, fs.*Ny];

fv_u = zeros(ng,nch,2,nch);
fv_u(:,:,1,1) = [zero, txx_r , txy_r , u_r.*txx  + u.*txx_r  + v_r.*txy  + v.*txy_r  + fc.*Tx_r,  fs.*Nx_r];
fv_u(:,:,1,2) = [zero, txx_ru, txy_ru, u_ru.*txx + u.*txx_ru             + v.*txy_ru + fc.*Tx_ru, fs.*Nx_ru];
fv_u(:,:,1,3) = [zero, txx_rv, txy_rv,             u.*txx_rv + v_rv.*txy + v.*txy_rv + fc.*Tx_rv, fs.*Nx_rv];
fv_u(:,:,1,4) = [zero, txx_rE, txy_rE,             u.*txx_rE             + v.*txy_rE + fc.*Tx_rE, fs.*Nx_rE];
fv_u(:,:,1,5) = [zero, txx_rN, txy_rN,             u.*txx_rN             + v.*txy_rN + fc_rN.*Tx, fs.*Nx_rN + fs_rN.*Nx];

fv_u(:,:,2,1) = [zero, txy_r , tyy_r , u_r.*txy  + u.*txy_r  + v_r.*tyy  + v.*tyy_r  + fc.*Ty_r,  fs.*Ny_r];
fv_u(:,:,2,2) = [zero, txy_ru, tyy_ru, u_ru.*txy + u.*txy_ru             + v.*tyy_ru + fc.*Ty_ru, fs.*Ny_ru];
fv_u(:,:,2,3) = [zero, txy_rv, tyy_rv,             u.*txy_rv + v_rv.*tyy + v.*tyy_rv + fc.*Ty_rv, fs.*Ny_rv];
fv_u(:,:,2,4) = [zero, txy_rE, tyy_rE,             u.*txy_rE             + v.*tyy_rE + fc.*Ty_rE, fs.*Ny_rE];
fv_u(:,:,2,5) = [zero, txy_rN, tyy_rN,             u.*txy_rN             + v.*tyy_rN + fc_rN.*Ty, fs.*Ny_rN + fs_rN.*Ny];

fv_q = zeros(ng,nch,2,nch,2);
fv_q(:,:,1,1,1) = [zero, txx_rx , txy_rx , u.*txx_rx  + v.*txy_rx  + fc.*Tx_rx,  fs.*Nx_rx];
fv_q(:,:,1,2,1) = [zero, txx_rux, txy_rux, u.*txx_rux + v.*txy_rux + fc.*Tx_rux, fs.*Nx_rux];
fv_q(:,:,1,3,1) = [zero, txx_rvx, txy_rvx, u.*txx_rvx + v.*txy_rvx + fc.*Tx_rvx, fs.*Nx_rvx];
fv_q(:,:,1,4,1) = [zero, txx_rEx, txy_rEx, u.*txx_rEx + v.*txy_rEx + fc.*Tx_rEx, fs.*Nx_rEx];
fv_q(:,:,1,5,1) = [zero,    zero,    zero,                    zero             , fs.*Nx_rNx];
fv_q(:,:,1,1,2) = [zero, txx_ry , txy_ry , u.*txx_ry  + v.*txy_ry,  zero             ];
fv_q(:,:,1,2,2) = [zero, txx_ruy, txy_ruy, u.*txx_ruy + v.*txy_ruy, zero             ];
fv_q(:,:,1,3,2) = [zero, txx_rvy, txy_rvy, u.*txx_rvy + v.*txy_rvy, zero             ];
fv_q(:,:,1,4,2) = [zero, txx_rEy, txy_rEy, u.*txx_rEy + v.*txy_rEy, zero             ];
fv_q(:,:,1,5,2) = [zero,    zero,    zero,                    zero, zero             ];

fv_q(:,:,2,1,1) = [zero, txy_rx , tyy_rx , u.*txy_rx  + v.*tyy_rx,  zero             ];
fv_q(:,:,2,2,1) = [zero, txy_rux, tyy_rux, u.*txy_rux + v.*tyy_rux, zero             ];
fv_q(:,:,2,3,1) = [zero, txy_rvx, tyy_rvx, u.*txy_rvx + v.*tyy_rvx, zero             ];
fv_q(:,:,2,4,1) = [zero, txy_rEx, tyy_rEx, u.*txy_rEx + v.*tyy_rEx, zero             ];
fv_q(:,:,2,5,1) = [zero,    zero,    zero,                    zero, zero             ];
fv_q(:,:,2,1,2) = [zero, txy_ry , tyy_ry , u.*txy_ry  + v.*tyy_ry  + fc.*Ty_ry,  fs.*Ny_ry];
fv_q(:,:,2,2,2) = [zero, txy_ruy, tyy_ruy, u.*txy_ruy + v.*tyy_ruy + fc.*Ty_ruy, fs.*Ny_ruy];
fv_q(:,:,2,3,2) = [zero, txy_rvy, tyy_rvy, u.*txy_rvy + v.*tyy_rvy + fc.*Ty_rvy, fs.*Ny_rvy];
fv_q(:,:,2,4,2) = [zero, txy_rEy, tyy_rEy, u.*txy_rEy + v.*tyy_rEy + fc.*Ty_rEy, fs.*Ny_rEy];
fv_q(:,:,2,5,2) = [zero,    zero,    zero,                                 zero, fs.*Ny_rNy];

f = f+fv;
f_udg = cat(4,f_u+fv_u,reshape(fv_q,ng,nch,2,2*nch));