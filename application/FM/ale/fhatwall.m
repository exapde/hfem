function [fh,fh_udg,fh_uh] = fhatwall(nl,p,udg,uh,param,time,signe)
%FHAT flux function
%   [fh,fhu,fhq,fhm] = fhat(nl,p,u,q,m,param)
%
%      NL(N,ND)              Normal N points
%      P(N,ND)               Coordinates for N points
%      U(N,NC)               Unknown vector for N points with NC components
%      Q(N,NC,ND)            Flux vector for N points with NC components in the
%                            coordinate directions
%      M(N,NC)               Hybrid unkowns for N points with NC components
%      PARAM                 Parameter list
%      FH(N,NC):              Volume flux at N points
%      FHU(N,NC,NC):         Jacobian of the flux flux vector w.r.t. U
%      FHQ(N,NC,NC,ND):      Jacobian of the flux flux vector w.r.t. Q
%      FHM(N,NC,NC):         Jacobian of the flux flux vector w.r.t. Lambda

[ng1,nch] = size(uh);

u = udg(:,1:nch);
umuh = (u-uh)';

[f,f_uhdg] = fluxwall(p,[uh,udg(:,nch+1:3*nch)],param,time);

[An,Anm] = getan(p,nl,uh,param,2,time);

n1d = spdiags(nl(:,1),0,size(nl,1),size(nl,1));
n2d = spdiags(nl(:,2),0,size(nl,1),size(nl,1));

%fh = multiprod(f(:,:,1),nl(:,1),2)+multiprod(f(:,:,2),nl(:,2),2) + multiprod(An,(u-uh),[2,3],2);
%fh = zeros(ng1,nch);

fh = n1d*f(:,:,1) + n2d*f(:,:,2);
temp = permute(An,[2,3,1]);
for ig = 1:ng1
    fh(ig,:) = fh(ig,:) + (temp(:,:,ig)*umuh(:,ig))';
end

% fh_uh = multiprod(permute(f_uhdg(:,:,1,1:nch),[1,2,4,3]),nl(:,1),[2,3],2) +   ...
%         multiprod(permute(f_uhdg(:,:,2,1:nch),[1,2,4,3]),nl(:,2),[2,3],2) + multiprod(Anm,(u-uh),[2,3],2) - An;

%n1d=full(n1d);n2d=full(n2d);
temp = reshape(n1d*reshape(f_uhdg(:,:,1,:),ng1,[]) + n2d*reshape(f_uhdg(:,:,2,:),ng1,[]),ng1,nch,[]);

fh_u = An;
fh_q = temp(:,:,nch+1:3*nch);
fh_uh = -An + temp(:,:,1:nch);

temp = reshape(permute(Anm,[2,4,3,1]),nch*nch,nch,ng1);
for ig = 1:ng1
    fh_uh(ig,:,:) = fh_uh(ig,:,:) + reshape(temp(:,:,ig)*umuh(:,ig),1,nch,nch);
end
    
fh_udg = cat(3,fh_u,fh_q);


function [f,f_UDG] = fluxwall(p,UDG,param,time)
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

[ng,nc] = size(UDG);
nch = nc/3;
nd  = 2;

p = mapping(p,time);

%V    = p(:,5:6);
G    = p(:,7:10);
g    = p(:,11);
gb   = p(:,12);
dgbdx = p(:,13:14);

g1  = 1./g;
gb1 = 1./gb;
gb2 = 1./(gb.*gb);

Ginv = [G(:,4).*g1 -G(:,2).*g1 -G(:,3).*g1 G(:,1).*g1];

qk = [UDG(:,5).*gb1   - UDG(:,1).*dgbdx(:,1).*gb2 ...
      UDG(:,6).*gb1   - UDG(:,2).*dgbdx(:,1).*gb2 ...
      UDG(:,7).*gb1   - UDG(:,3).*dgbdx(:,1).*gb2 ...
      UDG(:,8).*gb1   - UDG(:,4).*dgbdx(:,1).*gb2 ...
      UDG(:,9).*gb1   - UDG(:,1).*dgbdx(:,2).*gb2 ...
      UDG(:,10).*gb1  - UDG(:,2).*dgbdx(:,2).*gb2 ...
      UDG(:,11).*gb1  - UDG(:,3).*dgbdx(:,2).*gb2 ...
      UDG(:,12).*gb1  - UDG(:,4).*dgbdx(:,2).*gb2];

% udgm(:,1) = UDG(:,1).*gb1;
% udgm(:,2) = UDG(:,2).*gb1;
% udgm(:,3) = UDG(:,3).*gb1;
% udgm(:,4) = UDG(:,4).*gb1;
  
udgm = [UDG(:,1).*gb1 UDG(:,2).*gb1 UDG(:,3).*gb1 UDG(:,4).*gb1 ...
        qk(:,1).*Ginv(:,1)+qk(:,5).*Ginv(:,2) ...
        qk(:,2).*Ginv(:,1)+qk(:,6).*Ginv(:,2) ...
        qk(:,3).*Ginv(:,1)+qk(:,7).*Ginv(:,2) ...
        qk(:,4).*Ginv(:,1)+qk(:,8).*Ginv(:,2) ...
        qk(:,1).*Ginv(:,3)+qk(:,5).*Ginv(:,4) ...
        qk(:,2).*Ginv(:,3)+qk(:,6).*Ginv(:,4) ...
        qk(:,3).*Ginv(:,3)+qk(:,7).*Ginv(:,4) ...
        qk(:,4).*Ginv(:,3)+qk(:,8).*Ginv(:,4)];

udgm1_UDG1 = gb1;
% udgm2_UDG1 = 0;
% udgm3_UDG1 = 0;
% udgm4_UDG1 = 0;
udgm5_UDG1 = -(Ginv(:,1).*dgbdx(:,1).*gb2 + Ginv(:,2).*dgbdx(:,2).*gb2);
% udgm6_UDG1 = 0;
% udgm7_UDG1 = 0;
% udgm8_UDG1 = 0;
udgm9_UDG1 = -(Ginv(:,3).*dgbdx(:,1).*gb2 + Ginv(:,4).*dgbdx(:,2).*gb2);
% udgm10_UDG1 = 0;
% udgm11_UDG1 = 0;
% udgm12_UDG1 = 0;

% udgm1_UDG2 = 0;
udgm2_UDG2 = gb1;
% udgm3_UDG2 = 0;
% udgm4_UDG2 = 0;
% udgm5_UDG2 = 0; 
udgm6_UDG2 = -(Ginv(:,1).*dgbdx(:,1).*gb2 + Ginv(:,2).*dgbdx(:,2).*gb2);
% udgm7_UDG2 = 0;
% udgm8_UDG2 = 0;
% udgm9_UDG2 = 0;
udgm10_UDG2 = -(Ginv(:,3).*dgbdx(:,1).*gb2 + Ginv(:,4).*dgbdx(:,2).*gb2);
% udgm11_UDG2 = 0;
% udgm12_UDG2 = 0;

% udgm1_UDG3 = 0;
% udgm2_UDG3 = 0;
udgm3_UDG3 = gb1;
% udgm4_UDG3 = 0;
% udgm5_UDG3 = 0; 
% udgm6_UDG3 = 0;
udgm7_UDG3 = -(Ginv(:,1).*dgbdx(:,1).*gb2 + Ginv(:,2).*dgbdx(:,2).*gb2);
% udgm8_UDG3 = 0;
% udgm9_UDG3 = 0;
% udgm10_UDG3 = 0;
udgm11_UDG3 = -(Ginv(:,3).*dgbdx(:,1).*gb2 + Ginv(:,4).*dgbdx(:,2).*gb2);
% udgm12_UDG3 = 0;

% udgm1_UDG4 = 0;
% udgm2_UDG4 = 0;
% udgm3_UDG4 = 0;
udgm4_UDG4 = gb1;
% udgm5_UDG4 = 0; 
% udgm6_UDG4 = 0;
% udgm7_UDG4 = 0;
udgm8_UDG4 = -(Ginv(:,1).*dgbdx(:,1).*gb2 + Ginv(:,2).*dgbdx(:,2).*gb2);
% udgm9_UDG4 = 0;
% udgm10_UDG4 = 0;
% udgm11_UDG4 = 0;
udgm12_UDG4 = -(Ginv(:,3).*dgbdx(:,1).*gb2 + Ginv(:,4).*dgbdx(:,2).*gb2);

% udgm1_UDG5 = 0;
% udgm2_UDG5 = 0;
% udgm3_UDG5 = 0;
% udgm4_UDG5 = 0;
udgm5_UDG5 = Ginv(:,1).*gb1; 
% udgm6_UDG5 = 0;
% udgm7_UDG5 = 0;
% udgm8_UDG5 = 0;
udgm9_UDG5 = Ginv(:,3).*gb1;
% udgm10_UDG5 = 0;
% udgm11_UDG5 = 0;
% udgm12_UDG5 = 0;

% udgm1_UDG6 = 0;
% udgm2_UDG6 = 0;
% udgm3_UDG6 = 0;
% udgm4_UDG6 = 0;
% udgm5_UDG6 = 0;
udgm6_UDG6 = Ginv(:,1).*gb1; 
% udgm7_UDG6 = 0;
% udgm8_UDG6 = 0;
% udgm9_UDG6 = 0;
udgm10_UDG6 = Ginv(:,3).*gb1;
% udgm11_UDG6 = 0;
% udgm12_UDG6 = 0;

% udgm1_UDG7 = 0;
% udgm2_UDG7 = 0;
% udgm3_UDG7 = 0;
% udgm4_UDG7 = 0;
% udgm5_UDG7 = 0;
% udgm6_UDG7 = 0;
udgm7_UDG7 = Ginv(:,1).*gb1; 
% udgm8_UDG7 = 0;
% udgm9_UDG7 = 0;
% udgm10_UDG7 = 0;
udgm11_UDG7 = Ginv(:,3).*gb1;
% udgm12_UDG7 = 0;

% udgm1_UDG8 = 0;
% udgm2_UDG8 = 0;
% udgm3_UDG8 = 0;
% udgm4_UDG8 = 0;
% udgm5_UDG8 = 0;
% udgm6_UDG8 = 0;
% udgm7_UDG8 = 0;
udgm8_UDG8 = Ginv(:,1).*gb1; 
% udgm9_UDG8 = 0;
% udgm10_UDG8 = 0;
% udgm11_UDG8 = 0;
udgm12_UDG8 = Ginv(:,3).*gb1;

% udgm1_UDG9 = 0;
% udgm2_UDG9 = 0;
% udgm3_UDG9 = 0;
% udgm4_UDG9 = 0;
udgm5_UDG9 = Ginv(:,2).*gb1; 
% udgm6_UDG9 = 0;
% udgm7_UDG9 = 0;
% udgm8_UDG9 = 0;
udgm9_UDG9 = Ginv(:,4).*gb1;
% udgm10_UDG9 = 0;
% udgm11_UDG9 = 0;
% udgm12_UDG9 = 0;

% udgm1_UDG10 = 0;
% udgm2_UDG10 = 0;
% udgm3_UDG10 = 0;
% udgm4_UDG10 = 0;
% udgm5_UDG10 = 0;
udgm6_UDG10 = Ginv(:,2).*gb1; 
% udgm7_UDG10 = 0;
% udgm8_UDG10 = 0;
% udgm9_UDG10 = 0;
udgm10_UDG10 = Ginv(:,4).*gb1;
% udgm11_UDG10 = 0;
% udgm12_UDG10 = 0;

% udgm1_UDG11 = 0;
% udgm2_UDG11 = 0;
% udgm3_UDG11 = 0;
% udgm4_UDG11 = 0;
% udgm5_UDG11 = 0;
% udgm6_UDG11 = 0;
udgm7_UDG11 = Ginv(:,2).*gb1; 
% udgm8_UDG11 = 0;
% udgm9_UDG11 = 0;
% udgm10_UDG11 = 0;
udgm11_UDG11 = Ginv(:,4).*gb1;
% udgm12_UDG11 = 0;

% udgm1_UDG12 = 0;
% udgm2_UDG12 = 0;
% udgm3_UDG12 = 0;
% udgm4_UDG12 = 0;
% udgm5_UDG12 = 0;
% udgm6_UDG12 = 0;
% udgm7_UDG12 = 0;
udgm8_UDG12 = Ginv(:,2).*gb1; 
% udgm9_UDG12 = 0;
% udgm10_UDG12 = 0;
% udgm11_UDG12 = 0;
udgm12_UDG12 = Ginv(:,4).*gb1;

[fm,fm_udgm] = fluxwall_ns(p,udgm,param);

f(:,1,1) = g.*(fm(:,1,1).*Ginv(:,1) + fm(:,1,2).*Ginv(:,3));
f(:,2,1) = g.*(fm(:,2,1).*Ginv(:,1) + fm(:,2,2).*Ginv(:,3));
f(:,3,1) = g.*(fm(:,3,1).*Ginv(:,1) + fm(:,3,2).*Ginv(:,3));
f(:,4,1) = g.*(fm(:,4,1).*Ginv(:,1) + fm(:,4,2).*Ginv(:,3));

f(:,1,2) = g.*(fm(:,1,1).*Ginv(:,2) + fm(:,1,2).*Ginv(:,4));
f(:,2,2) = g.*(fm(:,2,1).*Ginv(:,2) + fm(:,2,2).*Ginv(:,4));
f(:,3,2) = g.*(fm(:,3,1).*Ginv(:,2) + fm(:,3,2).*Ginv(:,4));
f(:,4,2) = g.*(fm(:,4,1).*Ginv(:,2) + fm(:,4,2).*Ginv(:,4));

fn_UDG = zeros(ng,nch,nd,(nd+1)*nch);
for i1 = 1:nch
    for i2 = 1:nd
        fn_UDG(:,i1,i2,1) = fm_udgm(:,i1,i2,1).*udgm1_UDG1 + ...
                            fm_udgm(:,i1,i2,5).*udgm5_UDG1 + ...
                            fm_udgm(:,i1,i2,9).*udgm9_UDG1;

        fn_UDG(:,i1,i2,2) = fm_udgm(:,i1,i2,2).*udgm2_UDG2 + ...
                            fm_udgm(:,i1,i2,6).*udgm6_UDG2 + ...
                            fm_udgm(:,i1,i2,10).*udgm10_UDG2;

        fn_UDG(:,i1,i2,3) = fm_udgm(:,i1,i2,3).*udgm3_UDG3 + ...
                            fm_udgm(:,i1,i2,7).*udgm7_UDG3 + ...
                            fm_udgm(:,i1,i2,11).*udgm11_UDG3;

        fn_UDG(:,i1,i2,4) = fm_udgm(:,i1,i2,4).*udgm4_UDG4 + ...
                            fm_udgm(:,i1,i2,8).*udgm8_UDG4 + ...
                            fm_udgm(:,i1,i2,12).*udgm12_UDG4;

        fn_UDG(:,i1,i2,5) = fm_udgm(:,i1,i2,5).*udgm5_UDG5 + ...
                            fm_udgm(:,i1,i2,9).*udgm9_UDG5;

        fn_UDG(:,i1,i2,6) = fm_udgm(:,i1,i2,6).*udgm6_UDG6 + ...
                            fm_udgm(:,i1,i2,10).*udgm10_UDG6;

        fn_UDG(:,i1,i2,7) = fm_udgm(:,i1,i2,7).*udgm7_UDG7 + ...
                            fm_udgm(:,i1,i2,11).*udgm11_UDG7;
                        
        fn_UDG(:,i1,i2,8) = fm_udgm(:,i1,i2,8).*udgm8_UDG8 + ...
                            fm_udgm(:,i1,i2,12).*udgm12_UDG8;
                        
        fn_UDG(:,i1,i2,9) = fm_udgm(:,i1,i2,5).*udgm5_UDG9 + ...
                            fm_udgm(:,i1,i2,9).*udgm9_UDG9;                
                        
        fn_UDG(:,i1,i2,10) = fm_udgm(:,i1,i2,6).*udgm6_UDG10 + ...
                             fm_udgm(:,i1,i2,10).*udgm10_UDG10;                
                        
        fn_UDG(:,i1,i2,11) = fm_udgm(:,i1,i2,7).*udgm7_UDG11 + ...
                             fm_udgm(:,i1,i2,11).*udgm11_UDG11;                
                        
        fn_UDG(:,i1,i2,12) = fm_udgm(:,i1,i2,8).*udgm8_UDG12 + ...
                             fm_udgm(:,i1,i2,12).*udgm12_UDG12;                                
    end
end

f_UDG  = fn_UDG;
for i3=1:(nd+1)*nch    
    f_UDG(:,1,1,i3) = g.*(fn_UDG(:,1,1,i3).*Ginv(:,1) + fn_UDG(:,1,2,i3).*Ginv(:,3));
    f_UDG(:,2,1,i3) = g.*(fn_UDG(:,2,1,i3).*Ginv(:,1) + fn_UDG(:,2,2,i3).*Ginv(:,3));
    f_UDG(:,3,1,i3) = g.*(fn_UDG(:,3,1,i3).*Ginv(:,1) + fn_UDG(:,3,2,i3).*Ginv(:,3));
    f_UDG(:,4,1,i3) = g.*(fn_UDG(:,4,1,i3).*Ginv(:,1) + fn_UDG(:,4,2,i3).*Ginv(:,3));

    f_UDG(:,1,2,i3) = g.*(fn_UDG(:,1,1,i3).*Ginv(:,2) + fn_UDG(:,1,2,i3).*Ginv(:,4));
    f_UDG(:,2,2,i3) = g.*(fn_UDG(:,2,1,i3).*Ginv(:,2) + fn_UDG(:,2,2,i3).*Ginv(:,4));
    f_UDG(:,3,2,i3) = g.*(fn_UDG(:,3,1,i3).*Ginv(:,2) + fn_UDG(:,3,2,i3).*Ginv(:,4));
    f_UDG(:,4,2,i3) = g.*(fn_UDG(:,4,1,i3).*Ginv(:,2) + fn_UDG(:,4,2,i3).*Ginv(:,4));
end

function [f,f_udg] = fluxwall_ns(pg,udg,param)

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
fc = 1/(gam1*M2*Re*Pr);
           
V    = pg(:,5:6);
r    = udg(:,1);
ru   = udg(:,2);
rv   = udg(:,3);
rE   = udg(:,4);

rx   = udg(:,5);
rux  = udg(:,6);
rvx  = udg(:,7);
rEx  = udg(:,8);

ry   = udg(:,9);
ruy  = udg(:,10);
rvy  = udg(:,11);
rEy  = udg(:,12);

r1   = 1./r;
u    = ru.*r1;
v    = rv.*r1;
E    = rE.*r1;
q    = 0.5*(u.*u+v.*v);
p    = gam1*(rE-r.*q);
h    = E+p.*r1;

rh    = rE+p; %gam*rE-0.5*gam1*(ru^2+rv^2)/r 
rh_r  = 0.5*gam1*(u.^2+v.^2);
rh_ru = -gam1*u;
rh_rv = -gam1*v;
rh_rE = gam;
                                        
f = zeros(ng,nch,2);
f(:,:,1) = [ru-r.*V(:,1), ru.*u+p, rv.*u,   (ru-r.*V(:,1)).*h];
f(:,:,2) = [rv-r.*V(:,2), ru.*v,   rv.*v+p, (rv-r.*V(:,2)).*h];

f_u = zeros(ng,nch,2,nch);
f_u(:,:,1,1) = [zero-V(:,1), 0.5*((gam-3)*u.*u+gam1*v.*v), -u.*v,         2*gam1*u.*q-gam*E.*u-rh_r.*V(:,1)];
f_u(:,:,1,2) = [ one,                    (3-gam)*u,     v, gam*E-0.5*gam1*(3*u.*u+v.*v)-rh_ru.*V(:,1)];
f_u(:,:,1,3) = [zero,                      -gam1*v,     u,                   -gam1*u.*v-rh_rv.*V(:,1)];
f_u(:,:,1,4) = [zero,                     gam1*one,  zero,                        gam*u-rh_rE.*V(:,1)];

f_u(:,:,2,1) = [zero-V(:,2), -u.*v, 0.5*((gam-3)*v.*v+gam1*u.*u),         2*gam1*v.*q-gam*E.*v-rh_r.*V(:,2)];
f_u(:,:,2,2) = [zero,     v,                      -gam1*u,                   -gam1*u.*v-rh_ru.*V(:,2)];
f_u(:,:,2,3) = [ one,     u,                    (3-gam)*v, gam*E-0.5*gam1*(3*v.*v+u.*u)-rh_rv.*V(:,2)];
f_u(:,:,2,4) = [zero,  zero,                     gam1*one,                        gam*v-rh_rE.*V(:,2)];

u_r  = -u.*r1;
u_ru =  r1;
v_r  = -v.*r1;
v_rv =  r1;

ux  = (rux - rx.*u).*r1;
vx  = (rvx - rx.*v).*r1;
Ex  = (rEx - rx.*E).*r1;
qx  = u.*ux + v.*vx;
px  = gam1*(rEx - rx.*q - r.*qx);
Tx  = gam*M2*(px.*r - p.*rx).*r1.^2;

uy  = (ruy - ry.*u).*r1;
vy  = (rvy - ry.*v).*r1;
Ey  = (rEy - ry.*E).*r1;
qy  = u.*uy + v.*vy;
py  = gam1*(rEy - ry.*q - r.*qy);
Ty  = gam*M2*(py.*r - p.*ry).*r1.^2;

txx = Re1*c23*(2*ux - vy);
txy = Re1*(uy + vx);
tyy = Re1*c23*(2*vy - ux);

% Tx  = gam*M2*(px.*r - p.*rx).*r1.^2;
% Ty  = gam*M2*(py.*r - p.*ry).*r1.^2;

fv = zeros(ng,nch,2);
fv(:,:,1) = [zero, txx, txy, (u-V(:,1)).*txx + (v-V(:,2)).*txy + fc*Tx];
fv(:,:,2) = [zero, txy, tyy, (u-V(:,1)).*txy + (v-V(:,2)).*tyy + fc*Ty];

txx_r  =  Re1*c23*((4*ru.*rx-2*rv.*ry)-r.*(2*rux-rvy)).*r1.^3;
txx_ru = -Re1*c23*2*rx.*r1.^2;
txx_rv =  Re1*c23*ry.*r1.^2;
txx_rE =  zero;

txx_rx  = -Re1*c23*2*ru.*r1.^2;
txx_rux =  Re1*c23*2*r1;
txx_rvx =  zero;
txx_rEx =  zero;

txx_ry  =  Re1*c23*rv.*r1.^2;
txx_ruy =  zero;
txx_rvy = -Re1*c23.*r1;
txx_rEy =  zero;

txy_r  =  Re1*(2*(ru.*ry+rv.*rx)-r.*(ruy+rvx)).*r1.^3;
txy_ru = -Re1*ry.*r1.^2;
txy_rv = -Re1*rx.*r1.^2;
txy_rE =  zero;

txy_rx  = -Re1*rv.*r1.^2;
txy_rux =  zero;
txy_rvx =  Re1*r1;
txy_rEx =  zero;

txy_ry  = -Re1*ru.*r1.^2;
txy_ruy =  Re1*r1;
txy_rvy =  zero;
txy_rEy =  zero;

tyy_r  =  Re1*c23*((4*rv.*ry-2*ru.*rx)-r.*(2*rvy-rux)).*r1.^3;
tyy_ru =  Re1*c23*rx.*r1.^2;
tyy_rv = -Re1*c23*2*ry.*r1.^2;
tyy_rE =  zero;

tyy_rx  =  Re1*c23*ru.*r1.^2;
tyy_rux = -Re1*c23*r1;
tyy_rvx =  zero;
tyy_rEx =  zero;

tyy_ry  = -Re1*c23*2*rv.*r1.^2;
tyy_ruy =  zero;
tyy_rvy =  Re1*c23*2*r1;
tyy_rEy =  zero;

Tx_r  = -M2*gam*gam1*(rEx.*r.^2-2*rux.*r.*ru-2*rvx.*r.*rv-2*rE.*rx.*r+3*rx.*(ru.^2+rv.^2)).*r1.^4;
Tx_ru = -M2*gam*gam1*(r.*rux-2*ru.*rx).*r1.^3;
Tx_rv = -M2*gam*gam1*(r.*rvx-2*rv.*rx).*r1.^3;
Tx_rE = -M2*gam*gam1*rx.*r1.^2;

Tx_rx  =  M2*gam*gam1*(ru.^2+rv.^2-r.*rE).*r1.^3;
Tx_rux = -M2*gam*gam1*ru.*r1.^2;
Tx_rvx = -M2*gam*gam1*rv.*r1.^2;
Tx_rEx =  M2*gam*gam1*r1;

Ty_r  = -M2*gam*gam1*(rEy.*r.^2-2*ruy.*r.*ru-2*rvy.*r.*rv-2*rE.*ry.*r+3*ry.*(ru.^2+rv.^2)).*r1.^4;
Ty_ru = -M2*gam*gam1*(r.*ruy-2*ru.*ry).*r1.^3;
Ty_rv = -M2*gam*gam1*(r.*rvy-2*rv.*ry).*r1.^3;
Ty_rE = -M2*gam*gam1*ry.*r1.^2;

Ty_ry  = Tx_rx;
Ty_ruy = Tx_rux;
Ty_rvy = Tx_rvx;
Ty_rEy = Tx_rEx;

% fv(:,:,1) = [zero, txx, txy, u.*txx + v.*txy + fc*Tx];
% fv(:,:,2) = [zero, txy, tyy, u.*txy + v.*tyy + fc*Ty];

fv_u = zeros(ng,nch,2,nch);
fv_u(:,:,1,1) = [zero, txx_r , txy_r , u_r.*txx  + (u-V(:,1)).*txx_r  + v_r.*txy  + (v-V(:,2)).*txy_r  + fc*Tx_r ];
fv_u(:,:,1,2) = [zero, txx_ru, txy_ru, u_ru.*txx + (u-V(:,1)).*txx_ru             + (v-V(:,2)).*txy_ru + fc*Tx_ru];
fv_u(:,:,1,3) = [zero, txx_rv, txy_rv,             (u-V(:,1)).*txx_rv + v_rv.*txy + (v-V(:,2)).*txy_rv + fc*Tx_rv];
fv_u(:,:,1,4) = [zero, txx_rE, txy_rE,             (u-V(:,1)).*txx_rE             + (v-V(:,2)).*txy_rE + fc*Tx_rE];

fv_u(:,:,2,1) = [zero, txy_r , tyy_r , u_r.*txy  + (u-V(:,1)).*txy_r  + v_r.*tyy  + (v-V(:,2)).*tyy_r  + fc*Ty_r];
fv_u(:,:,2,2) = [zero, txy_ru, tyy_ru, u_ru.*txy + (u-V(:,1)).*txy_ru             + (v-V(:,2)).*tyy_ru + fc*Ty_ru];
fv_u(:,:,2,3) = [zero, txy_rv, tyy_rv,             (u-V(:,1)).*txy_rv + v_rv.*tyy + (v-V(:,2)).*tyy_rv + fc*Ty_rv];
fv_u(:,:,2,4) = [zero, txy_rE, tyy_rE,             (u-V(:,1)).*txy_rE             + (v-V(:,2)).*tyy_rE + fc*Ty_rE];

f_udg = zeros(ng,nch,2,nch,3);
f_udg(:,:,1,1,2) = [zero, txx_rx , txy_rx , (u-V(:,1)).*txx_rx  + (v-V(:,2)).*txy_rx  + fc*Tx_rx ];
f_udg(:,:,1,2,2) = [zero, txx_rux, txy_rux, (u-V(:,1)).*txx_rux + (v-V(:,2)).*txy_rux + fc*Tx_rux];
f_udg(:,:,1,3,2) = [zero, txx_rvx, txy_rvx, (u-V(:,1)).*txx_rvx + (v-V(:,2)).*txy_rvx + fc*Tx_rvx];
f_udg(:,:,1,4,2) = [zero, txx_rEx, txy_rEx, (u-V(:,1)).*txx_rEx + (v-V(:,2)).*txy_rEx + fc*Tx_rEx];
f_udg(:,:,1,1,3) = [zero, txx_ry , txy_ry , (u-V(:,1)).*txx_ry  + (v-V(:,2)).*txy_ry             ];
f_udg(:,:,1,2,3) = [zero, txx_ruy, txy_ruy, (u-V(:,1)).*txx_ruy + (v-V(:,2)).*txy_ruy            ];
f_udg(:,:,1,3,3) = [zero, txx_rvy, txy_rvy, (u-V(:,1)).*txx_rvy + (v-V(:,2)).*txy_rvy            ];
f_udg(:,:,1,4,3) = [zero, txx_rEy, txy_rEy, (u-V(:,1)).*txx_rEy + (v-V(:,2)).*txy_rEy            ];

f_udg(:,:,2,1,2) = [zero, txy_rx , tyy_rx , (u-V(:,1)).*txy_rx  + (v-V(:,2)).*tyy_rx             ];
f_udg(:,:,2,2,2) = [zero, txy_rux, tyy_rux, (u-V(:,1)).*txy_rux + (v-V(:,2)).*tyy_rux            ];
f_udg(:,:,2,3,2) = [zero, txy_rvx, tyy_rvx, (u-V(:,1)).*txy_rvx + (v-V(:,2)).*tyy_rvx            ];
f_udg(:,:,2,4,2) = [zero, txy_rEx, tyy_rEx, (u-V(:,1)).*txy_rEx + (v-V(:,2)).*tyy_rEx            ];
f_udg(:,:,2,1,3) = [zero, txy_ry , tyy_ry , (u-V(:,1)).*txy_ry  + (v-V(:,2)).*tyy_ry  + fc*Ty_ry ];
f_udg(:,:,2,2,3) = [zero, txy_ruy, tyy_ruy, (u-V(:,1)).*txy_ruy + (v-V(:,2)).*tyy_ruy + fc*Ty_ruy];
f_udg(:,:,2,3,3) = [zero, txy_rvy, tyy_rvy, (u-V(:,1)).*txy_rvy + (v-V(:,2)).*tyy_rvy + fc*Ty_rvy];
f_udg(:,:,2,4,3) = [zero, txy_rEy, tyy_rEy, (u-V(:,1)).*txy_rEy + (v-V(:,2)).*tyy_rEy + fc*Ty_rEy];

f = f+fv;
f_udg(:,:,:,:,1) = f_u+fv_u;
f_udg = reshape(f_udg,[ng nch 2 3*nch]);


