function [UDGn,UHn,UDG0,PDG0] = hdg_solve_bdf(master,mesh,app,UDG,UH,UDG0,PDG0,time,dt,torder)
%HDG_SOLVE_BDF Solve using the HDG method and Newton's iteraion - BDF timestepping
%   [UH,QH,UHAT] = HDG_SOLVE_BDF(MASTER,MESH,UH,QH,UHAT,APP,TIME,DT,TORDER)
%
%      MASTER:                  Master structure
%      MESH:                    Mesh structure
%      UH(NPL,NC,NE):           Vector of unknowns (initial guess)
%      QH(NPL,NC,2,NE):         Vector of gradients of U (initial guess)
%      UHAT(NC,3*NPS,NF):       Vector of U_hat's (initial guess)
%      APP:                     Application structure
%      TIME:                    Time
%      DT:                      Timestep
%      TORDER:                  Order of accuracy
%
%      UH(NPL,NC,NE):           Vector of unknowns
%      QH(NPL,NC,2,NE):         Vector of gradients of U
%      UHAT(NC,3*NPS,NF):       Vector of U_hat's
%
%      NPL:                     Number of DG nodes within an element
%      NC:                      Number of conservation equations solved (components)
%      NPS:                     Number of HDG nodes per edge (porder+1)
%      NE:                      Number of elements
%      NF:                      Number of faces
%

a = bdfcoeff(torder);

nc = app.nc;
nch = app.nch;

fc_u = a(1)/dt;
if app.wave
    fc_q = fc_u;
else
    fc_q = 1;
end

SH = a(2)*UDG0(:,:,:,1);
for i = 2:torder
    SH = SH + a(i+1)*UDG0(:,:,:,i);
end
SH = -SH/dt;

app.time = time;
app.fc_u = fc_u;
app.fc_q = fc_q;

[UDGn,UHn] = hdg_solve(master,mesh,app,UDG,UH,SH);

for i = torder:-1:2
    UDG0(:,:,:,i) = UDG0(:,:,:,i-1);         
end
UDG0(:,:,:,1) = UDGn(:,:,:);

if app.wave
    SHP = a(2)*PDG0(:,:,:,1);
    for i = 2:torder
        SHP = SHP + a(i+1)*PDG0(:,:,:,i);
    end
    SHP = -SHP/dt;
        
    for i = torder:-1:2
        PDG0(:,:,:,i) = PDG0(:,:,:,i-1);        
    end
    PDG0(:,:,:,1) = (UDGn(:,1:nch,:)+SHP)/fc_u;
end

function a = bdfcoeff(p)
% p - order

if p == 1 
    a =        [  1,  -1];
elseif p == 2
    a =  (1/2)*[  3,  -4,  1];
elseif p == 3
    a =  (1/6)*[ 11, -18,  9,  -2];
elseif p == 4
    a = (1/12)*[ 25, -48, 36, -16, 3];
else
    error('Invalid q combination');
end