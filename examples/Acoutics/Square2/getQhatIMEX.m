function [IMEXQhat] = getQhatIMEX(master,mesh,app,UDG)
%GETQHATIMEX computes the hybrid flux Qhat from explicit variables on the IMEX
%boundary. This flux will later be used by the Implicit resolution as a
%Neuman boundary condition.
%This code works only for the scalar wave equation in its present state.
%
% SYNTAX:  [IMEXQHAT] = GETQHATIMEX(MASTER,MESH,APP,UDG)
%
% INPUTS:
%    MASTER    - Master element structure
%    MESH      - Mesh structure
%    APP       - Application Strructure
%    UDG       - Vector of volumetric unknowns
%
% OUTPUTS:
%    IMEXQhat  - Hybrid flux Qhat = Q.n -tau.u computed on the IMEX boundary. 
%
% SEE ALSO: GETQHATIMEX
% 
% Author(s): Lauren Kolkman, Sebastien Terrana
% April 2018

% Get some dimensions
nd = mesh.nd;
nch = app.nch;
npf = master.npf;
ngf = master.ngf;
% IMEX faces
fimex = mesh.fimex;
nfimex = length(fimex);
% Initialize field IMEXVhat
QHAT = zeros(nch,npf,nfimex);

tau = app.arg{end};
shapfc = mkshape(mesh.porder,mesh.plocfc,mesh.plocfc,mesh.elemtype);
dshapft  = shapfc(:,:,2)';
perm = mesh.perm;

for i = 1:nfimex
    fi = fimex(i);
    e1 = mesh.f(fi,end-1);
    j = find(mesh.t2f(e1,:)==fi);
    u1 = UDG(perm(:,j),:,e1);
    pb = mesh.dgnodes(perm(:,j),:,e1);
    dpg = reshape(dshapft*pb,[npf nd]);
    jac = sqrt(dpg(:,1).^2+dpg(:,2).^2);
    nl   = [dpg(:,2)./jac,-dpg(:,1)./jac];
    qn = sum(u1(:,2:end).*nl,2);
    %exsol = exactsol(reshape(pb,[ngf nd 1]),app.time);
    %QHAT(1,:,i) = exsol(:,2,1).* nl(:,1) + exsol(:,3,1).* nl(:,2);
    QHAT(1,:,i) = qn - tau*u1(:,1);
    QHAT(1,:,i) = QHAT(1,end:-1:1,i);
end

% Get Qhat values on the Gauss points
shapft = master.shapft(:,:,1);
QHAT  = reshape(QHAT,[npf nfimex*nch]);
QHAT = shapft*QHAT;
IMEXQhat = reshape(QHAT,[ngf*nfimex nch]);

end
