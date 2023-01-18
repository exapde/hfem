function [UHAT,QHAT] = mkfhat(mesh,app,UDG)
%MKFHAT builds the hybrid UHAT and QHAT from the state UDG at each faces of
%the mesh in an explicit time integration context.
%This code works only for the scalar wave equation in its present state.
%
% SYNTAX:  [UHAT,QHAT] = MKFHAT(MESH,APP,UDG)
%
% INPUTS:
%    MESH - Mesh structure
%    APP  - Application Strructure
%    UDG  - Vector of volumetric unknowns
%
% OUTPUTS:
%    UHAT - Vector of hybrid unknowns on faces
%    QHAT - Vector of hybrid fluxes on faces
%
% NOTES: 
%    On the IMEX interface, UHAT is not computed but it is received from
%    the Implicit side, and then QHAT is computed.
%
% Author(s): Lauren Kolkman
% April 2018

global IMEXdata;

tau = app.arg{end};
shapfc = mkshape(mesh.porder,mesh.plocfc,mesh.plocfc,mesh.elemtype);
dshapft  = shapfc(:,:,2)';

nd = mesh.nd;
[npf,nfe] = size(mesh.perm);
[~,~,ne] = size(UDG);
perm = mesh.perm;
t2t = mkt2t(mesh.t,mesh.elemtype);
UHAT = zeros(npf,nfe,ne);
QHAT = zeros(npf,nfe,ne);
for i = 1:ne    
    e1 = i;    
    for j = 1:nfe % for each face of element i        
        e2 = t2t(i,j);
        u1 = UDG(perm(:,j),:,e1);        
        pb = mesh.dgnodes(perm(:,j),:,e1);
        dpg = reshape(dshapft*pb,[npf nd]);     
        jac = sqrt(dpg(:,1).^2+dpg(:,2).^2);
        nl   = [dpg(:,2)./jac,-dpg(:,1)./jac];
        if e2>0
            k  = mesh.t2f(e2,:)==mesh.t2f(e1,j);  % obtain the index of face i in the second element        
            u2 = UDG(perm(end:-1:1,k),:,e2);
            qn = sum((u1(:,2:end)-u2(:,2:end)).*nl,2);
            UHAT(:,j,i) = 0.5*(u1(:,1)+u2(:,1)) - (0.5/tau)*qn;            
            qn = sum(u1(:,2:end).*nl,2);
            QHAT(:,j,i) = qn - tau*(u1(:,1) - UHAT(:,j,i));
        else % Face on boundary
            faceid = mesh.t2f(e1,j);
            imexB = mesh.imexB; % Finds IMEX interface boundary number
            if mesh.f(faceid,end)==-imexB % The face is on the IMEX interface
                fid = find(mesh.fimex == faceid);
                UHAT(:,j,i) = IMEXdata(:,end:-1:1,fid);
                %Exact solution for debug purpose
                %exsol = exactsol(reshape(pb,[npf,nd,1]),app.time);
                %UHAT(:,j,i) = exsol(:,1,1);
            else % The face is on an external boundary
                UHAT(:,j,i) = 0;
            end
            qn = sum(u1(:,2:end).*nl,2);
            QHAT(:,j,i) = qn - tau*(u1(:,1) - UHAT(:,j,i));
        end
    end
end
UHAT = reshape(UHAT,[npf*nfe,ne]);
QHAT = reshape(QHAT,[npf*nfe,ne]);

