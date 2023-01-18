function udg = initL2(master,mesh,icase,gam)
% This routine was written by Cristian Ciuca (cciuca13@mit.edu).
% It initializes the HDG MHD computation by doing an L2 projection
% of the conserved variables at the element level. 
% Currently written for the Alfven polarised waves case (icase=4)
% and the iso-density magnetic vortex (icase=5).
% Input: master, mesh - required to compute the inverse of the 
% mass matrix and do the L2 projection;
%        icase - check which case to select proper IC 
%        gam - required by initial conditions (gamma)
% Output: the projected initial state of the vector of conserved 
% variables, i.e.
%   udg(:,1,:) :        rho
%   udg(:,2,:) :        rho*ux
%   udg(:,3,:) :        rho*uy
%   udg(:,4,:) :        rho*uz
%   udg(:,5,:) :        E
%   udg(:,6,:) :        Bx
%   udg(:,7,:) :        By
%   udg(:,8,:) :        By
%   udg(:,9,:) :        psi

npv = size(mesh.dgnodes,1);
ne = size(mesh.dgnodes,3);
nc = 4;
udg = zeros(npv,nc,ne);
if (icase == 4)
    udg = l2eprojection(mesh,master,@exactsolAlfven,[],[],9);  
elseif (icase ==5)
    udg = l2eprojection(mesh,master,@exactsolvortex,[],0,9); 
elseif (icase == 6)
    udg = l2eprojection(mesh,master,@initKH,[],[],4);
end
