function uh = initL2f(master,mesh,icase,gam)
% This routine was written by Cristian Ciuca (cciuca13@mit.edu).
% It initializes the HDG MHD computation by doing an L2 projection
% of the conserved variables at the face level (for uhat). 
% Currently written for the Alfven polarised waves case (icase=4)
% and the iso-density magnetic vortex (icase=5).
% Input: master, mesh - required to compute the inverse of the 
% mass matrix and do the L2 projection;
%        icase - check which case to select proper IC 
%        gam - required by initial conditions (gamma)
% Output: the projected initial state of the vector of conserved 
% variables, i.e.
%   uh(:,1,:) :        rho
%   uh(:,2,:) :        rho*ux
%   uh(:,3,:) :        rho*uy
%   uh(:,4,:) :        rho*uz
%   uh(:,5,:) :        E
%   uh(:,6,:) :        Bx
%   uh(:,7,:) :        By
%   uh(:,8,:) :        By
%   uh(:,9,:) :        psi

npf = master.npf;
nf =mesh.nf;
nc = 7;
uh = zeros(npf,nc,nf);
if (icase == 4)
    uh = l2fprojection(mesh,master,@exactsolAlfven,[],[],9);  
elseif (icase ==5)
    uh = l2fprojection(mesh,master,@exactsolvortex,[],0,7); 
elseif (icase == 6)
    uh = l2fprojection(mesh,master,@initKH,[],[],7);
elseif (icase == 2)
    uh = l2fprojection(mesh,master,@initOT,[],[],7);
elseif (icase == 3)
    uh = l2fprojection(mesh,master,@initRotor,[],[],7);
end