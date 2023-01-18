function [sr,dsr_dudg] = source(p,udg,param,time, app, master)
% function [sr,dsr_dudg] = source(p,udg,param,time, app, master)

[ng,nc] = size(udg);

if nc == 6
    nch = 2;
else 
    nch = 3;
end

sr = zeros(ng,nc);
dsr_dudg = zeros(ng,nch,nc);

return;

% ADDING POINT SOURCE CONTRIBUTION
% With a projection of the Dirac on the element's Legendre Basis.
% ngv = master.ngv;
% for i=1:size(app.ptsource.els,1)
%     n1 = app.ptsource.els(i);
%     pts = app.ptsource.pos(i,:);
%     [Lxo,fx,fy,fz] = tensorproduct(pts,master.porder);
%     [Lxp,fx,fy,fz] = tensorproduct(master.gpvl,master.porder);
%     ptforce = app.ptsource.ampl*Lxo*Lxp;
%     ptforce = ptforce' * app.ptsource.dir(i,:);
%     % TO DO  : EXPRIMER FORCE DANS REPERE GLOBAL
%     sr((n1-1)*ngv+1:n1*ngv,1:3) = ptforce;
%     %Ru(:,n1,1:3) = Ru(:,n1,1:3)+reshape(ptforce,[npv 1 ncu]);
% end
end

