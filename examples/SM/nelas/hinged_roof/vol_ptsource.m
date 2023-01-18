function [ Ru ] = vol_ptsource( Ru, app, master )
%VOL_PTSOURCE : Computes volumetric poit source using the Dirac property
%   Detailed explanation goes here

for i=1:size(app.ptsource.els,1)
    xe = app.ptsource.els(i);
    xn = app.ptsource.node(i);
    nodalforce = reshape(app.ptsource.ampl*app.ptsource.dir(i,:),[1 1 3]);
    Ru(xn,xe,1:3) = Ru(xn,xe,1:3) + nodalforce;

%   Generalized case when the Point Load is not located on a Node
%     pts = app.ptsource.locpos(i,:);
%     ptforce = mkshape(master.porder,master.plocvl,pts,1);
%     ptforce = app.ptsource.ampl*squeeze(ptforce(:,:,1))*app.ptsource.dir(i,:);
%     Ru(:,xe,1:3) = Ru(:,xe,1:3)+reshape(ptforce,[master.npv 1 3]);
end

end

