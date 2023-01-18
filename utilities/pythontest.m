function pythontest(master,mesh,app,UDG,UHAT,SH)


ng  = master.ngv;
nd  = master.nd;   
npf = size(mesh.perm,1);
nfe = size(mesh.perm,2);
ncu = app.ncu;
nc  = app.nc;
ngf = master.ngf;
ne  = size(mesh.dgnodes,3);
npv = size(mesh.dgnodes,1);

if isempty(SH)
    SH = zeros(npv,nc,ne);
end
%UHAT  = reshape(UHAT(:,mesh.elcon),[ncu npf nfe ne]);

[pge, minusJeInvDetJe, detJe] = volgeom(master,mesh.dgnodes);
[pgf, nlgf, detJf, XFace] = facegeom(master,mesh.dgnodes);
[F,F_udg,S,S_udg,udgg,pg] = volfluxes(master,app.arg,app.fc_u,app.time,app.tdep,UDG,SH(:,1:ncu,:),pge);
[FH, FH_udg, FH_uhat,udgf,uhatf] = facefluxes(master,app.arg,app.time,UDG,UHAT,pgf,nlgf);
[Rue, BDe] = volint(master,minusJeInvDetJe,detJe,F,F_udg,S,S_udg);
[Ru, BD, Rlu] = faceint(master,detJf,FH,FH_udg,FH_uhat,Rue,BDe);
[DUDG,UDG2,DUDG_DUHAT,Run,BDn,Rln] = elementsolve(master,app,mesh.dgnodes,SH,UHAT,UDG);
% [An,Anuh] = getan(reshape(nlgf,ngf*nfe*ne,[]),uhatf,app.arg,0);
% [Bn,Bnuh] = getan(reshape(nlgf,ngf*nfe*ne,[]),uhatf,app.arg,1);
[G, H, Rh, BH, BH_udg, BH_uhat, Rus] = uhatint(master,mesh.bf,app.bcs,app.bcm,app.arg,app.time,UHAT,UDG2-DUDG,mesh.dgnodes,Run,BDn,Rln);
[AE, FE] = elemmat(DUDG, DUDG_DUHAT, G, H, Rh);

F = reshape(F,[ng ne ncu nd]);
F_udg = reshape(F_udg,[ng ne ncu nd nc]);
S = reshape(S,[ng ne ncu]);
S_udg = reshape(S_udg,[ng ne ncu nc]);

FH = reshape(FH,[ngf nfe ne ncu]);
FH_udg = reshape(FH_udg,[ngf nfe ne ncu nc]);
FH_uhat = reshape(FH_uhat,[ngf nfe ne ncu ncu]);

% An = reshape(An,[ngf nfe ne ncu ncu]);
% Bn = reshape(Bn,[ngf nfe ne ncu ncu]);
% Anuh = reshape(Anuh,[ngf nfe ne ncu ncu ncu]);
% Bnuh = reshape(Bnuh,[ngf nfe ne ncu ncu ncu]);

BH = permute(reshape(BH,ngf,[],ncu),[2 1 3]);
BH_udg = permute(reshape(BH_udg,ngf,[],ncu,nc),[2 1 3 4]);
BH_uhat = permute(reshape(BH_uhat,ngf,[],ncu,ncu),[2 1 3 4]);

Rus = reshape(Rus,[npf,nfe,ne,ncu]);
% BH_udg = zeros(ngf*nbf,ncu,nc);
% BH_uh  = zeros(ngf*nbf,ncu,ncu);

save nstest.mat master mesh app UDG UHAT SH pge minusJeInvDetJe detJe pgf nlgf detJf XFace ...
     F F_udg S S_udg FH FH_udg FH_uhat Ru BD Rlu Run BDn Rln Rue BDe udgg pg udgf uhatf nlgf pgf ...
     DUDG DUDG_DUHAT G H Rh BH BH_udg BH_uhat Rus AE FE
  
% function master = pythonmaster(porder,pgauss,nd,elemtype,nodetype)
% 
% 
% master = masterelement(porder,pgauss,nd,elemtype,nodetype);
% dim    = master.nd;
% ngv    = master.ngv;
% ngf    = master.ngf;
% npv    = master.npv;
% npf    = master.npf;
% 
% master.shapvgdotshapvl  = zeros(npv*npv,ngv,dim+1);      
% for d=1:dim+1
%     master.shapvt(:,:,d) = master.shapvl(:,:,d)';
%     master.shapvg(:,:,d) = master.shapvl(:,:,d)*diag(master.gwvl);    
%     for ii=1:npv
%         for jj = 1:npv
%             master.shapvgdotshapvl((ii-1)*npv+jj,:,d) = master.shapvg(jj,:,d).*master.shapvl(ii,:,1);                    
%         end
%     end            
% end
% 
% % face shape functions and their derivatives 
% master.shapfgdotshapfc  = zeros(npf*npf,ngf,dim);      
% for d=1:dim
%     master.shapft(:,:,d) = master.shapfc(:,:,d)';
%     master.shapfg(:,:,d) = master.shapfc(:,:,d)*diag(master.gwfc);
%     for ii=1:npf
%         for jj = 1:npf
%             master.shapfgdotshapfc((ii-1)*npf+jj,:,d) = master.shapfg(jj,:,d).*master.shapfc(ii,:,1);                    
%         end
%     end            
% end
% 
% 
% 
