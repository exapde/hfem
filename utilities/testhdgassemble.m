function testhdgassemble(master,mesh,app,udg,uhat,uh,time,dt,nstage,torder)

fc_u = 1/(dt);
fc_q = 1;
app.itmax=4;
app.time = 0;
app.fc_u = fc_u;
app.fc_q = fc_q;
app.arg = {1.4,0.0,1.0};
sh = udg;    
for j=1:length(sh)
    sh{j} = sh{j}*fc_u;
end

[UDG,UH,UHAT] = hdg_solve_dirk(master,mesh,app,udg,uh,uhat,[],time,dt,nstage,torder);   

save testhdgsolvedirk.mat master mesh app udg uhat sh uh UDG UH UHAT time dt nstage torder 

% [UDG,UH,UHAT,AE,FE,DUDG_DUH,DUHAT_DUH,DUDGT,DUHAT,DUH,Un] = hdg_solve_dirk(master,mesh,app,udg,uh,uhat,[],time,dt,nstage,torder);   
% 
% save testhdgsolvedirk.mat master mesh app udg uhat sh uh UDG UH UHAT time dt nstage torder ...
%      AE FE DUDG_DUH DUHAT_DUH DUDGT DUHAT DUH Un

% [UDG,DUDG_DUH,AE,FE,UHAT,DUHAT_DUH] = hdg_assemble(master,mesh,app,udg,uh,sh,uhat);
% 
% ncu = size(uh,1);
% nc  = size(udg{1},2);
% ne  = size(mesh.dgnodes,3);
% nfe = size(mesh.perm,2);
% nqf = size(mesh.elcon,1)/nfe;
% is = find(mesh.subgrids==2);
% ut  = reshape(uh(:,mesh.elcon),[ncu nqf nfe ne]);
% dgnodes = mesh.dgnodes(:,:,is);
% ut = ut(:,:,:,is);
% [DUDG1,UDG1,DUDG1_DUH,DUHAT1,UHAT1,DUHAT1_DUH,G,H,R,BH,BH_udg,BH_uhat] = subgridsolve(master{2},app,dgnodes,sh{2},ut,udg{2},uhat{2});
% 
% ngf = size(BH,1);
% nbf = numel(BH)/(ngf*ncu);
% BH  = reshape(BH,ngf,nbf,ncu);
% BH_uhat = reshape(BH_uhat,ngf,nbf,ncu,ncu);
% BH_udg  = zeros(ngf,nbf,ncu,nc);
% 
% save testhdgassemble.mat master mesh app udg uhat sh uh ut UDG DUDG_DUH AE FE ...
%      UHAT DUHAT_DUH DUDG1 UDG1 DUDG1_DUH DUHAT1 UHAT1 DUHAT1_DUH G H R ...
%      BH BH_udg BH_uhat
% 

% ng  = master.ngv;
% nd  = master.nd;   
% npf = size(mesh.perm,1);
% nfe = size(mesh.perm,2);
% ncu = app.ncu;
% nc  = app.nc;
% ngf = master.ngf;
% ne  = size(mesh.dgnodes,3);
% npv = size(mesh.dgnodes,1);
% 
% if isempty(SH)
%     SH = zeros(npv,nc,ne);
% end
% %UHAT  = reshape(UHAT(:,mesh.elcon),[ncu npf nfe ne]);
% 
% [pge, minusJeInvDetJe, detJe] = volgeom(master,mesh.dgnodes);
% [pgf, nlgf, detJf, XFace] = facegeom(master,mesh.dgnodes);
% [F,F_udg,S,S_udg,udgg,pg] = volfluxes(master,app.arg,app.fc_u,app.time,app.tdep,UDG,SH(:,1:ncu,:),pge);
% [FH, FH_udg, FH_uhat,udgf,uhatf] = facefluxes(master,app.arg,app.time,UDG,UHAT,pgf,nlgf);
% [Rue, BDe] = volint(master,minusJeInvDetJe,detJe,F,F_udg,S,S_udg);
% [Ru, BD, Rlu] = faceint(master,detJf,FH,FH_udg,FH_uhat,Rue,BDe);
% [DUDG,UDG2,DUDG_DUHAT,Run,BDn,Rln] = elementsolve(master,app,mesh.dgnodes,SH,UHAT,UDG);
% % [An,Anuh] = getan(reshape(nlgf,ngf*nfe*ne,[]),uhatf,app.arg,0);
% % [Bn,Bnuh] = getan(reshape(nlgf,ngf*nfe*ne,[]),uhatf,app.arg,1);
% [G, H, Rh, BH, BH_udg, BH_uhat, Rus] = uhatint(master,mesh.bf,app.bcs,app.bcm,app.arg,app.time,UHAT,UDG2-DUDG,mesh.dgnodes,Run,BDn,Rln);
% [AE, FE] = elemmat(DUDG, DUDG_DUHAT, G, H, Rh);
% 
% F = reshape(F,[ng ne ncu nd]);
% F_udg = reshape(F_udg,[ng ne ncu nd nc]);
% S = reshape(S,[ng ne ncu]);
% S_udg = reshape(S_udg,[ng ne ncu nc]);
% 
% FH = reshape(FH,[ngf nfe ne ncu]);
% FH_udg = reshape(FH_udg,[ngf nfe ne ncu nc]);
% FH_uhat = reshape(FH_uhat,[ngf nfe ne ncu ncu]);
% 
% % An = reshape(An,[ngf nfe ne ncu ncu]);
% % Bn = reshape(Bn,[ngf nfe ne ncu ncu]);
% % Anuh = reshape(Anuh,[ngf nfe ne ncu ncu ncu]);
% % Bnuh = reshape(Bnuh,[ngf nfe ne ncu ncu ncu]);
% 
% BH = permute(reshape(BH,ngf,[],ncu),[2 1 3]);
% BH_udg = permute(reshape(BH_udg,ngf,[],ncu,nc),[2 1 3 4]);
% BH_uhat = permute(reshape(BH_uhat,ngf,[],ncu,ncu),[2 1 3 4]);
% 
% Rus = reshape(Rus,[npf,nfe,ne,ncu]);
% % BH_udg = zeros(ngf*nbf,ncu,nc);
% % BH_uh  = zeros(ngf*nbf,ncu,ncu);
% 
% save nstest.mat master mesh app UDG UHAT SH pge minusJeInvDetJe detJe pgf nlgf detJf XFace ...
%      F F_udg S S_udg FH FH_udg FH_uhat Ru BD Rlu Run BDn Rln Rue BDe udgg pg udgf uhatf nlgf pgf ...
%      DUDG DUDG_DUHAT G H Rh BH BH_udg BH_uhat Rus AE FE
%   
% function [UDGn,UHn,UHATn,AE,FE,DUDG_DUH,DUHAT_DUH,DUDGT,DUHAT,DUH,Un] = hdg_solve_dirk(master,mesh,app,UDG,UH,UHAT,PDG,time,dt,nstage,torder)
% 
% nsg = length(master);
% %nc = app.nc;
% ncu = app.ncu;
% 
% Un = UDG;
% for j=1:nsg    
%     [npv,nc,ne] = size(UDG{j});    
%      F{j}  = zeros(npv,nc,ne,nstage);    
% end
% 
% if app.wave   
%     for j=1:nsg                
%         [npv,nc,ne] = size(UDG{j});    
%         G{j}  = zeros(npv,nc,ne,nstage);        
%     end    
%     Pn  = PDG;
% end
% 
% [a,b,t] = dirkcoeff(nstage,torder);
% 
% for istage = 1:nstage
%     fprintf('DIRK stage :  %d\n', istage);
%     fc_u = 1/(dt*a(istage,istage));
%     if app.wave
%         fc_q = fc_u;
%     else
%         fc_q = 1;
%     end
% 
%     app.time = time+dt*t(istage);
%     app.fc_u = fc_u;
%     app.fc_q = fc_q;
%     
%     SH = UDG;    
%     for j=1:nsg
%         for jstage = 1:istage-1
%             SH{j} = SH{j} + dt*a(istage,jstage)*F{j}(:,:,:,jstage);
%         end
%         SH{j} = SH{j}*fc_u;
%     end
%             
%     [UDGn,UHn,UHATn,AE,FE,DUDG_DUH,DUHAT_DUH,DUDGT,DUHAT,DUH] = hdg_solve(master,mesh,app,UDG,UH,SH,UHAT);
%     
%     for j=1:nsg        
%         F{j}(:,:,:,istage) = UDGn{j}*fc_u - SH{j};
%         Un{j} = Un{j} + dt*b(istage)*F{j}(:,:,:,istage);
%     end
%     
%     if app.wave
%         SH = PDG;
%         for j=1:nsg
%             for jstage = 1:istage-1
%                 SH{j} = SH{j} + dt*a(istage,jstage)*G{j}(:,:,:,jstage);
%             end
%             SH{j} = SH{j}*fc_u;
%             
%             PDGn{j} = dt*a(istage,istage)*(UDGn{j}(:,1:ncu,:)+SH{j});
%             G{j}(:,:,:,istage) = PDGn{j}*fc_u-SH{j};
%             Pn{j} = Pn{j} + dt*b(istage)*G{j}(:,:,:,istage);            
%         end        
%     end
% end
% 
% UDGn = Un;
% 
% 
% function [a,b,t] = dirkcoeff(q,p)
% % q - number of stages
% % p - order of accuracy
% 
% % fprintf('DIRK (%d,%d) \n', p,q);
% 
% if q == 1 && p == 1
%     a = 1;
%     b = 1;
%     t = 1;
% elseif q == 1 && p == 2
%     % fprintf('Warning: Innacutate for non-wave formuations - use (2,2) or (3,3)');
%     a = 0.5;
%     b = 1;
%     t = 0.5;
% elseif q == 2 && p == 2
%     a = [1-0.5*sqrt(2), 0;
%            0.5*sqrt(2), 1-0.5*sqrt(2)];
%     b = [  0.5*sqrt(2), 1-0.5*sqrt(2)];
%     t = [1-0.5*sqrt(2), 1];
% elseif q == 2 && p == 3
%     % fprintf('Warning: Innacutate for non-wave formuations - use (2,2) or (3,3)');
%     a = [0.5+0.5/sqrt(3), 0;
%               -1/sqrt(3), 0.5 + 0.5/sqrt(3)];
%     b = [            0.5, 0.5];
%     t = [0.5+0.5/sqrt(3), 0.5-0.5/sqrt(3)];
% elseif q == 3 && p == 3
%     a1 = 0.4358665215;
%     t1 = (1+a1)/2;
%     b1 = -(6*a1^2-16*a1+1)/4;
%     b2 = (6*a1^2-20*a1+5)/4;    
%     a = [a1   ,  0, 0;
%          t1-a1, a1, 0;
%          b1   , b2, a1];
%     b = [b1, b2, a1];
%     t = [a1, t1, 1];    
% elseif q == 3 && p == 4
%     % fprintf('Warning: Innacutate for non-wave formuations - use (2,2) or (3,3)');
%     a1 = 2*cos(pi/18)/sqrt(3);
%     a = [ 0.5*(1+a1),            0,          0;
%              -0.5*a1,   0.5*(1+a1),          0;
%                 1+a1,    -(1+2*a1), 0.5*(1+a1)];
%     b = [ 1/(6*a1^2), 1-1/(3*a1^2), 1/(6*a1^2)];
%     t = [ 0.5*(1+a1),          0.5, 0.5*(1-a1)];
% else
%     error('Invalid (q,p) combination');
% end
% 
% function [UDG,UH,UHAT,AE,FE,DUDG_DUH,DUHAT_DUH,DUDGT,DUHAT,DUH] = hdg_solve(master,mesh,app,UDG,UH,SH,UHAT)
% 
% npf = size(mesh.perm,1);
% nfe = size(mesh.perm,2);
% ne  = mesh.ne;
% %nf  = mesh.nf;
% nc  = app.nc;
% ncu = app.ncu;
% %nd  = app.nd;
% nb  = npf*ncu;
% nsiz = mesh.nsiz;
% nsg = length(master);
% npv = size(mesh.dgnodes,1);
% 
% elcon   = mesh.elcon;
% il = zeros(nfe*npf*ncu,nfe*npf*ncu,ne);
% jl = zeros(nfe*npf*ncu,nfe*npf*ncu,ne);
% for i=1:ne
%     con = repmat((elcon(:,i)'-1)*ncu,ncu,1)+repmat((1:ncu)',1,nfe*npf);
%     con = reshape(con,nfe*npf*ncu,1);    
%     il(:,:,i) = repmat(con ,1,nfe*nb);
%     jl(:,:,i) = repmat(con',nfe*nb,1);        
% end
% 
% if isfield(app,'iterative') == 0
%     app.iterative = 0;
% elseif isempty(app.iterative)
%     app.iterative = 0;
% end
% 
% if iscell(master)==0
%     if isempty(SH)
%         SH = zeros(npv,nc,ne);
%     end
%     if app.wave==0
%         SH(:,ncu+1:end,:)=0;
%     end
% else
%     if isempty(SH)
%         SH = UDG;
%         for i=1:length(UDG);
%             SH{i}=0*SH{i};
%         end
%     end
%     if app.wave==0    
%         for i=1:length(UDG);
%             SH{i}(:,ncu+1:end,:,:)=0;
%         end        
%     end
% end
% 
% [UDG,DUDG_DUH,AE,FE,UHAT,DUHAT_DUH] = hdg_assemble(master,mesh,app,UDG,UH,SH,UHAT);  
% 
% F = sparse(reshape(il(:,1,:),nfe*nb*ne,1),ones(nfe*nb*ne,1),reshape(FE,nfe*nb*ne,1));                                                          
% KE = sparse(reshape(il,nfe*nfe*nb*nb*ne,1),reshape(jl,nfe*nfe*nb*nb*ne,1),reshape(AE,nfe*nfe*nb*nb*ne,1));        
% DUH = KE\F;
% DUH = -reshape(DUH,ncu,nsiz);        
% 
% DUHE = reshape(full(DUH(:,elcon(:))),[ncu,nfe*npf,ne]);    
% DUDGT = UDG;
% DUHAT = UHAT;
% for s=1:nsg
%     is = find(mesh.subgrids == s);
%     if (isempty(is) == 0) 
%         nref = master{s}.nref;
%         na = length(is);            
%         if nref==0
%             npv = size(UDG{s},1); 
%             for e=1:na
%                 tm = reshape(DUDG_DUH{s}(:,:,:,:,e),[npv*nc ncu*npf*nfe])*reshape(DUHE(:,:,is(e)),[npf*nfe*ncu 1]);
%                 DUDGT{s}(:,:,e) = reshape(tm,[npv nc]);
%             end                            
%         else                
%             npv = size(UDG{s},1); 
%             nes = size(master{s}.t,1); 
%             nf  = max(master{s}.elcon(:))/(npf);
%             for e=1:na
%                 tm = reshape(DUDG_DUH{s}(:,:,:,:,e),[npv*nc*nes ncu*npf*nfe])*reshape(DUHE(:,:,is(e)),[ncu*npf*nfe 1]);
%                 DUDGT{s}(:,:,(e-1)*nes+1:e*nes) = reshape(tm,[npv nc nes]);
%                 tm = reshape(DUHAT_DUH{s}(:,:,:,:,e),[ncu*npf*nf ncu*npf*nfe])*reshape(DUHE(:,:,is(e)),[ncu*npf*nfe 1]);
%                 DUHAT{s}(:,:,e) = reshape(tm,[ncu npf*nf]);
%             end            
%         end
%     end
% end        
% 
% 
% 
