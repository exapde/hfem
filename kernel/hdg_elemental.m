function [AE,FE,DUDG,DUDG_DUH] = hdg_elemental(master,app,dgnodes,bf,UDG,UH,SH)

npv = size(UDG,1);                
ne  = size(dgnodes,3);
npf = size(master.perm,1);
nfe = size(master.perm,2);
ncu = app.ncu;
nc  = app.nc;
ns  = 4200; %1050

nb = ceil(ne/ns);          
nk = 1:ns:ne;
nm = [nk(1:end); [nk(2:end)-1,ne]];    

UDG = permute(UDG,[1 3 2]);
SH  = permute(SH,[1 3 2]);
dgnodes = permute(dgnodes,[1 3 2]);
UH = permute(UH,[2 3 1]);

if app.getdqdg == 1
    nco = nc;
elseif app.getdqdg == 0
    nco = ncu;
else
    nco = 0;
end


AE  = zeros(ncu,nfe*npf,ncu,nfe*npf,ne);
FE  = zeros(ncu,nfe*npf,ne);
DUDG_DUH = zeros(npv,nco,ncu,npf*nfe,ne);
DUDG = zeros(npv,nco,ne);

% filename = 'Re';
% fileID = fopen([filename,'.bin'],'r');
% Re = fread(fileID,'double');
% fclose(fileID);
% Re = reshape(Re,size(FE));
% filename = 'He';
% fileID = fopen([filename,'.bin'],'r');
% He = fread(fileID,'double');
% fclose(fileID);
% He = reshape(He,size(AE));

for j=1:nb
    id = nm(1,j):nm(2,j);
    
    % compute volume integrals
    %BD = reshape(BD,[npv npv ne ncu nc]);
    [Ru, Rq, BD, M, C, L, Q, Ju, Jq, wrl] = volint(master,app,dgnodes(:,id,:),UDG(:,id,:),SH(:,id,:));
    
    %BDt = BD;
%     squeeze(Ru)
%     pause
    
%     [Ru2, Rq2, BD2, M2, C2] = mexvolint(master,app,permute(dgnodes(:,id,:),[1 3 2]),permute(UDG(:,id,:),[1 3 2]),permute(SH(:,id,:),[1 3 2]));
%     tm = permute(Ru2,[1 3 2]); max(abs(tm(:)-Ru(:)))
%     max(abs(Rq(:)-Rq2(:)))
%     tm = permute(BD2,[1 2 5 3 4]); max(abs(tm(:)-BD(:)))
%     max(abs(M(:)-M2(:)))
%     max(abs(C(:)-C2(:)))
%     pause
    
%     BD = permute(BD,[1 4 2 5 3]);
%     Ru = permute(Ru,[1 3 2]);
%     Ru(:,:,j)
%     pause
%     

    % compute face integrals        
    [Ru, Rq, Rh, BD, F, E, GK, H, Ju, Jq, Jh] = faceint(master,app,bf(:,id),dgnodes(:,id,:),UDG(:,id,:),UH(:,id,:),Ru,Rq,BD,Ju,Jq);            
%     squeeze(Rh)
%     pause
%    [BMiC,BMiE] = BMiCvol(master,app,BDt,M,C,E,wrl);    
%     M(:,:,j)
%     C(:,:,1,j)
%     C(:,:,2,j)
%     E(:,:,1,j)
%     E(:,:,2,j)
    %reshape(BD(:,:,:,:,j),[npv*ncu npv*nc])
    %reshape(F(:,:,:,:,j),[npv*ncu npf*nfe*ncu])
     %reshape(Ru(:,:,j),[1 npv*ncu])     

    %reshape(F(:,:,:,:,1),[npv*ncu npf*nfe*ncu])
    %reshape(Rh(:,:,1),[ncu, npf*nfe])
%     reshape(GK(:,:,:,:,1),[npf*nfe*ncu npv*nc]);
%     reshape(H(:,:,:,:,1),[npf*nfe*ncu npf*nfe*ncu])
%     pause
    
%     master.perm = master.perm-1;
%     master.permgeom = master.permgeom-1;
%     app.bcs = app.bcs';
%     [Ru2, Rq2, BD2, M2, C2, Juq2, Rh2, Jh2, E2, F2, G2, H2] = mexfaceint(master,app,permute(dgnodes(:,id,:),[1 3 2]),permute(UDG(:,id,:),[1 3 2]),permute(SH(:,id,:),[1 3 2]), permute(UH(:,id,:),[1 3 2]),bf(:,id));
%     max(abs(Ru(:)-Ru2(:)))
%     max(abs(Rq(:)-Rq2(:)))
%     max(abs(BD(:)-BD2(:)))
%     max(abs(M(:)-M2(:)))
%     max(abs(C(:)-C2(:)))  
%     max(abs(F(:)-F2(:)))
%     max(abs(Rh(:)-Rh2(:)))
%     max(abs(GK(:)-G2(:)))
%     max(abs(H(:)-H2(:)))
%     pause
    
    %Ru, *Rq, *BD, *M, *C, *Juq, *Rh, *Jh, *E, *F, *G, *H;        
    
    % perform schur compliment to obtain the elemental matrices and vectors
    if app.adjoint == 0
        [dudg, dudg_duh, ae, fe] = schur_primal(M, C, E, L, Q, BD, F, GK, H, Rq, Ru, Rh);
    else
        [dudg, dudg_duh, ae, fe] = schur_adjoint(M, C, E, L, Q, BD, F, GK, H, Jq, Ju, Jh);
    end
    
%    BMiCface(master,app,BMiC,BMiE,BD,M,C,E,wrb);    
    
%     nco=4;
%     dudg = dudg(:,1:nco,:);
%     dudg_duh = dudg_duh(:,1:nco,:,:,:);    
%     filename = 'elem';
%     fileID = fopen([filename,'.bin'],'r');
%     data = fread(fileID,'double');
%     fclose(fileID);
%     n2 = 0;
%     for i=1:200
%     n1 = n2+1; n2 = n1+npv*npv-1;
%     M1 = data(n1:n2);tm = M(:,:,i); err(i,1)=max(abs(tm(:)-M1(:)));
%     n1 = n2+1; n2=n1+npv*npv*2-1;    
%     C1 = data(n1:n2);tm = C(:,:,:,i); err(i,2)=max(abs(tm(:)-C1(:)));
%     n1 = n2+1; n2=n1+npv*npf*nfe*2-1;    
%     E1 = data(n1:n2);tm = E(:,:,:,i); err(i,3)=max(abs(tm(:)-E1(:)));
%     n1 = n2+1; n2=n1+npv*ncu*npv*nc-1;    
%     BD1 = data(n1:n2);tm = BD(:,:,:,:,i); err(i,4)=max(abs(tm(:)-BD1(:)));
%     n1 = n2+1; n2=n1+npv*ncu*npf*nfe*ncu-1;    
%     F1 = data(n1:n2);tm = F(:,:,:,:,i); err(i,5)=max(abs(tm(:)-F1(:)));
%     n1 = n2+1; n2=n1+npv*ncu-1;    
%     Ru1 = data(n1:n2);tm = Ru(:,:,i); err(i,6)=max(abs(tm(:)-Ru1(:)));
%     n1 = n2+1; n2=n1+npv*ncu*npf*nfe*nc-1;    
%     GK1 = data(n1:n2);tm = GK(:,:,:,:,i); err(i,7)=max(abs(tm(:)-GK1(:)));
%     n1 = n2+1; n2=n1+npf*nfe*ncu*npf*nfe*ncu-1;    
%     H1 = data(n1:n2);tm = H(:,:,:,:,i); err(i,8)=max(abs(tm(:)-H1(:)));
%     n1 = n2+1; n2=n1+npf*nfe*ncu-1;    
%     Rh1 = data(n1:n2);tm = Rh(:,:,i); err(i,9)=max(abs(tm(:)-Rh1(:)));    
%     n1 = n2+1; n2=n1+npf*nfe*ncu-1;    
%     Rh1 = data(n1:n2);tm = fe(:,:,i); err(i,10)=max(abs(tm(:)-Rh1(:)));    
%     n1 = n2+1; n2=n1+npf*nfe*ncu*npf*nfe*ncu-1;    
%     H1 = data(n1:n2);tm = ae(:,:,:,:,i); err(i,11)=max(abs(tm(:)-H1(:)));
%     n1 = n2+1; n2=n1+npv*ncu-1;    
%     Ru1 = data(n1:n2);tm = dudg(:,:,i); err(i,12)=max(abs(tm(:)-Ru1(:)));
%     n1 = n2+1; n2=n1+npv*ncu*npf*nfe*ncu-1;    
%     F1 = data(n1:n2);tm = dudg_duh(:,:,:,:,i); err(i,13)=max(abs(tm(:)+F1(:)));       
%     if max(abs(err(i,:)))>1e-9        
%         i 
%         err(i,:)
%         fe(:,:,i)
%         reshape(Rh1,size(fe(:,:,i)))
%         pause
%     end
%     end
%     pause
                
    DUDG(:,:,id)  = dudg(:,1:nco,:);
    DUDG_DUH(:,:,:,:,id) = dudg_duh(:,1:nco,:,:,:);
    AE(:,:,:,:,id) = ae;
    FE(:,:,id) = fe;        
end    


function [BMiC,BMiE] = BMiCvol(master,app,BD,M,C,E,wrl)

nc   = app.nc;
ncu  = app.ncu;
nd   = master.nd;
ngv  = master.ngv;
npv = size(M,1);
ne  = size(M,3);
npf = master.npf;
nfe = size(master.perm,2);
ndf = npf*nfe;

shapvt = master.shapvt(:,:,1);
shapvg = reshape(master.shapvg,[npv ngv*(nd+1)]);
shapvr = reshape(master.shapvg,[npv ngv (nd+1)]);
shapvr = shapvr(:,:,2:(nd+1));

%wrl = zeros(ngv*(nd+1),ne*ncu*nc);
wrl = reshape(wrl,[ngv,nd+1,ne,ncu,ncu,nd+1]);
wrq = permute(wrl(:,2:(nd+1),:,:,:,2:(nd+1)),[1 2 6 4 5 3]);
wrq = reshape(wrq,[ngv*nd*nd,ncu*ncu,ne]);
MiC = zeros(npv,npv,nd);
shapMiC = zeros(ngv,npv,nd);
shapshapMiC = zeros(npv,npv,ngv,nd,nd);
wrqshapMiC = zeros(ngv,nd,ncu,npv,ncu);
BMiC = zeros(npv,ncu,npv,ncu,ne);
MiE = zeros(npv,ndf,nd);
shapMiE = zeros(ngv,ndf,nd);
shapshapMiE = zeros(npv,ndf,ngv,nd,nd);
wrqshapMiE = zeros(ngv,nd,ncu,ncu,ndf);
BMiE = zeros(npv,ncu,ncu,ndf,ne);
err = zeros(ne,2);
for i = 1:ne
    Mi = inv(M(:,:,i));   
    for d = 1:nd
        MiC(:,:,d)  = Mi*C(:,:,d,i);    
        MiE(:,:,d)  = Mi*E(:,:,d,i);    
        shapMiC(:,:,d) = shapvt*MiC(:,:,d);
        shapMiE(:,:,d) = shapvt*MiE(:,:,d);
    end    
    for d = 1:nd
        for m = 1:nd
            for n = 1:ngv
                for j = 1:npv
                    for k = 1:npv                          
                        shapshapMiC(k,j,n,m,d) = shapvr(k,n,m)*shapMiC(n,j,d);
                    end
                end
            end
        end
    end
    BMiCt = reshape(shapshapMiC,[npv*npv ngv*nd*nd])*wrq(:,:,i);
    BMiCt = reshape(BMiCt,[npv npv ncu ncu]);
    BMiCt = permute(BMiCt,[1 3 2 4]);
    BMiC(:,:,:,:,i) = BMiCt;
        
    %wrqshapMiC = zeros(ngv,nd,ncu,npv,ncu);
    %wrq = reshape(wrl,[ngv,nd+1,ne,ncu,ncu,nd+1]);
    %shapMiC = zeros(ngv,npv,nd);
    for d = 1:ncu
        for m = 1:npv
            for n = 1:ncu
                for j = 1:nd
                    for k = 1:ngv                          
                        wrqshapMiC(k,j,n,m,d) = 0;
                        for l = 1:nd
                            wrqshapMiC(k,j,n,m,d) = wrqshapMiC(k,j,n,m,d) + ...
                                wrl(k,j+1,i,n,d,l+1)*shapMiC(k,m,l);
                        end
                    end
                end
            end
        end
    end    
    BMiCt = reshape(shapvr,[npv ngv*nd])*reshape(wrqshapMiC,[ngv*nd ncu*npv*ncu]);    
    
    %BD = reshape(BD,[npv npv ne ncu nc]);
    B    = reshape(BD(:,:,i,:,ncu+1:end),[npv npv ncu ncu nd]);
    B    = permute(B,[1 3 2 4 5]);        
    BMiCs = mapContractK(B,MiC,[1 2 4],[3 5],[],[1 3],2,[]);   
    BMiCs = permute(BMiCs,[1 2 4 3]);
    err(i,1) = max(abs(BMiCt(:)-BMiCs(:)));        
    
    for d = 1:nd
        for m = 1:nd
            for n = 1:ngv
                for j = 1:ndf
                    for k = 1:npv                          
                        shapshapMiE(k,j,n,m,d) = shapvr(k,n,m)*shapMiE(n,j,d);
                    end
                end
            end
        end
    end    
    BMiEt = reshape(shapshapMiE,[npv*ndf ngv*nd*nd])*wrq(:,:,i);
    BMiEt = reshape(BMiEt,[npv ndf ncu ncu]);
    BMiEt = permute(BMiEt,[1 3 4 2]);
    BMiE(:,:,:,:,i) = BMiEt;  
       
    for m = 1:ndf
        for d = 1:ncu        
            for n = 1:ncu
                for j = 1:nd
                    for k = 1:ngv                          
                        wrqshapMiE(k,j,n,d,m) = 0;
                        for l = 1:nd
                            wrqshapMiE(k,j,n,d,m) = wrqshapMiE(k,j,n,d,m) + ...
                                wrl(k,j+1,i,n,d,l+1)*shapMiE(k,m,l);
                        end
                    end
                end
            end
        end
    end    
    BMiEt = reshape(shapvr,[npv ngv*nd])*reshape(wrqshapMiE,[ngv*nd ncu*ncu*ndf]);    
        
    BMiEs = mapContractK(B,MiE,[1 2 4],[3 5],[],[1 3],2,[]);       
    err(i,2) = max(abs(BMiEt(:)-BMiEs(:)));
    
end
max(err)


function BMiC = BMiCface(master,app,BMiC,BMiE,BD,M,C,E,wrl)

nc   = app.nc;
ncu  = app.ncu;
nd   = master.nd;
ngv  = master.ngv;
npf  = master.npf;
ngf  = master.ngf;
nfe  = size(master.perm,2);
npv = size(M,1);
ne  = size(M,3);
ndf = npf*nfe;

% shapvt = master.shapvt(:,:,1);
% shapvg = reshape(master.shapvg,[npv ngv*(nd+1)]);
% shapvr = reshape(master.shapvg,[npv ngv (nd+1)]);
% shapvr = shapvr(:,:,2:(nd+1));

perm            = master.perm(:,:,1);
shapft          = master.shapft(:,:,1);
shapfg          = master.shapfg(:,:,1);
shapvt = zeros(ngf,nfe,npv);
shapvr = zeros(ngf,nfe,npv);
for j = 1:npv
    [indi,indj] = find(perm==j);
    for l = 1:length(indi)
        n = indi(l);
        k = indj(l);        
        for m = 1:ngf
            shapvt(m,k,j) = shapft(m,n);
            shapvr(m,k,j) = shapft(m,n)*master.gwfc(m);
        end        
    end
end
shapvt = reshape(shapvt,[ngf*nfe npv]);
shapvr = reshape(shapvr,[ngf*nfe npv]);

shapuhg = zeros(ngf,nfe,ndf);
for k = 1:nfe       
    for n = 1:npf
        j = (k-1)*npf+n;
        for m = 1:ngf       
            shapuhg(m,k,j) = shapft(m,n)*master.gwfc(m);
        end                
    end
end
shapuhg = reshape(shapvr,[ngf*nfe ndf]);

%wrk = reshape(bsxfun(@times,FH_udg,jac),[ngf nfe*ne*nch*nc]);
wrl = reshape(wrl,[ngf,nfe,ne,ncu,ncu,nd+1]);
wrq = permute(wrl(:,:,:,:,:,2:(nd+1)),[1 2 6 4 5 3]);
wrq = reshape(wrq,[ngf*nfe*nd,ncu*ncu,ne]);
MiC = zeros(npv,npv,nd);
%MiCf = zeros(npf,npf,nfe,nd);
MiCf = zeros(npf,nfe,npv,nd);
MiE = zeros(npv,ndf,nd);
shapMiC = zeros(ngf*nfe,npv,nd);
shapshapMiC = zeros(npv,npv,ngf*nfe,nd);
shapMiE = zeros(ngf*nfe,ndf,nd);
shapshapMiE = zeros(npv,ndf,ngf*nfe,nd);
shapshapGMiC = zeros(ndf,npv,ngf*nfe,nd);
shapshapGMiE = zeros(ndf,ndf,ngf*nfe,nd);
err = zeros(ne,3);
for i = 1:ne
    Mi = inv(M(:,:,i));   
    for d = 1:nd
        MiC(:,:,d)  = Mi*C(:,:,d,i);   
        MiE(:,:,d)  = Mi*E(:,:,d,i);
        shapMiC(:,:,d) = shapvt*MiC(:,:,d);
        shapMiE(:,:,d) = shapvt*MiE(:,:,d);
        MiCf(:,:,:,d) = reshape(MiC(perm(:),:,d),[npf nfe npv]);
%         for m = 1:nfe
%             MiCf(:,:,m,d) = MiC(perm(:,m),perm(:,m),d);            
% %             for k = 1:npf
% %                 for j = 1:npf
% %                     MiCf(j,k,m,d) = MiC(perm(j,m),perm(k,m),d);
% %                 end
% %             end
%         end
    end        
            
%     shapMiCf = reshape(shapft*reshape(MiCf,[npf npf*nfe*nd]),[ngf npf nfe nd]);
%     fhushapMiC = zeros(ngf,npf,nfe,ncu,ncu);   
%     for d = 1:nd
%         for k = 1:ncu
%             for j = 1:ncu
%                 for m = 1:nfe
%                     for n = 1:npf
%                         for l = 1:ngf
%                             fhushapMiC(l,n,m,j,k) = fhushapMiC(l,n,m,j,k) + ...
%                                 wrl(l,m,i,j,k,d+1)*shapMiCf(l,n,m,d);
%                         end
%                     end
%                 end
%             end
%         end
%     end                               
%     Df = reshape(shapfg*reshape(fhushapMiC,[ngf npf*nfe*ncu*ncu]),[npf npf nfe ncu ncu]);
%     BMiCf = zeros(npv,npv,1,ncu,ncu);
%     for is=1:nfe  % Adding face contributions - vector dependencies avoided
%         IP = perm(:,is);
%         BMiCf(IP,IP,1,:,:) = BMiCf(IP,IP,1,:,:) + Df(:,:,is,:,:);
%     end
        
    shapMiCf = reshape(shapft*reshape(MiCf,[npf nfe*npv*nd]),[ngf nfe npv nd]);    
    fhushapMiC = zeros(ngf,nfe,npv,ncu,ncu);     
    for d = 1:nd
        for k = 1:ncu
            for j = 1:ncu
                for n = 1:npv
                    for m = 1:nfe
                        for l = 1:ngf
                            fhushapMiC(l,m,n,j,k) = fhushapMiC(l,m,n,j,k) + ...
                                wrl(l,m,i,j,k,d+1)*shapMiCf(l,m,n,d);
                        end
                    end
                end
            end
        end
    end               
    Df = reshape(shapfg*reshape(fhushapMiC,[ngf nfe*npv*ncu*ncu]),[npf nfe npv ncu ncu]);
    BMiCf = zeros(npv,npv,ncu,ncu);
    for is = 1:nfe
        BMiCf(perm(:,is),:,:,:) = BMiCf(perm(:,is),:,:,:) + reshape(Df(:,is,:,:,:),[npf npv ncu ncu]);
    end
        
    % npf*ngf*nfe*npv*ncu*ncu
    % npv*npv*ngf*nfe*nd*ncu*ncu
    
    for d = 1:nd
        for m = 1:ngf*nfe            
            for j = 1:npv
                for k = 1:npv                          
                    shapshapMiC(k,j,m,d) = shapvr(m,k)*shapMiC(m,j,d);
                end
                for k = 1:ndf                          
                    shapshapGMiC(k,j,m,d) = shapuhg(m,k)*shapMiC(m,j,d);
                end
            end            
        end
    end    
    BMiCt = reshape(shapshapMiC,[npv*npv ngf*nfe*nd])*wrq(:,:,i);
    BMiCt = reshape(BMiCt,[npv npv ncu ncu]);  err(i,3) = max(abs(BMiCf(:)-BMiCt(:)));
    BMiCt = permute(BMiCt,[1 3 2 4]);
    BMiC(:,:,:,:,i)  = BMiC(:,:,:,:,i) + BMiCt;
            
    GMiCt = reshape(shapshapGMiC,[ndf*npv ngf*nfe*nd])*wrq(:,:,i);
    GMiCt = reshape(GMiCt,[ndf npv ncu ncu]);
    GMiCt = permute(GMiCt,[3 1 2 4]);    
    
    B    = reshape(BD(:,:,:,ncu+1:end,i),[npv,ncu,npv,ncu,nd]);               
    BMiCs = mapContractK(B,MiC,[1 2 4],[3 5],[],[1 3],2,[]);   
    BMiCs = permute(BMiCs,[1 2 4 3]);    
    BMiCt = BMiC(:,:,:,:,i);
    err(i,1) = max(abs(BMiCt(:)-BMiCs(:)));
    
    for d = 1:nd
        for m = 1:ngf*nfe            
            for j = 1:ndf
                for k = 1:npv                          
                    shapshapMiE(k,j,m,d) = shapvr(m,k)*shapMiE(m,j,d);
                end
                for k = 1:ndf                          
                    shapshapGMiE(k,j,m,d) = shapuhg(m,k)*shapMiE(m,j,d);
                end
            end            
        end
    end    
    BMiEt = reshape(shapshapMiE,[npv*ndf ngf*nfe*nd])*wrq(:,:,i);
    BMiEt = reshape(BMiEt,[npv ndf ncu ncu]);
    BMiEt = permute(BMiEt,[1 3 4 2]);
    BMiE(:,:,:,:,i)  = BMiE(:,:,:,:,i) + BMiEt;
        
    GMiEt = reshape(shapshapGMiE,[ndf*ndf ngf*nfe*nd])*wrq(:,:,i);
    GMiEt = reshape(GMiEt,[ndf ndf ncu ncu]);
    GMiEt = permute(GMiEt,[3 1 4 2]);
    
    BMiEs = mapContractK(B,MiE,[1 2 4],[3 5],[],[1 3],2,[]); 
    BMiEt = BMiE(:,:,:,:,i);
    err(i,2) = max(abs(BMiEt(:)-BMiEs(:)));    
end
max(err)











