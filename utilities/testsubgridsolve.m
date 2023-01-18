function testsubgridsolve(pgeom,porder,nref,ndim,elemtype,nodetype)

mastertest = mkmastersubgrid(pgeom,porder,2*pgeom,nref,ndim,elemtype,nodetype);

npm = size(mastertest.shapgeom,1);
npv = mastertest.npv;
ne = 10;
nes  = size(mastertest.t,1);
nd = ndim;
nc = 12;
ncu = 4;
nqf = size(mastertest.shapsf,2);
npf = size(mastertest.shapht,2);
nfe = size(mastertest.perm,2);
nf  = size(mastertest.f,1);

app.ncu = ncu;
app.nd = nd;
app.itmax = 1;
app.localsolve = 1;
app.fc_q = 1;
gam = 1.4;
epslm = 0.0;
Minf = 0.1;                  % Infinity conditions
pinf = 1/(gam*Minf^2);
alpha = 0*pi/180;
Re = 1000;
Pr = 0.72;
ui = [ 1, cos(alpha), sin(alpha), 0.5+pinf/(gam-1)];
app.iterative = 0;
app.localsolve = 1;
app.arg = {gam,epslm,Re,Pr,Minf,3};
app.bcm  = [2,1];  % 2: Wall, 1: Far-field
app.bcs  = [ui; ui];
app.wave = false;
app.tdep = false;
app.alag = false;
app.flg_q = 1;
app.flg_p = 0;
app.flg_g = 0;
app.fc_q = 1;
app.fc_u = 1;
app.fc_p = 0;
app.nd   = 2;
app.ncu  = 2+app.nd;         
app.nch  = app.ncu;          
app.nq   = app.nch*app.nd;   
app.nc   = app.nch+app.nq;   
app.time = [];
app.dtfc = [];
app.alpha = [];
app.itmax=1;      

dgnodes = rand(npm,ndim,ne);
udg = rand(npv,nc,nes*ne);
uh  = rand(ncu,nqf,nfe,ne);
uhat = rand(ncu,npf*nf,ne);
sh = 0*udg;

[DUDG,UDG,DUDG_DUH,DUHAT,UHAT,DUHAT_DUH,DUDG_DUHAT,Run,BDn,Rln,G,H,Rh,AE,FE,LE,K,F,L,dgx,UHATb] = subgridsolve(mastertest,app,dgnodes,sh,uh,udg,uhat);

bf = zeros(nfe,ne);
for k=1:ne                            
    for is=1:nfe  
        r = rand;
        if r>0.6
            bf(is,k) = -1;         
        end
        if r>0.8
            bf(is,k) = -2;         
        end
    end
end
[GA, HA, RA, FH, FH_udg, FH_uh,uda,uha,pga,nla,dga,Hh,Gh,FH_uha] = uhint(mastertest,app,udg,uh,dgnodes,bf);

save mastersubgridtest mastertest DUDG UDG DUDG_DUH DUDG_DUHAT DUHAT UHAT DUHAT_DUH ...
     Run BDn Rln G H Rh AE FE LE K F L app dgnodes sh uh udg uhat dgx UHATb ...
     GA HA RA bf FH FH_udg FH_uh uda uha pga nla dga Hh Gh FH_uha
 

% dgx = zeros(npm,nd,nes*ne);
% for i=1:ne
%     for k=1:nes
%         dgx(:,:,(i-1)*nes+k) = mastertest.shapgeom(:,:,k)*dgnodes(:,:,i);                           
%     end
% end

% [LE,Le,Lh,detJf] = sbouint(mastertest,app,dgx);
% 
% UDG = rand(npv,nc,nes);
% UH  = rand(ncu,nqf,nfe);
% UHAT = rand(ncu,npf,nfe,nes);
% 
% [BH, BH_udg, BH_uhat, an] = boufluxes(mastertest,UHAT,UH,UDG,dgx(:,:,1:nes));
% 
% Run = rand(npf, nfe, nes, ncu);
% BDn = rand(npf, npf, nfe, nes, ncu, nc);
% Rln = rand(npf, npf, nfe, nes, ncu, ncu);
% 
% [G, H, Rh, Gt] = uhatint(mastertest,UHAT,UH,UDG,dgx(:,:,1:nes),Run,BDn,Rln);
% %Gt = reshape(BDn,[npf npf*nfe*ne ncu nc]);
% 
% save mastersubgridtest mastertest Le Lh LE dgnodes dgx detJf ...
%      UDG UH UHAT BH BH_udg BH_uhat an Run BDn Rln G H Rh Gt

function [G, H, R, FH, FH_udg, FH_uh,uda,uha,pga,nla,dgx,Hh,Gh,wrn] = uhint(master,app,UDG,UH,dgnodes,bf)

shapht = master.shapht(:,:,1);
shapsg = master.shapsg;
perm = master.perm(:,:,1);

nsf = size(master.shapsf,3);
ngf = size(shapht,1);       
npf = size(shapht,2);
nes = size(master.t,1);
npv = size(UDG,1);
nc  = size(UDG,2);
ne  = size(UDG,3)/nes;    
nfe = size(perm,2);
ncu = size(UH,1);
nqf = size(UH,2);
npm = size(dgnodes,1);
nd  = size(dgnodes,2);

% UDG  : npv*nc*nes*ne  -> npv*nes*ne*nc -> npf*nsf*nfe*ne*nc -> ngf*nsf*nfe*ne*nc  
UDG  = reshape(permute(UDG,[1 3 2]),[npv*nes ne nc]);
udgg = reshape(UDG(master.ibdu(:),:,:),[npf nsf*nfe*ne*nc]);
udgg = reshape(shapht*udgg,[ngf*nsf*nfe*ne nc]);


% UH   : ncu*nqf*nfe*ne -> nqf*nfe*ne*nc -> ngf*nsf*nfe*ne*nc  
uhg = mapContractK(master.shapsf,UH,[1 3],2,[],2,[3 4 1],[]);
uhg  = reshape(uhg,[ngf*nsf*nfe*ne ncu]);

dgx  = zeros(npm,nd,nes,ne);
for i=1:ne
    for k=1:nes
        dgx(:,:,k,i) = master.shapgeom(:,:,k)*dgnodes(:,:,i);                           
    end
end

dgx = reshape(permute(dgx,[1 3 4 2]),[npm*nes ne nd]);
dgx = reshape(dgx(master.ibdgeom(:),:,:),[nqf nsf*nfe*ne nd]);
[pgf, nlgf, detJf] = bougeom(master,dgx); 
nlgf = reshape(nlgf,ngf*nsf*nfe*ne,nd);
pgf  = reshape(pgf,ngf*nsf*nfe*ne,nd);

nla=nlgf;pga= pgf;uda= udgg;uha= uhg;
[FH, FH_udg, FH_uh] = fhat(nlgf, pgf, udgg, uhg, app.arg, app.time); 

udgg = reshape(udgg,[ngf*nsf nfe ne nc]);
uhg  = reshape(uhg,[ngf*nsf nfe ne ncu]);
nlgf = reshape(nlgf,[ngf*nsf nfe ne nd]);
pgf  = reshape(pgf,[ngf*nsf nfe ne nd]);
for k=1:ne                            
    for is=1:nfe                            
        if bf(is,k) < 0                              
            ug = reshape(udgg(:,is,k,:),[ngf*nsf nc]);            
            uh = reshape(uhg(:,is,k,:),[ngf*nsf ncu]);
            nl = reshape(nlgf(:,is,k,:),[ngf*nsf nd]);
            pg = reshape(pgf(:,is,k,:),[ngf*nsf nd]);
            
            ib = app.bcm(-bf(is,k));        
            bv = repmat(app.bcs(-bf(is,k),:),[ngf*nsf 1]);            
            [fh, dfh_dudg, dfh_duh] = fbou(ib,bv,nl,pg,ug,uh,app.arg,app.time);                                                 

            im = (k-1)*ngf*nsf*nfe+(is-1)*ngf*nsf+1:(k-1)*ngf*nsf*nfe+is*ngf*nsf;            
            FH(im,:) = fh;
            FH_udg(im,:,:) = dfh_dudg;
            FH_uh(im,:,:) = dfh_duh;                   
        end                
    end
end
    

wrk = reshape(bsxfun(@times,FH,detJf(:)),[ngf nsf nfe*ne*ncu]);    
wrm = reshape(bsxfun(@times,FH_udg,detJf(:)),[ngf nsf nfe*ne*ncu*nc]);
wrn = reshape(bsxfun(@times,FH_uh,detJf(:)),[ngf nsf nfe*ne*ncu*ncu]);

R=0;
Hh=0;
for i=1:nsf
    R = R + reshape(shapsg(:,:,i)*reshape(wrk(:,i,:),[ngf nfe*ne*ncu]),[nqf nfe ne ncu]);    
    Hh = Hh + reshape(master.shapsgdotshapsf(:,:,i)*reshape(wrn(:,i,:),[ngf nfe*ne*ncu*ncu]),[nqf nqf nfe ne ncu ncu]);    
end

Gh = zeros(npf*nsf,nqf,nfe,ne,ncu,nc);
for i=1:nsf
    IP = (i-1)*npf+1:i*npf;
    Gh(IP,:,:,:,:,:) = reshape(master.shaphgdotshapsf(:,:,i)*reshape(wrm(:,i,:),[ngf nfe*ne*ncu*nc]),[npf nqf nfe ne ncu nc]);
end

ibdu = reshape(master.ibdu,[npf*nsf nfe]);
H = zeros(nqf*nfe,nqf*nfe,1,ne,ncu,ncu);
G = zeros(npv*nes,nqf,nfe,ne,ncu,nc); 
for is=1:nfe  
    IP = (is-1)*nqf+1:is*nqf;    
    H(IP,IP,1,:,:,:) = H(IP,IP,1,:,:,:) + Hh(:,:,is,:,:,:);
    G(ibdu(:,is),:,is,:,:,:) = G(ibdu(:,is),:,is,:,:,:) + Gh(:,:,is,:,:,:);    
end

R = permute(reshape(R,[nqf*nfe ne ncu]), [3 1 2]);
H = permute(reshape(H,[nqf*nfe nqf*nfe ne ncu ncu]),[4 1 5 2 3]);
G = permute(reshape(G,[npv nes nqf*nfe ne ncu nc]),[5 3 1 6 2 4]);

R = reshape(R,[ncu nqf*nfe ne]);
H = reshape(H,[ncu nqf*nfe ncu nqf*nfe ne]);
G = reshape(G,[ncu nqf*nfe npv nc*nes ne]);

function [DUDG,UDG,DUDG_DUH,DUHAT,UHAT,DUHAT_DUH,DUDG_DUHAT,Run,BDn,Rln,G,H,Rh,AE,FE,LE,K,F,L,dgx,UHATb] = subgridsolve(master,app,dgnodes,SH,UH,UDG,UHAT)

nes  = size(master.t,1);
npv  = size(UDG,1);
nc   = size(UDG,2);
ne   = size(UDG,3)/nes;
nqf  = size(UH,2);
npf  = size(master.perm,1);
nfe  = size(master.perm,2);
npm  = size(dgnodes,1);
ncu  = app.ncu;
nd   = app.nd;
ndf  = npf*nfe;
nbf  = ncu*ndf;
nf   = max(master.elcon(:))/(npf);
%nq   = nqf*nfe*ncu;            
            
il = zeros(ndf*ncu,ndf*ncu,nes);
jl = zeros(ndf*ncu,ndf*ncu,nes);
for i=1:nes
    con = repmat((master.elcon(:,i)'-1)*ncu,ncu,1)+repmat((1:ncu)',1,ndf);
    con = reshape(con,nfe*npf*ncu,1);    
    il(:,:,i) = repmat(con ,1,ndf*ncu);
    jl(:,:,i) = repmat(con',ndf*ncu,1);        
end                                        

dgx = zeros(npm,nd,nes*ne);
for i=1:ne
    for k=1:nes
        dgx(:,:,(i-1)*nes+k) = master.shapgeom(:,:,k)*dgnodes(:,:,i);                           
    end
end

DUHAT = UHAT;
DUDG_DUH = zeros(npv,nc*nes,ncu,nqf*nfe,ne);
DUHAT_DUH = zeros(ncu,npf*nf,ncu,nqf*nfe,ne);

eu = zeros(k,1);
for i=1:app.itmax
    UHATb = reshape(UHAT(:,master.elcon,:),[ncu npf nfe nes*ne]);    
    [DUDG,UDG,DUDG_DUHAT,Run,BDn,Rln] = elementsolve(master,app,dgx,SH,UHATb,UDG);                               
            
    for k=1:ne               
        m = (k-1)*nes+1:k*nes;                
        %[G, H, Rh] = uhatint(master,UHATb(:,:,:,m),UH(:,:,:,k),UDG(:,:,m),dgx(:,:,m),Run(:,:,m,:),BDn(:,:,:,m,:,:),Rln(:,:,:,m,:,:));            
        [G, H, Rh] = uhatint(master,UHATb(:,:,:,m),UH(:,:,:,k),UDG(:,:,m)-DUDG(:,:,m),dgx(:,:,m),Run(:,:,m,:),BDn(:,:,:,m,:,:),Rln(:,:,:,m,:,:));            
        [AE, FE] = elemmat(DUDG(:,:,m), DUDG_DUHAT(:,:,:,:,m), G, H, Rh);                            
                        
        K = sparse(reshape(il,nbf*nbf*nes,1),reshape(jl,nbf*nbf*nes,1),reshape(AE,nbf*nbf*nes,1));        
        F = sparse(reshape(il(:,1,:),nbf*nes,1),ones(nbf*nes,1),reshape(FE,nbf*nes,1));                                                          
        
        DUHAT(:,:,k) = -reshape(full(K\F),[ncu npf*nf]);
        UHAT(:,:,k) = UHAT(:,:,k) + DUHAT(:,:,k);         
        duhat = reshape(DUHAT(:,master.elcon,k),[ncu npf*nfe nes]);
        for e=1:nes            
            tm = mapContractK(DUDG_DUHAT(:,:,:,:,m(e)),duhat(:,:,e),[1 2],[3 4],[],[1 2],3,[]);   
            UDG(:,:,m(e)) = UDG(:,:,m(e)) + reshape(tm,[npv nc]);
            DUDG(:,:,m(e)) = DUDG(:,:,m(e)) + reshape(tm,[npv nc]);
        end        
        eu(k) = norm(duhat(:))/norm(reshape(UHAT(:,:,k),[ncu*npf*nf 1]));                 
    
        if eu(k)<=1e-8 || i==app.itmax
%             m = (k-1)*nes+1:k*nes;                
%             %[G, H, Rh] = uhatint(master,UHATb(:,:,:,m),UH(:,:,:,k),UDG(:,:,m),dgx(:,:,m),Run(:,:,m,:),BDn(:,:,:,m,:,:),Rln(:,:,:,m,:,:));            
%             [G, H, Rh] = uhatint(master,UHATb(:,:,:,m),UH(:,:,:,k),UDG(:,:,m)-DUDG(:,:,m),dgx(:,:,m),Run(:,:,m,:),BDn(:,:,:,m,:,:),Rln(:,:,:,m,:,:));            
%             [AE] = elemmat(DUDG(:,:,m), DUDG_DUHAT(:,:,:,:,m), G, H, Rh);                            
%             K = sparse(reshape(il,nbf*nbf*nes,1),reshape(jl,nbf*nbf*nes,1),reshape(AE,nbf*nbf*nes,1));                    
%             
            LE = sbouint(master,app,dgx(:,:,m));            
            L = zeros(ncu*npf*nf,ncu*nqf*nfe);
            for n=1:ncu*nqf*nfe
                L(:,n) = sparse(reshape(il(:,1,:),nbf*nes,1),ones(nbf*nes,1),reshape(LE(:,:,:,n),nbf*nes,1));                                                          
            end
            
            DUHAT_DUH(:,:,:,:,k) = reshape(K\L,[ncu npf*nf ncu nqf*nfe]);
            DUU = reshape(DUHAT_DUH(:,:,:,:,k), [ncu npf*nf ncu*nqf*nfe]);
            DUU = reshape(DUU(:,master.elcon,:),[ncu*npf*nfe nes ncu*nqf*nfe]);            
            tmp = zeros(npv*nc,nes,ncu*nqf*nfe);
            for e=1:nes                
                dudh = mapContractK(DUDG_DUHAT(:,:,:,:,m(e)),DUU(:,e,:),[1 2],[3 4],[],[1 2],3,[]);                   
                tmp(:,e,:) = reshape(dudh,[npv*nc,1,ncu*nqf*nfe]);
            end            
            DUDG_DUH(:,:,:,:,k) = reshape(tmp,[npv,nc*nes,ncu,nqf*nfe]);          
        end        
    end    
    if max(eu)<1e-8 || i==app.itmax
        break;
    end
end

function [G, H, Rh, Gt] = uhatint(master,UHAT,UH,UDG,dgnodes,Run,BDn,Rln)

perm = master.perm;
npv  = size(UDG,1);
ne   = size(UDG,3);
nc   = size(UDG,2);
npf  = size(perm,1);
nfe  = size(perm,2);
ncu  = numel(UHAT)/(npf*nfe*ne);
tm   = find(master.bf<0);
nbf  = length(tm(:));

[BH, BH_udg, BH_uhat, an] = boufluxes(master,UHAT,UH,UDG,dgnodes);

Run = reshape(Run,[npf*nfe*ne ncu]); 
Run(an,:) = reshape(master.shaphg(:,:,1)*BH,[npf*nbf ncu]);

Gt = reshape(BDn,[npf npf*nfe*ne ncu nc]);
Gt(:,an,:,:) = reshape(master.shaphgdotshaphc(:,:,1)*BH_udg,[npf npf*nbf ncu nc]);
Gt = reshape(Gt,[npf npf nfe ne ncu nc]); % [npf*nfs nqf nfe ne ncu nc]

Ht = reshape(-Rln,[npf npf*nfe*ne ncu ncu]); 
Ht(:,an,:,:) = reshape(master.shaphgdotshaphc(:,:,1)*BH_uhat,[npf npf*nbf ncu ncu]);
Ht  = reshape(Ht,[npf npf nfe ne ncu ncu]); % [nqf nqf nfe ne ncu ncu]

G   = zeros(npv,npf,nfe,ne,ncu,nc); % (npv*nes,nqf,nfe,ne,ncu,nc)
H   = zeros(npf*nfe,npf*nfe,1,ne,ncu,ncu);
for is=1:nfe  
    IP = (is-1)*npf+1:is*npf;
    G(perm(:,is),:,is,:,:,:) = G(perm(:,is),:,is,:,:,:) + Gt(:,:,is,:,:,:);
    H(IP,IP,1,:,:,:) = H(IP,IP,1,:,:,:) + Ht(:,:,is,:,:,:);
end

Rh  = permute(reshape(Run,[npf*nfe ne ncu]), [3 1 2]);
G  = permute(reshape(G,[npv npf*nfe ne ncu nc]),[4 2 1 5 3]);
H  = permute(reshape(H,[npf*nfe npf*nfe ne ncu ncu]),[4 1 5 2 3]);

function [BH, BH_udg, BH_uhat, an] = boufluxes(master,UHAT,UH,UDG,dgnodes)

shapht = master.shapht(:,:,1);
perm = master.perm(:,:,1);
an   = master.ibduhat;

nsf = size(master.shapsf,3);
ngf = size(shapht,1);       
npf = size(shapht,2);
nc  = size(UDG,2);
ne  = size(UDG,3);    
nfe = size(perm,2);
ncu = size(UH,1);
nqf = size(UH,2);
npm = size(dgnodes,1);
nd  = size(dgnodes,2);

% tm     = find(master.bf<0);
% nbf    = length(tm(:));
nbf = nsf*nfe;

dgx = reshape(permute(dgnodes,[1 3 2]),[npm*ne nd]);
dgx = reshape(dgx(master.ibdgeom,:),[nqf nsf*nfe nd]);
[~, ~, detJf] = bougeom(master,dgx); 

UHAT = reshape(UHAT,[ncu npf*nfe*ne]);
UHAT = reshape(UHAT(:,master.ibduhat),[ncu npf nsf nfe]);
UHAT = mapContractK(shapht,UHAT,1,2,[],2,[3 4 1],[]);   
UH   = mapContractK(master.shapsf,reshape(UH,[ncu nqf nfe]),[1 3],2,[],2,[3 1],[]);   
BH   = bsxfun(@times,reshape(UHAT,[ngf*nsf*nfe ncu])-reshape(UH,[ngf*nsf*nfe ncu]),detJf(:));
BH   = reshape(BH,ngf,nbf*ncu);

BH_uhat= zeros(ngf*nbf,ncu,ncu);
for i=1:ncu
    BH_uhat(:,i,i)=detJf(:);
end
BH_uhat = reshape(BH_uhat,ngf,nbf*ncu*ncu);

BH_udg  = zeros(ngf,nbf*ncu*nc);

% BH   = reshape(BH,ngf,nsf*nfe,ncu);
% BH   = permute(BH,[2,1,3]);
% BH_uhat = reshape(BH_uhat,ngf,nsf*nfe,ncu,ncu);
% BH_uhat = permute(BH_uhat,[2,1,3,4]);

function [LE,Le,Lh,detJf,dgx] = sbouint(master,app,dgnodes)

ne   = size(dgnodes,3);
npf  = size(master.perm,1);
nfe  = size(master.perm,2);
nqf  = size(master.permgeom,1);
ngf  = size(master.shapht,1);
ncu  = app.ncu;
npm  = size(dgnodes,1);
nd   = size(dgnodes,2);
nsf  = size(master.shapsf,3);

dgx = reshape(permute(dgnodes,[1 3 2]),[npm*ne nd]); 
dgx = reshape(dgx(master.ibdgeom,:),[nqf nsf*nfe nd]);
[~, ~, detJf] = bougeom(master,dgx); 

Le = zeros(npf,nsf,nqf,nfe);
detJf = reshape(detJf,[ngf nsf nfe]);
for i=1:nsf
    tm = master.shaphgdotshapsf(:,:,i)*reshape(detJf(:,i,:),[ngf nfe]);
    Le(:,i,:,:) = reshape(tm,[npf nqf nfe]);
end
%Le2 = reshape(Le,npf*nsf,[],nfe);

Lh = zeros(npf,nsf,nqf,ncu,nfe,ncu,nfe);
for i=1:nfe
    for j=1:ncu        
        Lh(:,:,:,j,i,j,i) = reshape(Le(:,:,:,i),[npf nsf nqf]);
    end
end

Lh = permute(Lh,[4 1 2 5 6 3 7]);
Lh = reshape(Lh,[ncu npf*nsf*nfe ncu*nqf*nfe]);
LE = zeros(ncu,npf*nfe*ne,ncu*nqf*nfe);
LE(:,master.ibduhat,:) = Lh;
LE = reshape(LE,[ncu npf*nfe ne ncu*nqf*nfe]);

%Le = permute(Le,[1 3 4 2]);



% ibduhat = reshape(master.ibduhat,[],nfe);
% LE2 = zeros(ncu,npf*nfe*ne,ncu,nqf,nfe);
% for i = 1:nfe
%     for j=1:ncu
%         LE2(j,ibduhat(:,i),j,:,i) = Le2(:,:,i);
%     end
% end
% max(abs(LE(:)-LE2(:)))
% %LE=LE2;
% LE = reshape(LE2,[ncu,npf*nfe,ne,ncu,nqf,nfe]);
% 
