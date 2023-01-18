function dgparallel(master,mesh,app,UDG,UH,SH,nproc)

% for i = 0:nproc-1
% filename = ['test' num2str(i)];
% fileID = fopen([filename,'.bin'],'r');
% data0 = fread(fileID,'double');
% fclose(fileID);
% 
% filename = ['out' num2str(i)];
% fileID = fopen([filename,'.bin'],'r');
% data1 = fread(fileID,'double');
% fclose(fileID);
% max(abs(data0(:)-data1(:)))
% end

ne = mesh.ne;
nf = mesh.nf;
npv = master.npv;
npf = master.npf;
nfe = size(master.perm,2);
ncu = app.ncu;
nc  = app.nc;

[dmd,meshp] = domaindecomposition2(mesh,nproc);

if isempty(SH)
    SH = zeros(npv,nc,ne);
end
for i = 1:nproc
    UDGp{i} = UDG(:,:,dmd{i}.elempart);
    SHp{i} = SH(:,:,dmd{i}.elempart);
    if strcmp(mesh.hybrid,'edg');
        UHp{i} = UH(:,dmd{i}.entpart);
    else
        UH = reshape(UH,[ncu npf nf]);
        UHp{i} = reshape(UH(:,:,dmd{i}.entpart),[ncu npf*length(dmd{i}.entpart)]);
    end
end

% [AE,FE] = hdg_elemental(master,app,mesh.dgnodes,mesh.bf,UDG,reshape(UH(:,mesh.elcon),[ncu npf*nfe ne]),SH);
% [Hg, Rg] = globalassembly(mesh, AE, FE);
% mesh = computeMDFordering(mesh, Hg);
% Pg = computeBILU0(mesh, Hg);
% restart = 20;
% tol = 1e-8;
% maxit = 1000;
% duh = gmres_bilu0(mesh, Pg, Hg, Rg, 0*Rg, restart, tol, maxit);
% duh = reshape(duh,size(Rg));

[H,R] = hdg_assemble(master,mesh,app,UDG,UH,0*UDG);
duh = H\R;
duh = reshape(duh,[],mesh.nf);
% max(abs(duh(:)-DUH(:)))
% pause

%Xg = matvec(mesh,Hg,Rg);

% [AE,FE] = hdg_elemental(master,app,mesh.dgnodes,mesh.bf,UDG,reshape(UH(:,mesh.elcon),[ncu ndf ne]),0*UDG);
% [Hg, Rg] = globalassembly(mesh, AE, FE);
% mesh = computeMDFordering(mesh, Hg, C1, w1);
% Mg = computeBILU0(mesh, Hg);
% restart = 20;
% tol = 1e-8;
% maxit = 1000;
% [x, flags, iter, rev] = gmres_bilu0(mesh, Mg, Hg, Rg, 0*Rg, restart, tol, maxit);

for i = 1:nproc      
    %uhp  = reshape(UH(:,mesh.elcon(:,dmd{i}.elempart)),[ncu npf*nfe size(meshp{i}.elcon,2)]);
    uhp = reshape(UHp{i}(:,meshp{i}.elcon),[ncu npf*nfe size(meshp{i}.elcon,2)]);    
    [ae{i},fe{i},dudg{i},dudg_duh{i}] = hdg_elemental(master,app,meshp{i}.dgnodes,meshp{i}.bf,UDGp{i},uhp,SHp{i});
%     tm = ae{i}-AE(:,:,:,:,dmd{i}.elempart); max(abs(tm(:)))
%     tm = fe{i}-FE(:,:,dmd{i}.elempart); max(abs(tm(:)))
%     tm = dudg{i}-DUDG(:,:,dmd{i}.elempart); max(abs(tm(:)))
%     tm = dudg_duh{i}-DUDG_DUH(:,:,:,:,dmd{i}.elempart); max(abs(tm(:)))
    
%     filename = ['out' num2str(i-1)];
%     fileID = fopen([filename,'.bin'],'r');
%     tm0 = fread(fileID,'double');
%     fclose(fileID);
%     max(abs(tm0(:)-ae{i}(:)))
%     pause
    
     [Hgp{i}, Rgp{i}, BJHg{i}, BKHg{i}] = globalassemblyp(meshp{i}, dmd{i}, ae{i}, fe{i});    
%     n = sum(dmd{i}.entpartpts(1:2));
%     max(max(abs(Rgp{i}-Rg(:,dmd{i}.entpart(1:n)))))
%     rp = meshp{i}.cbsr_rowpts;
%     ci = meshp{i}.cbsr_colind;
%     for j = 1:n
%         nj = rp(j+1)-rp(j);
%         for k = 1:nj
%             %Hgp{i}(:,:,rp(j)+k)
%             %[ci(rp(j)+1) ci(rp(j)+k)]
%             %[dmd{i}.entpart(ci(rp(j)+1)) dmd{i}.entpart(ci(rp(j)+k))]
%             m = mesh.cbsr_rowpts(dmd{i}.entpart(ci(rp(j)+1)));
%             q = find(mesh.cbsr_colind((m+1):(m+nj))==dmd{i}.entpart(ci(rp(j)+k)));
%             ta = Hgp{i}(:,:,rp(j)+k);
%             tb = Hg(:,:,m+q);
%             %max(abs(ta(:)-tb(:)))            
%         end
%     end    
%     
    filename = ['Rg' num2str(i-1)];
    fileID = fopen([filename,'.bin'],'r');
    Rg0 = fread(fileID,'double');
    fclose(fileID);
    max(abs(Rg0(1:numel(Rgp{i}(:)))-Rgp{i}(:)))    
    
%     max(abs(Rg0(:)))
%     max(abs(Rgp{i}(:)))
%     Rg0 = reshape(Rg0(1:numel(Rgp{i}(:))),size(Rgp{i}));
%     if i==2
%         for n = 1:size(Rg0,2)
%             err=Rgp{i}(:,n)-Rg0(:,n);
%             if max(abs(err(:)))>1e-8                
%                 n                
%             end
%         end
%     end    
%     pause

    filename = ['Hg' num2str(i-1)];
    fileID = fopen([filename,'.bin'],'r');
    Hg0 = fread(fileID,'double');
    fclose(fileID);
    max(abs(Hg0(1:numel(BJHg{i}(:)))-BJHg{i}(:)))
    
    filename = ['Kg' num2str(i-1)];
    fileID = fopen([filename,'.bin'],'r');
    Kg0 = fread(fileID,'double');
    fclose(fileID);
    max(abs(Kg0(1:numel(BKHg{i}(:)))-BKHg{i}(:)))
    
    
%     max(abs(Hg0(:)))
%     max(abs(BJHg{i}(:)))
%     Hg0 = reshape(Hg0,size(BJHg{i}));
%     if i==2
%         for n = 1:size(Hg0,3)
%             err=BJHg{i}(:,:,n)-Hg0(:,:,n);
%             if max(abs(err(:)))>1e-8
%                 BJHg{i}(:,:,n)
%                 Hg0(:,:,n)
%                 n
%                 pause
%             end
%         end
%     end    
end

% Simulate MPI send 
for i = 1:nproc      
    bsz = size(Rgp{i},1);
    nentsend = size(dmd{i}.entsend,1);
    buffsend{i} = zeros(bsz,nentsend);
    for j=1:nentsend
        buffsend{i}(:,j) = Rgp{i}(:,dmd{i}.entsend(j,2));
    end            
    
%     filename = ['buffsend' num2str(i-1)];
%     fileID = fopen([filename,'.bin'],'r');
%     buffsend0 = fread(fileID,'double');
%     fclose(fileID);
%     max(abs(buffsend0(:)-buffsend{i}(:)))    
%     pause
end

% Simulate matrix-vector product
for i = 1:nproc 
    Xgp{i} = matvecp1(meshp{i}.BJ_rowpts,meshp{i}.BJ_colind,BJHg{i},Rgp{i});        
end

% Simulate MPI receive 
for i = 1:nproc 
    [bsz,nrow] = size(Rgp{i});
    nentrecv = size(dmd{i}.entrecv,1);
    buffrecv{i} = zeros(bsz,nentrecv);        
    for j=1:length(dmd{i}.nbsd)
        % cpu i receive vectors from cpu k         
        k = dmd{i}.nbsd(j);        
        ik = find(dmd{i}.entrecv(:,1)==k);
        ki = find(dmd{k}.entsend(:,1)==i);        
        buffrecv{i}(:,ik) = buffsend{k}(:,ki);
    end
        
    filename = ['buffrecv' num2str(i-1)];
    fileID = fopen([filename,'.bin'],'r');
    buffrecv0 = fread(fileID,'double');
    fclose(fileID);
    max(abs(buffrecv0(:)-buffrecv{i}(:)))    
    %pause
    
    Rgp{i} = [Rgp{i} zeros(bsz,nentrecv)];  
    for j = 1:nentrecv        
        %k = dmd{i}.entrecv(j,2)-nrow;           
        %Rgp{i}(:,nrow+j) = buffrecv{i}(:,k);
        k = dmd{i}.entrecv(j,2);           
        Rgp{i}(:,k) = buffrecv{i}(:,j);
    end
           
    filename = ['Rgrecv' num2str(i-1)];
    fileID = fopen([filename,'.bin'],'r');
    Rg0 = fread(fileID,'double');
    fclose(fileID);
    Rg0 = reshape(Rg0,size(Rgp{i}));
    max(max(abs(Rg0(:)-Rgp{i}(:))))
%    max(max(abs(Rg0(:,(nrow+1):end)-Rgp{i}(:,(nrow+1):end))))
%     max(max(abs(Rgp{i}(:,(nrow+1):end))))
%     max(max(abs(Rg0(:,(nrow+1):end))))
%    Rgp{i}(:,(nrow+1):end)    
    
    Xgp{i} = matvecp2(meshp{i}.BK_rowpts,meshp{i}.BK_rowind,meshp{i}.BK_colind,Xgp{i},BKHg{i},Rgp{i}); 
    
    filename = ['x' num2str(i-1)];
    fileID = fopen([filename,'.bin'],'r');
    x0 = fread(fileID,'double');
    fclose(fileID);
    x0 = reshape(x0,size(Rgp{i}));
    max(max(abs(x0(:,1:nrow)-Xgp{i})))
    
    meshp{i}.BJ_nrows = length(meshp{i}.BJ_rowpts)-1;
    meshp{i}.BJ_nblks = length(meshp{i}.BJ_colind);
    meshp{i} = computeMDForderingp(meshp{i}, BJHg{i});
    filename = ['o2u' num2str(i-1)];
    fileID = fopen([filename,'.bin'],'r');
    o2u0 = fread(fileID,'int');
    fclose(fileID);    
    max(abs([meshp{i}.BILU0ordered2unordered-o2u0(1:meshp{i}.BJ_nrows)-1]))    
    
    Mg{i} = computeBILU0p(meshp{i}, BJHg{i});    
    filename = ['Mg' num2str(i-1)];
    fileID = fopen([filename,'.bin'],'r');
    Mg0 = fread(fileID,'double');
    fclose(fileID);
    max(abs(Mg0(1:numel(Mg{i}(:)))-Mg{i}(:)))./max(abs(Mg{i}(:)))    
    
    MgRg = applyBILU0p(meshp{i}, Mg{i}, Rgp{i});
    filename = ['MgRg' num2str(i-1)];
    fileID = fopen([filename,'.bin'],'r');
    MgRg0 = fread(fileID,'double');
    fclose(fileID);    
    MgRg0 = reshape(MgRg0,size(Rgp{i}));    
    max(max(abs(MgRg0(:)-MgRg(:))))    
    
    filename = ['duh' num2str(i-1)];
    fileID = fopen([filename,'.bin'],'r');
    duh0 = fread(fileID,'double');
    fclose(fileID);
    duh0 = reshape(duh0,size(Rgp{i}));
    e = dmd{i}.entpart(1:nrow);
    max(max(abs(duh0(:,1:nrow)-duh(:,e))))  
    i
    pause
end

% check matvec
% for i = 1:nproc
%     nrow = size(Xgp{i},2);
%     e = dmd{i}.entpart(1:nrow);
%     tm = Xg(:,e)-Xgp{i}; 
%     max(abs(tm(:)))
%     pause
% %     for j = 1:size(Xgp{i},2)
% %         e = dmd{i}.entpart(j);
% %         [Xg(:,e) Xgp{i}(:,j)]
% %         pause
% %     end
% end


function [Hg, Rg, BJHg, BKHg] = globalassemblyp(mesh, dmd, He, Re)

ne    = size(mesh.dgnodes,3);
nfe   = size(mesh.t2f,2); % number of faces per element    
%nf    = mesh.nf;
% nblks = mesh.cbsr_nblks;
% nrows = mesh.cbsr_nrows;
rowpts    = mesh.cbsr_rowpts;
colind    = mesh.cbsr_colind;
BJ_rowpts    = mesh.BJ_rowpts;
BJ_colind    = mesh.BJ_colind;
BK_rowpts    = mesh.BK_rowpts;
BK_colind    = mesh.BK_colind;
BK_rowind    = mesh.BK_rowind;

if strcmp(mesh.hybrid,'hdg') 
    ncu = size(Re,1);
    npf = size(Re,2)/nfe;
    bsz = ncu*npf;        
    He = reshape(He,[ncu npf nfe ncu npf nfe ne]);
    Re = reshape(Re,[ncu npf nfe ne]);    
    elcon = reshape(mesh.elcon,[npf nfe ne]);
    
    nrows = sum(dmd.entpartpts(1:2));
    nblks = mesh.cbsr_rowpts(nrows+1);    
    Hg = zeros(bsz,bsz,nblks);
    Rg = zeros(bsz,nrows);
    BJnblks = length(mesh.BJ_colind);
    BKnblks = length(mesh.BK_colind);
    BJHg = zeros(bsz,bsz,BJnblks);
    BKHg = zeros(bsz,bsz,BKnblks);
    for e = 1:ne
        t2f = mesh.t2f(e,:);                
        for a = 1:nfe
            i = t2f(a);  % block row i 
            if i<nrows+1
                ind = elcon(:,a,e) - (i-1)*npf;             
                tm = Re(:,ind,a,e);
                Rg(:,i) = Rg(:,i) + tm(:);            
            end
            for b = 1:nfe            
                j = t2f(b); % block column j                           
                for k = (rowpts(i)+1):rowpts(i+1)
                    if colind(k) == j, break; end
                end
                if (i<nrows+1) && (k<nblks+1)
                    jnd = elcon(:,b,e) - (j-1)*npf;  
                    Hg(:,:,k) = Hg(:,:,k) + reshape(He(:,ind,a,:,jnd,b,e),[bsz bsz]);                
                end                                
                if (i<nrows+1) && (j<nrows+1)
                    for k = (BJ_rowpts(i)+1):BJ_rowpts(i+1)
                        if BJ_colind(k) == j, break; end
                    end
                    jnd = elcon(:,b,e) - (j-1)*npf;
                    BJHg(:,:,k) = BJHg(:,:,k) + reshape(He(:,ind,a,:,jnd,b,e),[bsz bsz]);       
                end
                if (i<nrows+1) && (j>=nrows+1)                    
                    for m = 1:i
                        if BK_rowind(m) == i, break; end
                    end
                    for k = (BK_rowpts(m)+1):BK_rowpts(m+1)
                        if BK_colind(k) == j, break; end
                    end
                    jnd = elcon(:,b,e) - (j-1)*npf;
                    BKHg(:,:,k) = BKHg(:,:,k) + reshape(He(:,ind,a,:,jnd,b,e),[bsz bsz]);       
                end
            end
        end
    end    
elseif strcmp(mesh.hybrid,'edg') 
    nrows = sum(dmd.entpartpts(1:2));
    nblks = mesh.cbsr_rowpts(nrows+1);
    ncu = size(Re,1);
    bsz = ncu;    
    nn = size(mesh.elcon,1);
    He = reshape(He,[bsz nn bsz nn ne]);
    Re = reshape(Re,[bsz nn ne]);
    Hg = zeros(bsz,bsz,nblks);
    Rg = zeros(bsz,nrows);    
    BJnblks = length(mesh.BJ_colind);
    BKnblks = length(mesh.BK_colind);
    BJHg = zeros(bsz,bsz,BJnblks);
    BKHg = zeros(bsz,bsz,BKnblks);
    
    for e = 1:ne
        elc = mesh.elcon(:,e);
        for a = 1:nn
            i = elc(a);  % block row i
            if i<nrows+1
                Rg(:,i) = Rg(:,i) + Re(:,a,e);
            end
            for b = 1:nn            
                j = elc(b); % block column j
                for k = (rowpts(i)+1):rowpts(i+1)
                    if colind(k) == j, break; end
                end                
                if (i<nrows+1) && (k<nblks+1)
                    Hg(:,:,k) = Hg(:,:,k) + reshape(He(:,a,:,b,e),[bsz bsz]);                
                end
                if (i<nrows+1) && (j<nrows+1)
                    for k = (BJ_rowpts(i)+1):BJ_rowpts(i+1)
                        if BJ_colind(k) == j, break; end
                    end
                    BJHg(:,:,k) = BJHg(:,:,k) + reshape(He(:,a,:,b,e),[bsz bsz]);                                    
                end
                if (i<nrows+1) && (j>=nrows+1)                    
                    for m = 1:i
                        if BK_rowind(m) == i, break; end
                    end
                    for k = (BK_rowpts(m)+1):BK_rowpts(m+1)
                        if BK_colind(k) == j, break; end
                    end                    
                    BKHg(:,:,k) = BKHg(:,:,k) + reshape(He(:,a,:,b,e),[bsz bsz]);                                                        
                end
            end
        end
    end        
end


function mesh = computeMDForderingp(mesh, Hg)

nrows = mesh.BJ_nrows;
nblks = mesh.BJ_nblks;
maxcols = mesh.maxBlocksPerRow;
DeltaC = zeros(maxcols,maxcols);
rowpts    = mesh.BJ_rowpts;
colind    = mesh.BJ_colind;
w     = zeros(nrows,1);
%bsz  = size(Hg,1);
C  = zeros(nblks,1);
D  = zeros(nblks,1);
for i = 1:nrows
    for n = (rowpts(i)+1):(rowpts(i+1))
        tm = squeeze(Hg(:,:,rowpts(i)+1))\squeeze(Hg(:,:,n));
        C(n) = norm(tm(:));
        D(n) = norm(squeeze(Hg(:,:,n)),'fro');
    end
end

for k = 1:nrows
    n = (rowpts(k)+1):rowpts(k+1); % pointers to row k 
    c = colind(n);             % connectivity indices for row k    
    DeltaC = 0*DeltaC;
    for i = 2:length(n)
        ri  = c(i);              % row ri is a neighbor of row k  
        m = (rowpts(ri)+1):rowpts(ri+1); % pointers to row ri         
        a = find(colind(m)==k);      % find the position of row k in the neighboring list of row ri        
        ik = rowpts(ri)+a;  % find the connectivity index ik
        Cik = C(ik);             % get Cik
        Dik = D(ik);
        for j = 2:length(n)
            rj  = c(j);   % rj is a neighbor of row k  
            if j~=i                    % for each neigboring row rj that is different from row ri                
                kj = rowpts(k)+j; % find the connectivity index kj                                                             
                Ckj = C(kj);           % get Cik
                if isempty(find(colind(m)==rj, 1)) %If j does not neighbor i, then increase the discarded fill in
                    DeltaC(i,j) = Cik*Ckj;
                end
            end
        end
    end
    w(k) = norm(DeltaC(:));    
end


bign = 1e100;
p = zeros(nrows,1);
for l = 1:nrows
    [~,p(l)] = min(w);  
    w(p(l)) = bign;                % do not choose pl again
    rpl = (rowpts(p(l))+1):rowpts(p(l)+1); % pointers to row pl
    cpl = colind(rpl);                 % connectivity indices for row pl
    for e = 2:length(rpl)          % for each neighbor of row pl
        k = cpl(e);                % row k is a neighbor to row pl  
        % if row k is not numbered, then recompute the weight for row k
        if w(k) < bign                         
            n = (rowpts(k)+1):rowpts(k+1); % pointers to row k 
            c = colind(n);         % connectivity indices for row k
            nk = find(w(c)<bign);  % find neighbors that are not numbered
            ck = c(nk);            % unnumbered connectivity indices for row k
            DeltaC = 0*DeltaC;
            for i = 2:length(ck)
                ri  = ck(i);             % neighboring row ri to row k  
                m = (rowpts(ri)+1):rowpts(ri+1); % pointers to row ri         
                a = find(colind(m)==k);      % find the position of row k in the neighboring list of row ri                                
                ik = rowpts(ri)+a;  % find the connectivity index ik
                Cik = C(ik);             % get Cik
                Dik = D(ik);
                for j = 2:length(ck)
                    rj  = ck(j); % rj is a neighbor of row k  
                    if i~=j              % for each neigboring row rj that is different from row ri                                
                        kj = rowpts(k)+nk(j); 
                        Ckj = C(kj);
                        if isempty(find(colind(m)==rj, 1)) % If j does not neighbor i, then increase the discarded fill in
                            DeltaC(i,j) = Cik*Ckj;
                        end
                    end
                end
            end
            w(k) = norm(DeltaC(:));            
        end
    end
end

t = 0*p;
for i = 1:length(t)
    t(p(i)) = i;
end

mesh.BILU0ordered2unordered = p; % if p[i] = k then face k is pivoted at iteration i
mesh.BILU0unordered2ordered = t; % if t[k] = i then face k is pivoted at iteration i
% mesh.BILU0ordered2unordered = 1:mesh.nf;
% mesh.BILU0unordered2ordered = 1:mesh.nf;


function Mg = computeBILU0p(mesh, Hg)

nrows = mesh.BJ_nrows;
%nblks = mesh.cbsr_nblks;
rowpts    = mesh.BJ_rowpts;
colind    = mesh.BJ_colind;
ordered2unordered = mesh.BILU0ordered2unordered;
unordered2ordered = mesh.BILU0unordered2ordered;

Mg = Hg;
for j = 1:nrows
    rj = ordered2unordered(j);  % row rj    
    nj = (rowpts(rj)+1):rowpts(rj+1); % pointers to row rj 
    cj = colind(nj);               % connectivity indices for row rj        
    jj = rowpts(rj)+1;
    % inv(Ujj)
    Mg(:,:,jj) = inv(Mg(:,:,jj));
    for i = 2:length(nj)     % loop over each neighbor of row rj
        ri = cj(i);           % row ri is a neighbor of row rj 
        if unordered2ordered(ri) > j
           ni = (rowpts(ri)+1):rowpts(ri+1);
           ci = colind(ni);     % connectivity indices for row ri              
           aj = find(ci==rj);      % find the position of row rj in the neighboring list of row ri
           ij = rowpts(ri)+aj;
           % Lij = Uij*inv(Ujj)
           Mg(:,:,ij) = Mg(:,:,ij)*Mg(:,:,jj);
           ii = rowpts(ri)+1;
           ji = rowpts(rj)+i;
           % Uii = Uii - Lij*Uji
           Mg(:,:,ii) = Mg(:,:,ii) - Mg(:,:,ij)*Mg(:,:,ji);           
           for l = 2:length(ni)
               rij = ci(l); % row rij is a neighbor of both ri and rj
               if unordered2ordered(rij) > j
                   for m = 2:length(nj)
                       if rij == cj(m)       
                           il = rowpts(ri)+l;
                           jm = rowpts(rj)+m;
                           %  U_il <- U_il - L_ij*U_jm 
                           Mg(:,:,il) = Mg(:,:,il) - Mg(:,:,ij)*Mg(:,:,jm);
                       end
                   end
               end
           end
        end                    
    end
end


function v = applyBILU0p(mesh, Mg, x)

nrows = mesh.BJ_nrows;
%nblks = mesh.cbsr_nblks;
rowpts    = mesh.BJ_rowpts;
colind    = mesh.BJ_colind;
ordered2unordered = mesh.BILU0ordered2unordered;
unordered2ordered = mesh.BILU0unordered2ordered;

y = x;
% solve L y = x
for j = 1:nrows
    rj = ordered2unordered(j);  % row rj    
    nj = (rowpts(rj)+1):rowpts(rj+1); % pointers to row rj 
    cj = colind(nj);               % connectivity indices for row rj            
            
    % Compute y_j <- x_j - L_ji*y_i for all i<j neighboring j 
    y(:,rj) = x(:,rj);
    for i = 2:length(nj)  % loop over each neighbor of row rj
        ri = cj(i);           % row ri is a neighbor of row rj
        if unordered2ordered(ri) < j
            ji = rowpts(rj)+i;
            y(:,rj) = y(:,rj) - Mg(:,:,ji)*y(:,ri);            
        end
    end                
end

v = y;
% solve U v = y
for j = 1:nrows
    rj = ordered2unordered(nrows-j+1);  % row rj    
    nj = (rowpts(rj)+1):rowpts(rj+1); % pointers to row rj 
    cj = colind(nj);               % connectivity indices for row rj            
    jj = rowpts(rj)+1;
    
    % Compute v_j <- y_j - U_ji*y_i for all i>j neighboring j 
    v(:,rj) = y(:,rj);
    for i = 2:length(nj)  % loop over each neighbor of row rj
        ri = cj(i);           % row ri is a neighbor of row rj
        if unordered2ordered(ri) > unordered2ordered(rj)
            ji = rowpts(rj)+i;
            v(:,rj) = v(:,rj) - Mg(:,:,ji)*v(:,ri);            
        end
    end    
    
    % Compute v_j <- U_jj \ v_j 
    v(:,rj) = Mg(:,:,jj)*v(:,rj);
end

function Xg = matvecp1(rowpts,colind,Hg,Rg)

nrows = length(rowpts)-1;
bsz   = size(Rg,1);
Xg    = 0*Rg;
for i = 1:nrows
    k = (rowpts(i)+1):rowpts(i+1);
    j = colind(k);    
    Xg(:,i) = reshape(Hg(:,:,k),[bsz bsz*length(k)])*reshape(Rg(:,j),[bsz*length(k) 1]);
end

function Xg = matvecp2(rowpts,rowind,colind,Xg,Hg,Rg)

nrows = length(rowpts)-1;
bsz   = size(Rg,1);
for i = 1:nrows
    k = (rowpts(i)+1):rowpts(i+1);
    j = colind(k);    
    m = rowind(i);
    Xg(:,m) = Xg(:,m) + reshape(Hg(:,:,k),[bsz bsz*length(k)])*reshape(Rg(:,j),[bsz*length(k) 1]);
end

function [x, flags, iter, rev] = gmresp(mesh, M, A, b, x, restart, tol, maxit)

% check the number of input arguments
if nargin < 4
  error('Not enough input arguments.');
end

% get the dimension from b
N = size(b(:),1);

% Default parameters
if nargin < 5 || isempty(x),       x = zeros(N,1);  end;
if nargin < 6 || isempty(restart), restart = 20;    end;
if nargin < 7 || isempty(tol),     tol = 1e-6;      end;
if nargin < 8 || isempty(maxit),   maxit = N;       end;

b = applyBILU0(mesh, M, b);

% initialization
nrmb   = norm(b(:)); 
flags  = inf; 
iter   = 0; 
cycle  = 0;

% allocate memory
hh = zeros(restart+1,restart);
v  = zeros(N,restart+1);
e1 = zeros(restart+1,1); e1(1)=1;
rev = zeros(restart,1);

while (1) 
    % perform matrix-vector multiplication        
    r = matvec(mesh, A,x);
    
    % Apply preconditioner
    r = applyBILU0(mesh, M, r);
    
    r = b - r;
    
    % Normalize the residual
    beta = norm(r(:));
    %if beta 
    
    v(:,1) = r(:)/beta;    
    
    res  = beta;
%     iter = iter+1;
%     rev(iter) = res;    
    for j = 1:restart            
        
        % perform matrix-vector multiplication        
        r = matvec(mesh,A,reshape(v(:,j),size(x)));
        
        % Apply preconditioner
        r = applyBILU0(mesh, M, r);        
                        
        % Arnoldi process (i.e., GS orthogonalization on the Krylov supspace)        
        v(:,j+1) = r(:);
        for i = 1:j        
            hh(i,j) = v(:,i)'*v(:,j+1);
            v(:,j+1) = v(:,j+1) - hh(i,j)*v(:,i);        
        end
        hh(j+1,j) = norm(v(:,j+1));     
        if (hh(j+1,j) ~= 0.0)          
            v(:,j+1) = v(:,j+1)/hh(j+1,j);   
        else            
            break;
        end

        % solve the reduced system 
        y = hh(1:j+1,1:j)\(beta*e1(1:j+1));
        
        % compute the residual norm
        res = norm(beta*e1(1:j+1)-hh(1:j+1,1:j)*y);
        
        iter = iter + 1; 
        rev(iter) = res;
        
                % set flag=0 if convergence
        if res/nrmb <= tol
            flags = 0;
            break;
        end        
        
        % set flag=1 if reaching maximum iteration 
        if iter>=maxit
            flags = 1;
            break;
        end    
        
    end 
      
    % compute the solution    
    x(:) = x(:) + v(:,1:length(y))*y;    
    
    cycle = cycle + 1;

    % stop if converging or reaching the maximum iteration
    if flags < inf,      
        break; 
    end;     
end
fprintf('gmres(%d) converges at iteration %d to a solution with relative residual %g\n', [restart iter res/nrmb]);       
rev = rev/nrmb;


function [dmd,meshp] = domaindecomposition2(mesh,nproc)

if strcmp(mesh.hybrid,'hdg')
    meshp = domaindecomposition_hdg(mesh,nproc);
elseif strcmp(mesh.hybrid,'edg')    
    meshp = domaindecomposition_edg(mesh,nproc);    
end
dmd = meshp;

[npf,nfe] = size(mesh.perm);
for i=1:nproc % loop over each subdomain
    meshp{i}.perm = mesh.perm;
    meshp{i}.permgeom = mesh.permgeom;
    meshp{i}.dgnodes = mesh.dgnodes(:,:,meshp{i}.elempart);
    meshp{i}.bf = mesh.bf(:,meshp{i}.elempart);
    % compute local element-to-entity connectivities
    elcon = reshape(mesh.elcon,[npf nfe mesh.ne]);
    meshp{i}.elcon = elcon(:,:,meshp{i}.elempart);
    ne = length(meshp{i}.elempart);        
    %reshape(meshp{i}.elcon,[npf*nfe ne])
    if strcmp(mesh.hybrid,'hdg')
        t2f = mesh.t2f(meshp{i}.elempart,:);    
        for j = 1:ne
            for k = 1:nfe
                t2f(j,k) = find(t2f(j,k) == meshp{i}.entpart);
            end
        end
        for j = 1:ne
            for k = 1:nfe
                ii = elcon(:,k,meshp{i}.elempart(j)) - (mesh.t2f(meshp{i}.elempart(j),k)-1)*npf;                
                ind = ((t2f(j,k)-1)*npf+1):(t2f(j,k)*npf);
                meshp{i}.elcon(:,k,j) = ind(ii);
            end
        end
        meshp{i}.t2f = t2f;
    elseif strcmp(mesh.hybrid,'edg')        
        for j = 1:ne
            for k = 1:nfe        
                for m = 1:npf               
                    meshp{i}.elcon(m,k,j) = find(meshp{i}.elcon(m,k,j) == meshp{i}.entpart);
                end
            end
        end
        [~,meshp{i}.t2f] = mkt2f(mesh.t(meshp{i}.elempart,:),mesh.elemtype);
    end        
    meshp{i}.elcon = reshape(meshp{i}.elcon,[npf*nfe ne]);   
    meshp{i} = block_crs(meshp{i},mesh.hybrid);
    
    nownent = sum(meshp{i}.entpartpts(1:2));
    meshp{i}.BJ_rowpts = zeros(nownent+1,1);
    meshp{i}.BJ_colind = [];
    meshp{i}.BK_rowpts = 0;
    meshp{i}.BK_colind = [];
    meshp{i}.BK_rowind = 0;
    m = 1;
    for j = 1:nownent
        nj = (meshp{i}.cbsr_rowpts(j)+1):meshp{i}.cbsr_rowpts(j+1);
        cj = meshp{i}.cbsr_colind(nj);
        i1 = find(cj<=nownent);
        c1 = cj(i1);
        meshp{i}.BJ_rowpts(j+1) = meshp{i}.BJ_rowpts(j)+length(c1);
        meshp{i}.BJ_colind = [meshp{i}.BJ_colind; c1];
        i2 = find(cj>nownent);
        c2 = cj(i2);
        if ~isempty(c2)
            meshp{i}.BK_rowpts(m+1) = meshp{i}.BK_rowpts(m)+length(c2);
            meshp{i}.BK_colind = [meshp{i}.BK_colind; c2];
            meshp{i}.BK_rowind(m)  = j;
            m = m + 1;
        end
    end            
end

%[ae{i},fe{i},dudg{i},dudg_duh{i}] = hdg_elemental(masterp{i},appp{i},dgnodes{i},bf{i},udg{i},uh{i},sh{i});
    
function dmd = domaindecomposition_hdg(mesh,nproc)

[elempart,entpart] = meshpart(mesh.t2f,nproc);

% plot domain partition
% plotpart(mesh, nproc, elempart, entpart);

elempart = elempart+1;
entpart = entpart+1;

% loop over each subdomain
for i=1:nproc
    
    % list of elements in subdomain i                       
    indelem{i} = find(elempart==i); 
    
    % list of faces in subdomain i                       
    indent{i} = find(entpart==i); 
end

% loop over each subdomain
for i=1:nproc    
    % find all interface faces
    f2f = mesh.f2f(:,indent{i});
    i1 = find(f2f(:)>0);
    in = ones(size(f2f));
    in(i1) = ismember(f2f(i1),indent{i});
    [~,j1] = find(in==0);
    j1 = unique(j1);

    % list of all interface faces
    intf = indent{i}(j1);
    
    % find all faces that are connected to interface faces
    f2f = f2f(:,j1);
    i1 = find(f2f(:)>0);
    in = ones(size(f2f));
    in(i1) = ismember(f2f(i1),indent{i});
    f0 = unique(f2f(in==0));  % faces do not belong to subdomain i
    %f1 = unique(f2f(in==1));
    %f1 = f1(f1>0);            % faces belong to subdomain i    
    
    % [interior faces, own interface faces, other faces]   
    dmd{i}.entpart = [setdiff(indent{i},intf); intf; f0];  
    dmd{i}.entpartpts = [length(indent{i})-length(intf) length(intf) length(f0)];    
        
    % store faces received from neighboring subdomains to perform the matrix-vector product
    dmd{i}.entrecv = [0*f0 length(indent{i})+(1:1:length(f0))' f0];
    for k = 1:nproc
        if k ~= i
            in = ismember(f0,indent{k});
            dmd{i}.entrecv(in,1) = k;
        end
    end
    dmd{i}.entrecv = unique(dmd{i}.entrecv,'rows');    
    dmd{i}.nbsd = unique(dmd{i}.entrecv(:,1))';    
    
    % find all elements that are connected to interface faces
    elem = mesh.f(intf,end-1:end);
    i1 = find(elem(:)>0);
    in = ones(size(elem));
    in(i1) = ismember(elem(i1),indelem{i});
    elem0 = unique(elem(in==0)); % elements do not belong to subdomain i
    elem1 = unique(elem(in==1));
    elem1 = elem1(elem1>0);      % elements belong to subdomain i
    
    % [interior elements, interface elements, other elements]  
    dmd{i}.elempart = [setdiff(indelem{i},elem1); elem1; elem0];
    dmd{i}.elempartpts = [length(indelem{i})-length(elem1) length(elem1) length(elem0)];    
    
    % store elements received from neighboring subdomains to perform the matrix-vector assembly
    dmd{i}.elemrecv = [0*elem0 length(indelem{i})+(1:1:length(elem0))' elem0];
    for k = 1:nproc
        if k ~= i
            in = ismember(elem0,indelem{k});
            dmd{i}.elemrecv(in,1) = k;
        end
    end
    dmd{i}.elemrecv = unique(dmd{i}.elemrecv,'rows');        
end

% store faces sent to neighboring subdomains to perform the matrix-vector product
for k = 1:nproc
    dmd{k}.entsend = [];
end
for i = 1:nproc        
    for j = 1:length(dmd{i}.nbsd)
        % cpu k sends information to cpu i
        k = dmd{i}.nbsd(j);
        ii = dmd{i}.entrecv(:,1)==k;
        tm = dmd{i}.entrecv(ii,:);
        tm(:,1) = i;        
        for m = 1:size(tm,1)
            tm(m,2) = find(dmd{k}.entpart==tm(m,3));
        end
        dmd{k}.entsend = [dmd{k}.entsend; tm];        
    end    
end

% store elements sent to neighboring subdomains to assemble the linear system
for k = 1:nproc
    dmd{k}.elemsend = [];
end
for i = 1:nproc        
    for j = 1:length(dmd{i}.nbsd)
        % cpu k sends information to cpu i
        k = dmd{i}.nbsd(j);
        ii = dmd{i}.elemrecv(:,1)==k;
        tm = dmd{i}.elemrecv(ii,:);
        tm(:,1) = i;        
        for m = 1:size(tm,1)
            tm(m,2) = find(dmd{k}.elempart==tm(m,3));
        end
        dmd{k}.elemsend = [dmd{k}.elemsend; tm];        
    end    
end

for i = 1:nproc 
    dmd{i}.elempart'
    dmd{i}.elemsend 
    dmd{i}.elemrecv 
    dmd{i}.entpart'
    dmd{i}.entsend
    dmd{i}.entrecv
end

    
function dmd = domaindecomposition_edg(mesh,nproc)

[elempart,entpart] = meshpart(mesh.elcon',nproc);

% plot domain partition
% plotpart(mesh, nproc, elempart, entpart);

elempart = elempart+1;
entpart = entpart+1;

% loop over each subdomain
for i=1:nproc
    
    % list of elements in subdomain i                       
    indelem{i} = find(elempart==i); 
    
    % list of edg nodes in subdomain i                       
    indent{i} = find(entpart==i); 
end

% loop over each subdomain
for i=1:nproc
    % find all interface edgnodes    
    edgnumber = zeros(length(indent{i}),1); 
    for j = 1:length(indent{i})
        [~,j1] = find(mesh.elcon==indent{i}(j));
        j1 = unique(j1);
        edg1 = mesh.elcon(:,j1);
        edg1 = unique(edg1(:));
        if all(ismember(edg1,indent{i}))
            % edgnode indent{i}(j) is inside the subdomain i
            edgnumber(j) = 2;
        else
            % edgnode indent{i}(j) is on the interface between two subdomains
            edgnumber(j) = 1;            
        end
    end
        
    % list of all interface edgnodes
    intf = indent{i}(edgnumber==1);
    
    % find all neighboring edgnodes that are connected to interface edgnodes 
    edg0 = [];
    for j = 1:length(intf)
        [~,j1]=find(mesh.elcon==intf(j));
        j1 = unique(j1);
        edg1 = mesh.elcon(:,j1);
        edg1 = unique(edg1(:));
        in = ismember(edg1,indent{i});
        edg0 = [edg0; edg1(in==0)];
    end
    edg0 = unique(edg0);
    
    % [interior edgnodes, own interface edgnodes, other edgnodes]   
    dmd{i}.entpart = [setdiff(indent{i},intf); intf; edg0];  
    dmd{i}.entpartpts = [length(indent{i})-length(intf) length(intf) length(edg0)];    
    
    % store edgnodes received from neighboring subdomains to perform the matrix-vector product
    dmd{i}.entrecv = [0*edg0 length(indent{i})+(1:1:length(edg0))' edg0];
    for k = 1:nproc
        if k ~= i
            in = ismember(edg0,indent{k});
            dmd{i}.entrecv(in,1) = k;
        end
    end
    dmd{i}.entrecv = unique(dmd{i}.entrecv,'rows');    
    dmd{i}.nbsd = unique(dmd{i}.entrecv(:,1))';    
    
    % find all elements that are connected to interface edgnodes
    elem0 = [];
    elem1 = [];
    for j = 1:length(intf)
        [~,j1]=find(mesh.elcon==intf(j));
        j1 = unique(j1);        
        in = ismember(j1,indelem{i});
        elem0 = [elem0; j1(in==0)];
        elem1 = [elem1; j1(in==1)];
    end
    elem0 = unique(elem0);
    elem1 = unique(elem1);
    in = zeros(length(elem1),1);        
    for j = 1:length(elem1)        
        if all(ismember(mesh.elcon(:,elem1(j)),indent{i}))            
            in(j) = 1;
        end
    end  
    elem1 = elem1(in==0);    
    
    % [interior elements, interface elements, other elements]  
    dmd{i}.elempart = [setdiff(indelem{i},elem1); elem1; elem0];
    dmd{i}.elempartpts = [length(indelem{i})-length(elem1) length(elem1) length(elem0)];    
    
    % store elements received from neighboring subdomains to perform the matrix-vector assembly
    dmd{i}.elemrecv = [0*elem0 length(indelem{i})+(1:1:length(elem0))' elem0];
    for k = 1:nproc
        if k ~= i
            in = ismember(elem0,indelem{k});
            dmd{i}.elemrecv(in,1) = k;
        end
    end
    dmd{i}.elemrecv = unique(dmd{i}.elemrecv,'rows');        
end

% store edgnodes sent to neighboring subdomains to perform the matrix-vector product
for k = 1:nproc
    dmd{k}.entsend = [];
end
for i = 1:nproc        
    for j = 1:length(dmd{i}.nbsd)
        % cpu k sends information to cpu i
        k = dmd{i}.nbsd(j);
        ii = dmd{i}.entrecv(:,1)==k;
        tm = dmd{i}.entrecv(ii,:);
        tm(:,1) = i;        
        for m = 1:size(tm,1)
            tm(m,2) = find(dmd{k}.entpart==tm(m,3));
        end
        dmd{k}.entsend = [dmd{k}.entsend; tm];        
    end    
end

% store elements sent to neighboring subdomains to assemble the linear system
for k = 1:nproc
    dmd{k}.elemsend = [];
end
for i = 1:nproc        
    for j = 1:length(dmd{i}.nbsd)
        % cpu k sends information to cpu i
        k = dmd{i}.nbsd(j);
        ii = dmd{i}.elemrecv(:,1)==k;
        tm = dmd{i}.elemrecv(ii,:);
        tm(:,1) = i;        
        for m = 1:size(tm,1)            
            tm(m,2) = find(dmd{k}.elempart==tm(m,3));
        end
        dmd{k}.elemsend = [dmd{k}.elemsend; tm];        
    end    
end

for i = 1:nproc 
    dmd{i}.elempart'
    dmd{i}.elemsend 
    dmd{i}.elemrecv 
    dmd{i}.entpart'
    dmd{i}.entsend
    dmd{i}.entrecv
end

