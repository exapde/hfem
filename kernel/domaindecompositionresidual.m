function meshp = domaindecompositionresidual(mesh,nproc)

[elempart,~] = meshpart(mesh,nproc,mesh.p,mesh.t);
elempart = elempart+1;
%entpart = entpart+1;

[npf,nfe] = size(mesh.perm);
elcon = reshape(mesh.elcon,[npf nfe mesh.ne]);

mesh.t2t = mkt2t(mesh.t,mesh.elemtype);

% list of elements in subdomain i                       
for i=1:nproc        
    
    meshp{i}.elempart = find(elempart==i);     
    
    % find all interface elements
    t2t = mesh.t2t(meshp{i}.elempart,:);
    i1 = find(t2t(:)>0);
    in = ones(size(t2t));
    in(i1) = ismember(t2t(i1),meshp{i}.elempart);    
    [j1,~] = find(in==0);
    j1 = unique(j1);     
    inte = meshp{i}.elempart(j1);
    % [interface elements, interior elements]
    meshp{i}.elempart = [inte; setdiff(meshp{i}.elempart,inte)];
    meshp{i}.elempartpts = [length(inte) length(meshp{i}.elempart)-length(inte)];
        
    meshp{i}.perm = mesh.perm;
    meshp{i}.permgeom = mesh.permgeom;
    meshp{i}.dgnodes = mesh.dgnodes(:,:,meshp{i}.elempart);    
    meshp{i}.elcon = elcon(:,:,meshp{i}.elempart);
    meshp{i}.bf = mesh.bf(:,meshp{i}.elempart);    
            
    % interfaces have bf=0
%     for j = 1:length(inte)
%         t2t = mesh.t2t(inte(j),:);
%         i1 = find(t2t(:)>0);
%         in = ones(size(t2t));
%         in(i1) = ismember(t2t(i1),meshp{i}.elempart);    
%         j1 = find(in==0);
%         meshp{i}.bf(j1,j) = 0;
%     end
    
    ne = length(meshp{i}.elempart);
    t2f = mesh.t2f(meshp{i}.elempart,:);            
    meshp{i}.entpart = unique(t2f(:));
        
    % find all interface faces
    f2f = mesh.f2f(:,meshp{i}.entpart);
    i1 = find(f2f(:)>0);
    in = ones(size(f2f));
    in(i1) = ismember(f2f(i1),meshp{i}.entpart);
    [~,j1] = find(in==0);
    j1 = unique(j1);        
    intf = meshp{i}.entpart(j1);
    % [interface faces, interior faces]
    meshp{i}.entpart = [intf; setdiff(meshp{i}.entpart,intf)];
    meshp{i}.entpartpts = [length(intf) length(meshp{i}.entpart)-length(intf)];
    
    % fix t2f and elcon
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
    
    % fix bf so that interface faces have bf=0
    ind = meshp{i}.t2f'<=length(intf);
    meshp{i}.bf(ind) = 0;            
    
    % [neighboring subdomains, local indices, global indices ]
    meshp{i}.intf = [0*intf (1:length(intf))' intf];
    
    % find elcon for neighboring elements
    meshp{i}.elconintf = zeros(npf,length(intf));
    for j = 1:length(intf)
        fj = intf(j);
        ej = mesh.f(fj,end-1:end);
        
        kf = mesh.t2f(ej,:);      % obtain neighboring faces 
        i1 = (kf(1,:)==fj);  % obtain the index of face i in the 1st element
        i2 = (kf(2,:)==fj);  % obtain the index of face i in the 2nd element            
                        
        j1 = elcon(:,i1,ej(1)) - (fj-1)*npf;        
        j2 = elcon(:,i2,ej(2)) - (fj-1)*npf;        
        
        if ismember(ej(1),meshp{i}.elempart)            
            meshp{i}.elconintf(:,j) = j2;
            %[j1 meshp{i}.elcon(:,i1,ej(1)==meshp{i}.elempart) j2]
        else            
            meshp{i}.elconintf(:,j) = j1;            
            %[j2 meshp{i}.elcon(:,i2,ej(2)==meshp{i}.elempart) j1]
        end                
    end
    
    meshp{i}.elcon = reshape(meshp{i}.elcon,[npf*nfe ne]);    
end

% determine neighboring subdomains
for i = 1:nproc        
    for j = 1:nproc
        if j~= i
            in = ismember(meshp{i}.intf(:,end), meshp{j}.intf(:,end));
            meshp{i}.intf(in,1) = j;
        end
    end
end

% entsend, entsendpts, nbsd
% entrecv, entrecvpts, nbsd
for i = 1:nproc        
    meshp{i}.nbsd = unique(meshp{i}.intf(:,1));
    [~,in] = sortrows(meshp{i}.intf(:,1:2));
    meshp{i}.intf = meshp{i}.intf(in,:);
    
    meshp{i}.entsend = meshp{i}.intf(:,2);
    for j = 1:length(meshp{i}.nbsd)
        meshp{i}.entsendpts(j) = length(find(meshp{i}.intf(:,1)==meshp{i}.nbsd(j)));    
    end
    meshp{i}.entrecv = meshp{i}.entsend;
    meshp{i}.entrecvpts = meshp{i}.entsendpts;
end



%     for j = 1:length(meshp{i}.nbsd)
%         meshp{i}.entsendpts(j) = length(find(meshp{i}.entsend(:,1)==meshp{i}.nbsd(j)));
%         meshp{i}.entrecvpts(j) = length(find(meshp{i}.entrecv(:,1)==meshp{i}.nbsd(j)));
%         meshp{i}.elemsendpts(j) = length(find(meshp{i}.elemsend(:,1)==meshp{i}.nbsd(j)));
%         meshp{i}.elemrecvpts(j) = length(find(meshp{i}.elemrecv(:,1)==meshp{i}.nbsd(j)));
%     end

%     /* copy some portion of sys.x to buffsend */
%     for (i=0; i<sys.nentsend; i++) {
%         ri = sys.entsend[i];
%         DCOPY(&bsz, &sys.x[bsz*ri], &inc, &sys.buffsend[bsz*i], &inc);
%     }
% 
%     /* non-blocking send */
%     Int nsend, psend = 0;
%     for (n=0; n<sys.nnbsd; n++) {
%         neighbor = sys.nbsd[n];
%         nsend = sys.entsendpts[n]*bsz;
%         if (nsend>0) {
%             MPI_Isend(&sys.buffsend[psend], nsend, MPI_DOUBLE, neighbor, 0,
%                    MPI_COMM_WORLD, &sys.requests[request_counter]);
%             psend += nsend;
%             request_counter += 1;
%         }
%     }
% 
%     /* non-blocking receive */
%     Int nrecv, precv = 0;
%     for (n=0; n<sys.nnbsd; n++) {
%         neighbor = sys.nbsd[n];
%         nrecv = sys.entrecvpts[n]*bsz;
%         if (nrecv>0) {
%             MPI_Irecv(&sys.buffrecv[precv], nrecv, MPI_DOUBLE, neighbor, 0,
%                    MPI_COMM_WORLD, &sys.requests[request_counter]);
%             precv += nrecv;
%             request_counter += 1;
%         }
%     }

%     /* copy buffrecv to sys.x */
%     for (i=0; i<sys.nentrecv; i++) {
%         ri = sys.entrecv[i];
%         DCOPY(&bsz, &sys.buffrecv[bsz*i], &inc, &sys.x[bsz*ri], &inc);
%     }

%      n = meshp{i}.entpartpts(1);     
%      for j = 1:n
%          fj = meshp{i}.intf(j,end);
%          for k = 1:nproc
%              if k~=i
%                  
%              end
%          end
%      end
     
%     intf = meshp{i}.entpart(1:n);    
%     meshp{i}.entsend = [0*intf (1:n)' intf];
    
%         % find all faces that are connected to interface faces
%     f2f = f2f(:,j1);
%     i1 = find(f2f(:)>0);
%     in = ones(size(f2f));
%     in(i1) = ismember(f2f(i1),meshp{i}.entpart);
%     f0 = unique(f2f(in==0));  % faces do not belong to subdomain i
%     f0
%     pause
%     
%     % store faces received from neighboring subdomains to perform the matrix-vector product
%     meshp{i}.entrecv = [0*f0 length(meshp{i}.entpart)+(1:1:length(f0))' f0];
%     for k = 1:nproc
%         if k ~= i
%             in = ismember(f0,meshp{k}.entpart);
%             meshp{i}.entrecv(in,1) = k;
%         end
%     end
%     meshp{i}.entrecv
%     meshp{i}.entrecv = unique(meshp{i}.entrecv,'rows');    
%     meshp{i}.nbsd = unique(meshp{i}.entrecv(:,1))';        
%     if meshp{i}.nbsd(1)==0
%         error('something wrong');
%     end



