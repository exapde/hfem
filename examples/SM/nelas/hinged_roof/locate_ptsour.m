function [ app ] = locate_ptsour( mesh, app )
%LOCATE_PTSOUR Locate the face to which each point source belongs, and,
% on that face, the closest node to the point source.

npts = size(app.ptsource.globpos,1);
app.ptsource.els = zeros(npts,1);
app.ptsource.nf = zeros(npts,1);
app.ptsource.node = zeros(npts,1);

% Case of a point source located at the surface 
if app.ptsource.surf==1
    for m=1:npts
        bm = app.ptsource.ibpl(m);
        pos= app.ptsource.globpos(m,:);
        ib = find(mesh.f(:,end)==-bm);
        for i = 1:length(ib)
            ni = ib(i);
            fi = mesh.f(ni,:);
            pi = mesh.p(fi(1:end-2),:);
            isin = 0;
            for j=1:size(pi,1)
                if norm(pi(j,:) - pos)<1e-7
                    isin = 1;
                    break;
                end
            end
            if isin==1
                xe = fi(end-1);
                xf = find(mesh.t2f(xe,:)==ni);
                pf = mesh.dgnodes(mesh.perm(:,xf),:,xe);
                for j=1:size(pf,1)
                    if norm(pf(j,:) - pos)<1e-7
                        jf = j;
                        break;
                    end
                end
                break;
            end
        end
        app.ptsource.els(m) = xe;
        app.ptsource.nf(m)  = xf;
        app.ptsource.node(m) = jf;
    end
else % Case of a Volumetric point Source
    for m=1:npts
        pos= app.ptsource.globpos(m,:);
        for nel=1:mesh.ne
            for j=1:size(mesh.dgnodes,1)
                pe =mesh.dgnodes(j,:,nel);
                if norm(pe - pos)<1e-7
                    jf = j;
                    xe = nel;
                    break;
                end
            end
        end
        app.ptsource.els(m) = xe;
        app.ptsource.node(m) = jf;
    end
end


