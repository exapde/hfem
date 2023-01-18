function hh=meshplot(mesh,opts,expr)
%MESHPLOT  Plot Mesh Structure (with straight edges)
%    MESHPLOT(MESH,[OPTS])
%
%    MESH:       Mesh structure
%    OPTS:       (logical)
%      OPTS(1):  Plot dgnodes (default=0)
%      OPTS(2):  Plot triangle numbers (default=0)
%      OPTS(3):  Plot node numbers (default=0)


if nargin<2 | isempty(opts), opts=[0]; end
if length(opts)<2, opts=[opts,0]; end
if length(opts)<3, opts=[opts,0]; end
if nargin<3, expr=[]; end

p=mesh.p;
t=mesh.t;

dim=size(p,2);

if dim < 2 | dim > 3,
    error('Only can handle dim=2 or dim=3');
end

hh=[];
if dim == 2
    pars={'facecolor',[.8,1,.8],'edgecolor','k','Linew',1};
    clf,hh=[hh;patch('faces',t,'vertices',p,pars{:})];
    view(2),axis equal
    
elseif dim ==3
    f=mesh.f;
    if ~isempty(expr) & ~ischar(expr)
        tri1=[];
        for it=1:size(f,1)
            bndnbr=-f(it,5);
            if bndnbr>0
                if ismember(bndnbr,expr)
                    tri1=[tri1;f(it,1:3)];
                end
            end
        end
    else
        tri1=surftri(p,t);
        if ~isempty(expr)
            incl=find(eval(expr));
            tincl=any(ismember(t,incl),2);
            t1=t(tincl,:);
            tri1=tri1(any(ismember(tri1,incl),2),:);
            tri2=surftri(p,t1);
            tri2=setdiff(tri2,tri1,'rows');
            h=trimesh(tri2,p(:,1),p(:,2),p(:,3));
            hh=[hh;h];
            set(h,'facecolor',[.6,.8,1],'edgecolor','k');
            hold on
        end
    end
    h=trimesh(tri1,p(:,1),p(:,2),p(:,3));
    hh=[hh;h];
    hold off
    set(h,'facecolor',[.8,1,.8],'edgecolor','k');
    axis equal
    fancycamera
end

if opts(1)
    if dim == 2       
            xx=squeeze(mesh.nodgeom(:,1,:));
            yy=squeeze(mesh.nodgeom(:,2,:));
            zz=0*xx;
    elseif dim ==3
            xx=squeeze(mesh.nodgeom(:,1,:));
            yy=squeeze(mesh.nodgeom(:,2,:));
            zz=squeeze(mesh.nodgeom(:,3,:));
    end
    line(xx(:),yy(:),zz(:),'lines','n','marker','.','markersi',16,'col','b');
end


if opts(2)
    if dim == 2    
        for it=1:size(t,1)
            pmid=mean(p(t(it,:),:),1);
            txtpars={'fontname','times','fontsize',20,'horizontala','center','BackgroundColor',[1,1,1]};
            text(pmid(1),pmid(2),pmid(3),num2str(it),txtpars{:});
        end
    end
end

if opts(3)
    if dim == 2
            for i=1:size(p,1)
                txtpars={'fontname','times','fontsize',20,'fontweight','bold', ...
                    'horizontala','center','col','w', 'BackgroundColor',[0.5,0.5,0.5]};
                text(p(i,1),p(i,2),num2str(i),txtpars{:});
            end
    end
end
        
if nargout<1, clear hh; end
