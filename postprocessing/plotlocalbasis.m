function plotlocalbasis(porder, dim, elemtype, nodetype)


[plc,tlc] = masternodes(porder,dim,elemtype,nodetype);

if elemtype==0
    figure(1); clf;
    plot([0 1 0 0],[0 0 1 0],'k-','LineWidth',1);
    hold on;
    plot(plc(:,1),plc(:,2),'o','LineWidth',2,'MarkerSize',8);
    set(gca,'FontSize',18);
    axis equal;
    axis tight;
    xlabel('x','FontSize',20); ylabel('y','FontSize',20);
    set(gca, 'LooseInset', get(gca, 'TightInset'));
    fn = ['trinode' num2str(porder)];
    print('-dpng',fn);

%     [p,t]=trianglemesh(9);
%     mesh = mkmesh(p,t,porder,{'true'},elemtype,nodetype);
%     gpvl = reshape(permute(mesh.dgnodes,[1 3 2]),[],2);
%     shapvl = mkshape(porder,plc,gpvl,elemtype);
%     
%     for i = 1:size(plc,1)
%         u = reshape(shapvl(i,:,1),[size(plc,1) 1 mesh.ne]);
%         figure(2);clf; 
%         scaplot(mesh,u,[],2);        
%         axis off; colormap jet; colorbar('FontSize',12);      
%         set(gca, 'LooseInset', get(gca, 'TightInset'));
%         fn = ['trishape' num2str(porder) num2str(i)];
%         print('-dpng',fn);
%     end    

    [p,t]=squaremesh(9);
    mesh = mkmesh(p,t,porder,{'true'},elemtype,nodetype);    
    figure(3);clf; 
    meshplot(mesh);
    hold on;
    for i = 1:size(t,1)
        ti = t(i,:);
        pi = p(ti,:);
        in = 0;
        for j = 1:3
            if norm(pi(j,:)-[0.5 0.5+1/8])<1e-8
                in = 1;
                break;
            end
        end
        if in==1
            patch('vertices',p,'faces',ti,'cdata',ones(3,1), ...
                       'facecol','r','edgec','k');            
        end
    end
    gpvl = reshape(permute(mesh.dgnodes,[1 3 2]),[],2);
    plot(gpvl(:,1),gpvl(:,2),'ob');
    axis off; axis tight;
    fn = ['trisupport2' num2str(porder)];
    print('-dpng',fn);
end


xi = linspace(0,1,201)';
[plc,tlc] = masternodes(porder,1,0,1);
shapvl = mkshape(porder,plc,xi,0);
figure(1);clf; 
plot(xi,shapvl(:,:,1),'LineWidth',2);
set(gca,'FontSize',18);
axis tight;
xlabel('\xi','FontSize',22); 
set(gca, 'LooseInset', get(gca, 'TightInset'));
fn = ['lineshape' num2str(porder)];
print('-dpng',fn);


 