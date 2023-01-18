function [pb,tb] = mkbeamfrom2plates(p1, t1, th1, p2, t2, th2)

in = xiny(p1, p2);
i1 = find(in>0);
% figure(1);clf; hold on;
% boundaryplot(p1,t1,{'facecolor','g','edgecolor','k','Linew',1}); 
% boundaryplot(p2,t2,{'facecolor','g','edgecolor','k','Linew',1}); 
% q = p1(i1,:);
% plot3(q(:,1),q(:,2),q(:,3),'or','MarkerSize',6,'LineWidth',2);
% view(3);

e1 = beamline(t1, i1);
nd = size(p1,2);
n1 = length(e1);
pb = zeros(4,nd,n1);
tb = zeros(n1-1,8);
for i = 1:n1
    q = p1(e1(i),:);
    d = sqrt((p1(:,1)-q(1)).^2+(p1(:,2)-q(2)).^2+(p1(:,3)-q(3)).^2);
    [~,j1] = min(abs(d-th1));    
    d = sqrt((p2(:,1)-q(1)).^2+(p2(:,2)-q(2)).^2+(p2(:,3)-q(3)).^2);
    [~,j2] = min(abs(d-th2));
    q1 = p1(j1,:);
    q2 = p2(j2,:);
    q3 = q1 + q2 - q;    
    pb(:,:,i) = [q; q1; q3; q2];
    if i<n1
        tb(i,:)  = (i-1)*4+[1 2 3 4 5 6 7 8];
    end
end
pb = reshape(permute(pb,[1 3 2]),[4*n1 nd]);

% figure(1);clf; hold on;
% boundaryplot(p1,t1,{'facecolor','g','edgecolor','k','Linew',1}); 
% boundaryplot(p2,t2,{'facecolor','g','edgecolor','k','Linew',1}); 
% boundaryplot(pb,tb,{'facecolor','r','edgecolor','k','Linew',1}); 
% axis equal; axis tight;
% view([75 27]);        

function e1 = beamline(t1, i1) 

n1 = length(i1);
e1 = zeros(n1,3);
for i = 1:n1
    k = i1(i);  % node k
    [I,~] = find(t1==k);
    in = t1(I,:);
    in = unique(in(:));
    im = setdiff(intersect(in,i1),k); % find neighboring nodes that belong to i1
    if length(im)==1
        e1(i,:) = [k im(1) 0];
    elseif length(im)==2
        e1(i,:) = [k im(1) im(2)];
    else
        error('something wrong');
    end    
end

c1 = zeros(n1,1);
if n1>2
    in = find(e1(:,3)==0);
    c1(1) = e1(in(1),1);
    c1(n1) = e1(in(2),1);    
    for i = 2:n1-1
        k = c1(i-1);
        m = e1(:,1)==k;
        if ismember(e1(m,2),c1(1:i-1))==0
            c1(i) = e1(m,2);
        else
            c1(i) = e1(m,3);
        end
    end    
else
    c1 = e1(:,1);
end
e1 = c1;

