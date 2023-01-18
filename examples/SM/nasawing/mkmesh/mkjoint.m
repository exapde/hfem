function [joinp, joint, edge] = mkjoint(center, tangent, thickness, isplot)

normal = [-tangent(2,:); tangent(1,:)];

nn = length(thickness);
if nn==2    
    [joinp, joint, edge] =  case2(center, normal, tangent, thickness);    
elseif nn==3                
    [joinp, joint, edge] =  case3(center, normal, tangent, thickness);       
elseif nn==4
    [joinp, joint, edge] =  case4(center, normal, tangent, thickness);                
end

if nargin<4
    isplot=0;
end

if isplot==1
figure(1); clf;    
hold on;
plot(center(1),center(2),'*','MarkerSize',8);
for i = 1:size(joint,1)            
    k = joint(i,:)>0;
    p = joinp(:,joint(i,k));
    x = [p p(:,1)]; x= x';
    plot(x(:,1),x(:,2),'-r','linewidth',2.5);            
end
for i = 1:nn
    x =[center center+5*tangent(:,i)*thickness(i)]; x=x';
    plot(x(:,1),x(:,2),'-g','linewidth',0.5);                    
    x = [edge(:,1,i) edge(:,1,i)+ 5*tangent(:,i)*thickness(i)];  x=x';
    plot(x(:,1),x(:,2),'-b','linewidth',2.5);
    x = [edge(:,2,i) edge(:,2,i)+ 5*tangent(:,i)*thickness(i)];  x=x';
    plot(x(:,1),x(:,2),'-b','linewidth',2.5);
end                
axis equal;
end

function [joinp, joint, edge] =  case2(center, normal, tangent, thickness)

nd = length(center);
nn = length(thickness);
q = zeros(nd,2,nn);
p = zeros(nd,2,nn);    
for i = 1:nn
    q(:,1,i) = center + 0.5*normal(:,i)*thickness(i);
    q(:,2,i) = center - 0.5*normal(:,i)*thickness(i); 
end    

tn = dot(tangent(:,1),tangent(:,2));    
if tn<=-0.9999
    th = (thickness(1)+thickness(2))/4;
    for i = 1:nn    
        p(:,1,i) = q(:,1,i) + tangent(:,i)*th;
        p(:,2,i) = q(:,2,i) + tangent(:,i)*th;    
        joinp = [p(:,1,1) p(:,2,1) p(:,1,2) p(:,2,2)];
        edge1 = [p(:,1,1) p(:,2,1)];
        edge2 = [p(:,1,2) p(:,2,2)];           
    end
    joint = [1 2 3 4];
elseif tn<=-0.8
    p1 = lineintersection(q(:,1,1), tangent(:,1), q(:,1,2), tangent(:,2));
    p2 = lineintersection(q(:,1,1), tangent(:,1), q(:,2,2), tangent(:,2));
    p3 = lineintersection(q(:,2,1), tangent(:,1), q(:,1,2), tangent(:,2));
    p4 = lineintersection(q(:,2,1), tangent(:,1), q(:,2,2), tangent(:,2));        
    
    p0 = segmentintersection(p1, normal(:,2),  center, tangent(:,2));
    if p0(2)>0
        pa = lineintersection(p1, normal(:,2), p2, tangent(:,2));
        pb = lineintersection(p4, normal(:,1), p2, tangent(:,1));
        joinp = [p1 pa p2 p3 pb p4];    
        joint = [1 2 3 4; 3 5 6 4];      
        edge1 = [pb p4];            
        edge2 = [p1 pa];        
    else
        pa = lineintersection(p1, normal(:,1), p3, tangent(:,1));
        pb = lineintersection(p4, normal(:,2), p3, tangent(:,2));
        joinp = [p1 pa p2 p3 pb p4];    
        joint = [1 2 4 3; 4 5 6 3];        
        edge1 = [p1 pa];
        edge2 = [pb p4];            
    end        
else 
    p1 = lineintersection(q(:,1,1), tangent(:,1), q(:,1,2), tangent(:,2));
    p2 = lineintersection(q(:,1,1), tangent(:,1), q(:,2,2), tangent(:,2));
    p3 = lineintersection(q(:,2,1), tangent(:,1), q(:,1,2), tangent(:,2));
    p4 = lineintersection(q(:,2,1), tangent(:,1), q(:,2,2), tangent(:,2));        
    joinp = [p1 p2 p3 p4];        
    
    dist = zeros(4,1);
    for i = 1:4
        dist(i) = norm(joinp(:,i)-center);
    end    
    tt = center+max(dist)*0.5*(tangent(:,1)+tangent(:,2));    
    for i = 1:4
        dist(i) = norm(joinp(:,i)-tt);
    end    
    [~,j] = min(dist);
    if j==1
        pa = lineintersection(p1, normal(:,1), q(:,2,1), tangent(:,1));
        pb = lineintersection(p1, normal(:,2), q(:,2,2), tangent(:,2));
        joinp = [pa p1 pb p4];
    elseif j==2
        pa = lineintersection(p2, normal(:,1), q(:,2,1), tangent(:,1));
        pb = lineintersection(p2, normal(:,2), q(:,1,2), tangent(:,2));
        joinp = [pa p2 pb p3];
        % Test if swapped/degenerated quad
        if dot(normal(:,2),p3-pa)<0.
            joinp = [p3 p2 pb pa];
        end
        if norm(p3-pa)<0.25*thickness(2)
            joinp = [p4 p2 pb pa];
        end
    elseif j==3
        pa = lineintersection(p3, normal(:,1), q(:,1,1), tangent(:,1));
        pb = lineintersection(p3, normal(:,2), q(:,2,2), tangent(:,2));
        joinp = [pa p3 pb p2];
    else
        pa = lineintersection(p4, normal(:,1), q(:,1,1), tangent(:,1));
        pb = lineintersection(p4, normal(:,2), q(:,1,2), tangent(:,2));
        joinp = [pa p4 pb p1];
    end       
    edge1 = joinp(:,[1 2]);
    edge2 = joinp(:,[2 3]);
    joint = [1 2 3 4];
end 

edge = zeros(nd,2,2);
edge(:,:,1) = edge1;
edge(:,:,2) = edge2;

function joinp = joinm(center, tangent, normal, thickness, p1, p2, m)

if norm(p1+tangent(:,m)-center)>norm(p2+tangent(:,m)-center)
    pa = p1;
    pb = p1+thickness(m)*tangent(:,m);
    pc = lineintersection(pb, normal(:,m), p2, tangent(:,m));
    pd = p2;    
    joinp = [pa pd pc pb];
else
    pa = p2;
    pb = p2+thickness(m)*tangent(:,m);
    pc = lineintersection(pb, normal(:,m), p1, tangent(:,m));
    pd = p1;
    joinp = [pd pa pb pc];
end

function [joinp, joint, edge] =  case3(center, normal, tangent, thickness)

nd = length(center);
nn = length(thickness);

% First step : ordering tangent vectors with respect to their angle with
% the first tangent vector.
tn = zeros(1,nn);
tn(1) = 0.;
for i = 2:nn
    % Relative angle to first tangent vector = [0,1] in local coordinates
    tn(i) = acos(tangent(2,i));
    if tangent(1,i)>0
        tn(i) = 2*pi-tn(i);
    end
    % Previous ordering...
    %tn(i) = acos(dot([1; 0],tangent(:,i))/norm(tangent(:,i)));       
    %if tangent(2,i)<0
    %    tn(i) = 2*pi-tn(i);
    %end
end
[~,ind] = sort(tn);
thickness = thickness(ind);
normal = normal(:,ind);
tangent = tangent(:,ind);

q = zeros(nd,2,nn);    
for i = 1:nn
    q(:,1,i) = center + 0.5*normal(:,i)*thickness(i);
    q(:,2,i) = center - 0.5*normal(:,i)*thickness(i);
end    

tn = zeros(3,1);
tn(1) = dot(tangent(:,1),tangent(:,2));    
tn(2) = dot(tangent(:,1),tangent(:,3));   
tn(3) = dot(tangent(:,2),tangent(:,3));   
[tm,jm] = min(tn);
if tm<=-0.99
    if jm==1
%         p2 = lineintersection(q(:,1,2), tangent(:,2), q(:,2,3), tangent(:,3));
%         p3 = lineintersection(q(:,1,3), tangent(:,3), q(:,2,1), tangent(:,1));
%         joinp3 = joinm(center, tangent, normal, thickness, p2, p3, 3);                
%         p1 = lineintersection(p2, normal(:,2), q(:,2,2), tangent(:,2));
%         p4 = lineintersection(p3, normal(:,1), q(:,1,1), tangent(:,1));
%         
%         joinp = [p1 p2 p3 p4 joinp3(:,3:4)];
%         joint = [1 2 3 4; 2 3 5 6];
%         edge = zeros(nd,2,nn);
%         edge(:,:,1) = joinp(:,[3 4]);
%         edge(:,:,2) = joinp(:,[1 2]);
%         edge(:,:,3) = joinp(:,[5 6]);
%         [~,ind] = sort(ind);
%         edge = edge(:,:,ind);        
        p2 = lineintersection(q(:,1,2), tangent(:,2), q(:,2,3), tangent(:,3));
        p3 = lineintersection(q(:,1,3), tangent(:,3), q(:,2,1), tangent(:,1));
        p1 = lineintersection(p2, normal(:,2), q(:,2,2), tangent(:,2));
        p4 = lineintersection(p3, normal(:,1), q(:,1,1), tangent(:,1));
        
        joinp = [p1 p2 p3 p4];
        joint = [1 2 3 4];
        edge = zeros(nd,2,nn);
        edge(:,:,1) = joinp(:,[3 4]);
        edge(:,:,2) = joinp(:,[1 2]);
        edge(:,:,3) = joinp(:,[2 3]);
        [~,ind] = sort(ind);
        edge = edge(:,:,ind);        
    end
    if jm==2
%         p1 = lineintersection(q(:,1,1), tangent(:,1), q(:,2,2), tangent(:,2));
%         p2 = lineintersection(q(:,1,2), tangent(:,2), q(:,2,3), tangent(:,3));
%         joinp2 = joinm(center, tangent, normal, thickness, p1, p2, 2);        
%         p3 = lineintersection(p1, normal(:,1), q(:,2,1), tangent(:,1));
%         p4 = lineintersection(p2, normal(:,3), q(:,1,3), tangent(:,3));
%         
%         joinp = [p1 p2 p4 p3 joinp2(:,3:4)];
%         joint = [1 2 3 4; 1 2 5 6];
%         edge = zeros(nd,2,nn);
%         edge(:,:,1) = joinp(:,[1 4]);
%         edge(:,:,2) = joinp(:,[5 6]);
%         edge(:,:,3) = joinp(:,[2 3]);
%         [~,ind] = sort(ind);
%         edge = edge(:,:,ind);        
        p1 = lineintersection(q(:,1,1), tangent(:,1), q(:,2,2), tangent(:,2));
        p2 = lineintersection(q(:,1,2), tangent(:,2), q(:,2,3), tangent(:,3));       
        p3 = lineintersection(p1, normal(:,1), q(:,2,1), tangent(:,1));
        p4 = lineintersection(p2, normal(:,3), q(:,1,3), tangent(:,3));
        
        joinp = [p1 p2 p4 p3];
        joint = [1 2 3 4];
        edge = zeros(nd,2,nn);
        edge(:,:,1) = joinp(:,[1 4]);
        edge(:,:,2) = joinp(:,[1 2]);
        edge(:,:,3) = joinp(:,[2 3]);
        [~,ind] = sort(ind);
        edge = edge(:,:,ind);        
    end
    if jm==3
%         p1 = lineintersection(q(:,1,1), tangent(:,1), q(:,2,2), tangent(:,2));
%         p3 = lineintersection(q(:,1,3), tangent(:,3), q(:,2,1), tangent(:,1));
%         joinp1 = joinm(center, tangent, normal, thickness, p1, p3, 1);        
%         p2 = lineintersection(p1, normal(:,2), q(:,1,2), tangent(:,2));
%         p4 = lineintersection(p3, normal(:,3), q(:,2,3), tangent(:,3));
%         
%         joinp = [p1 p2 p4 p3 joinp1(:,3:4)];
%         joint = [1 2 3 4; 1 4 5 6];
%         edge = zeros(nd,2,nn);
%         edge(:,:,1) = joinp(:,[5 6]);
%         edge(:,:,2) = joinp(:,[1 2]);
%         edge(:,:,3) = joinp(:,[3 4]);
%         [~,ind] = sort(ind);
%         edge = edge(:,:,ind);        
        p1 = lineintersection(q(:,1,1), tangent(:,1), q(:,2,2), tangent(:,2));
        p3 = lineintersection(q(:,1,3), tangent(:,3), q(:,2,1), tangent(:,1));
        p2 = lineintersection(p1, normal(:,2), q(:,1,2), tangent(:,2));
        p4 = lineintersection(p3, normal(:,3), q(:,2,3), tangent(:,3));
        
        joinp = [p1 p2 p4 p3];
        joint = [1 2 3 4];
        edge = zeros(nd,2,nn);
        edge(:,:,1) = joinp(:,[1 4]);
        edge(:,:,2) = joinp(:,[1 2]);
        edge(:,:,3) = joinp(:,[3 4]);
        [~,ind] = sort(ind);
        edge = edge(:,:,ind);        
    end
else
%     p1 = lineintersection(q(:,1,1), tangent(:,1), q(:,2,2), tangent(:,2));
%     p2 = lineintersection(q(:,1,2), tangent(:,2), q(:,2,3), tangent(:,3));
%     p3 = lineintersection(q(:,1,3), tangent(:,3), q(:,2,1), tangent(:,1));
%     joinp1 = joinm(center, tangent, normal, thickness, p1, p3, 1);
%     joinp2 = joinm(center, tangent, normal, thickness, p1, p2, 2);
%     joinp3 = joinm(center, tangent, normal, thickness, p2, p3, 3);
% 
%     joinp = [p1 p2 p3 joinp1(:,3:4) joinp2(:,3:4) joinp3(:,3:4)];
%     joint = [1 2 3 0; 1 3 4 5; 1 2 6 7; 2 3 8 9];
% 
%     edge = zeros(nd,2,nn);
%     edge(:,:,1) = joinp1(:,3:4);
%     edge(:,:,2) = joinp2(:,3:4);
%     edge(:,:,3) = joinp3(:,3:4);
%     [~,ind] = sort(ind);
%     edge = edge(:,:,ind);
    p1 = lineintersection(q(:,1,1), tangent(:,1), q(:,2,2), tangent(:,2));
    p2 = lineintersection(q(:,1,2), tangent(:,2), q(:,2,3), tangent(:,3));
    p3 = lineintersection(q(:,1,3), tangent(:,3), q(:,2,1), tangent(:,1));

    edge = zeros(nd,2,nn);
    if jm==1
        if norm(p2-center)>norm(p3-center)
            p4 = lineintersection(p2, normal(:,3), q(:,1,3), tangent(:,3));            
            joinp = [p2 p1 p3 p4];
            edge(:,:,3) = joinp(:,[1 4]);
        else
            p4 = lineintersection(p3, normal(:,3), q(:,2,3), tangent(:,3));            
            joinp = [p2 p1 p3 p4];
            edge(:,:,3) = joinp(:,[3 4]);
        end        
        edge(:,:,1) = joinp(:,[2 3]);
        edge(:,:,2) = joinp(:,[1 2]);             
    elseif jm==2
        if norm(p1-center)>norm(p2-center)
            p4 = lineintersection(p1, normal(:,2), q(:,1,2), tangent(:,2));            
            joinp = [p1 p3 p2 p4];
            edge(:,:,2) = joinp(:,[1 4]);
        else
            p4 = lineintersection(p2, normal(:,2), q(:,2,2), tangent(:,2));            
            joinp = [p1 p3 p2 p4];
            edge(:,:,2) = joinp(:,[3 4]);
        end        
        edge(:,:,1) = joinp(:,[1 2]);
        edge(:,:,3) = joinp(:,[2 3]);    
    elseif jm==3
        if norm(p3-center)>norm(p1-center)
            p4 = lineintersection(p3, normal(:,1), q(:,1,1), tangent(:,1));            
            joinp = [p3 p2 p1 p4];
            edge(:,:,1) = joinp(:,[1 4]);
        else
            p4 = lineintersection(p1, normal(:,1), q(:,2,1), tangent(:,1));            
            joinp = [p3 p2 p1 p4];
            edge(:,:,1) = joinp(:,[3 4]);
        end        
        edge(:,:,2) = joinp(:,[2 3]);
        edge(:,:,3) = joinp(:,[1 2]);                     
    end        
    
    joint = [1 2 3 4];
    [~,ind] = sort(ind);
    edge = edge(:,:,ind);
end

function [joinp, joint, edge] =  case4(center, normal, tangent, thickness)

nd = length(center);
nn = length(thickness);

tn = zeros(1,4);
for i = 1:4
    tn(i) = acos(dot([1; 0],tangent(:,i))/norm(tangent(:,i)));       
    if tangent(2,i)<0
        tn(i) = 2*pi-tn(i);
    end
end
[~,ind] = sort(tn);
thickness = thickness(ind);
normal = normal(:,ind);
tangent = tangent(:,ind);        

q = zeros(nd,2,nn);    
for i = 1:nn
    q(:,1,i) = center + 0.5*normal(:,i)*thickness(i);
    q(:,2,i) = center - 0.5*normal(:,i)*thickness(i);
end    

tn = zeros(4,1);
tn(1) = dot(tangent(:,1),tangent(:,2));   
tn(2) = dot(tangent(:,2),tangent(:,3));   
tn(3) = dot(tangent(:,3),tangent(:,4));   
tn(4) = dot(tangent(:,4),tangent(:,1));   
%t13 = dot(tangent(:,1),tangent(:,3));   
%t24 = dot(tangent(:,4),tangent(:,2)); 
[tm,jm] = min(tn);
if (tm<=-0.99) %&& (tm<min([t13 t24]))
    if jm==1
        p2 = lineintersection(q(:,1,2), tangent(:,2), q(:,2,3), tangent(:,3));
        p3 = lineintersection(q(:,1,3), tangent(:,3), q(:,2,4), tangent(:,4));
        p4 = lineintersection(q(:,1,4), tangent(:,4), q(:,2,1), tangent(:,1));
        p1 = lineintersection(p2, normal(:,2), q(:,2,2), tangent(:,2));
        p5 = lineintersection(p4, normal(:,1), q(:,1,1), tangent(:,1));
        p6 = 0.5*(p1+p5);
        joinp3 = joinm(center, tangent, normal, thickness, p2, p3, 3);
        joinp4 = joinm(center, tangent, normal, thickness, p3, p4, 4);
        
        joinp = [p1 p2 p3 p4 p5 p6 joinp3(:,3:4) joinp4(:,3:4)];
        joint = [6 5 4 3; 1 6 3 2; 2 3 7 8; 3 4 9 10];

        edge = zeros(nd,2,nn);
        edge(:,:,1) = [p4 p5];
        edge(:,:,2) = [p1 p2];
        edge(:,:,3) = joinp3(:,3:4);
        edge(:,:,4) = joinp4(:,3:4);
        [~,ind] = sort(ind);
        edge = edge(:,:,ind);        
    end
    if jm==2
        p1 = lineintersection(q(:,1,1), tangent(:,1), q(:,2,2), tangent(:,2));    
        p3 = lineintersection(q(:,1,3), tangent(:,3), q(:,2,4), tangent(:,4));
        p4 = lineintersection(q(:,1,4), tangent(:,4), q(:,2,1), tangent(:,1));    
        p2 = lineintersection(p1, normal(:,2), q(:,1,2), tangent(:,2));
        p5 = lineintersection(p3, normal(:,3), q(:,2,3), tangent(:,3));
        p6 = 0.5*(p2+p5);
        joinp1 = joinm(center, tangent, normal, thickness, p1, p4, 1);
        joinp4 = joinm(center, tangent, normal, thickness, p3, p4, 4);
        
        joinp = [p1 p2 p3 p4 p5 p6 joinp1(:,3:4) joinp4(:,3:4)];
        joint = [1 4 7 8; 4 1 2 6; 3 4 6 5; 3 4 9 10];

        edge = zeros(nd,2,nn);
        edge(:,:,1) = joinp1(:,3:4);
        edge(:,:,2) = [p1 p2];
        edge(:,:,3) = [p5 p3];
        edge(:,:,4) = joinp4(:,3:4);
        [~,ind] = sort(ind);
        edge = edge(:,:,ind);        
    end
    if jm==3
        p1 = lineintersection(q(:,1,1), tangent(:,1), q(:,2,2), tangent(:,2));
        p2 = lineintersection(q(:,1,2), tangent(:,2), q(:,2,3), tangent(:,3));
        p4 = lineintersection(q(:,1,4), tangent(:,4), q(:,2,1), tangent(:,1));        
        p3 = lineintersection(p2, normal(:,3), q(:,1,3), tangent(:,3));
        p5 = lineintersection(p4, normal(:,4), q(:,2,4), tangent(:,4));
        p6 = 0.5*(p3+p5);
        joinp1 = joinm(center, tangent, normal, thickness, p1, p4, 1);
        joinp2 = joinm(center, tangent, normal, thickness, p1, p2, 2);
        
        joinp = [p1 p2 p3 p4 p5 p6 joinp1(:,3:4) joinp2(:,3:4)];
        joint = [1 4 7 8; 1 2 9 10; 2 1 6 3; 1 4 5 6];

        edge = zeros(nd,2,nn);
        edge(:,:,1) = joinp1(:,3:4);
        edge(:,:,2) = joinp2(:,3:4);
        edge(:,:,3) = [p2 p3];
        edge(:,:,4) = [p4 p5];
        [~,ind] = sort(ind);
        edge = edge(:,:,ind);        
    end
    if jm==4
        p1 = lineintersection(q(:,1,1), tangent(:,1), q(:,2,2), tangent(:,2));
        p2 = lineintersection(q(:,1,2), tangent(:,2), q(:,2,3), tangent(:,3));
        p3 = lineintersection(q(:,1,3), tangent(:,3), q(:,2,4), tangent(:,4));        
        p4 = lineintersection(p1, normal(:,1), q(:,2,1), tangent(:,1));
        p5 = lineintersection(p3, normal(:,4), q(:,1,4), tangent(:,4));
        p6 = 0.5*(p4+p5);        
        joinp2 = joinm(center, tangent, normal, thickness, p1, p2, 2);
        joinp3 = joinm(center, tangent, normal, thickness, p2, p3, 3);
        
        joinp = [p1 p2 p3 p4 p5 p6 joinp2(:,3:4) joinp3(:,3:4)];
        joint = [4 6 2 1; 1 2 7 8; 2 3 9 10; 5 3 2 6];

        edge = zeros(nd,2,nn);
        edge(:,:,1) = [p1 p4];
        edge(:,:,2) = joinp2(:,3:4);
        edge(:,:,3) = joinp3(:,3:4);
        edge(:,:,4) = [p3 p5];
        [~,ind] = sort(ind);
        edge = edge(:,:,ind);                
    end
else
%     p1 = lineintersection(q(:,1,1), tangent(:,1), q(:,2,2), tangent(:,2));
%     p2 = lineintersection(q(:,1,2), tangent(:,2), q(:,2,3), tangent(:,3));
%     p3 = lineintersection(q(:,1,3), tangent(:,3), q(:,2,4), tangent(:,4));
%     p4 = lineintersection(q(:,1,4), tangent(:,4), q(:,2,1), tangent(:,1));
%     joinp1 = joinm(center, tangent, normal, thickness, p1, p4, 1);
%     joinp2 = joinm(center, tangent, normal, thickness, p1, p2, 2);
%     joinp3 = joinm(center, tangent, normal, thickness, p2, p3, 3);
%     joinp4 = joinm(center, tangent, normal, thickness, p3, p4, 4);
% 
%     joinp = [p1 p2 p3 p4 joinp1(:,3:4) joinp2(:,3:4) joinp3(:,3:4) joinp4(:,3:4)];
%     joint = [1 2 3 4; 1 4 5 6; 1 2 7 8; 2 3 9 10; 3 4 11 12];
% 
%     edge = zeros(nd,2,nn);
%     edge(:,:,1) = joinp1(:,3:4);
%     edge(:,:,2) = joinp2(:,3:4);
%     edge(:,:,3) = joinp3(:,3:4);
%     edge(:,:,4) = joinp4(:,3:4);
%     [~,ind] = sort(ind);
%     edge = edge(:,:,ind);

    p1 = lineintersection(q(:,1,1), tangent(:,1), q(:,2,2), tangent(:,2));
    p2 = lineintersection(q(:,1,2), tangent(:,2), q(:,2,3), tangent(:,3));
    p3 = lineintersection(q(:,1,3), tangent(:,3), q(:,2,4), tangent(:,4));
    p4 = lineintersection(q(:,1,4), tangent(:,4), q(:,2,1), tangent(:,1));

    joinp = [p1 p2 p3 p4];
    joint = [1 2 3 4];

    edge = zeros(nd,2,nn);
    edge(:,:,1) = joinp(:,[1 4]);
    edge(:,:,2) = joinp(:,[1 2]);
    edge(:,:,3) = joinp(:,[2 3]);
    edge(:,:,4) = joinp(:,[3 4]);
    [~,ind] = sort(ind);
    edge = edge(:,:,ind);
end

return;

nn = 2;
nd = 2;
center = zeros(nd,1);
tangent = -1 + 2*rand(nd,nn);
for i = 1:nn
    tangent(:,i) = tangent(:,i)/norm(tangent(:,i));
end
%tangent = [1 -1; 0 0];
thickness = 0.1*(0.2+0.6*rand(1,nn));

%tangent = [0 0 -1; 1 -1 0];

if nd == 2
    normal = [-tangent(2,:); tangent(1,:)];
else
    normal = zeros(nn,nd);
    for i = 1:nn
        [~,j] = min(abs(tangent(:,i)));
        if j==1
            normal(2,i) = -tangent(3,i);
            normal(3,i) =  tangent(2,i);
        elseif i==2
            normal(1,i) =  tangent(3,i);
            normal(3,i) = -tangent(1,i);            
        elseif i==3
            normal(1,i) = -tangent(2,i);
            normal(2,i) =  tangent(1,i);                
        end                
        normal(:,i) = normal(:,i)/norm(normal(:,i));
    end    
end
joinp = mkjoint(center, tangent, thickness);



% function [joinp1,x1] = joinmn(q, tangent, m, n)
% 
% nd = size(tangent,1);
% t1 = segmentintersection(q(:,1,m), tangent(:,m), q(:,1,n), tangent(:,n));
% t2 = segmentintersection(q(:,1,m), tangent(:,m), q(:,2,n), tangent(:,n));
% t3 = segmentintersection(q(:,2,m), tangent(:,m), q(:,1,n), tangent(:,n));
% t4 = segmentintersection(q(:,2,m), tangent(:,m), q(:,2,n), tangent(:,n));                
% tt = [t1 t2 t3 t4];
% ind = find(tt(1,:)>0 & tt(2,:)>0);
% if isempty(ind)==0
%     pp = zeros(nd,4);
%     [~,j] = max(sum(tt(:,ind).^2+tt(:,ind).^2,1));
%     i = ind(j);    
%     th = tt(1,i);
%     pp(:,1) = q(:,1,m) + tangent(:,m)*th;
%     pp(:,2) = q(:,2,m) + tangent(:,m)*th;
%     th = tt(2,i);
%     pp(:,3) = q(:,1,n) + tangent(:,n)*th;
%     pp(:,4) = q(:,2,n) + tangent(:,n)*th;    
%     if norm(pp(:,1)-pp(:,3))<1e-12
%         joinp1 = pp(:,[2 1 4]);
%         x1 = pp(:,1);
%     elseif norm(pp(:,1)-pp(:,4))<1e-12
%         joinp1 = pp(:,[2 1 3]);
%         x1 = pp(:,1);
%     elseif norm(pp(:,2)-pp(:,3))<1e-12
%         joinp1 = pp(:,[1 2 4]);
%         x1 = pp(:,2);
%     elseif norm(pp(:,2)-pp(:,4))<1e-12
%         joinp1 = pp(:,[1 2 3]);
%         x1 = pp(:,2);    
%     end
% else
%     x1 = [];
%     joinp1 = [];
% end
% 

% function [joinp, joint, edge] =  case3(center, normal, tangent, thickness)
% 
% 
% nd = length(center);
% nn = length(thickness);
% tn = zeros(3,1);
% tn(1) = dot(tangent(:,1),tangent(:,2));    
% tn(2) = dot(tangent(:,1),tangent(:,3));   
% tn(3) = dot(tangent(:,2),tangent(:,3));   
% [~,j] = min(tn);
% if j==1
%     if thickness(1)<thickness(2)
%         ind = [1 2 3];
%     else
%         ind = [2 1 3];
%     end
% elseif j==2
%     if thickness(1)<thickness(3)
%         ind = [1 3 2];
%     else
%         ind = [3 1 2];
%     end
% elseif j==3
%     if thickness(2)<thickness(3)
%         ind = [2 3 1];
%     else
%         ind = [3 2 1];
%     end
% end        
% thickness = thickness(ind);
% normal = normal(:,ind);
% tangent = tangent(:,ind);        
% 
% q = zeros(nd,2,nn);    
% for i = 1:nn
%     q(:,1,i) = center + 0.5*normal(:,i)*thickness(i);
%     q(:,2,i) = center - 0.5*normal(:,i)*thickness(i);
% end    
% 
% % determine the nodes of the joint
% [joinp1,x1] = joinmn(q, tangent, 1, 3);
% [joinp2,x2] = joinmn(q, tangent, 2, 3);
% th = max(thickness);
% if norm(x1)>norm(x2)        
%     % determine 
%     p1 = lineintersection(joinp1(:,1), tangent(:,1), q(:,1,2), tangent(:,2));
%     p2 = lineintersection(joinp1(:,1), tangent(:,1), q(:,2,2), tangent(:,2));
%     p3 = lineintersection(joinp1(:,3), tangent(:,3), q(:,1,2), tangent(:,2));
%     p4 = lineintersection(joinp1(:,3), tangent(:,3), q(:,2,2), tangent(:,2));    
%     dist = zeros(4,1);
%     if isempty(p1)==0
%         same1 = sameside(q(:,1,2), q(:,2,2), center+tangent(:,2), p1);
%         same2 = sameside(center, center+tangent(:,2), joinp1(:,1), p1);
%         if (same1==1) && (same2==1)
%             dist(1) = norm(p1-q(:,1,2));
%         end
%     end
%     if isempty(p2)==0
%         same1 = sameside(q(:,1,2), q(:,2,2), center+tangent(:,2), p2);
%         same2 = sameside(center, center+tangent(:,2), joinp1(:,1), p2);    
%         if (same1==1) && (same2==1)
%             dist(2) = norm(p2-q(:,2,2));
%         end
%     end
%     if isempty(p3)==0
%         same1 = sameside(q(:,1,2), q(:,2,2), center+tangent(:,2), p3);
%         same2 = sameside(center, center+tangent(:,2), joinp1(:,3), p3);
%         if (same1==1) && (same2==1)
%             dist(3) = norm(p3-q(:,1,2));
%         end
%     end
%     if isempty(p4)==0
%         same1 = sameside(q(:,1,2), q(:,2,2), center+tangent(:,2), p4);
%         same2 = sameside(center, center+tangent(:,2), joinp1(:,3), p4);
%         if (same1==1) && (same2==1)
%             dist(4) = norm(p4-q(:,2,2));
%         end
%     end
%     th = max(dist);    
%         
%     joinp = [joinp1 q(:,1,2)+tangent(:,2)*th q(:,2,2)+tangent(:,2)*th];
%     ind = ind([1 3 2]);
% else
%     p1 = lineintersection(joinp2(:,1), tangent(:,2), q(:,1,1), tangent(:,1));
%     p2 = lineintersection(joinp2(:,1), tangent(:,2), q(:,2,1), tangent(:,1));
%     p3 = lineintersection(joinp2(:,3), tangent(:,3), q(:,1,1), tangent(:,1));
%     p4 = lineintersection(joinp2(:,3), tangent(:,3), q(:,2,1), tangent(:,1));    
%     dist = zeros(4,1);
%     if isempty(p1)==0
%         same1 = sameside(q(:,1,1), q(:,2,1), center+tangent(:,1), p1);
%         same2 = sameside(center, center+tangent(:,1), joinp2(:,1), p1);
%         if (same1==1) && (same2==1)
%             dist(1) = norm(p1-q(:,1,1));
%         end
%     end
%     if isempty(p2)==0
%         same1 = sameside(q(:,1,1), q(:,2,1), center+tangent(:,1), p2);
%         same2 = sameside(center, center+tangent(:,1), joinp2(:,1), p2);    
%         if (same1==1) && (same2==1)
%             dist(2) = norm(p2-q(:,2,1));
%         end
%     end
%     if isempty(p3)==0
%         same1 = sameside(q(:,1,1), q(:,2,1), center+tangent(:,1), p3);
%         same2 = sameside(center, center+tangent(:,1), joinp2(:,3), p3);
%         if (same1==1) && (same2==1)
%             dist(3) = norm(p3-q(:,1,1));
%         end
%     end
%     if isempty(p4)==0
%         same1 = sameside(q(:,1,1), q(:,2,1), center+tangent(:,1), p4);
%         same2 = sameside(center, center+tangent(:,1), joinp2(:,3), p4);
%         if (same1==1) && (same2==1)
%             dist(4) = norm(p4-q(:,2,1));
%         end
%     end
%     th = max(dist);    
%     
%     joinp = [joinp2 q(:,1,1)+tangent(:,1)*th q(:,2,1)+tangent(:,1)*th];
%     ind = ind([2 3 1]);
% end  
% 
% % reorder the nodes if necessary
% quad = joinp(:,[1 3 4 5]);    
% p1 = lineintersection(quad(:,1), quad(:,3)-quad(:,1),  quad(:,2), quad(:,4)-quad(:,2));
% quad = joinp(:,[1 3 5 4]);
% p2 = lineintersection(quad(:,1), quad(:,3)-quad(:,1),  quad(:,2), quad(:,4)-quad(:,2));
% c = mean(quad,2);
% if norm(p1-c)>norm(p2-c)
%     joinp = joinp(:,[1 2 3 5 4]);
% end
% 
% % determine the three edges of the joint
% edge = zeros(nd,2,3);
% edge(:,:,1) = joinp(:,1:2);
% edge(:,:,2) = joinp(:,2:3);
% 
% dist4 = zeros(3,1);
% dist5 = zeros(3,1);
% for i = 1:3
%     dist4(i) = norm(joinp(:,4)-joinp(:,i));
%     dist5(i) = norm(joinp(:,5)-joinp(:,i));
% end
% if min(dist4)<1e-6
%     [~,j] = min(dist4);
%     edge(:,:,3) = joinp(:,[j 5]);
%     joinp(:,4) = [];
%     joint = [1 2 3 4];
% elseif min(dist5)<1e-6
%     [~,j] = min(dist5);
%     edge(:,:,3) = joinp(:,[j 4]);
%     joinp(:,5) = [];
%     joint = [1 2 3 4];
% else
%     edge(:,:,3) = joinp(:,4:5);    
%     % insert an additional node to form two quadrilaterals
%     p1 = 0.5*(joinp(:,1)+joinp(:,5));
%     p2 = 0.5*(joinp(:,3)+joinp(:,4));
%     if norm(p1-joinp(:,3))<norm(p2-joinp(:,1))
%         joinp = [joinp p1];    
%         joint = [1 2 3 6; 3 4 5 6];
%     else
%         joinp = [joinp(:,1:3) p2 joinp(:,4:5)];    
%         joint = [1 2 3 4; 4 5 6 1];
%     end        
% end
% 
% [~,ind] = sort(ind);
% edge = edge(:,:,ind);


