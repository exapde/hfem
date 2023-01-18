function [p, t, e] = mkbeam(beamline, normal, thickness, shellcenter, cs)

[np, nd, ns] = size(normal);
joinp = zeros(4, nd, np);
joint = zeros(np, 4);
edge = zeros(2, nd, ns, np);
if ns==2
    for i=1:np
        nmv = reshape(normal(i,:,:),[nd ns])';
        [joinp(:,:,i), joint(i,:), edge(:,:,:,i)] = case2(beamline(i,:), nmv, thickness, shellcenter, cs);
        joint(i,:) = joint(i,:) + (i-1)*4;
    end
elseif ns==3
    for i=1:np
        nmv = reshape(normal(i,:,:),[nd ns])';
        [joinp(:,:,i), joint(i,:), edge(:,:,:,i)] = case3(beamline(i,:), nmv, thickness, shellcenter, cs);
        joint(i,:) = joint(i,:) + (i-1)*4;
    end    
end

p = reshape(permute(joinp, [1 3 2]),[4*np nd]);
t = zeros(np-1,8);
for i = 1:np-1
    t(i,:) = [joint(i,:) joint(i+1,:)];
end
e = reshape(permute(edge, [1 4 2 3]),[2*np nd ns]);

function d = distancetoplane(p, v, x)

e = -dot(v,p);
n = size(x,1);
d = zeros(n,1);
for i = 1:n
    d(i) = sum(v.*x(i,:))+e;
end

function [joinp, joint, edge] = case2(beampoint, normal, thickness, shellcenter, cs)

nd = length(beampoint);
edge = zeros(2,nd,2);
if cs == 0
    p1 =  beampoint + 0.5*thickness(1)*normal(1,:);
    p2 =  beampoint - 0.5*thickness(1)*normal(1,:);
    p3 =  p2 + thickness(2)*normal(2,:);
    p4 =  p1 + thickness(2)*normal(2,:);    
    edge(:,:,1) = [p1; p2];
    if norm(shellcenter(2,:)-0.5*(p1+p4))<norm(shellcenter(2,:)-0.5*(p2+p3))
        edge(:,:,1) = [p1; p4];
    else
        edge(:,:,1) = [p2; p3];
    end                
%     cp = (p1+p2+p3+p4)/4;    
%     d = distancetoplane(p1, normal(1,:), [cp; shellcenter(2,:)]);
%     if d(1)*d(2)<0
%         edge(:,:,2) = [p1; p4];
%     else
%         edge(:,:,2) = [p2; p3];
%     end    
elseif cs == 1
    p1 =  beampoint + 0.5*thickness(1)*normal(1,:);
    p2 =  beampoint - 0.5*thickness(1)*normal(1,:);
    p3 =  p2 - thickness(2)*normal(2,:);
    p4 =  p1 - thickness(2)*normal(2,:);    
    edge(:,:,1) = [p1; p2];
    if norm(shellcenter(2,:)-0.5*(p1+p4))<norm(shellcenter(2,:)-0.5*(p2+p3))
        edge(:,:,1) = [p1; p4];
    else
        edge(:,:,1) = [p2; p3];
    end            
%     cp = (p1+p2+p3+p4)/4;    
%     d = distancetoplane(p1, normal(1,:), [cp; shellcenter(2,:)]);
%     if d(1)*d(2)<0
%         edge(:,:,2) = [p1; p4];
%     else
%         edge(:,:,2) = [p2; p3];
%     end    
elseif cs == 2
    p1 =  beampoint + 0.5*thickness(2)*normal(2,:);
    p2 =  beampoint - 0.5*thickness(2)*normal(2,:);
    p3 =  p2 + thickness(1)*normal(1,:);
    p4 =  p1 + thickness(1)*normal(1,:);    
    edge(:,:,2) = [p1; p2];
    if norm(shellcenter(1,:)-0.5*(p1+p4))<norm(shellcenter(1,:)-0.5*(p2+p3))
        edge(:,:,1) = [p1; p4];
    else
        edge(:,:,1) = [p2; p3];
    end            
%     cp = (p1+p2+p3+p4)/4;    
%     d = distancetoplane(p1, normal(2,:), [cp; shellcenter(1,:)]);
%     if d(1)*d(2)<0
%         edge(:,:,1) = [p1; p4];
%     else
%         edge(:,:,1) = [p2; p3];
%     end    
elseif cs == 3    
    p1 =  beampoint + 0.5*thickness(2)*normal(2,:);
    p2 =  beampoint - 0.5*thickness(2)*normal(2,:);
    p3 =  p2 - thickness(1)*normal(1,:);
    p4 =  p1 - thickness(1)*normal(1,:);    
    edge(:,:,2) = [p1; p2];
    if norm(shellcenter(1,:)-0.5*(p1+p4))<norm(shellcenter(1,:)-0.5*(p2+p3))
        edge(:,:,1) = [p1; p4];
    else
        edge(:,:,1) = [p2; p3];
    end        
%     cp = (p1+p2+p3+p4)/4;    
%     d = distancetoplane(p1, normal(2,:), [cp; shellcenter(1,:)]);    
%     if d(1)*d(2)<0
%         edge(:,:,1) = [p1; p4];
%     else
%         edge(:,:,1) = [p2; p3];
%     end    
end    
joinp = [p1; p2; p3; p4];
joint = [1 2 3 4];


function [joinp, joint, edge] = case3(beampoint, normal, thickness, shellcenter, cs)

nd = length(beampoint);
edge = zeros(2,nd,3);

p1 = beampoint - 0.5*thickness(3)*normal(3,:);
p2 = beampoint + 0.5*thickness(3)*normal(3,:);
edge(:,:,3) = [p1; p2];  

if cs==0
    p3 = p2 + thickness(1)*normal(1,:);
    p4 = p1 + thickness(1)*normal(1,:);    
    %cp = (p1+p2+p3+p4)/4;  
    %d = distancetoplane(p1, normal(3,:), [cp; shellcenter(1,:)]);
    %if d(1)*d(2)<0    
    if norm(shellcenter(1,:)-0.5*(p1+p4))<norm(shellcenter(1,:)-0.5*(p2+p3))
        p3 = p2 + thickness(2)*normal(2,:);
        edge(:,:,1) = [p1; p4];
        edge(:,:,2) = [p2; p3];
    else
        p4 = p1 + thickness(2)*normal(2,:);    
        edge(:,:,1) = [p2; p3];
        edge(:,:,2) = [p1; p4];
    end    
elseif cs==1
    p3 = p2 - thickness(1)*normal(1,:);
    p4 = p1 - thickness(1)*normal(1,:);    
    %cp = (p1+p2+p3+p4)/4;  
    %d = distancetoplane(p1, normal(3,:), [cp; shellcenter(1,:)]);
    %if d(1)*d(2)<0    
    if norm(shellcenter(1,:)-0.5*(p1+p4))<norm(shellcenter(1,:)-0.5*(p2+p3))
        p3 = p2 - thickness(2)*normal(2,:);
        edge(:,:,1) = [p1; p4];
        edge(:,:,2) = [p2; p3];
    else
        p4 = p1 - thickness(2)*normal(2,:);    
        edge(:,:,1) = [p2; p3];
        edge(:,:,2) = [p1; p4];
    end        
end

joinp = [p1; p2; p3; p4];
joint = [1 2 3 4];




