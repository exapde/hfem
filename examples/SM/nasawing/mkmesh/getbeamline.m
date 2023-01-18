function [beamline,normal,tangent,shellcenter] = getbeamline(p, ts, ind, tol)

n = length(ind);
if n==2
    [f1, ~, n1, n2, c1, c2] = intersectshells(p,ts{ind(1)},p,ts{ind(2)},tol);    
    normal = n1;
    normal(:,:,2) = n2;    
    shellcenter = [c1; c2];
elseif n==3
    [f1, ~, n1, n3, c1, c3] = intersectshells(p,ts{ind(1)},p,ts{ind(3)},tol);
    [~, ~, n2, ~, c2, ~] = intersectshells(p,ts{ind(2)},p,ts{ind(3)},tol);    
    normal = n1;
    normal(:,:,2) = n2;
    normal(:,:,3) = n3;    
    shellcenter = [c1; c2; c3];    
elseif n==4
    [f1, ~, n1, n4, c1, c4] = intersectshells(p,ts{ind(1)},p,ts{ind(4)},tol);
    [~, ~, n2, n3, c2, c3] = intersectshells(p,ts{ind(2)},p,ts{ind(3)},tol);    
    normal = n1;
    normal(:,:,2) = n2;
    normal(:,:,3) = n3;    
    normal(:,:,4) = n4;    
    shellcenter = [c1; c2; c3; c4];    
end
beamline = p([f1(:,1); f1(end,2)],:);

m = size(beamline,1);
tangent = normal;
for i = 1:n
    for j = 1:m
        a1 = beamline(j,:) + [-normal(j,2,i) normal(j,1,i) normal(j,3,i)];
        a2 = beamline(j,:) + [normal(j,2,i) -normal(j,1,i) normal(j,3,i)];
        if norm(shellcenter(i,:)-a1)<norm(shellcenter(i,:)-a2)
            tangent(j,:,i) = [-normal(j,2,i) normal(j,1,i) normal(j,3,i)];
        else
            tangent(j,:,i) = [normal(j,2,i) -normal(j,1,i) normal(j,3,i)];
        end
    end    
end




