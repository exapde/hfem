function [pb, tb, eb, pp] = mkbeamplate(p, ts, pp, th, ind, cb, tol)

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
end

beamline = p([f1(:,1); f1(end,2)],:);
thn = th(ind);
[pb, tb, eb] = mkbeam(beamline, normal, thn, shellcenter, cb);

for i = 1:n
    pp{ind(i)} = fixp(pp{ind(i)}, eb(:,:,i));
end   
