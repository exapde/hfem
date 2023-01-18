function [e1, e2] = intersectmesh(p1,t1,p2,t2,tol)
% intersect (p1,t1) and (p2,t2) 

if nargin<5
    tol = 1e-10;
end

nt1 = size(t1,1);
nt2 = size(t2,1);
e1 = [];
e2 = [];
for i = 1:nt1
    pi = p1(t1(i,:),:);
    for j = 1:nt2
        pj = p2(t2(j,:),:);
        in = xiny(pi,pj, tol);
        if min(in)>0
            e1 = [e1 i];
            e2 = [e2 j];
        end
    end
end

function in = xiny(x,y,tol)
% Determine if each row of x is a member of y
% If row j of x is a member of y and x(j,:) = y(k,:) then in(j) = k
% Else in(j) = 0

[m,dim] = size(x);

in = zeros(m,1);
if dim==2
    for j=1:m
        d2 = (y(:,1)-x(j,1)).^2 + (y(:,2)-x(j,2)).^2;
        [md,id] = min(d2);
        if md<tol, in(j)=id; end
    end
else
    for j=1:m
        d2 = (y(:,1)-x(j,1)).^2 + (y(:,2)-x(j,2)).^2 + (y(:,3)-x(j,3)).^2;
        [md,id] = min(d2);
        if md<tol, in(j)=id; end
    end
end


