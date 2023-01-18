function [p3,t3] = QtoT( p2, q, thickness, nlayers)

dz = thickness/nlayers;

lnd = true(size(p2,1),1);
lnd(q) = 0;
np = size(p2,1)-sum(lnd);

p2(lnd,1) = Inf;
[Y,I] = sort(p2(:,1));
J(I) = 1:size(p2,1);
q = J(q);

p2 = p2(I(1:np),:);
t2 = reshape([q(:,1:3)'; q(:,1)'; q(:,3:4)'],3,2*size(q,1))';

% Re-order element nodes

ne = size(t2,1);
[Y,I] = min(t2,[],2);
t = [t2,t2];
for i = 1:ne
    t2(i,:) = [t(i,I(i)),t(i,I(i)+1),t(i,I(i)+2)];
end

% Create 3D mesh

p3 = zeros(np*(nlayers+1),3);
zc = ones(np,1)*(0:dz:thickness);
p3 = [repmat(p2,nlayers+1,1), zc(:)];

t3 = zeros(3*ne*nlayers,4);

% First layer
for i = 1:ne
    in = t2(i,:);
    jn = in + np;
    t3(3*(i-1)+1,:) = [jn(1), jn(3), jn(2), in(1)];
    if in(2) < in(3)
        t3(3*(i-1)+2,:) = [in(1), in(2), in(3), jn(3)];
        t3(3*(i-1)+3,:) = [in(1), in(2), jn(3), jn(2)];
    else
        t3(3*(i-1)+2,:) = [in(1), in(2), in(3), jn(2)];
        t3(3*(i-1)+3,:) = [in(3), in(1), jn(2), jn(3)];
    end
end

for i = 2:nlayers
    t3(3*ne*(i-1)+1:3*ne*i,:) = t3(3*ne*(i-2)+1:3*ne*(i-1),:) + np;
end

