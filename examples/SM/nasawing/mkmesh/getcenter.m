function [c1,c2] = getcenter(p1,t1)

npe = size(t1,2);
nd = size(p1,2);
ne = size(t1,1);
pn = zeros(npe,ne,nd);
for i = 1:nd
    pn(:,:,i) = reshape(p1(t1,i),[ne npe])';
end
c2 = reshape(mean(pn,1),[ne nd]);
pn = reshape(pn,[npe*ne nd]);
c1 = mean(pn,1);


