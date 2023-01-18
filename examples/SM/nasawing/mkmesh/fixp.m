function pnew = fixp(pold, e)

n = size(e,1);
pnew = pold;
for i = 1:n
    d = (pold(:,1)-e(i,1)).^2 + (pold(:,2)-e(i,2)).^2 + (pold(:,3)-e(i,3)).^2; 
    [~,j] = min(d);
    pnew(j,:) = e(i,:);
end
