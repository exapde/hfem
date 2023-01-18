function [f1, f2, n1, n2, c1, c2] = intersectshells(p1, t1, p2, t2, tol)

npe = size(t1,2);
if npe == 3
    plocfc = [0  0; 1  0; 0  1];    
    shapfc = mkshape(1,plocfc,plocfc,0);
    shapfc  = reshape(permute(shapfc(:,:,2:end),[2 3 1]),[npe*2 npe]);

elseif npe == 4
    plocfc = [0  0; 1  0; 1  1; 0  1];
    shapfc = mkshape(1,plocfc,plocfc,1);
    shapfc  = reshape(permute(shapfc(:,:,2:end),[2 3 1]),[npe*2 npe]);
end

[f1, f2, n1, n2] = interedge(p1,t1,p2,t2,shapfc,tol);

nd = size(p1,2);
ne = size(t1,1);
pn = zeros(npe,ne,nd);
for i = 1:nd
    pn(:,:,i) = reshape(p1(t1,i),[ne npe])';
end
pn = reshape(pn,[npe*ne nd]);
c1 = mean(pn,1);
ne = size(t2,1);
pn = zeros(npe,ne,nd);
for i = 1:nd
    pn(:,:,i) = reshape(p2(t2,i),[ne npe])';
end
pn = reshape(pn,[npe*ne nd]);
c2 = mean(pn,1);

isplot=0;
if isplot==1
    
figure(1); clf;
patch('vertices',p1,'faces',t1,'facecol','g','edgec','k');
hold on;
patch('vertices',p2,'faces',t2,'facecol','g','edgec','k');
for i = 1:size(f1,1)
    pi = p1(f1(i,1:2),:);
    plot3(pi(:,1),pi(:,2),pi(:,3),'-r','LineWidth',1.5);
    
    pa = pi;
    pb = pa - n1(i:i+1,:);
    pi = [pa(1,:); pb(1,:)];
    plot3(pi(:,1),pi(:,2),pi(:,3),'-b','LineWidth',1.5);
    pi = [pa(2,:); pb(2,:)];
    plot3(pi(:,1),pi(:,2),pi(:,3),'-b','LineWidth',1.5);
    
    pb = pa + n2(i:i+1,:);
    pi = [pa(1,:); pb(1,:)];
    plot3(pi(:,1),pi(:,2),pi(:,3),'-b','LineWidth',1.5);
    pi = [pa(2,:); pb(2,:)];
    plot3(pi(:,1),pi(:,2),pi(:,3),'-b','LineWidth',1.5);
end
axis equal
axis tight
colormap(jet);
colorbar
view([75 27])

end

% tol = 1e-8;
% n = length(ind);
% for i = 2:n
%     [f{1}, f{i}, n{1}, n{i}] = interedge(p{1},t{1},p{i},t{i},shapfq,tol);
% end


function nlg = normalvector(shapfc,pn)

[npe,nd] = size(pn);
dpg = shapfc*reshape(pn,[npe nd]);
dpg = permute(reshape(dpg,[npe nd-1 nd]), [1 3 2]);    

nlg = zeros(npe,nd);
nlg(:,1) = dpg(:,2,1).*dpg(:,3,2) - dpg(:,3,1).*dpg(:,2,2);
nlg(:,2) = dpg(:,3,1).*dpg(:,1,2) - dpg(:,1,1).*dpg(:,3,2);
nlg(:,3) = dpg(:,1,1).*dpg(:,2,2) - dpg(:,2,1).*dpg(:,1,2);
jac = sqrt(nlg(:,1).^2+nlg(:,2).^2+nlg(:,3).^2);
nlg   = bsxfun(@rdivide, nlg, jac);
nlg = reshape(nlg, [npe nd]);

function [k,l] = reorderf(f)

f = f(:,1:2);
g = unique(f(:));
count = histc(f(:), g);
root = g(count==1);

m = size(f,1);
n = length(g);
k = zeros(m,1);
l = zeros(n,1);
l(1) = root(1);
for i = 2:n
    for j = 1:m
        if (f(j,1)==l(i-1)) || (f(j,2)==l(i-1))
            a = setdiff(f(j,:),l(i-1));
            if ismember(a, l(1:i-1))==0
                l(i) = a;
                k(i-1) = j;
                break;
            end
        end
    end
end

if max(abs((unique(k)-(1:m)')))>0
    error('something wrong');
end

function [f1, f2, n1, n2] = interedge(p1,t1,p2,t2,shapfq,tol)

if nargin<5
    tol = 1e-10;
end

elemtype1 = size(t1,2)-3;
elemtype2 = size(t2,2)-3;
    
[f1,~,~] = mkt2f(t1,elemtype1);
[f2,~,~] = mkt2f(t2,elemtype2);

nf1 = size(f1,1);
nf2 = size(f2,1);
e1 = [];
e2 = [];
for i = 1:nf1    
    pi = p1(f1(i,1:2),:);
    for j = 1:nf2
        pj = p2(f2(j,1:2),:);        
        pk = pi-pj;% [pi(1,1)-pj(1,1) pi(2,1)-pj(2,1) pi(1,2)-pj(1,2) pi(2,2)-pj(2,2)];
        d1 = max(abs(pk(:)));
        pk = pi([2 1],:)-pj; %[pi(2,1)-pj(1,1) pi(1,1)-pj(2,1) pi(2,2)-pj(1,2) pi(1,2)-pj(2,2)];
        d2 = max(abs(pk(:)));        
        if min(d1,d2)<tol
            e1 = [e1 i];
            e2 = [e2 j];            
        end
    end
end
f1 = f1(e1,:);
f2 = f2(e2,:);

[k,l] = reorderf(f1);
f1 = f1(k,:);
f2 = f2(k,:);
if f1(1,1)~=l(1)
    f1 = f1(:,[2 1 3 4]);
end
if f2(1,1)~=l(1)
    f2 = f2(:,[2 1 3 4]);
end

nd = size(p1,2);
m = size(f1,1);
n1 = zeros(m+1,nd);
n2 = zeros(m+1,nd);
for i = 1:m
    ti = t1(f1(i,3),:);
    pi = p1(ti,:);
    ni = normalvector(shapfq,pi);
    na = ni(ti==f1(i,1),:);  
    nb = ni(ti==f1(i,2),:);  
    n1(i,:) = n1(i,:) + na;
    n1(i+1,:) = n1(i+1,:) + nb;            
    
    ti = t2(f2(i,3),:);
    pi = p2(ti,:);
    ni = normalvector(shapfq,pi);
    na = ni(ti==f2(i,1),:);  
    nb = ni(ti==f2(i,2),:);  
    n2(i,:) = n2(i,:) + na;
    n2(i+1,:) = n2(i+1,:) + nb;            
end
n1(2:m,:) = n1(2:m,:)/2;
n2(2:m,:) = n2(2:m,:)/2;





% 
% function nt = normaltangentvectors(pn, p1, p2, shapft, shapfq)
% 
% [npe,nd] = size(pn);
% if npe==3
%     nl = normalvector(shapft,pn);
% else
%     nl = normalvector(shapfq,pn);
% end
% 
% % tangent vectors
% tl = p2-p1;
% tl = tl/norm(tl);
% 
% nt = zeros(3,nd);
% pm = mean(pn,1);
% for i = 1:npe
%     if norm(pn(i,:)-p1)<1e-12
%         nt(1,:) = nl(i,:);
%         nt(2,:) = tl;
%         t1 = cross(nl(i,:),tl);
%         t1 = t1/norm(1);
%         t2 = -t1;        
%         if norm(p1+t1-pm)<(p1+t2-pm)            
%             nt(3,:) = t1;
%         else
%             nt(3,:) = t2;
%         end
%     end        
% end


