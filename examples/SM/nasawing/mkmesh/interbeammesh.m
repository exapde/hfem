function [f1, f2, l1, l2] = interbeammesh(p1,t1,p2,t2,tol)

nd = size(p1,2);
npf = size(t1,2);
porder = 1;
if npf==4
    elemtype=1;
    plocfc = [0  0; 1  0; 1  1; 0  1];
else
    elemtype=0;
    plocfc = [0  0; 1  0; 0  1];    
end
gpfc = plocfc;
ngf = size(gpfc,1);
shapfc = mkshape(porder,plocfc,gpfc,elemtype);
dshapft  = reshape(permute(shapfc(:,:,2:nd),[2 3 1]),[ngf*(nd-1) npf]);


nl = length(l1);
for i = 1:nl
    ni = length(l1{i});
    for j = 1:ni
        k1 = l1{i}(j);
        e1 = f1(k1,:);  
        % [a1, b1] thickness
        a1 = thickness1(e1(3)); 
        m1 = normalvector(dshapft, p1(t1(e1(3),:),:));
        if e1(4)>0
            b1 = thickness1(e1(4));
            n1 = normalvector(dshapft, p1(t1(e1(4),:),:));
        end                        
                
        k2 = l2{i}(j);
        e2 = f2(k2,:);  
        a2 = thickness2(e2(3)); 
        m2 = normalvector(dshapft, p2(t2(e2(3),:),:));
        if e2(4)>0
            b2 = thickness2(e2(4));
            n2 = normalvector(dshapft, p2(t2(e2(4),:),:));
        end                
        
    end
end

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
l1 = interlines(f1(:,1:2));
l2 = interlines(f2(:,1:2));

nd = size(p1,2);
figure(1); clf;
patch('vertices',p1,'faces',t1,'facecol','g','edgec','k');
hold on;
patch('vertices',p2,'faces',t2,'facecol','g','edgec','k');
hold on;
if nd==2
    for i = 1:size(f1,1)
        pi = p1(f1(i,1:2),:);
        plot(pi(:,1),pi(:,2),'-r');
    end
else
    for i = 1:size(f1,1)
        pi = p1(f1(i,1:2),:);
        plot3(pi(:,1),pi(:,2),pi(:,3),'-r','LineWidth',1.5);
    end
end
axis equal
axis tight
colormap(jet);
colorbar
view([75 27])


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


function listlines = interlines(f)

nf = size(f,1);
listnb = cell(nf,1);
roots = zeros(nf,1);
for i=1:nf
    fi = f(i,1:2);
    in = find(f(:,1)==fi(1));
    in = unique([in; find(f(:,2)==fi(1))]);
    in = unique([in; find(f(:,1)==fi(2))]);
    in = unique([in; find(f(:,2)==fi(2))]);
    in = setdiff(in,i); 
    listnb{i} = in;
    roots(i) = length(in);
end

in = find(roots==1);
j = 1;
for i = 1:length(in)
    k = in(i);
    root = 1;
    for m = 1:j-1
        if ismember(k,listlines{m})
            root = 0;
            break;
        end
    end
    if root==1 
        listlines{j} = k;        
        while (1)
            newnb = setdiff(listnb{k},listlines{j});        
            if isempty(newnb)==1
                break;
            else
                listlines{j} = [listlines{j}; newnb];
                k = newnb;
            end
        end
        j = j + 1;
    end
end

function nlg = normalvector(dshapft,pn)

nd = size(pn,2);
ngf = size(dshapft,1)/(nd-1);
dpg = dshapft*pn;
dpg = reshape(dpg,[ngf,nd-1,nd]);    

nlg = zeros(ngf,nd);
nlg(:,1) = dpg(:,1,2).*dpg(:,2,3) - dpg(:,1,3).*dpg(:,2,2);
nlg(:,2) = dpg(:,1,3).*dpg(:,2,1) - dpg(:,1,1).*dpg(:,2,3);
nlg(:,3) = dpg(:,1,1).*dpg(:,2,2) - dpg(:,1,2).*dpg(:,2,1);
jac = sqrt(nlg(:,1).^2+nlg(:,2).^2+nlg(:,3).^2);
nlg   = bsxfun(@rdivide, nlg, jac);

return;

