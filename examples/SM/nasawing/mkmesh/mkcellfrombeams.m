function [pc,tc,cp] = mkcellfrombeams(pb,tb,ind)

nd = 3;
n = length(ind);
p = zeros(2,nd,n);
for i = 1:n
    j = ind(i);
    ti = tb{j}(1,1:4);
    p(1,:,i) = mean(pb{j}(ti,:),1);
    ti = tb{j}(end,5:8);
    p(2,:,i) = mean(pb{j}(ti,:),1);
end

ct = zeros(4,n);
cp = zeros(4,3,n);
for i = 1:n
    pi = p(:,:,i);
    j = 1;
    if i == 1
        j = n;    
    end
    pj = p(:,:,j);
    d(1) = norm(pi(1,:)-pj(1,:));
    d(2) = norm(pi(1,:)-pj(2,:));
    d(3) = norm(pi(2,:)-pj(1,:));
    d(4) = norm(pi(2,:)-pj(2,:));
    [~,k] = min(d);
    if k<=2
        ct(:,i) = tb{ind(i)}(1,1:4);
    else
        ct(:,i) = tb{ind(i)}(end,5:8);
    end    
    cp(:,:,i) = pb{ind(i)}(ct(:,i),:);
end
tc = [1 2 3 4 5 6 7 8];
 
p1 = [];
p2 = [];
for i = 2:n
    in = xiny(cp(:,:,i), cp(:,:,1));  
    kk = find(in>0);
    nk = length(kk);
    if nk==0
        d = zeros(4,4);
        for m = 1:4
            for k = 1:4
                d(k,m) = norm(cp(k,:,1)-cp(m,:,i));
            end
        end
        [a,b] = min(d);
        [~,m] = min(a);
        k = b(m);
        in1 = [k:1:4 1:1:(k-1)];
        in2 = [m:1:4 1:1:(m-1)];
        pc = [cp(in1,:,1); cp(in2,:,i)];        
        figure(1);clf;hold on;
        for i = 1:length(ind)
            boundaryplot(pb{ind(i)},tb{ind(i)},{'facecolor','r','edgecolor','k','Linew',1}); 
        end
        boundaryplot(pc,tc,{'facecolor','g','edgecolor','k','Linew',1}); 
        % tm=reshape(permute(cp,[1 3 2]),[12 3]);
        % plot3(tm(:,1),tm(:,2),tm(:,3),'bo','MarkerSize',10,'LineWidth',1.5);
        axis equal
        axis tight
        colormap(jet);
        colorbar
        view([75 27]);                 
        return;
    elseif nk==2
         if (kk(1)-1==0) && (kk(2)-2==0)
            mm = [4; 3];
        elseif (kk(1)-2==0) && (kk(2)-3==0)
            mm = [1; 4];
        elseif (kk(1)-3==0) && (kk(2)-4==0)
            mm = [2; 1];
        elseif (kk(1)-1==0) && (kk(2)-4==0)
            mm = [2; 3];
        end    
        p1 = [p1; cp(kk,:,i)];
        p2 = [p2; cp(mm,:,i)];
    else
        error('Something Wrong');
    end
end

[p1,ia] = unique(p1,'rows');
p2 = p2(ia,:);

in = xiny(cp(:,:,1), p1);  
kk = find(in>0);
nk = length(kk);
if nk==2
    im = in(kk);
    if (kk(1)-1==0) && (kk(2)-2==0)
        pa = p2(im(1),:) + cp(3,:,1) - cp(2,:,1);
        pb = p2(im(2),:) + cp(4,:,1) - cp(1,:,1);
        pc = [cp(:,:,1); p2(im(1),:); p2(im(2),:); pa; pb];
    elseif (kk(1)-2==0) && (kk(2)-3==0)
        pa = p2(im(1),:) + cp(1,:,1) - cp(2,:,1);
        pb = p2(im(2),:) + cp(4,:,1) - cp(3,:,1);
        pc = [cp(:,:,1); pa; p2(im(1),:); p2(im(2),:); pb];    
    elseif (kk(1)-3==0) && (kk(2)-4==0)
        pa = p2(im(1),:) + cp(1,:,1) - cp(4,:,1);
        pb = p2(im(2),:) + cp(2,:,1) - cp(3,:,1);
        pc = [cp(:,:,1); pa; pb; p2(im(1),:); p2(im(2),:)];        
    elseif (kk(1)-1==0) && (kk(2)-4==0)
        pa = p2(im(1),:) + cp(2,:,1) - cp(1,:,1);
        pb = p2(im(2),:) + cp(3,:,1) - cp(4,:,1);
        pc = [cp(:,:,1); p2(im(1),:); pa; pb; p2(im(2),:)];            
    end    
elseif nk==3
    j = find(in==0);
    im = in(kk);
    if j==1 
        p = p2(im(1),:) + p2(im(3),:)-p2(im(2),:);
        pc = [cp(:,:,1); p; p2(im,:)];
    elseif j==2
        p = p2(im(1),:) + p2(im(2),:)-p2(im(3),:);
        pc = [cp(:,:,1); p2(im(1),:); p; p2(im(2),:); p2(im(3),:)];
    elseif j==3
        p = p2(im(2),:) + p2(im(3),:)-p2(im(1),:);
        pc = [cp(:,:,1); p2(im(1),:); p2(im(2),:); p; p2(im(3),:)];
    elseif j==4
        p = p2(im(1),:) + p2(im(3),:)-p2(im(2),:);
        pc = [cp(:,:,1); p2(im(1),:); p2(im(2),:); p2(im(3),:); p];
    end
elseif nk==4 
    pc = [cp(:,:,1); p2(in,:)];
end

% figure(1);clf;hold on;
% for i = 1:length(ind)
%     boundaryplot(pb{ind(i)},tb{ind(i)},{'facecolor','r','edgecolor','k','Linew',1}); 
% end
% boundaryplot(pc,tc,{'facecolor','g','edgecolor','k','Linew',1}); 
% % tm=reshape(permute(cp,[1 3 2]),[12 3]);
% % plot3(tm(:,1),tm(:,2),tm(:,3),'bo','MarkerSize',10,'LineWidth',1.5);
% axis equal
% axis tight
% colormap(jet);
% colorbar
% view([75 27]);        
% 


% pc = zeros(8,3);
% pc(1:4,:) = cp(:,:,1);
% tc = [1 2 3 4 5 6 7 8];
% 
% 
% cp = reshape(permute(cp,[1 3 2]),[4*np nd]);
% in = xiny(cp, pc(1:4,:));   
% po = cp(in==0,:);
% m = size(po,1);
% if m>4 || m<=1
%     error('something wrong');
% end
% 
% if m==2
%     
% end
% 
% 
% 
