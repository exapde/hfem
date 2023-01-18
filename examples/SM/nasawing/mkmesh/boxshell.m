p = round([X Y Z]/1e-8)*1e-8;
comp = 2;  t2 = CQUAD(ID==comp,:); th2 = thickness(ID==comp);
comp = 3;  t3 = CQUAD(ID==comp,:); th3 = thickness(ID==comp);

[shell2, ns2, es2] = breakshells(p, t2, th2);
[shell3, ns3, es3] = breakshells(p, t3, th3);

ts2 = shell2{3};
ts3{1} = shell3{1}{1};
for i=2:1:(length(shell3{1})-1)
    ts3{i}=shell3{1}{i+1};
end
ind = [1 2 3 6 5 4 7 8];
for i=1:8
    ts{i} = t2(ts2{ind(i)},:);
    th{i} = th2(ts2{ind(i)});
end
ind = [1 2 3 4 5];
for i=1:5
    ts{8+i} = t3(ts3{ind(i)},:);
    th{8+i} = th3(ts3{ind(i)});
end

npf = 4; nd = 3;
figure(1); clf;
hold on;
for m = 1:length(ts)      
    ne = size(ts{m},1);    
    pn = zeros(npf,ne,nd);
    for k = 1:nd
        pn(:,:,k) = reshape(p(ts{m},k),[ne npf])';
    end        
    patch(pn(:,:,1),pn(:,:,2),pn(:,:,3),th{m}','edgecolor','k');
    view([75 27])        
end
axis equal
axis tight
colormap(jet);
colorbar
view([75 27])        

[f1, f2, n1, n2, c1, c2] = intersectshells(p,ts{1},p,ts{9},1e-6);
beamline = p([f1(:,1); f1(end,2)],:);
normal = n1;
normal(:,:,2) = n2;
thn = 10*[th{1}(1) th{9}(1)];
shellcenter = [c1; c2];

[pb{1}, tb{1}, eb{1}] = mkbeam(beamline, normal, thn, shellcenter, 2);
[pp{1}, tp{1}] = mkplate(p, ts{1}, thn(1), 1);
[pp{9}, tp{9}] = mkplate(p, ts{9}, thn(2), 0);

pp{1} = fixp(pp{1}, eb{1}(:,:,1));
pp{9} = fixp(pp{9}, eb{1}(:,:,2));

figure(2); clf;
pars={'facecolor','r','edgecolor','k','Linew',1};
boundaryplot(pb{1},tb{1},pars); 
hold on;
pars={'facecolor','g','edgecolor','k','Linew',1};
boundaryplot(pp{1},tp{1},pars); 
pars={'facecolor','b','edgecolor','k','Linew',1};
boundaryplot(pp{9},tp{9},pars); 
view(3);


[f1, f10, n1, n10, c1, c10] = intersectshells(p,ts{1},p,ts{10},1e-6);
[f2, ~, n2, ~, c2, ~] = intersectshells(p,ts{2},p,ts{10},1e-6);
beamline = p([f1(:,1); f1(end,2)],:);
normal = n1;
normal(:,:,2) = n2;
normal(:,:,3) = n10;
thn = 10*[th{1}(1) th{2}(1) th{10}(1)];
shellcenter = [c1; c2; c10];

[pb{2}, tb{2}, eb{2}] = mkbeam(beamline, normal, thn, shellcenter, 0);

%[pp{1}, tp{1}] = mkplate(p, ts{1}, thn(1), 1);
[pp{2}, tp{2}] = mkplate(p, ts{2}, thn(2), 1);
[pp{10}, tp{10}] = mkplate(p, ts{10}, thn(3), 0);

pp{1} = fixp(pp{1}, eb{2}(:,:,1));
pp{2} = fixp(pp{2}, eb{2}(:,:,2));
pp{10} = fixp(pp{10}, eb{2}(:,:,3));

figure(3); clf;
pars={'facecolor','r','edgecolor','k','Linew',1};
boundaryplot(pb{1},tb{1},pars); 
hold on;
pars={'facecolor','r','edgecolor','k','Linew',1};
boundaryplot(pb{2},tb{2},pars); 
pars={'facecolor','g','edgecolor','k','Linew',1};
boundaryplot(pp{1},tp{1},pars); 
pars={'facecolor','g','edgecolor','k','Linew',1};
boundaryplot(pp{2},tp{2},pars); 
pars={'facecolor','g','edgecolor','k','Linew',1};
boundaryplot(pp{9},tp{9},pars); 
pars={'facecolor','g','edgecolor','k','Linew',1};
boundaryplot(pp{10},tp{10},pars); 
view(3);


[f2, f11, n2, n11, c2, c11] = intersectshells(p,ts{2},p,ts{11},1e-6);
[f3, ~, n3, ~, c3, ~] = intersectshells(p,ts{3},p,ts{11},1e-6);
beamline = p([f2(:,1); f2(end,2)],:);
normal = n2;
normal(:,:,2) = n3;
normal(:,:,3) = n11;
thn = 10*[th{2}(1) th{3}(1) th{11}(1)];
shellcenter = [c2; c3; c11];

[pb{3}, tb{3}, eb{3}] = mkbeam(beamline, normal, thn, shellcenter, 0);

[pp{3}, tp{3}] = mkplate(p, ts{3}, thn(2), 1);
[pp{11}, tp{11}] = mkplate(p, ts{11}, thn(3), 0);

pp{2} = fixp(pp{2}, eb{3}(:,:,1));
pp{3} = fixp(pp{3}, eb{3}(:,:,2));
pp{11} = fixp(pp{11}, eb{3}(:,:,3));

figure(4); clf;
pars={'facecolor','r','edgecolor','k','Linew',1};
boundaryplot(pb{1},tb{1},pars); 
hold on;
pars={'facecolor','r','edgecolor','k','Linew',1};
boundaryplot(pb{2},tb{2},pars); 
pars={'facecolor','r','edgecolor','k','Linew',1};
boundaryplot(pb{3},tb{3},pars); 
pars={'facecolor','g','edgecolor','k','Linew',1};
boundaryplot(pp{1},tp{1},pars); 
pars={'facecolor','g','edgecolor','k','Linew',1};
boundaryplot(pp{2},tp{2},pars); 
pars={'facecolor','g','edgecolor','k','Linew',1};
boundaryplot(pp{3},tp{3},pars); 
pars={'facecolor','g','edgecolor','k','Linew',1};
boundaryplot(pp{9},tp{9},pars); 
pars={'facecolor','g','edgecolor','k','Linew',1};
boundaryplot(pp{10},tp{10},pars); 
pars={'facecolor','g','edgecolor','k','Linew',1};
boundaryplot(pp{11},tp{11},pars); 
view(3);


[f3, f12, n3, n12, c3, c12] = intersectshells(p,ts{3},p,ts{12},1e-6);
[f4, ~, n4, ~, c4, ~] = intersectshells(p,ts{4},p,ts{12},1e-6);
beamline = p([f3(:,1); f3(end,2)],:);
normal = n3;
normal(:,:,2) = n4;
normal(:,:,3) = n12;
thn = 10*[th{3}(1) th{4}(1) th{12}(1)];
shellcenter = [c3; c4; c12];

[pb{4}, tb{4}, eb{4}] = mkbeam(beamline, normal, thn, shellcenter, 0);

[pp{4}, tp{4}] = mkplate(p, ts{4}, thn(2), 1);
[pp{12}, tp{12}] = mkplate(p, ts{12}, thn(3), 0);

pp{3} = fixp(pp{3}, eb{4}(:,:,1));
pp{4} = fixp(pp{4}, eb{4}(:,:,2));
pp{12} = fixp(pp{12}, eb{4}(:,:,3));

figure(5); clf;
hold on;
pars={'facecolor','r','edgecolor','k','Linew',1};
ind = [1 2 3 4];
for i = 1:length(ind)
    boundaryplot(pb{ind(i)},tb{ind(i)},pars); 
end
pars={'facecolor','g','edgecolor','k','Linew',1};
ind = [1 2 3 4 9 10 11 12];
for i = 1:length(ind)
    boundaryplot(pp{ind(i)},tp{ind(i)},pars); 
end
view(3);

[f4, f13, n4, n13, c4, c13] = intersectshells(p,ts{4},p,ts{13},1e-6);
beamline = p([f4(:,1); f4(end,2)],:);
normal = n4;
normal(:,:,2) = n13;
thn = 10*[th{4}(1) th{13}(1)];
shellcenter = [c4; c13];

[pb{5}, tb{5}, eb{5}] = mkbeam(beamline, normal, thn, shellcenter, 2);

[pp{13}, tp{13}] = mkplate(p, ts{13}, thn(2), 0);

pp{4} = fixp(pp{4}, eb{5}(:,:,1));
pp{13} = fixp(pp{13}, eb{5}(:,:,2));

figure(6); clf;
hold on;
pars={'facecolor','r','edgecolor','k','Linew',1};
ind = [1 2 3 4 5];
for i = 1:length(ind)
    boundaryplot(pb{ind(i)},tb{ind(i)},pars); 
end
pars={'facecolor','g','edgecolor','k','Linew',1};
ind = [1 2 3 4 9 10 11 12 13];
for i = 1:length(ind)
    boundaryplot(pp{ind(i)},tp{ind(i)},pars); 
end
view(3);

[f5, f13, n5, n13, c5, c13] = intersectshells(p,ts{5},p,ts{13},1e-6);
beamline = p([f5(:,1); f5(end,2)],:);
normal = n5;
normal(:,:,2) = n13;
thn = 10*[th{5}(1) th{13}(1)];
shellcenter = [c5; c13];

[pb{6}, tb{6}, eb{6}] = mkbeam(beamline, normal, thn, shellcenter, 3);

[pp{5}, tp{5}] = mkplate(p, ts{5}, thn(1), 2);

pp{5} = fixp(pp{5}, eb{6}(:,:,1));
pp{13} = fixp(pp{13}, eb{6}(:,:,2));

figure(7); clf;
hold on;
pars={'facecolor','r','edgecolor','k','Linew',1};
ind = [1 2 3 4 5 6];
for i = 1:length(ind)
    boundaryplot(pb{ind(i)},tb{ind(i)},pars); 
end
pars={'facecolor','g','edgecolor','k','Linew',1};
ind = [1 2 3 4 9 10 11 12 13 5];
for i = 1:length(ind)
    boundaryplot(pp{ind(i)},tp{ind(i)},pars); 
end
view(3);


[f5, f12, n5, n12, c5, c12] = intersectshells(p,ts{5},p,ts{12},1e-6);
[f6, ~, n6, ~, c6, ~] = intersectshells(p,ts{6},p,ts{12},1e-6);
beamline = p([f5(:,1); f5(end,2)],:);
normal = n5;
normal(:,:,2) = n6;
normal(:,:,3) = n12;
thn = 10*[th{5}(1) th{6}(1) th{12}(1)];
shellcenter = [c5; c6; c12];

[pb{7}, tb{7}, eb{7}] = mkbeam(beamline, normal, thn, shellcenter, 1);
[pp{6}, tp{6}] = mkplate(p, ts{6}, thn(2), 2);
pp{5} = fixp(pp{5}, eb{7}(:,:,1));
pp{6} = fixp(pp{6}, eb{7}(:,:,2));
pp{12} = fixp(pp{12}, eb{7}(:,:,3));

figure(8); clf;
hold on;
pars={'facecolor','r','edgecolor','k','Linew',1};
ind = [1 2 3 4 5 6 7];
for i = 1:length(ind)
    boundaryplot(pb{ind(i)},tb{ind(i)},pars); 
end
pars={'facecolor','g','edgecolor','k','Linew',1};
ind = [1 2 3 4 9 10 11 12 13 5 6];
for i = 1:length(ind)
    boundaryplot(pp{ind(i)},tp{ind(i)},pars); 
end
view(3);

[f6, f11, n6, n11, c6, c11] = intersectshells(p,ts{6},p,ts{11},1e-6);
[f7, ~, n7, ~, c7, ~] = intersectshells(p,ts{7},p,ts{11},1e-6);
beamline = p([f6(:,1); f6(end,2)],:);
normal = n6;
normal(:,:,2) = n7;
normal(:,:,3) = n11;
thn = 10*[th{6}(1) th{7}(1) th{11}(1)];
shellcenter = [c6; c7; c11];

[pb{8}, tb{8}, eb{8}] = mkbeam(beamline, normal, thn, shellcenter, 1);
[pp{7}, tp{7}] = mkplate(p, ts{7}, thn(2), 2);
pp{6} = fixp(pp{6}, eb{8}(:,:,1));
pp{7} = fixp(pp{7}, eb{8}(:,:,2));
pp{11} = fixp(pp{11}, eb{8}(:,:,3));

figure(8); clf;
hold on;
pars={'facecolor','r','edgecolor','k','Linew',1};
ind = [1 2 3 4 5 6 7 8];
for i = 1:length(ind)
    boundaryplot(pb{ind(i)},tb{ind(i)},pars); 
end
pars={'facecolor','g','edgecolor','k','Linew',1};
ind = [1 2 3 4 9 10 11 12 13 5 6 7];
for i = 1:length(ind)
    boundaryplot(pp{ind(i)},tp{ind(i)},pars); 
end
view(3);


[f7, f10, n7, n10, c7, c10] = intersectshells(p,ts{7},p,ts{10},1e-6);
[f8, ~, n8, ~, c8, ~] = intersectshells(p,ts{8},p,ts{10},1e-6);
beamline = p([f7(:,1); f7(end,2)],:);
normal = n7;
normal(:,:,2) = n8;
normal(:,:,3) = n10;
thn = 10*[th{7}(1) th{8}(1) th{10}(1)];
shellcenter = [c7; c8; c10];

[pb{9}, tb{9}, eb{9}] = mkbeam(beamline, normal, thn, shellcenter, 1);
[pp{8}, tp{8}] = mkplate(p, ts{8}, thn(2), 2);
pp{7} = fixp(pp{7}, eb{9}(:,:,1));
pp{8} = fixp(pp{8}, eb{9}(:,:,2));
pp{10} = fixp(pp{10}, eb{9}(:,:,3));

figure(9); clf;
hold on;
pars={'facecolor','r','edgecolor','k','Linew',1};
ind = [1 2 3 4 5 6 7 8 9];
for i = 1:length(ind)
    boundaryplot(pb{ind(i)},tb{ind(i)},pars); 
end
pars={'facecolor','g','edgecolor','k','Linew',1};
ind = [1 2 3 4 9 10 11 12 13 5 6 7 8];
for i = 1:length(ind)
    boundaryplot(pp{ind(i)},tp{ind(i)},pars); 
end
view(3);

[f8, f9, n8, n9, c8, c9] = intersectshells(p,ts{8},p,ts{9},1e-6);
beamline = p([f8(:,1); f8(end,2)],:);
normal = n8;
normal(:,:,2) = n9;
thn = 10*[th{8}(1) th{9}(1)];
shellcenter = [c8; c9];

[pb{10}, tb{10}, eb{10}] = mkbeam(beamline, normal, thn, shellcenter, 3);
pp{8} = fixp(pp{8}, eb{10}(:,:,1));
pp{9} = fixp(pp{9}, eb{10}(:,:,2));

figure(10); clf;
hold on;
pars={'facecolor','r','edgecolor','k','Linew',1};
ind = [1 2 3 4 5 6 7 8 9 10];
for i = 1:length(ind)
    boundaryplot(pb{ind(i)},tb{ind(i)},pars); 
end
pars={'facecolor','g','edgecolor','k','Linew',1};
ind = [1 2 3 4 9 10 11 12 13 5 6 7 8];
for i = 1:length(ind)
    boundaryplot(pp{ind(i)},tp{ind(i)},pars); 
end
view(3);

% 
% % make the shell panels
% m = 11;
% n = 11;
% [p{1},t{1}]=squaremesh(m,n,0,1);
% p{1}(:,3) = 0;
% 
% [p{2},t{2}]=squaremesh(m,n,0,1);
% p{2}(:,3) = 1;
% 
% [p{3},t{3}]=squaremesh(m,n,0,1);
% p{3} = [0*p{3}(:,1) p{3}(:,2) p{3}(:,1)];
% 
% [p{4},t{4}]=squaremesh(m,n,0,1);
% p{4} = [ones(length(p{4}(:,1)),1) p{4}(:,2) p{4}(:,1)];
% 
% [p{5},t{5}]=squaremesh(m,n,0,1);
% p{5} = [p{5}(:,1) 0*p{5}(:,1) p{5}(:,2)];
% 
% [p{6},t{6}]=squaremesh(m,n,0,1);
% p{6} = [p{6}(:,1) ones(length(p{6}(:,1)),1) p{6}(:,2)];
% 
% figure(1); clf;
% patch('vertices',p{1},'faces',t{1},'facecol','r','edgec','k');
% hold on;
% patch('vertices',p{2},'faces',t{2},'facecol','g','edgec','k');
% patch('vertices',p{3},'faces',t{3},'facecol','b','edgec','k');
% patch('vertices',p{4},'faces',t{4},'facecol','y','edgec','k');
% patch('vertices',p{5},'faces',t{5},'facecol','m','edgec','k');
% patch('vertices',p{6},'faces',t{6},'facecol','c','edgec','k');
% view(3);
% 
% % thickness of each shell panel
% tmin = 0.01;
% tmax = 0.04;
% thickness = tmin + rand(1,length(p))*(tmax-tmin);
% 
% % connection information
% con = [3 5; 5 4; 4 6; 6 3];
% 
% for n = 1:size(con)
%     [f{1}{i}, f{i}{1}, l{1}{i}, l{i}{1}] = interedgemesh(p{1},t{1},p{i},t{i},tol);
% end
% 
% 
% 
% 
% 
% 
