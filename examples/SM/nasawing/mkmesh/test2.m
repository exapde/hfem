
m = 11;
n = 11;

[p{1},t{1}]=squaremesh(m,n,0,1);
p{1}(:,3) = 0;

[p{2},t{2}]=squaremesh(m,n,0,1);
p{2}(:,3) = 1;

[p{3},t{3}]=squaremesh(m,n,0,1);
p{3} = [0*p{3}(:,1) p{3}(:,2) p{3}(:,1)];

[p{4},t{4}]=squaremesh(m,n,0,1);
p{4} = [ones(length(p{4}(:,1)),1) p{4}(:,2) p{4}(:,1)];

[p{5},t{5}]=squaremesh(m,n,0,1);
p{5} = [p{5}(:,1) 0*p{5}(:,1) p{5}(:,2)];

[p{6},t{6}]=squaremesh(m,n,0,1);
p{6} = [p{6}(:,1) ones(length(p{6}(:,1)),1) p{6}(:,2)];

nd = size(p{1},2);
figure(1); clf;
patch('vertices',p{1},'faces',t{1},'facecol','r','edgec','k');
hold on;
patch('vertices',p{2},'faces',t{2},'facecol','g','edgec','k');
patch('vertices',p{3},'faces',t{3},'facecol','b','edgec','k');
patch('vertices',p{4},'faces',t{4},'facecol','y','edgec','k');
patch('vertices',p{5},'faces',t{5},'facecol','m','edgec','k');
patch('vertices',p{6},'faces',t{6},'facecol','c','edgec','k');
view(3);

ind = [6 5 3 4];
tol=1e-8;
for i = 3:6
    [f{1}{i}, f{i}{1}, l{1}{i}, l{i}{1}] = interedgemesh(p{1},t{1},p{i},t{i},tol);
    [f{2}{i}, f{i}{2}, l{2}{i}, l{i}{2}] = interedgemesh(p{2},t{2},p{i},t{i},tol);
    k = ind(i-2);
    [f{k}{i}, f{i}{k}, l{k}{i}, l{i}{k}] = interedgemesh(p{k},t{k},p{i},t{i},tol);
end




