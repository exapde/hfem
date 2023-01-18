function mesh = mkmesh_L_imex(porder)

% xmin = -10;
% xmax =  11;
% ymin = -10;
% ymax =  10;

xmin = 0;
xmax =  1;
ymin = 0;
ymax =  1;
xmid = 0.5;
ymid = 0.5;

% rectangle
% pv1 = [xmin ymin; xmax ymin; xmax ymax; xmid ymax; xmid ymid; xmin ymid];
%pv1 = [xmin ymin; xmax ymin; xmax ymax; xmin ymax];
m  = 5;
a  = 2;
nodes1 = logdec(linspace(0,0.5,m),a)';
nodes2 = loginc(linspace(0.5,1,m),a)';
nodes3 = linspace(0,0.5,m)';
nodes4 = linspace(0.5,1,m)';

nodes_skewed  = [nodes1; nodes2(2:end)];
nodes_even = [nodes3; nodes4(2:end)];

mid = zeros(m,1);

for i = 1:m
    mid(i) = 0.5;
end


n  = length(nodes_skewed);
n2 = length(nodes1);

pv1 = [nodes_even zeros(n,1); ones(n-1,1) nodes_even(2:end); nodes4(end-1:-1:1) ones(n2-1,1); mid(2:end) nodes2(end-1:-1:1); nodes1(end-1:-1:1) mid(2:end) ;zeros(n2-2,1) nodes3(end-1:-1:2)];
% pv1 = [s zeros(n,1); ones(n-1,1) s(2:end); s(end-1:-1:1) ones(n-1,1); zeros(n-2,1) s(end-1:-1:2)];

H = 0.5; %dictates element size
[p,t]=polymesh({pv1},[1],[1,0],[H,1.75]);
[p,t] = fixmesh(p,t);

s1 = strcat('all(p(:,2)<',int2str(eval('ymin')),'+1e-6)');
s2 = strcat('all(p(:,1)>',int2str(eval('xmax')),'-1e-6)');

s3 = strcat('all(p(:,2)>',int2str(eval('ymax')),'-1e-6) & all(p(:,1)> 0.5 -1e-6)');
s4 = strcat('all(p(:,2)> 0.5 -1e-6) & all(p(:,1)> 0.5 -1e-6) & all(p(:,1)< 0.5 +1e-6)');
s5 = strcat('all(p(:,1)< 0.5 +1e-6) & all(p(:,2)> 0.5 -1e-6) & all(p(:,2)< 0.5 +1e-6)');

s6 = strcat('all(p(:,2)< 0.5 +1e-6) & all(p(:,1)< ',int2str(eval('xmin')),'+1e-6)');

bndexpr={s1,s2,s3,s4,s5,s6};   

mesh = mkmesh(p,t,porder,bndexpr,0,1);

function [x,y] = linenodes(x1,y1,x2,y2,N)

xi = loginc(linspace(0,0.5,N),3);
xi = [xi 1-xi(end-1:-1:1)]';
x = x1 + xi*(x2-x1);
y = y1 + xi*(y2-y1);

function y=loginc(x,alpha)

a=min(x(:));b=max(x(:));
y = (a + (b-a)*(exp(alpha*(x-a)/(b-a))-1)/(exp(alpha)-1));
