function [p,t]=prismmesh(n)

if nargin<1, n=2; end

[p,t]=cubemesh(n);
px=p(:,1); py=p(:,2); 
ix=find(all(px(t)+py(t)<=1,2));
t=t(ix,:);
[p,t]=fixmesh(p,t);
