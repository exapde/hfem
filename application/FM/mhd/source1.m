function [s,s_udg] = source(pg,udg,param,time)
%SOURCE.M
%    [S,S_UDG] = SOURCE(PG,UDG,PARAM,TIME)
[ng,nc] = size(udg);
nch = 7;
alpha1 = param{2};
alpha2 = sqrt(0.18 * alpha1);
alpha = alpha1^2/alpha2^2;

s = zeros(ng,nch);
s(:,7) = -udg(:,7)*alpha;

s_udg = zeros(ng,nch,nc);
s_udg(:,7,7) = -alpha; 

