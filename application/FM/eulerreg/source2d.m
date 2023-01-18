function [s,s_udg] = source2d(pg,udg,param,time)
%SOURCE2D.M
%    [S,S_UDG] = SOURCE2D(PG,UDG,PARAM,TIME)
[ng,nc] = size(udg);
nch = 4;
s = zeros(ng,nch);
s_udg = zeros(ng,nch,nc);
