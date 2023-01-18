function [s,s_udg] = source2(pg,udg,param,time)
%SOURCE2.M
%    [S,S_UDG] = SOURCE2(PG,UDG,PARAM,TIME)
[ng,nc] = size(udg);
nch = 6;
s = zeros(ng,nch);
s_udg = zeros(ng,nch,nc);
