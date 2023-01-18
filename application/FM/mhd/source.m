function [s,s_udg] = source(pg,udg,param,time)
%SOURCE.M
%    [S,S_UDG] = SOURCE(PG,UDG,PARAM,TIME)
[ng,nc] = size(udg);
nch = 7;
s = zeros(ng,nch);
s_udg = zeros(ng,nch,nc);
