function [s,s_udg] = source3d(pg,udg,param,time)
%SOURCE3D.M
%    [S,S_UDG] = SOURCE3D(PG,UDG,PARAM,TIME)
[ng,nc] = size(udg);
nch = 5;
s = zeros(ng,nch);
s_udg = zeros(ng,nch,nc);
