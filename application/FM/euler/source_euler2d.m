function [s,s_udg] = source_euler2d(pg,udg,param,time)
%SOURCE_EULER2D.M
%    [S,S_UDG] = SOURCE_EULER2D(PG,UDG,PARAM,TIME)
[ng,nc] = size(udg);
nch = 4;
s = zeros(ng,nch);
s_udg = zeros(ng,nch,nc);
