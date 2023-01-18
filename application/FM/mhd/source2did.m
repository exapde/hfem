function [s,s_udg] = source2did(pg,udg,param,time)
%SOURCE2DID.M
%    [S,S_UDG] = SOURCE2DID(PG,UDG,PARAM,TIME)
[ng,nc] = size(udg);
nch = 8;
s = zeros(ng,nch);
s_udg = zeros(ng,nch,nc);
