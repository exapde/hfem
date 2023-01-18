function [s,s_udg] = source2d(pg,udg,param,time)
%SOURCE2D
%    [S,S_UDG] = SOURCE2D(PG,UDG,PARAM,TIME)

%    This function was generated by the Symbolic Math Toolbox version 7.1.
%    16-Apr-2018 10:15:52
[ng,nc] = size(udg);
nch = 7;
nd = 2;
one = ones(ng,1);
param4 = param{4};
u7 = udg(:,7);
zero = zeros(ng,1);
s = [zero;zero;zero;zero;zero;zero;-param4.*u7];
if nargout > 1
    s_udg = [zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;-one.*param4];
end
s = reshape(s,ng,nch);
s_udg = reshape(s_udg,ng,nch,nc);