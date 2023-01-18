function [fb,fb_udg,fb_uh] = fbou2d(ib,uinf,nl,pg,udg,uh,param,time)
%FBOU2D.M
%    [FB,FB_UDG,FB_UH] = FBOU2D(IB,UINF,NL,PG,UDG,UH,PARAM,TIME)
[ng,nc] = size(udg);
nch = 4;
switch (ib)
    case 1
        fb = uinf-uh;
        fb_udg = zeros(ng,nch,nc);
        fb_uh = -ones(ng,nch,nch);
    case 2
        [fb,fb_udg,fb_uh] = fhat(nl,pg,udg,uh,param,time);
        fb = fb - uinf;
    otherwise
         error('unknown boundary type');
end
