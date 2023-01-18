function [fb,fb_udg,fb_uh] = fbou(ib,uinf,nl,pg,udg,uh,param,time)
%FBOU.M
%    [FB,FB_UDG,FB_UH] = FBOU(IB,UINF,NL,PG,UDG,UH,PARAM,TIME)
[ng,nc] = size(udg);
nch = 7;
switch (ib)
    case 1
        fb = uinf-uh;
        fb_udg = zeros(ng,nch,nc);
        fb_uh = -ones(ng,nch,nch);
    case 2
        [fb,fb_udg,fb_uh] = fhat(nl,pg,udg,uh,param,time);
        fb = fb - uinf;
    case 3 % periodic boundary conditions
        [fb,fb_udg,fb_uh] = fhat(nl,pg,udg,uh,param,time); 
    otherwise
         error('unknown boundary type');
end
