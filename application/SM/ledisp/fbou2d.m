function [fb,fb_udg,fb_uh] = fbou2d(ib,uinf,nl,pg,udg,uh,param,time)
%FBOU2D.M
%    [FB,FB_UDG,FB_UH] = FBOU2D(IB,UINF,NL,PG,UDG,UH,PARAM,TIME)
[ng,nc] = size(udg);
nch = 2;
nd = 2;
switch (ib)
    case 1
        one = ones(ng,1);
        uh1 = uh(:,1);
        uh2 = uh(:,2);
        uinf1 = uinf(:,1);
        uinf2 = uinf(:,2);
        zero = zeros(ng,1);
        fb = [-one.*(uh1-uinf1);-uh2+uinf2];
        if nargout > 1
            fb_udg = [zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero];
        end
        if nargout > 2
            fb_uh = [-one;zero;zero;-one];
        end
        fb = reshape(fb,ng,nch);
        fb_udg = reshape(fb_udg,ng,nch,nc);
        fb_uh = reshape(fb_uh,ng,nch,nch);
    case 2
        nl1 = nl(:,1);
        nl2 = nl(:,2);
        one = ones(ng,1);
        param1 = param{1};
        param2 = param{2};
        param3 = param{3};
        u1 = udg(:,1);
        u2 = udg(:,2);
        u3 = udg(:,3);
        u4 = udg(:,4);
        u5 = udg(:,5);
        u6 = udg(:,6);
        uh1 = uh(:,1);
        uh2 = uh(:,2);
        uinf1 = uinf(:,1);
        uinf2 = uinf(:,2);
        zero = zeros(ng,1);
        t2 = u3+u6;
        t3 = param2.*t2;
        t4 = u4+u5;
        fb = [-uinf1+param3.*(u1-uh1)+nl1.*(t3+param1.*u3.*2.0)+nl2.*param1.*t4;-uinf2+param3.*(u2-uh2)+nl2.*(t3+param1.*u6.*2.0)+nl1.*param1.*t4];
        if nargout > 1
            t5 = one.*param3;
            t6 = nl2.*one.*param1;
            t7 = nl1.*one.*param1;
            t8 = param1.*2.0;
            t9 = param2+t8;
            fb_udg = [t5;zero;zero;t5;nl1.*one.*t9;nl2.*one.*param2;t6;t7;t6;t7;nl1.*one.*param2;nl2.*one.*t9];
        end
        if nargout > 2
            fb_uh = [-t5;zero;zero;-t5];
        end
        fb = reshape(fb,ng,nch);
        fb_udg = reshape(fb_udg,ng,nch,nc);
        fb_uh = reshape(fb_uh,ng,nch,nch);
    case 3
        nl1 = nl(:,1);
        nl2 = nl(:,2);
        one = ones(ng,1);
        param1 = param{1};
        param2 = param{2};
        param3 = param{3};
        u2 = udg(:,2);
        u3 = udg(:,3);
        u4 = udg(:,4);
        u5 = udg(:,5);
        u6 = udg(:,6);
        uh1 = uh(:,1);
        uh2 = uh(:,2);
        uinf1 = uinf(:,1);
        uinf2 = uinf(:,2);
        zero = zeros(ng,1);
        fb = [-one.*(uh1-uinf1);-uinf2+param3.*(u2-uh2)+nl2.*(param2.*(u3+u6)+param1.*u6.*2.0)+nl1.*param1.*(u4+u5)];
        if nargout > 1
            t2 = nl1.*one.*param1;
            t3 = one.*param3;
            fb_udg = [zero;zero;zero;t3;zero;nl2.*one.*param2;zero;t2;zero;t2;zero;nl2.*one.*(param1.*2.0+param2)];
        end
        if nargout > 2
            fb_uh = [-one;zero;zero;-t3];
        end
        fb = reshape(fb,ng,nch);
        fb_udg = reshape(fb_udg,ng,nch,nc);
        fb_uh = reshape(fb_uh,ng,nch,nch);
    case 4
        nl1 = nl(:,1);
        nl2 = nl(:,2);
        one = ones(ng,1);
        param1 = param{1};
        param2 = param{2};
        param3 = param{3};
        u1 = udg(:,1);
        u3 = udg(:,3);
        u4 = udg(:,4);
        u5 = udg(:,5);
        u6 = udg(:,6);
        uh1 = uh(:,1);
        uh2 = uh(:,2);
        uinf1 = uinf(:,1);
        uinf2 = uinf(:,2);
        zero = zeros(ng,1);
        fb = [-uinf1+param3.*(u1-uh1)+nl1.*(param2.*(u3+u6)+param1.*u3.*2.0)+nl2.*param1.*(u4+u5);-uh2+uinf2];
        if nargout > 1
            t2 = nl2.*one.*param1;
            t3 = one.*param3;
            fb_udg = [t3;zero;zero;zero;nl1.*one.*(param1.*2.0+param2);zero;t2;zero;t2;zero;nl1.*one.*param2;zero];
        end
        if nargout > 2
            fb_uh = [-t3;zero;zero;-one];
        end
        fb = reshape(fb,ng,nch);
        fb_udg = reshape(fb_udg,ng,nch,nc);
        fb_uh = reshape(fb_uh,ng,nch,nch);
    otherwise
         error('unknown boundary type');
end
