function [fb,fb_udg,fb_uh] = fbou3d(ib,uinf,nl,pg,udg,uh,param,time)
%FBOU3D.M
%    [FB,FB_UDG,FB_UH] = FBOU3D(IB,UINF,NL,PG,UDG,UH,PARAM,TIME)
[ng,nc] = size(udg);
nch = 3;
nd = 3;
switch (ib)
    case 1
        one = ones(ng,1);
        uh1 = uh(:,1);
        uh2 = uh(:,2);
        uh3 = uh(:,3);
        uinf1 = uinf(:,1);
        uinf2 = uinf(:,2);
        uinf3 = uinf(:,3);
        zero = zeros(ng,1);
        fb = [-one.*(uh1-uinf1);-uh2+uinf2;-uh3+uinf3];
        if nargout > 1
            fb_udg = [zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero];
        end
        if nargout > 2
            fb_uh = [-one;zero;zero;zero;-one;zero;zero;zero;-one];
        end
        fb = reshape(fb,ng,nch);
        fb_udg = reshape(fb_udg,ng,nch,nc);
        fb_uh = reshape(fb_uh,ng,nch,nch);
    case 2
        nl1 = nl(:,1);
        nl2 = nl(:,2);
        nl3 = nl(:,3);
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
        u7 = udg(:,7);
        u8 = udg(:,8);
        u9 = udg(:,9);
        u10 = udg(:,10);
        u11 = udg(:,11);
        u12 = udg(:,12);
        uh1 = uh(:,1);
        uh2 = uh(:,2);
        uh3 = uh(:,3);
        uinf1 = uinf(:,1);
        uinf2 = uinf(:,2);
        uinf3 = uinf(:,3);
        zero = zeros(ng,1);
        t2 = u4+u8+u12;
        t3 = param2.*t2;
        t4 = u5+u7;
        t5 = u6+u10;
        t6 = u9+u11;
        fb = [-uinf1+param3.*(u1-uh1)+nl1.*(t3+param1.*u4.*2.0)+nl2.*param1.*t4+nl3.*param1.*t5;-uinf2+param3.*(u2-uh2)+nl2.*(t3+param1.*u8.*2.0)+nl1.*param1.*t4+nl3.*param1.*t6;-uinf3+param3.*(u3-uh3)+nl3.*(t3+param1.*u12.*2.0)+nl1.*param1.*t5+nl2.*param1.*t6];
        if nargout > 1
            t7 = one.*param3;
            t8 = nl1.*one.*param1;
            t9 = nl2.*one.*param1;
            t10 = param1.*2.0;
            t11 = param2+t10;
            t12 = nl3.*one.*param2;
            t13 = nl3.*one.*param1;
            t14 = nl1.*one.*param2;
            t15 = nl2.*one.*param2;
            fb_udg = [t7;zero;zero;zero;t7;zero;zero;zero;t7;nl1.*one.*t11;t15;t12;t9;t8;zero;t13;zero;t8;t9;t8;zero;t14;nl2.*one.*t11;t12;zero;t13;t9;t13;zero;t8;zero;t13;t9;t14;t15;nl3.*one.*t11];
        end
        if nargout > 2
            fb_uh = [-t7;zero;zero;zero;-t7;zero;zero;zero;-t7];
        end
        fb = reshape(fb,ng,nch);
        fb_udg = reshape(fb_udg,ng,nch,nc);
        fb_uh = reshape(fb_uh,ng,nch,nch);
    case 3
        nl1 = nl(:,1);
        nl2 = nl(:,2);
        nl3 = nl(:,3);
        one = ones(ng,1);
        param1 = param{1};
        param2 = param{2};
        param3 = param{3};
        u2 = udg(:,2);
        u3 = udg(:,3);
        u4 = udg(:,4);
        u5 = udg(:,5);
        u6 = udg(:,6);
        u7 = udg(:,7);
        u8 = udg(:,8);
        u9 = udg(:,9);
        u10 = udg(:,10);
        u11 = udg(:,11);
        u12 = udg(:,12);
        uh1 = uh(:,1);
        uh2 = uh(:,2);
        uh3 = uh(:,3);
        uinf1 = uinf(:,1);
        uinf2 = uinf(:,2);
        uinf3 = uinf(:,3);
        zero = zeros(ng,1);
        t2 = u4+u8+u12;
        t3 = param2.*t2;
        t4 = u9+u11;
        fb = [-one.*(uh1-uinf1);-uinf2+param3.*(u2-uh2)+nl2.*(t3+param1.*u8.*2.0)+nl1.*param1.*(u5+u7)+nl3.*param1.*t4;-uinf3+param3.*(u3-uh3)+nl3.*(t3+param1.*u12.*2.0)+nl1.*param1.*(u6+u10)+nl2.*param1.*t4];
        if nargout > 1
            t5 = one.*param3;
            t6 = nl1.*one.*param1;
            t7 = nl3.*one.*param2;
            t8 = nl3.*one.*param1;
            t9 = nl2.*one.*param1;
            t10 = nl2.*one.*param2;
            t11 = param1.*2.0;
            t12 = param2+t11;
            fb_udg = [zero;zero;zero;zero;t5;zero;zero;zero;t5;zero;t10;t7;zero;t6;zero;zero;zero;t6;zero;t6;zero;zero;nl2.*one.*t12;t7;zero;t8;t9;zero;zero;t6;zero;t8;t9;zero;t10;nl3.*one.*t12];
        end
        if nargout > 2
            fb_uh = [-one;zero;zero;zero;-t5;zero;zero;zero;-t5];
        end
        fb = reshape(fb,ng,nch);
        fb_udg = reshape(fb_udg,ng,nch,nc);
        fb_uh = reshape(fb_uh,ng,nch,nch);
    case 4
        nl1 = nl(:,1);
        nl2 = nl(:,2);
        nl3 = nl(:,3);
        one = ones(ng,1);
        param1 = param{1};
        param2 = param{2};
        param3 = param{3};
        u1 = udg(:,1);
        u3 = udg(:,3);
        u4 = udg(:,4);
        u5 = udg(:,5);
        u6 = udg(:,6);
        u7 = udg(:,7);
        u8 = udg(:,8);
        u9 = udg(:,9);
        u10 = udg(:,10);
        u11 = udg(:,11);
        u12 = udg(:,12);
        uh1 = uh(:,1);
        uh2 = uh(:,2);
        uh3 = uh(:,3);
        uinf1 = uinf(:,1);
        uinf2 = uinf(:,2);
        uinf3 = uinf(:,3);
        zero = zeros(ng,1);
        t2 = u4+u8+u12;
        t3 = param2.*t2;
        t4 = u6+u10;
        fb = [-uinf1+param3.*(u1-uh1)+nl1.*(t3+param1.*u4.*2.0)+nl2.*param1.*(u5+u7)+nl3.*param1.*t4;-uh2+uinf2;-uinf3+param3.*(u3-uh3)+nl3.*(t3+param1.*u12.*2.0)+nl2.*param1.*(u9+u11)+nl1.*param1.*t4];
        if nargout > 1
            t5 = one.*param3;
            t6 = nl2.*one.*param1;
            t7 = nl3.*one.*param2;
            t8 = nl3.*one.*param1;
            t9 = nl1.*one.*param1;
            t10 = nl1.*one.*param2;
            t11 = param1.*2.0;
            t12 = param2+t11;
            fb_udg = [t5;zero;zero;zero;zero;zero;zero;zero;t5;nl1.*one.*t12;zero;t7;t6;zero;zero;t8;zero;t9;t6;zero;zero;t10;zero;t7;zero;zero;t6;t8;zero;t9;zero;zero;t6;t10;zero;nl3.*one.*t12];
        end
        if nargout > 2
            fb_uh = [-t5;zero;zero;zero;-one;zero;zero;zero;-t5];
        end
        fb = reshape(fb,ng,nch);
        fb_udg = reshape(fb_udg,ng,nch,nc);
        fb_uh = reshape(fb_uh,ng,nch,nch);
    case 5
        nl1 = nl(:,1);
        nl2 = nl(:,2);
        nl3 = nl(:,3);
        one = ones(ng,1);
        param1 = param{1};
        param2 = param{2};
        param3 = param{3};
        u1 = udg(:,1);
        u2 = udg(:,2);
        u4 = udg(:,4);
        u5 = udg(:,5);
        u6 = udg(:,6);
        u7 = udg(:,7);
        u8 = udg(:,8);
        u9 = udg(:,9);
        u10 = udg(:,10);
        u11 = udg(:,11);
        u12 = udg(:,12);
        uh1 = uh(:,1);
        uh2 = uh(:,2);
        uh3 = uh(:,3);
        uinf1 = uinf(:,1);
        uinf2 = uinf(:,2);
        uinf3 = uinf(:,3);
        zero = zeros(ng,1);
        t2 = u4+u8+u12;
        t3 = param2.*t2;
        t4 = u5+u7;
        fb = [-uinf1+param3.*(u1-uh1)+nl1.*(t3+param1.*u4.*2.0)+nl3.*param1.*(u6+u10)+nl2.*param1.*t4;-uinf2+param3.*(u2-uh2)+nl2.*(t3+param1.*u8.*2.0)+nl3.*param1.*(u9+u11)+nl1.*param1.*t4;-uh3+uinf3];
        if nargout > 1
            t5 = one.*param3;
            t6 = nl2.*one.*param1;
            t7 = nl1.*one.*param1;
            t8 = param1.*2.0;
            t9 = param2+t8;
            t10 = nl3.*one.*param1;
            t11 = nl1.*one.*param2;
            t12 = nl2.*one.*param2;
            fb_udg = [t5;zero;zero;zero;t5;zero;zero;zero;zero;nl1.*one.*t9;t12;zero;t6;t7;zero;t10;zero;zero;t6;t7;zero;t11;nl2.*one.*t9;zero;zero;t10;zero;t10;zero;zero;zero;t10;zero;t11;t12;zero];
        end
        if nargout > 2
            fb_uh = [-t5;zero;zero;zero;-t5;zero;zero;zero;-one];
        end
        fb = reshape(fb,ng,nch);
        fb_udg = reshape(fb_udg,ng,nch,nc);
        fb_uh = reshape(fb_uh,ng,nch,nch);
    case 6
        nl1 = nl(:,1);
        nl2 = nl(:,2);
        nl3 = nl(:,3);
        one = ones(ng,1);
        param1 = param{1};
        param2 = param{2};
        param3 = param{3};
        u3 = udg(:,3);
        u4 = udg(:,4);
        u6 = udg(:,6);
        u8 = udg(:,8);
        u9 = udg(:,9);
        u10 = udg(:,10);
        u11 = udg(:,11);
        u12 = udg(:,12);
        uh1 = uh(:,1);
        uh2 = uh(:,2);
        uh3 = uh(:,3);
        uinf1 = uinf(:,1);
        uinf2 = uinf(:,2);
        uinf3 = uinf(:,3);
        zero = zeros(ng,1);
        fb = [-one.*(uh1-uinf1);-uh2+uinf2;-uinf3+nl3.*(param1.*u12.*2.0+param2.*(u4+u8+u12))+param3.*(u3-uh3)+nl1.*param1.*(u6+u10)+nl2.*param1.*(u9+u11)];
        if nargout > 1
            t2 = nl3.*one.*param2;
            t3 = nl1.*one.*param1;
            t4 = nl2.*one.*param1;
            t5 = one.*param3;
            fb_udg = [zero;zero;zero;zero;zero;zero;zero;zero;t5;zero;zero;t2;zero;zero;zero;zero;zero;t3;zero;zero;zero;zero;zero;t2;zero;zero;t4;zero;zero;t3;zero;zero;t4;zero;zero;nl3.*one.*(param1.*2.0+param2)];
        end
        if nargout > 2
            fb_uh = [-one;zero;zero;zero;-one;zero;zero;zero;-t5];
        end
        fb = reshape(fb,ng,nch);
        fb_udg = reshape(fb_udg,ng,nch,nc);
        fb_uh = reshape(fb_uh,ng,nch,nch);
    case 7
        nl1 = nl(:,1);
        nl2 = nl(:,2);
        nl3 = nl(:,3);
        one = ones(ng,1);
        param1 = param{1};
        param2 = param{2};
        param3 = param{3};
        u1 = udg(:,1);
        u4 = udg(:,4);
        u5 = udg(:,5);
        u6 = udg(:,6);
        u7 = udg(:,7);
        u8 = udg(:,8);
        u10 = udg(:,10);
        u12 = udg(:,12);
        uh1 = uh(:,1);
        uh2 = uh(:,2);
        uh3 = uh(:,3);
        uinf1 = uinf(:,1);
        uinf2 = uinf(:,2);
        uinf3 = uinf(:,3);
        zero = zeros(ng,1);
        fb = [-uinf1+nl1.*(param1.*u4.*2.0+param2.*(u4+u8+u12))+param3.*(u1-uh1)+nl2.*param1.*(u5+u7)+nl3.*param1.*(u6+u10);-uh2+uinf2;-uh3+uinf3];
        if nargout > 1
            t2 = nl2.*one.*param1;
            t3 = nl3.*one.*param1;
            t4 = nl1.*one.*param2;
            t5 = one.*param3;
            fb_udg = [t5;zero;zero;zero;zero;zero;zero;zero;zero;nl1.*one.*(param1.*2.0+param2);zero;zero;t2;zero;zero;t3;zero;zero;t2;zero;zero;t4;zero;zero;zero;zero;zero;t3;zero;zero;zero;zero;zero;t4;zero;zero];
        end
        if nargout > 2
            fb_uh = [-t5;zero;zero;zero;-one;zero;zero;zero;-one];
        end
        fb = reshape(fb,ng,nch);
        fb_udg = reshape(fb_udg,ng,nch,nc);
        fb_uh = reshape(fb_uh,ng,nch,nch);
    case 8
        nl1 = nl(:,1);
        nl2 = nl(:,2);
        nl3 = nl(:,3);
        one = ones(ng,1);
        param1 = param{1};
        param2 = param{2};
        param3 = param{3};
        u2 = udg(:,2);
        u4 = udg(:,4);
        u5 = udg(:,5);
        u7 = udg(:,7);
        u8 = udg(:,8);
        u9 = udg(:,9);
        u11 = udg(:,11);
        u12 = udg(:,12);
        uh1 = uh(:,1);
        uh2 = uh(:,2);
        uh3 = uh(:,3);
        uinf1 = uinf(:,1);
        uinf2 = uinf(:,2);
        uinf3 = uinf(:,3);
        zero = zeros(ng,1);
        fb = [-one.*(uh1-uinf1);-uinf2+nl2.*(param1.*u8.*2.0+param2.*(u4+u8+u12))+param3.*(u2-uh2)+nl1.*param1.*(u5+u7)+nl3.*param1.*(u9+u11);-uh3+uinf3];
        if nargout > 1
            t2 = nl1.*one.*param1;
            t3 = nl3.*one.*param1;
            t4 = nl2.*one.*param2;
            t5 = one.*param3;
            fb_udg = [zero;zero;zero;zero;t5;zero;zero;zero;zero;zero;t4;zero;zero;t2;zero;zero;zero;zero;zero;t2;zero;zero;nl2.*one.*(param1.*2.0+param2);zero;zero;t3;zero;zero;zero;zero;zero;t3;zero;zero;t4;zero];
        end
        if nargout > 2
            fb_uh = [-one;zero;zero;zero;-t5;zero;zero;zero;-one];
        end
        fb = reshape(fb,ng,nch);
        fb_udg = reshape(fb_udg,ng,nch,nc);
        fb_uh = reshape(fb_uh,ng,nch,nch);
    otherwise
         error('unknown boundary type');
end
