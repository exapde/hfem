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
        x1 = pg(:,1);
        x2 = pg(:,2);
        zero = zeros(ng,1);
        fb = [one.*(-uh1+uinf1+x1);-uh2+uinf2+x2];
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
        lambda = param{2};
        mu = param{1};
        nl1 = nl(:,1);
        nl2 = nl(:,2);
        one = ones(ng,1);
        q11 = udg(:,3);
        q12 = udg(:,5);
        q21 = udg(:,4);
        q22 = udg(:,6);
        tau = param{3};
        u1 = udg(:,1);
        u2 = udg(:,2);
        uh1 = uh(:,1);
        uh2 = uh(:,2);
        uinf1 = uinf(:,1);
        uinf2 = uinf(:,2);
        zero = zeros(ng,1);
        t2 = q11.*q22;
        t4 = q12.*q21;
        t3 = t2-t4;
        t5 = 1.0./t3;
        t6 = log(t3);
        t8 = lambda.*t6;
        t7 = mu-t8;
        fb = [-uinf1+nl1.*(mu.*q11-q22.*t5.*t7)+nl2.*(mu.*q12+q21.*t5.*t7)+tau.*(u1-uh1);-uinf2+nl1.*(mu.*q21+q12.*t5.*t7)+nl2.*(mu.*q22-q11.*t5.*t7)+tau.*(u2-uh2)];
        if nargout > 1
            t9 = one.*tau;
            t10 = 1.0./t3.^2;
            t11 = q22.^2;
            t12 = q12.*q22.*t7.*t10;
            t13 = lambda.*q12.*q22.*t10;
            t14 = t12+t13;
            t15 = q12.^2;
            t16 = q21.*q22.*t7.*t10;
            t17 = lambda.*q21.*q22.*t10;
            t18 = t16+t17;
            t19 = q21.^2;
            t20 = t5.*t7;
            t21 = q12.*q21.*t7.*t10;
            t22 = lambda.*q12.*q21.*t10;
            t23 = t20+t21+t22;
            t24 = q11.*q22.*t7.*t10;
            t25 = lambda.*q11.*q22.*t10;
            t26 = q11.*q21.*t7.*t10;
            t27 = lambda.*q11.*q21.*t10;
            t28 = t26+t27;
            t29 = q11.*q12.*t7.*t10;
            t30 = lambda.*q11.*q12.*t10;
            t31 = t29+t30;
            t32 = q11.^2;
            fb_udg = [t9;zero;zero;t9;nl1.*(mu+lambda.*t10.*t11+t7.*t10.*t11)-nl2.*t18;-nl1.*t14+nl2.*(t24+t25-t5.*t7);-nl1.*t14+nl2.*t23;nl1.*(mu+lambda.*t10.*t15+t7.*t10.*t15)-nl2.*t31;nl2.*(mu+lambda.*t10.*t19+t7.*t10.*t19)-nl1.*t18;nl1.*t23-nl2.*t28;-nl2.*t28+nl1.*(-t20+t24+t25);nl2.*(mu+lambda.*t10.*t32+t7.*t10.*t32)-nl1.*t31];
        end
        if nargout > 2
            fb_uh = [-t9;zero;zero;-t9];
        end
        fb = reshape(fb,ng,nch);
        fb_udg = reshape(fb_udg,ng,nch,nc);
        fb_uh = reshape(fb_uh,ng,nch,nch);
    case 3
        lambda = param{2};
        mu = param{1};
        nl1 = nl(:,1);
        nl2 = nl(:,2);
        one = ones(ng,1);
        q11 = udg(:,3);
        q12 = udg(:,5);
        q21 = udg(:,4);
        q22 = udg(:,6);
        tau = param{3};
        u2 = udg(:,2);
        uh1 = uh(:,1);
        uh2 = uh(:,2);
        uinf1 = uinf(:,1);
        uinf2 = uinf(:,2);
        x1 = pg(:,1);
        zero = zeros(ng,1);
        t2 = q11.*q22;
        t4 = q12.*q21;
        t3 = t2-t4;
        t5 = 1.0./t3;
        t6 = log(t3);
        t8 = lambda.*t6;
        t7 = mu-t8;
        fb = [one.*(-uh1+uinf1+x1);-uinf2+nl1.*(mu.*q21+q12.*t5.*t7)+nl2.*(mu.*q22-q11.*t5.*t7)+tau.*(u2-uh2)];
        if nargout > 1
            t9 = 1.0./t3.^2;
            t10 = q12.^2;
            t11 = q11.*q12.*t7.*t9;
            t12 = lambda.*q11.*q12.*t9;
            t13 = t11+t12;
            t14 = q11.^2;
            t15 = one.*tau;
            fb_udg = [zero;zero;zero;t15;zero;-nl1.*(lambda.*q12.*q22.*t9+q12.*q22.*t7.*t9)+nl2.*(-t5.*t7+lambda.*q11.*q22.*t9+q11.*q22.*t7.*t9);zero;nl1.*(mu+lambda.*t9.*t10+t7.*t9.*t10)-nl2.*t13;zero;-nl2.*(lambda.*q11.*q21.*t9+q11.*q21.*t7.*t9)+nl1.*(t5.*t7+lambda.*q12.*q21.*t9+q12.*q21.*t7.*t9);zero;nl2.*(mu+lambda.*t9.*t14+t7.*t9.*t14)-nl1.*t13];
        end
        if nargout > 2
            fb_uh = [-one;zero;zero;-t15];
        end
        fb = reshape(fb,ng,nch);
        fb_udg = reshape(fb_udg,ng,nch,nc);
        fb_uh = reshape(fb_uh,ng,nch,nch);
    case 4
        lambda = param{2};
        mu = param{1};
        nl1 = nl(:,1);
        nl2 = nl(:,2);
        one = ones(ng,1);
        q11 = udg(:,3);
        q12 = udg(:,5);
        q21 = udg(:,4);
        q22 = udg(:,6);
        tau = param{3};
        u1 = udg(:,1);
        uh1 = uh(:,1);
        uh2 = uh(:,2);
        uinf1 = uinf(:,1);
        uinf2 = uinf(:,2);
        x2 = pg(:,2);
        zero = zeros(ng,1);
        t2 = q11.*q22;
        t4 = q12.*q21;
        t3 = t2-t4;
        t5 = 1.0./t3;
        t6 = log(t3);
        t8 = lambda.*t6;
        t7 = mu-t8;
        fb = [-uinf1+nl1.*(mu.*q11-q22.*t5.*t7)+nl2.*(mu.*q12+q21.*t5.*t7)+tau.*(u1-uh1);-uh2+uinf2+x2];
        if nargout > 1
            t9 = 1.0./t3.^2;
            t10 = q22.^2;
            t11 = q21.*q22.*t7.*t9;
            t12 = lambda.*q21.*q22.*t9;
            t13 = t11+t12;
            t14 = q21.^2;
            t15 = t5.*t7;
            t16 = one.*tau;
            fb_udg = [t16;zero;zero;zero;nl1.*(mu+lambda.*t9.*t10+t7.*t9.*t10)-nl2.*t13;zero;-nl1.*(lambda.*q12.*q22.*t9+q12.*q22.*t7.*t9)+nl2.*(t15+lambda.*q12.*q21.*t9+q12.*q21.*t7.*t9);zero;nl2.*(mu+lambda.*t9.*t14+t7.*t9.*t14)-nl1.*t13;zero;nl1.*(-t15+lambda.*q11.*q22.*t9+q11.*q22.*t7.*t9)-nl2.*(lambda.*q11.*q21.*t9+q11.*q21.*t7.*t9);zero];
        end
        if nargout > 2
            fb_uh = [-t16;zero;zero;-one];
        end
        fb = reshape(fb,ng,nch);
        fb_udg = reshape(fb_udg,ng,nch,nc);
        fb_uh = reshape(fb_uh,ng,nch,nch);
    otherwise
         error('unknown boundary type');
end
