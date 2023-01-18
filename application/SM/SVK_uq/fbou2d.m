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
        t2 = q11.^2;
        t3 = t2.*(1.0./2.0);
        t4 = q21.^2;
        t5 = t4.*(1.0./2.0);
        t6 = q12.^2;
        t7 = t6.*(1.0./2.0);
        t8 = q22.^2;
        t9 = t8.*(1.0./2.0);
        t10 = t3+t5+t7+t9-1.0;
        t11 = lambda.*t10;
        t12 = q11.*q12.*(1.0./2.0);
        t13 = q21.*q22.*(1.0./2.0);
        t14 = t12+t13;
        t15 = t3+t5-1.0./2.0;
        t16 = mu.*t15.*2.0;
        t17 = t11+t16;
        t18 = t7+t9-1.0./2.0;
        t19 = mu.*t18.*2.0;
        t20 = t11+t19;
        fb = [-uinf1+nl1.*(q11.*t17+mu.*q12.*t14.*2.0)+nl2.*(q12.*t20+mu.*q11.*t14.*2.0)+tau.*(u1-uh1);-uinf2+nl1.*(q21.*t17+mu.*q22.*t14.*2.0)+nl2.*(q22.*t20+mu.*q21.*t14.*2.0)+tau.*(u2-uh2)];
        if nargout > 1
            t21 = one.*tau;
            t22 = lambda.*q11;
            t23 = mu.*q11.*2.0;
            t24 = t22+t23;
            t25 = mu.*q12.*q22;
            t26 = mu.*t14.*2.0;
            t27 = lambda.*q21;
            t28 = mu.*q21.*2.0;
            t29 = t27+t28;
            t30 = lambda.*q11.*q12;
            t31 = mu.*q11.*q12;
            t32 = t26+t30+t31;
            t33 = lambda.*q12;
            t34 = mu.*q12.*2.0;
            t35 = t33+t34;
            t36 = lambda.*q12.*q21;
            t37 = mu.*q11.*q22;
            t38 = t36+t37;
            t39 = mu.*q11.*q21;
            t40 = lambda.*q11.*q22;
            t41 = mu.*q12.*q21;
            t42 = t40+t41;
            t43 = lambda.*q21.*q22;
            t44 = mu.*q21.*q22;
            t45 = t26+t43+t44;
            t46 = lambda.*q22;
            t47 = mu.*q22.*2.0;
            t48 = t46+t47;
            fb_udg = [t21;zero;zero;t21;nl2.*t32+nl1.*(t11+t16+mu.*t6+q11.*t24);nl2.*t42+nl1.*(t25+q21.*t24);nl2.*t38+nl1.*(t25+q11.*t29);nl2.*t45+nl1.*(t11+t16+mu.*t8+q21.*t29);nl1.*t32+nl2.*(t11+t19+mu.*t2+q12.*t35);nl1.*t38+nl2.*(t39+q22.*t35);nl1.*t42+nl2.*(t39+q12.*t48);nl1.*t45+nl2.*(t11+t19+mu.*t4+q22.*t48)];
        end
        if nargout > 2
            fb_uh = [-t21;zero;zero;-t21];
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
        t2 = q11.^2;
        t3 = t2.*(1.0./2.0);
        t4 = q21.^2;
        t5 = t4.*(1.0./2.0);
        t6 = q12.^2;
        t7 = t6.*(1.0./2.0);
        t8 = q22.^2;
        t9 = t8.*(1.0./2.0);
        t10 = t3+t5+t7+t9-1.0;
        t11 = lambda.*t10;
        t12 = q11.*q12.*(1.0./2.0);
        t13 = q21.*q22.*(1.0./2.0);
        t14 = t12+t13;
        t15 = t3+t5-1.0./2.0;
        t16 = mu.*t15.*2.0;
        t17 = mu.*t14.*2.0;
        t18 = lambda.*q21.*q22;
        t19 = mu.*q21.*q22;
        t20 = t17+t18+t19;
        t21 = t7+t9-1.0./2.0;
        t22 = mu.*t21.*2.0;
        fb = [one.*(-uh1+uinf1+x1);-uinf2+tau.*(u2-uh2)+nl1.*(q21.*(t11+t16)+mu.*q22.*t14.*2.0)+nl2.*(q22.*(t11+t22)+mu.*q21.*t14.*2.0)];
        if nargout > 1
            t23 = one.*tau;
            fb_udg = [zero;zero;zero;t23;zero;nl1.*(q21.*(lambda.*q11+mu.*q11.*2.0)+mu.*q12.*q22)+nl2.*(lambda.*q11.*q22+mu.*q12.*q21);zero;nl1.*(t11+t16+mu.*t8+q21.*(lambda.*q21+mu.*q21.*2.0))+nl2.*t20;zero;nl2.*(q22.*(lambda.*q12+mu.*q12.*2.0)+mu.*q11.*q21)+nl1.*(lambda.*q12.*q21+mu.*q11.*q22);zero;nl2.*(t11+t22+mu.*t4+q22.*(lambda.*q22+mu.*q22.*2.0))+nl1.*t20];
        end
        if nargout > 2
            fb_uh = [-one;zero;zero;-t23];
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
        t2 = q11.^2;
        t3 = t2.*(1.0./2.0);
        t4 = q21.^2;
        t5 = t4.*(1.0./2.0);
        t6 = q12.^2;
        t7 = t6.*(1.0./2.0);
        t8 = q22.^2;
        t9 = t8.*(1.0./2.0);
        t10 = t3+t5+t7+t9-1.0;
        t11 = lambda.*t10;
        t12 = q11.*q12.*(1.0./2.0);
        t13 = q21.*q22.*(1.0./2.0);
        t14 = t12+t13;
        t15 = t3+t5-1.0./2.0;
        t16 = mu.*t15.*2.0;
        t17 = mu.*t14.*2.0;
        t18 = lambda.*q11.*q12;
        t19 = mu.*q11.*q12;
        t20 = t17+t18+t19;
        t21 = t7+t9-1.0./2.0;
        t22 = mu.*t21.*2.0;
        fb = [-uinf1+tau.*(u1-uh1)+nl1.*(q11.*(t11+t16)+mu.*q12.*t14.*2.0)+nl2.*(q12.*(t11+t22)+mu.*q11.*t14.*2.0);-uh2+uinf2+x2];
        if nargout > 1
            t23 = one.*tau;
            fb_udg = [t23;zero;zero;zero;nl1.*(t11+t16+mu.*t6+q11.*(lambda.*q11+mu.*q11.*2.0))+nl2.*t20;zero;nl1.*(q11.*(lambda.*q21+mu.*q21.*2.0)+mu.*q12.*q22)+nl2.*(lambda.*q12.*q21+mu.*q11.*q22);zero;nl2.*(t11+t22+mu.*t2+q12.*(lambda.*q12+mu.*q12.*2.0))+nl1.*t20;zero;nl2.*(q12.*(lambda.*q22+mu.*q22.*2.0)+mu.*q11.*q21)+nl1.*(lambda.*q11.*q22+mu.*q12.*q21);zero];
        end
        if nargout > 2
            fb_uh = [-t23;zero;zero;-one];
        end
        fb = reshape(fb,ng,nch);
        fb_udg = reshape(fb_udg,ng,nch,nc);
        fb_uh = reshape(fb_uh,ng,nch,nch);
    otherwise
         error('unknown boundary type');
end
