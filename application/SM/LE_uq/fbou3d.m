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
        x1 = pg(:,1);
        x2 = pg(:,2);
        x3 = pg(:,3);
        zero = zeros(ng,1);
        fb = [one.*(-uh1+uinf1+x1);-uh2+uinf2+x2;-uh3+uinf3+x3];
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
        lambda = param{2};
        mu = param{1};
        nl1 = nl(:,1);
        nl2 = nl(:,2);
        nl3 = nl(:,3);
        one = ones(ng,1);
        q11 = udg(:,4);
        q12 = udg(:,7);
        q13 = udg(:,10);
        q21 = udg(:,5);
        q22 = udg(:,8);
        q23 = udg(:,11);
        q31 = udg(:,6);
        q32 = udg(:,9);
        q33 = udg(:,12);
        tau = param{3};
        u1 = udg(:,1);
        u2 = udg(:,2);
        u3 = udg(:,3);
        uh1 = uh(:,1);
        uh2 = uh(:,2);
        uh3 = uh(:,3);
        uinf1 = uinf(:,1);
        uinf2 = uinf(:,2);
        uinf3 = uinf(:,3);
        zero = zeros(ng,1);
        t2 = q11+q22+q33+3.0;
        t3 = lambda.*t2;
        t4 = q12+q21;
        t5 = q13+q31;
        t6 = q23+q32;
        fb = [-uinf1+nl1.*(t3+mu.*(q11+1.0).*2.0)+tau.*(u1-uh1)+mu.*nl2.*t4+mu.*nl3.*t5;-uinf2+nl2.*(t3+mu.*(q22+1.0).*2.0)+tau.*(u2-uh2)+mu.*nl1.*t4+mu.*nl3.*t6;-uinf3+nl3.*(t3+mu.*(q33+1.0).*2.0)+tau.*(u3-uh3)+mu.*nl1.*t5+mu.*nl2.*t6];
        if nargout > 1
            t7 = one.*tau;
            t8 = mu.*nl1.*one;
            t9 = mu.*nl2.*one;
            t10 = mu.*2.0;
            t11 = lambda+t10;
            t12 = lambda.*nl3.*one;
            t13 = mu.*nl3.*one;
            t14 = lambda.*nl1.*one;
            t15 = lambda.*nl2.*one;
            fb_udg = [t7;zero;zero;zero;t7;zero;zero;zero;t7;nl1.*one.*t11;t15;t12;t9;t8;zero;t13;zero;t8;t9;t8;zero;t14;nl2.*one.*t11;t12;zero;t13;t9;t13;zero;t8;zero;t13;t9;t14;t15;nl3.*one.*t11];
        end
        if nargout > 2
            fb_uh = [-t7;zero;zero;zero;-t7;zero;zero;zero;-t7];
        end
        fb = reshape(fb,ng,nch);
        fb_udg = reshape(fb_udg,ng,nch,nc);
        fb_uh = reshape(fb_uh,ng,nch,nch);
    case 3
        lambda = param{2};
        mu = param{1};
        nl1 = nl(:,1);
        nl2 = nl(:,2);
        nl3 = nl(:,3);
        one = ones(ng,1);
        q11 = udg(:,4);
        q12 = udg(:,7);
        q13 = udg(:,10);
        q21 = udg(:,5);
        q22 = udg(:,8);
        q23 = udg(:,11);
        q31 = udg(:,6);
        q32 = udg(:,9);
        q33 = udg(:,12);
        tau = param{3};
        u2 = udg(:,2);
        u3 = udg(:,3);
        uh1 = uh(:,1);
        uh2 = uh(:,2);
        uh3 = uh(:,3);
        uinf1 = uinf(:,1);
        uinf2 = uinf(:,2);
        uinf3 = uinf(:,3);
        x1 = pg(:,1);
        zero = zeros(ng,1);
        t2 = q11+q22+q33+3.0;
        t3 = lambda.*t2;
        t4 = q23+q32;
        fb = [one.*(-uh1+uinf1+x1);-uinf2+nl2.*(t3+mu.*(q22+1.0).*2.0)+tau.*(u2-uh2)+mu.*nl1.*(q12+q21)+mu.*nl3.*t4;-uinf3+nl3.*(t3+mu.*(q33+1.0).*2.0)+tau.*(u3-uh3)+mu.*nl1.*(q13+q31)+mu.*nl2.*t4];
        if nargout > 1
            t5 = one.*tau;
            t6 = mu.*nl1.*one;
            t7 = lambda.*nl3.*one;
            t8 = mu.*nl3.*one;
            t9 = mu.*nl2.*one;
            t10 = lambda.*nl2.*one;
            t11 = mu.*2.0;
            t12 = lambda+t11;
            fb_udg = [zero;zero;zero;zero;t5;zero;zero;zero;t5;zero;t10;t7;zero;t6;zero;zero;zero;t6;zero;t6;zero;zero;nl2.*one.*t12;t7;zero;t8;t9;zero;zero;t6;zero;t8;t9;zero;t10;nl3.*one.*t12];
        end
        if nargout > 2
            fb_uh = [-one;zero;zero;zero;-t5;zero;zero;zero;-t5];
        end
        fb = reshape(fb,ng,nch);
        fb_udg = reshape(fb_udg,ng,nch,nc);
        fb_uh = reshape(fb_uh,ng,nch,nch);
    case 4
        lambda = param{2};
        mu = param{1};
        nl1 = nl(:,1);
        nl2 = nl(:,2);
        nl3 = nl(:,3);
        one = ones(ng,1);
        q11 = udg(:,4);
        q12 = udg(:,7);
        q13 = udg(:,10);
        q21 = udg(:,5);
        q22 = udg(:,8);
        q23 = udg(:,11);
        q31 = udg(:,6);
        q32 = udg(:,9);
        q33 = udg(:,12);
        tau = param{3};
        u1 = udg(:,1);
        u3 = udg(:,3);
        uh1 = uh(:,1);
        uh2 = uh(:,2);
        uh3 = uh(:,3);
        uinf1 = uinf(:,1);
        uinf2 = uinf(:,2);
        uinf3 = uinf(:,3);
        x2 = pg(:,2);
        zero = zeros(ng,1);
        t2 = q11+q22+q33+3.0;
        t3 = lambda.*t2;
        t4 = q13+q31;
        fb = [-uinf1+nl1.*(t3+mu.*(q11+1.0).*2.0)+tau.*(u1-uh1)+mu.*nl2.*(q12+q21)+mu.*nl3.*t4;-uh2+uinf2+x2;-uinf3+nl3.*(t3+mu.*(q33+1.0).*2.0)+tau.*(u3-uh3)+mu.*nl2.*(q23+q32)+mu.*nl1.*t4];
        if nargout > 1
            t5 = one.*tau;
            t6 = mu.*nl2.*one;
            t7 = lambda.*nl3.*one;
            t8 = mu.*nl3.*one;
            t9 = mu.*nl1.*one;
            t10 = lambda.*nl1.*one;
            t11 = mu.*2.0;
            t12 = lambda+t11;
            fb_udg = [t5;zero;zero;zero;zero;zero;zero;zero;t5;nl1.*one.*t12;zero;t7;t6;zero;zero;t8;zero;t9;t6;zero;zero;t10;zero;t7;zero;zero;t6;t8;zero;t9;zero;zero;t6;t10;zero;nl3.*one.*t12];
        end
        if nargout > 2
            fb_uh = [-t5;zero;zero;zero;-one;zero;zero;zero;-t5];
        end
        fb = reshape(fb,ng,nch);
        fb_udg = reshape(fb_udg,ng,nch,nc);
        fb_uh = reshape(fb_uh,ng,nch,nch);
    case 5
        lambda = param{2};
        mu = param{1};
        nl1 = nl(:,1);
        nl2 = nl(:,2);
        nl3 = nl(:,3);
        one = ones(ng,1);
        q11 = udg(:,4);
        q12 = udg(:,7);
        q13 = udg(:,10);
        q21 = udg(:,5);
        q22 = udg(:,8);
        q23 = udg(:,11);
        q31 = udg(:,6);
        q32 = udg(:,9);
        q33 = udg(:,12);
        tau = param{3};
        u1 = udg(:,1);
        u2 = udg(:,2);
        uh1 = uh(:,1);
        uh2 = uh(:,2);
        uh3 = uh(:,3);
        uinf1 = uinf(:,1);
        uinf2 = uinf(:,2);
        uinf3 = uinf(:,3);
        x3 = pg(:,3);
        zero = zeros(ng,1);
        t2 = q11+q22+q33+3.0;
        t3 = lambda.*t2;
        t4 = q12+q21;
        fb = [-uinf1+nl1.*(t3+mu.*(q11+1.0).*2.0)+tau.*(u1-uh1)+mu.*nl3.*(q13+q31)+mu.*nl2.*t4;-uinf2+nl2.*(t3+mu.*(q22+1.0).*2.0)+tau.*(u2-uh2)+mu.*nl3.*(q23+q32)+mu.*nl1.*t4;-uh3+uinf3+x3];
        if nargout > 1
            t5 = one.*tau;
            t6 = mu.*nl2.*one;
            t7 = mu.*nl1.*one;
            t8 = mu.*2.0;
            t9 = lambda+t8;
            t10 = mu.*nl3.*one;
            t11 = lambda.*nl1.*one;
            t12 = lambda.*nl2.*one;
            fb_udg = [t5;zero;zero;zero;t5;zero;zero;zero;zero;nl1.*one.*t9;t12;zero;t6;t7;zero;t10;zero;zero;t6;t7;zero;t11;nl2.*one.*t9;zero;zero;t10;zero;t10;zero;zero;zero;t10;zero;t11;t12;zero];
        end
        if nargout > 2
            fb_uh = [-t5;zero;zero;zero;-t5;zero;zero;zero;-one];
        end
        fb = reshape(fb,ng,nch);
        fb_udg = reshape(fb_udg,ng,nch,nc);
        fb_uh = reshape(fb_uh,ng,nch,nch);
    case 6
        lambda = param{2};
        mu = param{1};
        nl1 = nl(:,1);
        nl2 = nl(:,2);
        nl3 = nl(:,3);
        one = ones(ng,1);
        q11 = udg(:,4);
        q13 = udg(:,10);
        q22 = udg(:,8);
        q23 = udg(:,11);
        q31 = udg(:,6);
        q32 = udg(:,9);
        q33 = udg(:,12);
        tau = param{3};
        u3 = udg(:,3);
        uh1 = uh(:,1);
        uh2 = uh(:,2);
        uh3 = uh(:,3);
        uinf1 = uinf(:,1);
        uinf2 = uinf(:,2);
        uinf3 = uinf(:,3);
        x1 = pg(:,1);
        x2 = pg(:,2);
        zero = zeros(ng,1);
        fb = [one.*(-uh1+uinf1+x1);-uh2+uinf2+x2;-uinf3+nl3.*(lambda.*(q11+q22+q33+3.0)+mu.*(q33+1.0).*2.0)+tau.*(u3-uh3)+mu.*nl1.*(q13+q31)+mu.*nl2.*(q23+q32)];
        if nargout > 1
            t2 = lambda.*nl3.*one;
            t3 = mu.*nl1.*one;
            t4 = mu.*nl2.*one;
            t5 = one.*tau;
            fb_udg = [zero;zero;zero;zero;zero;zero;zero;zero;t5;zero;zero;t2;zero;zero;zero;zero;zero;t3;zero;zero;zero;zero;zero;t2;zero;zero;t4;zero;zero;t3;zero;zero;t4;zero;zero;nl3.*one.*(lambda+mu.*2.0)];
        end
        if nargout > 2
            fb_uh = [-one;zero;zero;zero;-one;zero;zero;zero;-t5];
        end
        fb = reshape(fb,ng,nch);
        fb_udg = reshape(fb_udg,ng,nch,nc);
        fb_uh = reshape(fb_uh,ng,nch,nch);
    case 7
        lambda = param{2};
        mu = param{1};
        nl1 = nl(:,1);
        nl2 = nl(:,2);
        nl3 = nl(:,3);
        one = ones(ng,1);
        q11 = udg(:,4);
        q12 = udg(:,7);
        q13 = udg(:,10);
        q21 = udg(:,5);
        q22 = udg(:,8);
        q31 = udg(:,6);
        q33 = udg(:,12);
        tau = param{3};
        u1 = udg(:,1);
        uh1 = uh(:,1);
        uh2 = uh(:,2);
        uh3 = uh(:,3);
        uinf1 = uinf(:,1);
        uinf2 = uinf(:,2);
        uinf3 = uinf(:,3);
        x2 = pg(:,2);
        x3 = pg(:,3);
        zero = zeros(ng,1);
        fb = [-uinf1+nl1.*(lambda.*(q11+q22+q33+3.0)+mu.*(q11+1.0).*2.0)+tau.*(u1-uh1)+mu.*nl2.*(q12+q21)+mu.*nl3.*(q13+q31);-uh2+uinf2+x2;-uh3+uinf3+x3];
        if nargout > 1
            t2 = mu.*nl2.*one;
            t3 = mu.*nl3.*one;
            t4 = lambda.*nl1.*one;
            t5 = one.*tau;
            fb_udg = [t5;zero;zero;zero;zero;zero;zero;zero;zero;nl1.*one.*(lambda+mu.*2.0);zero;zero;t2;zero;zero;t3;zero;zero;t2;zero;zero;t4;zero;zero;zero;zero;zero;t3;zero;zero;zero;zero;zero;t4;zero;zero];
        end
        if nargout > 2
            fb_uh = [-t5;zero;zero;zero;-one;zero;zero;zero;-one];
        end
        fb = reshape(fb,ng,nch);
        fb_udg = reshape(fb_udg,ng,nch,nc);
        fb_uh = reshape(fb_uh,ng,nch,nch);
    case 8
        lambda = param{2};
        mu = param{1};
        nl1 = nl(:,1);
        nl2 = nl(:,2);
        nl3 = nl(:,3);
        one = ones(ng,1);
        q11 = udg(:,4);
        q12 = udg(:,7);
        q21 = udg(:,5);
        q22 = udg(:,8);
        q23 = udg(:,11);
        q32 = udg(:,9);
        q33 = udg(:,12);
        tau = param{3};
        u2 = udg(:,2);
        uh1 = uh(:,1);
        uh2 = uh(:,2);
        uh3 = uh(:,3);
        uinf1 = uinf(:,1);
        uinf2 = uinf(:,2);
        uinf3 = uinf(:,3);
        x1 = pg(:,1);
        x3 = pg(:,3);
        zero = zeros(ng,1);
        fb = [one.*(-uh1+uinf1+x1);-uinf2+nl2.*(lambda.*(q11+q22+q33+3.0)+mu.*(q22+1.0).*2.0)+tau.*(u2-uh2)+mu.*nl1.*(q12+q21)+mu.*nl3.*(q23+q32);-uh3+uinf3+x3];
        if nargout > 1
            t2 = mu.*nl1.*one;
            t3 = mu.*nl3.*one;
            t4 = lambda.*nl2.*one;
            t5 = one.*tau;
            fb_udg = [zero;zero;zero;zero;t5;zero;zero;zero;zero;zero;t4;zero;zero;t2;zero;zero;zero;zero;zero;t2;zero;zero;nl2.*one.*(lambda+mu.*2.0);zero;zero;t3;zero;zero;zero;zero;zero;t3;zero;zero;t4;zero];
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
