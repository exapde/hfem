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
        t2 = q11+q22+2.0;
        t3 = lambda.*t2;
        t4 = q12+q21;
        fb = [-uinf1+nl1.*(t3+mu.*(q11+1.0).*2.0)+tau.*(u1-uh1)+mu.*nl2.*t4;-uinf2+nl2.*(t3+mu.*(q22+1.0).*2.0)+tau.*(u2-uh2)+mu.*nl1.*t4];
        if nargout > 1
            t5 = one.*tau;
            t6 = mu.*nl2.*one;
            t7 = mu.*nl1.*one;
            t8 = mu.*2.0;
            t9 = lambda+t8;
            fb_udg = [t5;zero;zero;t5;nl1.*one.*t9;lambda.*nl2.*one;t6;t7;t6;t7;lambda.*nl1.*one;nl2.*one.*t9];
        end
        if nargout > 2
            fb_uh = [-t5;zero;zero;-t5];
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
        fb = [one.*(-uh1+uinf1+x1);-uinf2+tau.*(u2-uh2)+nl2.*(lambda.*(q11+q22+2.0)+mu.*(q22+1.0).*2.0)+mu.*nl1.*(q12+q21)];
        if nargout > 1
            t2 = mu.*nl1.*one;
            t3 = one.*tau;
            fb_udg = [zero;zero;zero;t3;zero;lambda.*nl2.*one;zero;t2;zero;t2;zero;nl2.*one.*(lambda+mu.*2.0)];
        end
        if nargout > 2
            fb_uh = [-one;zero;zero;-t3];
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
        fb = [-uinf1+tau.*(u1-uh1)+nl1.*(lambda.*(q11+q22+2.0)+mu.*(q11+1.0).*2.0)+mu.*nl2.*(q12+q21);-uh2+uinf2+x2];
        if nargout > 1
            t2 = mu.*nl2.*one;
            t3 = one.*tau;
            fb_udg = [t3;zero;zero;zero;nl1.*one.*(lambda+mu.*2.0);zero;t2;zero;t2;zero;lambda.*nl1.*one;zero];
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
