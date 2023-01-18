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
        [fb,fb_udg,fb_uh] = fhat(nl,pg,udg,uh,param,time);
        fb = fb - uinf;
    case 3
        kappa = param{2};
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
        t2 = q12.*q21;
        t4 = q11.*q22;
        t3 = t2-t4+1.0;
        t5 = t2-t4;
        t6 = 1.0./t5;
        fb = [one.*(-uh1+uinf1+x1);uinf2+tau.*(u2-uh2)+nl2.*(mu.*q22-kappa.*q11.*t3+mu.*q11.*t6)+nl1.*(mu.*q21+kappa.*q12.*t3-(mu.*q12)./(t2-q11.*q22))];
        if nargout > 1
            t7 = 1.0./t5.^2;
            t8 = q12.^2;
            t9 = 1.0./(t2-t4).^2;
            t10 = mu.*t6;
            t11 = q11.^2;
            t12 = kappa.*q11.*q12;
            t13 = mu.*q11.*q12.*t9;
            t14 = t12+t13;
            t15 = one.*tau;
            fb_udg = [zero;zero;zero;t15;zero;nl2.*(t10-kappa.*t3+kappa.*q11.*q22+mu.*q11.*q22.*t7)-nl1.*(kappa.*q12.*q22+mu.*q12.*q22.*t7);zero;-nl2.*t14+nl1.*(mu+kappa.*t8+mu.*t8.*t9);zero;nl1.*(-t10+kappa.*t3+kappa.*q12.*q21+mu.*q12.*q21.*t9)-nl2.*(kappa.*q11.*q21+mu.*q11.*q21.*t9);zero;-nl1.*t14+nl2.*(mu+kappa.*t11+mu.*t9.*t11)];
        end
        if nargout > 2
            fb_uh = [-one;zero;zero;-t15];
        end
        fb = reshape(fb,ng,nch);
        fb_udg = reshape(fb_udg,ng,nch,nc);
        fb_uh = reshape(fb_uh,ng,nch,nch);
    case 4
        kappa = param{2};
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
        t2 = q12.*q21;
        t4 = q11.*q22;
        t3 = t2-t4+1.0;
        fb = [uinf1+nl1.*(mu.*q11+(mu.*q22)./(t2-t4)-kappa.*q22.*t3)+tau.*(u1-uh1)+nl2.*(mu.*q12+kappa.*q21.*t3-(mu.*q21)./(t2-q11.*q22));-uh2+uinf2+x2];
        if nargout > 1
            t5 = q22.^2;
            t7 = t2-t4;
            t6 = 1.0./t7.^2;
            t8 = 1.0./(t2-t4).^2;
            t9 = q21.^2;
            t10 = kappa.*q21.*q22;
            t11 = 1.0./t7;
            t12 = kappa.*t3;
            t13 = one.*tau;
            fb_udg = [t13;zero;zero;zero;nl1.*(mu+kappa.*t5+mu.*t5.*t6)-nl2.*(t10+mu.*q21.*q22.*t6);zero;nl2.*(t12-mu.*t11+kappa.*q12.*q21+mu.*q12.*q21.*t8)-nl1.*(kappa.*q12.*q22+mu.*q12.*q22.*t8);zero;nl2.*(mu+kappa.*t9+mu.*t8.*t9)-nl1.*(t10+mu.*q21.*q22.*t8);zero;nl1.*(-t12+mu.*t11+kappa.*q11.*q22+mu.*q11.*q22.*t8)-nl2.*(kappa.*q11.*q21+mu.*q11.*q21.*t8);zero];
        end
        if nargout > 2
            fb_uh = [-t13;zero;zero;-one];
        end
        fb = reshape(fb,ng,nch);
        fb_udg = reshape(fb_udg,ng,nch,nc);
        fb_uh = reshape(fb_uh,ng,nch,nch);
    otherwise
         error('unknown boundary type');
end
