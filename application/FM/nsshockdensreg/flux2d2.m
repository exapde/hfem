function [f,f_udg] = flux2d2(pg,udg,param,time)
%FLUX2D2
%    [F,F_UDG] = FLUX2D2(PG,UDG,PARAM,TIME)

%    This function was generated by the Symbolic Math Toolbox version 5.8.
%    30-Mar-2013 23:56:39
[ng,nc] = size(udg);
nch = 4;
nd = 2;
one = ones(ng,1);
param1 = param{1};
param3 = param{3};
param4 = param{4};
param8 = param{8};
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
zero = zeros(ng,1);
t2 = 1.0./param8;
t3 = t2.*u1.*2.0e1;
t4 = exp(t3);
t5 = t4+2.0e1;
t6 = log(t5);
t7 = 1.0./t6;
t8 = 1.0./t5;
t9 = 1.0./param8.^2;
t10 = 1.0./t6.^2;
t11 = u2.^2;
t12 = 1.0./param3;
t13 = t9.*t10.*t11.*2.0e2;
t14 = u3.^2;
t15 = t9.*t10.*t14.*2.0e2;
t16 = t13+t15;
t29 = param8.*t6.*t16.*(1.0./2.0e1);
t17 = -t29+u4;
t18 = param1-1.0;
t28 = t2.*t4.*t7.*t8.*u3.*u5.*2.0e1;
t19 = -t28+u7;
t20 = t2.*t7.*t19.*2.0e1;
t30 = t2.*t4.*t7.*t8.*u2.*u9.*2.0e1;
t21 = -t30+u10;
t22 = t2.*t7.*t21.*2.0e1;
t23 = t20+t22;
t27 = t2.*t4.*t7.*t8.*u2.*u5.*2.0e1;
t24 = -t27+u6;
t34 = t2.*t4.*t7.*t8.*u3.*u9.*2.0e1;
t25 = -t34+u11;
t71 = t2.*t7.*t24.*4.0e1;
t72 = t2.*t7.*t25.*2.0e1;
t26 = t71-t72;
t31 = t12.*t23;
t32 = t2.*t7.*u2.*u3.*2.0e1;
t33 = t31+t32;
t35 = t17.*t18;
t36 = t2.*t7.*u4.*2.0e1;
t37 = t2.*t7.*t17.*t18.*2.0e1;
t38 = t36+t37;
t39 = t2.*t7.*t24.*2.0e1;
t98 = t2.*t7.*t25.*4.0e1;
t40 = t39-t98;
t41 = 1.0./param4;
t42 = 1.0./t18;
t43 = t2.*u1.*4.0e1;
t44 = exp(t43);
t45 = 1.0./t5.^2;
t46 = 1.0./param8.^3;
t47 = 1.0./t6.^3;
t48 = t4.*t8.*t16;
t49 = t4.*t8.*t11.*t46.*t47.*8.0e3;
t50 = t4.*t8.*t14.*t46.*t47.*8.0e3;
t51 = t49+t50;
t81 = param8.*t6.*t51.*(1.0./2.0e1);
t52 = t48-t81;
t53 = t7.*t9.*t44.*t45.*u3.*u5.*4.0e2;
t54 = t9.*t10.*t44.*t45.*u3.*u5.*4.0e2;
t77 = t4.*t7.*t8.*t9.*u3.*u5.*4.0e2;
t55 = t53+t54-t77;
t56 = t2.*t7.*t55.*2.0e1;
t57 = t7.*t9.*t44.*t45.*u2.*u9.*4.0e2;
t58 = t9.*t10.*t44.*t45.*u2.*u9.*4.0e2;
t85 = t4.*t7.*t8.*t9.*u2.*u9.*4.0e2;
t59 = t57+t58-t85;
t60 = t2.*t7.*t59.*2.0e1;
t86 = t4.*t8.*t9.*t10.*t19.*4.0e2;
t87 = t4.*t8.*t9.*t10.*t21.*4.0e2;
t61 = t56+t60-t86-t87;
t62 = t7.*t9.*t44.*t45.*u2.*u5.*4.0e2;
t63 = t9.*t10.*t44.*t45.*u2.*u5.*4.0e2;
t76 = t4.*t7.*t8.*t9.*u2.*u5.*4.0e2;
t64 = t62+t63-t76;
t65 = t2.*t7.*t64.*4.0e1;
t66 = t7.*t9.*t44.*t45.*u3.*u9.*4.0e2;
t67 = t9.*t10.*t44.*t45.*u3.*u9.*4.0e2;
t90 = t4.*t7.*t8.*t9.*u3.*u9.*4.0e2;
t68 = t66+t67-t90;
t69 = t4.*t8.*t9.*t10.*t25.*4.0e2;
t70 = t65+t69-t2.*t7.*t68.*2.0e1-t4.*t8.*t9.*t10.*t24.*8.0e2;
t73 = t9.*t10.*t24.*u2.*4.0e2;
t74 = t9.*t10.*t19.*u3.*4.0e2;
t75 = t73+t74;
t78 = param8.*t6.*t75.*(1.0./2.0e1);
t79 = t4.*t8.*t16.*u5;
t80 = t78+t79-u8;
t82 = param8.*t6.*t18.*t80.*(1.0./2.0e1);
t83 = t4.*t8.*t17.*t18.*u5;
t84 = t82+t83;
t88 = t12.*t61;
t89 = t88-t4.*t8.*t9.*t10.*u2.*u3.*4.0e2;
t91 = t2.*t7.*t18.*t52.*2.0e1;
t92 = t4.*t8.*t9.*t10.*u4.*4.0e2;
t93 = t4.*t8.*t9.*t10.*t17.*t18.*4.0e2;
t94 = t91+t92+t93;
t95 = t2.*t7.*t64.*2.0e1;
t96 = t4.*t8.*t9.*t10.*t25.*8.0e2;
t97 = t95+t96-t2.*t7.*t68.*4.0e1-t4.*t8.*t9.*t10.*t24.*4.0e2;
t99 = t9.*t10.*t21.*u2.*4.0e2;
t100 = t9.*t10.*t25.*u3.*4.0e2;
t101 = t99+t100;
t102 = param8.*t6.*t101.*(1.0./2.0e1);
t103 = t4.*t8.*t16.*u9;
t104 = t102+t103-u12;
t105 = param8.*t6.*t18.*t104.*(1.0./2.0e1);
t106 = t4.*t8.*t17.*t18.*u9;
t107 = t105+t106;
f = [u2;t35+t12.*t26.*(2.0./3.0)+t2.*t7.*t11.*2.0e1;t33;t38.*u2+t2.*t7.*t12.*t23.*u3.*2.0e1+t2.*t7.*t12.*t26.*u2.*(4.0e1./3.0)-param1.*t9.*t10.*t12.*t41.*t42.*t84.*4.0e2;u3;t33;t35-t12.*t40.*(2.0./3.0)+t2.*t7.*t14.*2.0e1;t38.*u3+t2.*t7.*t12.*t23.*u2.*2.0e1-t2.*t7.*t12.*t40.*u3.*(4.0e1./3.0)-param1.*t9.*t10.*t12.*t41.*t42.*t107.*4.0e2];
if nargout > 1
    t108 = t2.*t7.*u3.*2.0e1;
    t109 = t108-t4.*t8.*t9.*t10.*t12.*u9.*4.0e2;
    t110 = t2.*t7.*t12.*t23.*2.0e1;
    t111 = t2.*t7.*u2.*2.0e1;
    t112 = t111-t4.*t8.*t9.*t10.*t12.*u5.*4.0e2;
    t113 = one.*t18;
    t114 = t2.*t7.*2.0e1;
    t115 = t2.*t7.*t18.*2.0e1;
    t116 = t114+t115;
    t117 = t2.*t7.*t12.*2.0e1;
    t118 = t4.*t8.*t17.*t18;
    t119 = param8.*t6.*t18.*t52.*(1.0./2.0e1);
    t120 = t118+t119;
    t121 = t9.*t10.*t12.*u3.*4.0e2;
    t122 = t9.*t10.*t12.*u2.*4.0e2;
    t123 = t2.*t7.*t12.*(8.0e1./3.0);
    t124 = param1.*t2.*t7.*t12.*t41.*2.0e1;
    f_udg = [zero;-t18.*t52+t12.*t70.*(2.0./3.0)-t4.*t8.*t9.*t10.*t11.*4.0e2;t89;-t94.*u2+t2.*t7.*t12.*t61.*u3.*2.0e1+t2.*t7.*t12.*t70.*u2.*(4.0e1./3.0)-t4.*t8.*t9.*t10.*t12.*t23.*u3.*4.0e2-t4.*t8.*t9.*t10.*t12.*t26.*u2.*(8.0e2./3.0)-param1.*t9.*t10.*t12.*t41.*t42.*(t4.*t8.*t18.*t80+param8.*t6.*t18.*(param8.*t6.*(t9.*t10.*t55.*u3.*4.0e2+t9.*t10.*t64.*u2.*4.0e2-t4.*t8.*t19.*t46.*t47.*u3.*1.6e4-t4.*t8.*t24.*t46.*t47.*u2.*1.6e4).*(1.0./2.0e1)+t4.*t8.*t75-t4.*t8.*t51.*u5+t2.*t4.*t8.*t16.*u5.*2.0e1-t2.*t16.*t44.*t45.*u5.*2.0e1).*(1.0./2.0e1)-t4.*t8.*t18.*t52.*u5+t2.*t4.*t8.*t17.*t18.*u5.*2.0e1-t2.*t17.*t18.*t44.*t45.*u5.*2.0e1).*4.0e2+param1.*t4.*t8.*t12.*t41.*t42.*t46.*t47.*t84.*1.6e4;zero;t89;-t18.*t52-t12.*t97.*(2.0./3.0)-t4.*t8.*t9.*t10.*t14.*4.0e2;-t94.*u3+t2.*t7.*t12.*t61.*u2.*2.0e1-t2.*t7.*t12.*t97.*u3.*(4.0e1./3.0)-t4.*t8.*t9.*t10.*t12.*t23.*u2.*4.0e2+t4.*t8.*t9.*t10.*t12.*t40.*u3.*(8.0e2./3.0)-param1.*t9.*t10.*t12.*t41.*t42.*(t4.*t8.*t18.*t104+param8.*t6.*t18.*(param8.*t6.*(t9.*t10.*t59.*u2.*4.0e2+t9.*t10.*t68.*u3.*4.0e2-t4.*t8.*t21.*t46.*t47.*u2.*1.6e4-t4.*t8.*t25.*t46.*t47.*u3.*1.6e4).*(1.0./2.0e1)+t4.*t8.*t101-t4.*t8.*t51.*u9+t2.*t4.*t8.*t16.*u9.*2.0e1-t2.*t16.*t44.*t45.*u9.*2.0e1).*(1.0./2.0e1)-t4.*t8.*t18.*t52.*u9+t2.*t4.*t8.*t17.*t18.*u9.*2.0e1-t2.*t17.*t18.*t44.*t45.*u9.*2.0e1).*4.0e2+param1.*t4.*t8.*t12.*t41.*t42.*t46.*t47.*t107.*1.6e4;one;t2.*t7.*u2.*4.0e1-t2.*t7.*t18.*u2.*2.0e1-t4.*t8.*t9.*t10.*t12.*u5.*(1.6e3./3.0);t109;t36+t37+t2.*t7.*t12.*t26.*(4.0e1./3.0)-t9.*t10.*t11.*t18.*4.0e2-t4.*t8.*t12.*t46.*t47.*u2.*u5.*(3.2e4./3.0)-t4.*t8.*t12.*t46.*t47.*u3.*u9.*8.0e3-param1.*t9.*t10.*t12.*t41.*t42.*(param8.*t6.*t18.*(param8.*t6.*(t9.*t10.*t24.*4.0e2-t4.*t8.*t46.*t47.*u2.*u5.*8.0e3).*(1.0./2.0e1)+t4.*t8.*t9.*t10.*u2.*u5.*4.0e2).*(1.0./2.0e1)-t2.*t4.*t7.*t8.*t18.*u2.*u5.*2.0e1).*4.0e2;zero;t109;t2.*t7.*t18.*u2.*-2.0e1+t4.*t8.*t9.*t10.*t12.*u5.*(8.0e2./3.0);t110-t9.*t10.*t18.*u2.*u3.*4.0e2+t4.*t8.*t12.*t46.*t47.*u3.*u5.*(1.6e4./3.0)-t4.*t8.*t12.*t46.*t47.*u2.*u9.*8.0e3-param1.*t9.*t10.*t12.*t41.*t42.*(param8.*t6.*t18.*(param8.*t6.*(t9.*t10.*t21.*4.0e2-t4.*t8.*t46.*t47.*u2.*u9.*8.0e3).*(1.0./2.0e1)+t4.*t8.*t9.*t10.*u2.*u9.*4.0e2).*(1.0./2.0e1)-t2.*t4.*t7.*t8.*t18.*u2.*u9.*2.0e1).*4.0e2;zero;t2.*t7.*t18.*u3.*-2.0e1+t4.*t8.*t9.*t10.*t12.*u9.*(8.0e2./3.0);t112;t110-t9.*t10.*t18.*u2.*u3.*4.0e2-t4.*t8.*t12.*t46.*t47.*u3.*u5.*8.0e3+t4.*t8.*t12.*t46.*t47.*u2.*u9.*(1.6e4./3.0)-param1.*t9.*t10.*t12.*t41.*t42.*(param8.*t6.*t18.*(param8.*t6.*(t9.*t10.*t19.*4.0e2-t4.*t8.*t46.*t47.*u3.*u5.*8.0e3).*(1.0./2.0e1)+t4.*t8.*t9.*t10.*u3.*u5.*4.0e2).*(1.0./2.0e1)-t2.*t4.*t7.*t8.*t18.*u3.*u5.*2.0e1).*4.0e2;one;t112;t2.*t7.*u3.*4.0e1-t2.*t7.*t18.*u3.*2.0e1-t4.*t8.*t9.*t10.*t12.*u9.*(1.6e3./3.0);t36+t37-t9.*t10.*t14.*t18.*4.0e2-t2.*t7.*t12.*t40.*(4.0e1./3.0)-t4.*t8.*t12.*t46.*t47.*u2.*u5.*8.0e3-t4.*t8.*t12.*t46.*t47.*u3.*u9.*(3.2e4./3.0)-param1.*t9.*t10.*t12.*t41.*t42.*(param8.*t6.*t18.*(param8.*t6.*(t9.*t10.*t25.*4.0e2-t4.*t8.*t46.*t47.*u3.*u9.*8.0e3).*(1.0./2.0e1)+t4.*t8.*t9.*t10.*u3.*u9.*4.0e2).*(1.0./2.0e1)-t2.*t4.*t7.*t8.*t18.*u3.*u9.*2.0e1).*4.0e2;zero;t113;zero;t116.*u2-param1.*t4.*t8.*t9.*t10.*t12.*t41.*u5.*4.0e2;zero;zero;t113;t116.*u3-param1.*t4.*t8.*t9.*t10.*t12.*t41.*u9.*4.0e2;zero;t4.*t8.*t9.*t10.*t12.*u2.*(-1.6e3./3.0);t4.*t8.*t9.*t10.*t12.*u3.*-4.0e2;t4.*t8.*t11.*t12.*t46.*t47.*(-3.2e4./3.0)-t4.*t8.*t12.*t14.*t46.*t47.*8.0e3-param1.*t9.*t10.*t12.*t41.*t42.*t120.*4.0e2;zero;t4.*t8.*t9.*t10.*t12.*u3.*-4.0e2;t4.*t8.*t9.*t10.*t12.*u2.*(8.0e2./3.0);t4.*t8.*t12.*t46.*t47.*u2.*u3.*(-8.0e3./3.0);zero;t123;zero;t9.*t10.*t12.*u2.*(1.6e3./3.0)-param1.*t9.*t10.*t12.*t41.*u2.*4.0e2;zero;zero;t2.*t7.*t12.*(-4.0e1./3.0);t9.*t10.*t12.*u3.*(-8.0e2./3.0);zero;zero;t117;t121-param1.*t9.*t10.*t12.*t41.*u3.*4.0e2;zero;t117;zero;t122;zero;zero;zero;t124;zero;zero;zero;zero;zero;t4.*t8.*t9.*t10.*t12.*u3.*(8.0e2./3.0);t4.*t8.*t9.*t10.*t12.*u2.*-4.0e2;t4.*t8.*t12.*t46.*t47.*u2.*u3.*(-8.0e3./3.0);zero;t4.*t8.*t9.*t10.*t12.*u2.*-4.0e2;t4.*t8.*t9.*t10.*t12.*u3.*(-1.6e3./3.0);t4.*t8.*t11.*t12.*t46.*t47.*-8.0e3-t4.*t8.*t12.*t14.*t46.*t47.*(3.2e4./3.0)-param1.*t9.*t10.*t12.*t41.*t42.*t120.*4.0e2;zero;zero;t117;t121;zero;t117;zero;t122-param1.*t9.*t10.*t12.*t41.*u2.*4.0e2;zero;t2.*t7.*t12.*(-4.0e1./3.0);zero;t9.*t10.*t12.*u2.*(-8.0e2./3.0);zero;zero;t123;t9.*t10.*t12.*u3.*(1.6e3./3.0)-param1.*t9.*t10.*t12.*t41.*u3.*4.0e2;zero;zero;zero;zero;zero;zero;zero;t124];
end
f = reshape(f,ng,nch,nd);
f_udg = reshape(f_udg,ng,nch,nd,nc);
