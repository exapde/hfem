function [fh,fh_udg,fh_uh] = fhat3d(nl,pg,udg,uh,param,time)
%FHAT3D
%    [FH,FH_UDG,FH_UH] = FHAT3D(NL,PG,UDG,UH,PARAM,TIME)

%    This function was generated by the Symbolic Math Toolbox version 6.3.
%    17-Nov-2017 09:38:32
[ng,nc] = size(udg);
nch = 9;
nd = 3;
nl1 = nl(:,1);
nl2 = nl(:,2);
nl3 = nl(:,3);
one = ones(ng,1);
param1 = param{1};
param2 = param{2};
u1 = udg(:,1);
u2 = udg(:,2);
u3 = udg(:,3);
u4 = udg(:,4);
u5 = udg(:,5);
u6 = udg(:,6);
u7 = udg(:,7);
u8 = udg(:,8);
u9 = udg(:,9);
uh1 = uh(:,1);
uh2 = uh(:,2);
uh3 = uh(:,3);
uh4 = uh(:,4);
uh5 = uh(:,5);
uh6 = uh(:,6);
uh7 = uh(:,7);
uh8 = uh(:,8);
uh9 = uh(:,9);
zero = zeros(ng,1);
t2 = uh2.^2;
t3 = uh3.^2;
t4 = uh4.^2;
t5 = uh7.^2;
t6 = uh8.^2;
t7 = 1.0./uh1;
t8 = uh6.^2;
t9 = 1.0./uh1.^2;
t10 = t7.*uh2.*uh7;
t11 = t10-t7.*uh3.*uh6;
t12 = t7.*uh2.*uh8;
t13 = t12-t7.*uh4.*uh6;
t14 = t7.*uh3.*uh8;
t15 = t14-t7.*uh4.*uh7;
t16 = param2.^2;
t17 = one.*param2;
t18 = nl1.*param1.*t2;
t19 = nl1.*param1.*t3;
t20 = nl1.*param1.*t4;
t21 = nl1.*uh1.*uh5.*2.0;
t22 = param2.*uh1.*uh2.*2.0;
t23 = nl2.*uh1.*uh6.*uh7.*2.0;
t24 = nl3.*uh1.*uh6.*uh8.*2.0;
t25 = nl1.*param1.*t8.*uh1;
t26 = nl1.*param1.*t5.*uh1;
t27 = nl1.*param1.*t6.*uh1;
t28 = t18+t19+t20+t21+t22+t23+t24+t25+t26+t27-nl1.*t2.*3.0-nl1.*t3-nl1.*t4-nl1.*t5.*uh1.*2.0-nl1.*t6.*uh1.*2.0-nl2.*uh2.*uh3.*2.0-nl3.*uh2.*uh4.*2.0-param2.*u2.*uh1.*2.0-nl1.*param1.*uh1.*uh5.*2.0;
t29 = nl2.*param1.*t2;
t30 = nl2.*param1.*t3;
t31 = nl2.*param1.*t4;
t32 = nl2.*uh1.*uh5.*2.0;
t33 = param2.*uh1.*uh3.*2.0;
t34 = nl1.*uh1.*uh6.*uh7.*2.0;
t35 = nl3.*uh1.*uh7.*uh8.*2.0;
t36 = nl2.*param1.*t8.*uh1;
t37 = nl2.*param1.*t5.*uh1;
t38 = nl2.*param1.*t6.*uh1;
t39 = t29+t30+t31+t32+t33+t34+t35+t36+t37+t38-nl2.*t2-nl2.*t3.*3.0-nl2.*t4-nl2.*t6.*uh1.*2.0-nl2.*t8.*uh1.*2.0-nl1.*uh2.*uh3.*2.0-nl3.*uh3.*uh4.*2.0-param2.*u3.*uh1.*2.0-nl2.*param1.*uh1.*uh5.*2.0;
t40 = nl3.*param1.*t2;
t41 = nl3.*param1.*t3;
t42 = nl3.*param1.*t4;
t43 = nl3.*uh1.*uh5.*2.0;
t44 = param2.*uh1.*uh4.*2.0;
t45 = nl1.*uh1.*uh6.*uh8.*2.0;
t46 = nl2.*uh1.*uh7.*uh8.*2.0;
t47 = nl3.*param1.*t8.*uh1;
t48 = nl3.*param1.*t5.*uh1;
t49 = nl3.*param1.*t6.*uh1;
t50 = t40+t41+t42+t43+t44+t45+t46+t47+t48+t49-nl3.*t2-nl3.*t3-nl3.*t4.*3.0-nl3.*t5.*uh1.*2.0-nl3.*t8.*uh1.*2.0-nl1.*uh2.*uh4.*2.0-nl2.*uh3.*uh4.*2.0-param2.*u4.*uh1.*2.0-nl3.*param1.*uh1.*uh5.*2.0;
t51 = param1.*t2.*uh2;
t52 = param1.*t3.*uh2;
t53 = param1.*t4.*uh2;
t54 = uh1.*uh3.*uh6.*uh7.*2.0;
t55 = uh1.*uh4.*uh6.*uh8.*2.0;
t56 = param1.*t8.*uh1.*uh2;
t57 = param1.*t5.*uh1.*uh2;
t58 = param1.*t6.*uh1.*uh2;
t59 = t51+t52+t53+t54+t55+t56+t57+t58-t2.*uh2-t3.*uh2-t4.*uh2-t5.*uh1.*uh2.*2.0-t6.*uh1.*uh2.*2.0-param1.*uh1.*uh2.*uh5.*2.0;
t60 = 1.0./uh1.^3;
t61 = param1.*t3.*uh3;
t62 = param1.*t2.*uh3;
t63 = param1.*t4.*uh3;
t64 = uh1.*uh2.*uh6.*uh7.*2.0;
t65 = uh1.*uh4.*uh7.*uh8.*2.0;
t66 = param1.*t8.*uh1.*uh3;
t67 = param1.*t5.*uh1.*uh3;
t68 = param1.*t6.*uh1.*uh3;
t69 = t61+t62+t63+t64+t65+t66+t67+t68-t2.*uh3-t3.*uh3-t4.*uh3-t6.*uh1.*uh3.*2.0-t8.*uh1.*uh3.*2.0-param1.*uh1.*uh3.*uh5.*2.0;
t70 = param1.*t4.*uh4;
t71 = param1.*t2.*uh4;
t72 = param1.*t3.*uh4;
t73 = uh1.*uh2.*uh6.*uh8.*2.0;
t74 = uh1.*uh3.*uh7.*uh8.*2.0;
t75 = param1.*t8.*uh1.*uh4;
t76 = param1.*t5.*uh1.*uh4;
t77 = param1.*t6.*uh1.*uh4;
t78 = t70+t71+t72+t73+t74+t75+t76+t77-t2.*uh4-t3.*uh4-t4.*uh4-t5.*uh1.*uh4.*2.0-t8.*uh1.*uh4.*2.0-param1.*uh1.*uh4.*uh5.*2.0;
fh = [nl1.*uh2+nl2.*uh3+nl3.*uh4+param2.*(u1-uh1);t7.*t28.*(-1.0./2.0);t7.*t39.*(-1.0./2.0);t7.*t50.*(-1.0./2.0);param2.*(u5-uh5)-nl1.*t9.*t59.*(1.0./2.0)-nl2.*t9.*t69.*(1.0./2.0)-nl3.*t9.*t78.*(1.0./2.0);nl2.*t11+nl3.*t13+nl1.*uh9+param2.*(u6-uh6);-nl1.*t11+nl3.*t15+nl2.*uh9+param2.*(u7-uh7);-nl1.*t13-nl2.*t15+nl3.*uh9+param2.*(u8-uh8);param2.*(u9-uh9)+nl1.*t16.*uh6+nl2.*t16.*uh7+nl3.*t16.*uh8];
if nargout > 1
    fh_udg = [t17;zero;zero;zero;zero;zero;zero;zero;zero;zero;t17;zero;zero;zero;zero;zero;zero;zero;zero;zero;t17;zero;zero;zero;zero;zero;zero;zero;zero;zero;t17;zero;zero;zero;zero;zero;zero;zero;zero;zero;t17;zero;zero;zero;zero;zero;zero;zero;zero;zero;t17;zero;zero;zero;zero;zero;zero;zero;zero;zero;t17;zero;zero;zero;zero;zero;zero;zero;zero;zero;t17;zero;zero;zero;zero;zero;zero;zero;zero;zero;t17];
end
if nargout > 2
    t79 = t9.*uh2.*uh7;
    t80 = t79-t9.*uh3.*uh6;
    t81 = t9.*uh2.*uh8;
    t82 = t81-t9.*uh4.*uh6;
    t83 = t9.*uh3.*uh8;
    t84 = t83-t9.*uh4.*uh7;
    t85 = nl1.*uh3.*2.0;
    t86 = nl2.*uh2.*2.0;
    t87 = nl3.*uh4.*2.0;
    t88 = param1.*uh2.*uh3.*2.0;
    t89 = uh1.*uh6.*uh7.*2.0;
    t90 = t88+t89-uh2.*uh3.*2.0;
    t91 = t6.*uh1.*2.0;
    t92 = param1.*uh1.*uh5.*2.0;
    t93 = nl3.*t7.*uh8;
    t94 = nl1.*uh4.*2.0;
    t95 = nl3.*uh2.*2.0;
    t96 = nl2.*uh4.*2.0;
    t97 = nl3.*uh3.*2.0;
    t98 = nl1.*uh2.*2.0;
    t99 = nl2.*uh3.*2.0;
    t100 = param1.*uh2.*uh4.*2.0;
    t101 = uh1.*uh6.*uh8.*2.0;
    t102 = t100+t101-uh2.*uh4.*2.0;
    t103 = param1.*uh3.*uh4.*2.0;
    t104 = uh1.*uh7.*uh8.*2.0;
    t105 = t103+t104-uh3.*uh4.*2.0;
    t106 = t8.*uh1.*2.0;
    t107 = t5.*uh1.*2.0;
    t108 = nl1.*t7.*uh6;
    t109 = nl2.*t7.*uh7;
    t110 = nl3.*uh1.*uh8.*2.0;
    t111 = uh1.*uh4.*uh8.*2.0;
    t112 = nl3.*t7.*uh4;
    t113 = nl1.*uh1.*uh6.*2.0;
    t114 = nl2.*uh1.*uh7.*2.0;
    t115 = uh1.*uh2.*uh6.*2.0;
    t116 = uh1.*uh3.*uh7.*2.0;
    t117 = nl1.*t7.*uh2;
    t118 = nl2.*t7.*uh3;
    t119 = nl1.*one;
    t120 = nl2.*one;
    t121 = nl3.*one;
    fh_uh = [-t17;t9.*t28.*(1.0./2.0)-t7.*(nl1.*t5.*-2.0-nl1.*t6.*2.0+nl1.*uh5.*2.0-param2.*u2.*2.0+param2.*uh2.*2.0+nl1.*param1.*t5+nl1.*param1.*t6+nl1.*param1.*t8-nl1.*param1.*uh5.*2.0+nl2.*uh6.*uh7.*2.0+nl3.*uh6.*uh8.*2.0).*(1.0./2.0);t9.*t39.*(1.0./2.0)-t7.*(nl2.*t6.*-2.0-nl2.*t8.*2.0+nl2.*uh5.*2.0-param2.*u3.*2.0+param2.*uh3.*2.0+nl2.*param1.*t5+nl2.*param1.*t6+nl2.*param1.*t8-nl2.*param1.*uh5.*2.0+nl1.*uh6.*uh7.*2.0+nl3.*uh7.*uh8.*2.0).*(1.0./2.0);t9.*t50.*(1.0./2.0)-t7.*(nl3.*t5.*-2.0-nl3.*t8.*2.0+nl3.*uh5.*2.0-param2.*u4.*2.0+param2.*uh4.*2.0+nl3.*param1.*t5+nl3.*param1.*t6+nl3.*param1.*t8-nl3.*param1.*uh5.*2.0+nl1.*uh6.*uh8.*2.0+nl2.*uh7.*uh8.*2.0).*(1.0./2.0);-one.*(nl1.*t9.*(t5.*uh2.*-2.0-t6.*uh2.*2.0+param1.*t5.*uh2+param1.*t6.*uh2+param1.*t8.*uh2-param1.*uh2.*uh5.*2.0+uh3.*uh6.*uh7.*2.0+uh4.*uh6.*uh8.*2.0).*(1.0./2.0)+nl2.*t9.*(t6.*uh3.*-2.0-t8.*uh3.*2.0+param1.*t5.*uh3+param1.*t6.*uh3+param1.*t8.*uh3-param1.*uh3.*uh5.*2.0+uh2.*uh6.*uh7.*2.0+uh4.*uh7.*uh8.*2.0).*(1.0./2.0)+nl3.*t9.*(t5.*uh4.*-2.0-t8.*uh4.*2.0+param1.*t5.*uh4+param1.*t6.*uh4+param1.*t8.*uh4-param1.*uh4.*uh5.*2.0+uh2.*uh6.*uh8.*2.0+uh3.*uh7.*uh8.*2.0).*(1.0./2.0)-nl1.*t59.*t60-nl2.*t60.*t69-nl3.*t60.*t78);-one.*(nl2.*t80+nl3.*t82);one.*(nl1.*t80-nl3.*t84);one.*(nl1.*t82+nl2.*t84);zero;t119;one.*t7.*(t87+t99+nl1.*uh2.*6.0-param2.*uh1.*2.0-nl1.*param1.*uh2.*2.0).*(1.0./2.0);one.*t7.*(t85+t86-nl2.*param1.*uh2.*2.0).*(1.0./2.0);one.*t7.*(t94+t95-nl3.*param1.*uh2.*2.0).*(1.0./2.0);-one.*(nl1.*t9.*(t2.*3.0+t3+t4+t91+t92+t107-param1.*t2.*3.0-param1.*t3-param1.*t4-param1.*t5.*uh1-param1.*t6.*uh1-param1.*t8.*uh1).*(-1.0./2.0)+nl2.*t9.*t90.*(1.0./2.0)+nl3.*t9.*t102.*(1.0./2.0));one.*(t93+t109);-nl1.*one.*t7.*uh7;-nl1.*one.*t7.*uh8;zero;t120;one.*t7.*(t85+t86-nl1.*param1.*uh3.*2.0).*(1.0./2.0);one.*t7.*(t87+t98+nl2.*uh3.*6.0-param2.*uh1.*2.0-nl2.*param1.*uh3.*2.0).*(1.0./2.0);one.*t7.*(t96+t97-nl3.*param1.*uh3.*2.0).*(1.0./2.0);-one.*(nl2.*t9.*(t2+t3.*3.0+t4+t91+t92+t106-param1.*t2-param1.*t3.*3.0-param1.*t4-param1.*t5.*uh1-param1.*t6.*uh1-param1.*t8.*uh1).*(-1.0./2.0)+nl1.*t9.*t90.*(1.0./2.0)+nl3.*t9.*t105.*(1.0./2.0));-nl2.*one.*t7.*uh6;one.*(t93+t108);-nl2.*one.*t7.*uh8;zero;t121;one.*t7.*(t94+t95-nl1.*param1.*uh4.*2.0).*(1.0./2.0);one.*t7.*(t96+t97-nl2.*param1.*uh4.*2.0).*(1.0./2.0);one.*t7.*(t98+t99+nl3.*uh4.*6.0-param2.*uh1.*2.0-nl3.*param1.*uh4.*2.0).*(1.0./2.0);-one.*(nl3.*t9.*(t2+t3+t4.*3.0+t92+t106+t107-param1.*t2-param1.*t3-param1.*t4.*3.0-param1.*t5.*uh1-param1.*t6.*uh1-param1.*t8.*uh1).*(-1.0./2.0)+nl1.*t9.*t102.*(1.0./2.0)+nl2.*t9.*t105.*(1.0./2.0));-nl3.*one.*t7.*uh6;-nl3.*one.*t7.*uh7;one.*(t108+t109);zero;zero;one.*t7.*(nl1.*uh1.*2.0-nl1.*param1.*uh1.*2.0).*(-1.0./2.0);one.*t7.*(nl2.*uh1.*2.0-nl2.*param1.*uh1.*2.0).*(-1.0./2.0);one.*t7.*(nl3.*uh1.*2.0-nl3.*param1.*uh1.*2.0).*(-1.0./2.0);one.*(-param2+nl1.*param1.*t7.*uh2+nl2.*param1.*t7.*uh3+nl3.*param1.*t7.*uh4);zero;zero;zero;zero;zero;one.*t7.*(t110+t114+nl1.*param1.*uh1.*uh6.*2.0).*(-1.0./2.0);one.*t7.*(nl1.*uh1.*uh7.*2.0-nl2.*uh1.*uh6.*4.0+nl2.*param1.*uh1.*uh6.*2.0).*(-1.0./2.0);one.*t7.*(nl1.*uh1.*uh8.*2.0-nl3.*uh1.*uh6.*4.0+nl3.*param1.*uh1.*uh6.*2.0).*(-1.0./2.0);-one.*(nl2.*t9.*(uh1.*uh2.*uh7.*2.0-uh1.*uh3.*uh6.*4.0+param1.*uh1.*uh3.*uh6.*2.0).*(1.0./2.0)+nl3.*t9.*(uh1.*uh2.*uh8.*2.0-uh1.*uh4.*uh6.*4.0+param1.*uh1.*uh4.*uh6.*2.0).*(1.0./2.0)+nl1.*t9.*(t111+t116+param1.*uh1.*uh2.*uh6.*2.0).*(1.0./2.0));-one.*(param2+t112+t118);nl1.*one.*t7.*uh3;nl1.*one.*t7.*uh4;nl1.*one.*t16;zero;one.*t7.*(nl1.*uh1.*uh7.*-4.0+nl2.*uh1.*uh6.*2.0+nl1.*param1.*uh1.*uh7.*2.0).*(-1.0./2.0);one.*t7.*(t110+t113+nl2.*param1.*uh1.*uh7.*2.0).*(-1.0./2.0);one.*t7.*(nl2.*uh1.*uh8.*2.0-nl3.*uh1.*uh7.*4.0+nl3.*param1.*uh1.*uh7.*2.0).*(-1.0./2.0);-one.*(nl1.*t9.*(uh1.*uh2.*uh7.*-4.0+uh1.*uh3.*uh6.*2.0+param1.*uh1.*uh2.*uh7.*2.0).*(1.0./2.0)+nl3.*t9.*(uh1.*uh3.*uh8.*2.0-uh1.*uh4.*uh7.*4.0+param1.*uh1.*uh4.*uh7.*2.0).*(1.0./2.0)+nl2.*t9.*(t111+t115+param1.*uh1.*uh3.*uh7.*2.0).*(1.0./2.0));nl2.*one.*t7.*uh2;-one.*(param2+t112+t117);nl2.*one.*t7.*uh4;nl2.*one.*t16;zero;one.*t7.*(nl1.*uh1.*uh8.*-4.0+nl3.*uh1.*uh6.*2.0+nl1.*param1.*uh1.*uh8.*2.0).*(-1.0./2.0);one.*t7.*(nl2.*uh1.*uh8.*-4.0+nl3.*uh1.*uh7.*2.0+nl2.*param1.*uh1.*uh8.*2.0).*(-1.0./2.0);one.*t7.*(t113+t114+nl3.*param1.*uh1.*uh8.*2.0).*(-1.0./2.0);-one.*(nl1.*t9.*(uh1.*uh2.*uh8.*-4.0+uh1.*uh4.*uh6.*2.0+param1.*uh1.*uh2.*uh8.*2.0).*(1.0./2.0)+nl2.*t9.*(uh1.*uh3.*uh8.*-4.0+uh1.*uh4.*uh7.*2.0+param1.*uh1.*uh3.*uh8.*2.0).*(1.0./2.0)+nl3.*t9.*(t115+t116+param1.*uh1.*uh4.*uh8.*2.0).*(1.0./2.0));nl3.*one.*t7.*uh2;nl3.*one.*t7.*uh3;-one.*(param2+t117+t118);nl3.*one.*t16;zero;zero;zero;zero;zero;t119;t120;t121;-t17];
end
fh = reshape(fh,ng,nch);
fh_udg = reshape(fh_udg,ng,nch,nc);
fh_uh = reshape(fh_uh,ng,nch,nch);
