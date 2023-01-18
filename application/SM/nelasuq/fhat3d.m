function [fh,fh_udg,fh_uh] = fhat3d(nl,pg,udg,uh,param,time)
%FHAT3D
%    [FH,FH_UDG,FH_UH] = FHAT3D(NL,PG,UDG,UH,PARAM,TIME)

%    This function was generated by the Symbolic Math Toolbox version 7.1.
%    18-Mar-2017 21:08:12
[ng,nc] = size(udg);
nch = 3;
nd = 3;
kappa = param{2};
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
zero = zeros(ng,1);
t2 = q21.*q32;
t30 = q22.*q31;
t3 = t2-t30;
t4 = q11.*q23.*q32;
t5 = q12.*q21.*q33;
t6 = q13.*q22.*q31;
t10 = q11.*q22.*q33;
t11 = q12.*q23.*q31;
t12 = q13.*q21.*q32;
t7 = t4+t5+t6-t10-t11-t12+1.0;
t8 = q21.*q33;
t32 = q23.*q31;
t9 = t8-t32;
t13 = q22.*q33;
t31 = q23.*q32;
t14 = t13-t31;
t15 = t4+t5+t6-t10-t11-t12;
t16 = 1.0./t15;
t17 = q11.*q32;
t35 = q12.*q31;
t18 = t17-t35;
t19 = q11.*q33;
t36 = q13.*q31;
t20 = t19-t36;
t21 = q12.*q33;
t37 = q13.*q32;
t22 = t21-t37;
t23 = q11.*q22;
t38 = q12.*q21;
t24 = t23-t38;
t25 = q11.*q23;
t39 = q13.*q21;
t26 = t25-t39;
t27 = q12.*q23;
t40 = q13.*q22;
t28 = t27-t40;
fh = [tau.*(u1-uh1)+nl3.*(mu.*q13+(mu.*t3)./(t4+t5+t6-q11.*q22.*q33-q12.*q23.*q31-q13.*q21.*q32)-kappa.*t3.*t7)+nl2.*(mu.*q12+kappa.*t7.*t9-mu.*t9.*t16)+nl1.*(mu.*q11-kappa.*t7.*t14+mu.*t14.*t16);tau.*(u2-uh2)+nl3.*(mu.*q23+kappa.*t7.*t18-mu.*t16.*t18)+nl2.*(mu.*q22-kappa.*t7.*t20+mu.*t16.*t20)+nl1.*(mu.*q21+kappa.*t7.*t22-mu.*t16.*t22);tau.*(u3-uh3)+nl3.*(mu.*q33-kappa.*t7.*t24+mu.*t16.*t24)+nl2.*(mu.*q32+kappa.*t7.*t26-mu.*t16.*t26)+nl1.*(mu.*q31-kappa.*t7.*t28+mu.*t16.*t28)];
if nargout > 1
    t29 = one.*tau;
    t33 = 1.0./(t4+t5+t6-t10-t11-t12).^2;
    t34 = t14.^2;
    t41 = mu.*q32.*t16;
    t42 = mu.*q33.*t16;
    t43 = kappa.*t14.*t22;
    t44 = mu.*t14.*t22.*t33;
    t45 = t43+t44;
    t46 = t22.^2;
    t47 = mu.*q22.*t16;
    t48 = mu.*q23.*t16;
    t49 = kappa.*t14.*t28;
    t50 = mu.*t14.*t28.*t33;
    t51 = t49+t50;
    t52 = nl1.*t51;
    t53 = mu.*q12.*t16;
    t54 = mu.*q13.*t16;
    t55 = kappa.*t22.*t28;
    t56 = mu.*t22.*t28.*t33;
    t57 = t55+t56;
    t58 = t28.^2;
    t59 = kappa.*t9.*t14;
    t60 = mu.*t9.*t14.*t33;
    t61 = t59+t60;
    t62 = t9.^2;
    t63 = kappa.*t9.*t22;
    t64 = kappa.*q33.*t7;
    t65 = mu.*t9.*t22.*t33;
    t66 = -t42+t63+t64+t65;
    t67 = kappa.*t9.*t28;
    t68 = kappa.*q23.*t7;
    t69 = mu.*t9.*t28.*t33;
    t70 = -t48+t67+t68+t69;
    t71 = mu.*q31.*t16;
    t72 = kappa.*t14.*t20;
    t73 = mu.*t14.*t20.*t33;
    t74 = kappa.*t9.*t20;
    t75 = mu.*t9.*t20.*t33;
    t76 = t74+t75;
    t77 = kappa.*t20.*t22;
    t78 = mu.*t20.*t22.*t33;
    t79 = t77+t78;
    t80 = t20.^2;
    t81 = kappa.*t20.*t28;
    t82 = kappa.*q13.*t7;
    t83 = mu.*t20.*t28.*t33;
    t84 = -t54+t81+t82+t83;
    t85 = mu.*q21.*t16;
    t86 = kappa.*t14.*t26;
    t87 = mu.*t14.*t26.*t33;
    t88 = kappa.*t9.*t26;
    t89 = mu.*t9.*t26.*t33;
    t90 = t88+t89;
    t91 = nl2.*t90;
    t92 = mu.*q11.*t16;
    t93 = kappa.*t22.*t26;
    t94 = mu.*t22.*t26.*t33;
    t95 = kappa.*t20.*t26;
    t96 = mu.*t20.*t26.*t33;
    t97 = t95+t96;
    t98 = kappa.*t26.*t28;
    t99 = mu.*t26.*t28.*t33;
    t100 = t98+t99;
    t101 = t26.^2;
    t102 = kappa.*t3.*t9;
    t103 = mu.*t3.*t9.*t33;
    t104 = t102+t103;
    t105 = kappa.*t3.*t14;
    t106 = mu.*t3.*t14.*t33;
    t107 = t105+t106;
    t108 = t3.^2;
    t109 = kappa.*t3.*t20;
    t110 = kappa.*q31.*t7;
    t111 = mu.*t3.*t20.*t33;
    t112 = -t71+t109+t110+t111;
    t113 = kappa.*t3.*t22;
    t114 = kappa.*q32.*t7;
    t115 = mu.*t3.*t22.*t33;
    t116 = -t41+t113+t114+t115;
    t117 = kappa.*t3.*t26;
    t118 = kappa.*q21.*t7;
    t119 = mu.*t3.*t26.*t33;
    t120 = -t85+t117+t118+t119;
    t121 = kappa.*t3.*t28;
    t122 = kappa.*q22.*t7;
    t123 = mu.*t3.*t28.*t33;
    t124 = -t47+t121+t122+t123;
    t125 = kappa.*t9.*t18;
    t126 = mu.*t9.*t18.*t33;
    t127 = kappa.*t14.*t18;
    t128 = mu.*t14.*t18.*t33;
    t129 = kappa.*t3.*t18;
    t130 = mu.*t3.*t18.*t33;
    t131 = t129+t130;
    t132 = kappa.*t18.*t20;
    t133 = mu.*t18.*t20.*t33;
    t134 = t132+t133;
    t135 = kappa.*t18.*t22;
    t136 = mu.*t18.*t22.*t33;
    t137 = t135+t136;
    t138 = t18.^2;
    t139 = kappa.*t18.*t26;
    t140 = kappa.*q11.*t7;
    t141 = mu.*t18.*t26.*t33;
    t142 = -t92+t139+t140+t141;
    t143 = kappa.*t18.*t28;
    t144 = kappa.*q12.*t7;
    t145 = mu.*t18.*t28.*t33;
    t146 = -t53+t143+t144+t145;
    t147 = kappa.*t9.*t24;
    t148 = mu.*t9.*t24.*t33;
    t149 = kappa.*t14.*t24;
    t150 = mu.*t14.*t24.*t33;
    t151 = kappa.*t3.*t24;
    t152 = mu.*t3.*t24.*t33;
    t153 = t151+t152;
    t154 = nl3.*t153;
    t155 = kappa.*t20.*t24;
    t156 = mu.*t20.*t24.*t33;
    t157 = kappa.*t22.*t24;
    t158 = mu.*t22.*t24.*t33;
    t159 = kappa.*t18.*t24;
    t160 = mu.*t18.*t24.*t33;
    t161 = t159+t160;
    t162 = kappa.*t24.*t26;
    t163 = mu.*t24.*t26.*t33;
    t164 = t162+t163;
    t165 = kappa.*t24.*t28;
    t166 = mu.*t24.*t28.*t33;
    t167 = t165+t166;
    t168 = t24.^2;
    fh_udg = [t29;zero;zero;zero;t29;zero;zero;zero;t29;-nl2.*t61+nl3.*t107+nl1.*(mu+kappa.*t34+mu.*t33.*t34);nl2.*(t42+t72+t73-kappa.*q33.*t7)-nl3.*(t41+t127+t128-kappa.*q32.*t7)-nl1.*t45;t52-nl2.*(t48+t86+t87-kappa.*q23.*t7)+nl3.*(t47+t149+t150-kappa.*q22.*t7);-nl1.*t45+nl2.*t66-nl3.*t116;-nl2.*t79+nl3.*t137+nl1.*(mu+kappa.*t46+mu.*t33.*t46);nl2.*(t54+t93+t94-kappa.*q13.*t7)-nl3.*(t53+t157+t158-kappa.*q12.*t7)-nl1.*t57;t52-nl2.*t70+nl3.*t124;-nl1.*t57+nl2.*t84-nl3.*t146;-nl2.*t100+nl3.*t167+nl1.*(mu+kappa.*t58+mu.*t33.*t58);-nl1.*t61-nl3.*t104+nl2.*(mu+kappa.*t62+mu.*t33.*t62);nl3.*(t71+t125+t126-kappa.*q31.*t7)+nl1.*t66-nl2.*t76;t91-nl3.*(t85+t147+t148-kappa.*q21.*t7)-nl1.*t70;-nl2.*t76+nl3.*t112+nl1.*(t42-t64+t72+t73);-nl1.*t79-nl3.*t134+nl2.*(mu+kappa.*t80+mu.*t33.*t80);nl3.*(t92+t155+t156-kappa.*q11.*t7)+nl1.*t84-nl2.*t97;t91-nl3.*t120-nl1.*(t48-t68+t86+t87);-nl2.*t97+nl3.*t142+nl1.*(t54-t82+t93+t94);-nl1.*t100-nl3.*t164+nl2.*(mu+kappa.*t101+mu.*t33.*t101);-nl2.*t104+nl1.*t107+nl3.*(mu+kappa.*t108+mu.*t33.*t108);nl2.*t112-nl1.*t116-nl3.*t131;t154-nl2.*t120+nl1.*t124;-nl3.*t131-nl1.*(t41-t114+t127+t128)+nl2.*(t71-t110+t125+t126);-nl2.*t134+nl1.*t137+nl3.*(mu+kappa.*t138+mu.*t33.*t138);nl2.*t142-nl1.*t146-nl3.*t161;t154+nl1.*(t47-t122+t149+t150)-nl2.*(t85-t118+t147+t148);-nl3.*t161-nl1.*(t53-t144+t157+t158)+nl2.*(t92-t140+t155+t156);-nl2.*t164+nl1.*t167+nl3.*(mu+kappa.*t168+mu.*t33.*t168)];
end
if nargout > 2
    fh_uh = [-t29;zero;zero;zero;-t29;zero;zero;zero;-t29];
end
fh = reshape(fh,ng,nch);
fh_udg = reshape(fh_udg,ng,nch,nc);
fh_uh = reshape(fh_uh,ng,nch,nch);
