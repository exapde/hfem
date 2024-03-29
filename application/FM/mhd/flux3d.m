function [f,f_udg] = flux3d(pg,udg,param,time)
%FLUX3D
%    [F,F_UDG] = FLUX3D(PG,UDG,PARAM,TIME)

%    This function was generated by the Symbolic Math Toolbox version 6.3.
%    17-Nov-2017 09:38:18
[ng,nc] = size(udg);
nch = 9;
nd = 3;
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
zero = zeros(ng,1);
t2 = u2.^2;
t3 = 1.0./u1.^2;
t4 = u6.^2;
t5 = t4.*(1.0./2.0);
t6 = u7.^2;
t7 = t6.*(1.0./2.0);
t8 = u8.^2;
t9 = t8.*(1.0./2.0);
t10 = 1.0./u1;
t11 = param1-1.0;
t12 = t2.*t3.*(1.0./2.0);
t13 = u3.^2;
t14 = t3.*t13.*(1.0./2.0);
t15 = u4.^2;
t16 = t3.*t15.*(1.0./2.0);
t17 = t12+t14+t16;
t18 = t17.*u1;
t19 = t5+t7+t9+t18-u5;
t20 = t10.*u2.*u3;
t21 = t20-u6.*u7;
t22 = t5+t7+t9;
t23 = t10.*t22;
t24 = t10.*u5;
t36 = t10.*t11.*t19;
t25 = t23+t24-t36;
t26 = t10.*u2.*u6;
t27 = t10.*u3.*u7;
t28 = t10.*u4.*u8;
t29 = t26+t27+t28;
t30 = t10.*u3.*u6;
t31 = param2.^2;
t32 = t10.*u2.*u4;
t33 = t32-u6.*u8;
t34 = t10.*u3.*u4;
t35 = t34-u7.*u8;
t37 = t10.*u4.*u6;
t38 = t10.*u4.*u7;
t39 = 1.0./u1.^3;
t40 = t2.*t39;
t41 = t13.*t39;
t42 = t15.*t39;
t43 = t40+t41+t42;
t45 = t43.*u1;
t44 = t12+t14+t16-t45;
t46 = t3.*u2.*u6;
t47 = t3.*u3.*u7;
t48 = t3.*u4.*u8;
t49 = t46+t47+t48;
t50 = t3.*t22;
t51 = t3.*u5;
t52 = t10.*t11.*t44;
t55 = t3.*t11.*t19;
t53 = t50+t51+t52-t55;
t54 = t3.*u2.*u7;
t56 = t3.*u2.*u8;
t57 = t3.*u3.*u8;
t58 = t10.*u3;
t59 = t10.*u4;
t60 = -t10.*u6.*u7-t3.*t11.*u2.*u3;
t61 = t10.*u2;
t62 = t10.*u6;
t63 = t10.*u8;
t64 = -t10.*u6.*u8-t3.*t11.*u2.*u4;
t65 = -t10.*u7.*u8-t3.*t11.*u3.*u4;
t66 = t10.*u7;
t67 = one.*t11;
t68 = t10.*t11;
t69 = t10+t68;
t73 = t10.*t11.*u6;
t70 = t62-t73;
t71 = t10.*u2.*u7;
t72 = u6-t11.*u6;
t74 = t10.*u2.*u8;
t78 = t10.*t11.*u7;
t75 = t66-t78;
t76 = one.*t31;
t77 = u7-t11.*u7;
t79 = t10.*u3.*u8;
f = [u2;-t5+t7+t9+t2.*t10-t11.*t19;t21;t33;t25.*u2-t29.*u6;u9;t30-t10.*u2.*u7;t37-t10.*u2.*u8;t31.*u6;u3;t21;t5-t7+t9+t10.*t13-t11.*t19;t35;t25.*u3-t29.*u7;-t30+t71;u9;t38-t10.*u3.*u8;t31.*u7;u4;t33;t35;t5+t7-t9+t10.*t15-t11.*t19;t25.*u4-t29.*u8;-t37+t74;-t38+t79;u9;t31.*u8];
if nargout > 1
    t82 = t11.*u8;
    t80 = -t82+u8;
    t83 = t10.*t11.*u8;
    t81 = t63-t83;
    f_udg = [zero;-t2.*t3-t11.*t44;-t3.*u2.*u3;-t3.*u2.*u4;t49.*u6-t53.*u2;zero;t54-t3.*u3.*u6;t56-t3.*u4.*u6;zero;zero;-t3.*u2.*u3;-t3.*t13-t11.*t44;-t3.*u3.*u4;t49.*u7-t53.*u3;-t54+t3.*u3.*u6;zero;t57-t3.*u4.*u7;zero;zero;-t3.*u2.*u4;-t3.*u3.*u4;-t3.*t15-t11.*t44;t49.*u8-t53.*u4;-t56+t3.*u4.*u6;-t57+t3.*u4.*u7;zero;zero;one;t10.*u2.*2.0-t10.*t11.*u2;t58;t59;t23+t24-t36-t4.*t10-t2.*t3.*t11;zero;-t10.*u7;-t10.*u8;zero;zero;t58;-t10.*t11.*u2;zero;t60;t66;zero;zero;zero;zero;t59;zero;-t10.*t11.*u2;t64;t63;zero;zero;zero;zero;-t10.*t11.*u3;t61;zero;t60;zero;t62;zero;zero;one;t61;t10.*u3.*2.0-t10.*t11.*u3;t59;t23+t24-t36-t6.*t10-t3.*t11.*t13;-t62;zero;-t63;zero;zero;zero;t59;-t10.*t11.*u3;t65;zero;t63;zero;zero;zero;-t10.*t11.*u4;zero;t61;t64;zero;zero;t62;zero;zero;zero;-t10.*t11.*u4;t58;t65;zero;zero;t66;zero;one;t61;t58;t10.*u4.*2.0-t10.*t11.*u4;t23+t24-t36-t8.*t10-t3.*t11.*t15;-t62;-t66;zero;zero;zero;t67;zero;zero;t69.*u2;zero;zero;zero;zero;zero;zero;t67;zero;t69.*u3;zero;zero;zero;zero;zero;zero;zero;t67;t69.*u4;zero;zero;zero;zero;zero;-u6-t11.*u6;-u7;-u8;-t27-t28+t70.*u2-t10.*u2.*u6.*2.0;zero;t58;t59;t76;zero;-u7;t72;zero;-t71+t70.*u3;-t58;zero;zero;zero;zero;-u8;zero;t72;-t74+t70.*u4;-t59;zero;zero;zero;zero;t77;-u6;zero;-t30+t75.*u2;zero;-t61;zero;zero;zero;-u6;-u7-t11.*u7;-u8;-t26-t28+t75.*u3-t10.*u3.*u7.*2.0;t61;zero;t59;t76;zero;zero;-u8;t77;-t79+t75.*u4;zero;-t59;zero;zero;zero;t80;zero;-u6;-t37+t81.*u2;zero;zero;-t61;zero;zero;zero;t80;-u7;-t38+t81.*u3;zero;zero;-t58;zero;zero;-u6;-u7;-t82-u8;-t26-t27+t81.*u4-t10.*u4.*u8.*2.0;t61;t58;zero;t76;zero;zero;zero;zero;zero;one;zero;zero;zero;zero;zero;zero;zero;zero;zero;one;zero;zero;zero;zero;zero;zero;zero;zero;zero;one;zero];
end
f = reshape(f,ng,nch,nd);
f_udg = reshape(f_udg,ng,nch,nd,nc);
