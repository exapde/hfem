function [f,f_udg] = flux2did(pg,udg,param,time)
%FLUX2DID
%    [F,F_UDG] = FLUX2DID(PG,UDG,PARAM,TIME)

%    This function was generated by the Symbolic Math Toolbox version 6.3.
%    06-Dec-2017 14:46:32
[ng,nc] = size(udg);
nch = 8;
nd = 2;
one = ones(ng,1);
param1 = param{1};
u1 = udg(:,1);
u2 = udg(:,2);
u3 = udg(:,3);
u4 = udg(:,4);
u5 = udg(:,5);
u6 = udg(:,6);
u7 = udg(:,7);
u8 = udg(:,8);
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
t47 = t10.*t11.*t19;
t25 = t23+t24-t47;
t26 = t10.*u2.*u6;
t27 = t10.*u3.*u7;
t28 = t10.*u4.*u8;
t29 = t26+t27+t28;
t30 = t10.*u2.*u7;
t31 = 1.0./u1.^3;
t32 = t2.*t31;
t33 = t13.*t31;
t34 = t15.*t31;
t35 = t32+t33+t34;
t37 = t35.*u1;
t36 = t12+t14+t16-t37;
t38 = t3.*u2.*u6;
t39 = t3.*u3.*u7;
t40 = t3.*u4.*u8;
t41 = t38+t39+t40;
t42 = t3.*t22;
t43 = t3.*u5;
t44 = t10.*t11.*t36;
t45 = t42+t43+t44-t3.*t11.*t19;
t46 = t3.*u3.*u6;
t48 = t10.*u3;
t49 = t10.*u7;
t50 = -t10.*u6.*u7-t3.*t11.*u2.*u3;
t51 = t10.*u2;
t52 = t10.*u4;
t53 = t10.*u8;
t54 = t10.*u6;
t55 = one.*t11;
t56 = t10.*t11;
t57 = t10+t56;
t58 = t54-t10.*t11.*u6;
t59 = t10.*u3.*u6;
f = [u2;-t5+t7+t9+t2.*t10-t11.*t19;t21;-u6.*u8+t10.*u2.*u4;t25.*u2-t29.*u6;zero;t30-t10.*u3.*u6;t10.*u2.*u8-t10.*u4.*u6;u3;t21;t5-t7+t9+t10.*t13-t11.*t19;-u7.*u8+t10.*u3.*u4;t25.*u3-t29.*u7;-t30+t59;zero;t10.*u3.*u8-t10.*u4.*u7];
if nargout > 1
    t60 = t49-t10.*t11.*u7;
    t61 = u8-t11.*u8;
    t62 = t53-t10.*t11.*u8;
    f_udg = [zero;-t2.*t3-t11.*t36;-t3.*u2.*u3;-t3.*u2.*u4;t41.*u6-t45.*u2;zero;t46-t3.*u2.*u7;-t3.*u2.*u8+t3.*u4.*u6;zero;-t3.*u2.*u3;-t3.*t13-t11.*t36;-t3.*u3.*u4;t41.*u7-t45.*u3;-t46+t3.*u2.*u7;zero;-t3.*u3.*u8+t3.*u4.*u7;one;t10.*u2.*2.0-t10.*t11.*u2;t48;t52;t23+t24-t47-t4.*t10-t2.*t3.*t11;zero;t49;t53;zero;t48;-t10.*t11.*u2;zero;t50;-t49;zero;zero;zero;-t10.*t11.*u3;t51;zero;t50;zero;-t10.*u6;zero;one;t51;t10.*u3.*2.0-t10.*t11.*u3;t52;t23+t24-t47-t6.*t10-t3.*t11.*t13;t54;zero;t53;zero;-t10.*t11.*u4;zero;t51;-t10.*u6.*u8-t3.*t11.*u2.*u4;zero;zero;-t54;zero;zero;-t10.*t11.*u4;t48;-t10.*u7.*u8-t3.*t11.*u3.*u4;zero;zero;-t49;zero;t55;zero;zero;t57.*u2;zero;zero;zero;zero;zero;t55;zero;t57.*u3;zero;zero;zero;zero;-u6-t11.*u6;-u7;-u8;-t27-t28+t58.*u2-t10.*u2.*u6.*2.0;zero;-t48;-t52;zero;-u7;u6-t11.*u6;zero;-t30+t58.*u3;t48;zero;zero;zero;u7-t11.*u7;-u6;zero;-t59+t60.*u2;zero;t51;zero;zero;-u6;-u7-t11.*u7;-u8;-t26-t28+t60.*u3-t10.*u3.*u7.*2.0;-t51;zero;-t52;zero;t61;zero;-u6;t62.*u2-t10.*u4.*u6;zero;zero;t51;zero;zero;t61;-u7;t62.*u3-t10.*u4.*u7;zero;zero;t48];
end
f = reshape(f,ng,nch,nd);
f_udg = reshape(f_udg,ng,nch,nd,nc);
