function [f,f_udg] = sensor2d2(pg,udg,param,time)
%SENSOR2D2
%    [F,F_UDG] = SENSOR2D2(PG,UDG,PARAM,TIME)

%    This function was generated by the Symbolic Math Toolbox version 5.8.
%    27-Mar-2013 09:10:51
[ng,nc] = size(udg);
nch = 4;
nd = 2;
param1 = param{1};
param8 = param{8};
u1 = udg(:,1);
u2 = udg(:,2);
u3 = udg(:,3);
u4 = udg(:,4);
u5 = udg(:,5);
u6 = udg(:,6);
u9 = udg(:,9);
u11 = udg(:,11);
zero = zeros(ng,1);
t2 = 1.0./param8;
t3 = t2.*u1.*2.0e1;
t4 = exp(t3);
t5 = t4+2.0e1;
t6 = log(t5);
t7 = 1.0./t6;
t8 = 1.0./param8.^2;
t9 = 1.0./t6.^2;
t10 = u2.^2;
t11 = u3.^2;
t12 = param1-1.0;
t13 = t8.*t9.*t10.*4.0e2;
t14 = t8.*t9.*t11.*4.0e2;
t15 = t13+t14;
t16 = t12.*t15;
t17 = t8.*t9.*t10.*2.0e2;
t18 = t8.*t9.*t11.*2.0e2;
t19 = t17+t18;
t37 = param8.*t6.*t19.*(1.0./2.0e1);
t20 = -t37+u4;
t21 = param1.*t2.*t7.*t12.*t20.*4.0e1;
t22 = t16+t21;
t23 = param1+1.0;
t24 = 1.0./t23;
t25 = t22.*t24;
t26 = 1.0./sqrt(t25);
t32 = t2.*t7.*u2.*u5.*2.0e1;
t27 = -t32+u6;
t34 = t2.*t7.*u3.*u9.*2.0e1;
t28 = -t34+u11;
t29 = 1.0./t5;
t30 = 1.0./param8.^3;
t31 = 1.0./t6.^3;
t33 = t2.*t7.*t27.*2.0e1;
t35 = t2.*t7.*t28.*2.0e1;
t36 = t33+t35;
f = t26.*t36.*1.0e1-5.0;
if nargout > 1
    t38 = 1.0./t25.^(3.0./2.0);
    t39 = t2.*t7.*t26.*2.0e2;
    f_udg = [t26.*(t4.*t8.*t9.*t27.*t29.*4.0e2+t4.*t8.*t9.*t28.*t29.*4.0e2-t4.*t29.*t30.*t31.*u2.*u5.*8.0e3-t4.*t29.*t30.*t31.*u3.*u9.*8.0e3).*-1.0e1+t24.*t36.*t38.*(t12.*(t4.*t10.*t29.*t30.*t31.*1.6e4+t4.*t11.*t29.*t30.*t31.*1.6e4)+param1.*t2.*t7.*t12.*(t4.*t19.*t29-param8.*t6.*(t4.*t10.*t29.*t30.*t31.*8.0e3+t4.*t11.*t29.*t30.*t31.*8.0e3).*(1.0./2.0e1)).*4.0e1+param1.*t4.*t8.*t9.*t12.*t20.*t29.*8.0e2).*5.0;t8.*t9.*t26.*u5.*-4.0e3-t24.*t36.*t38.*(t8.*t9.*t12.*u2.*8.0e2-param1.*t8.*t9.*t12.*u2.*8.0e2).*5.0;t8.*t9.*t26.*u9.*-4.0e3-t24.*t36.*t38.*(t8.*t9.*t12.*u3.*8.0e2-param1.*t8.*t9.*t12.*u3.*8.0e2).*5.0;param1.*t2.*t7.*t12.*t24.*t36.*t38.*-2.0e2;t8.*t9.*t26.*u2.*-4.0e3;t39;zero;zero;t8.*t9.*t26.*u3.*-4.0e3;zero;t39;zero];
end
% f = reshape(f,ng,nch,nd);
% f_udg = reshape(f_udg,ng,nch,nd,nc);
