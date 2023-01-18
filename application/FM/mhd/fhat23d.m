function [fh,fh_udg,fh_uh] = fhat23d(nl,pg,udg,uh,param,time)
%FHAT23D
%    [FH,FH_UDG,FH_UH] = FHAT23D(NL,PG,UDG,UH,PARAM,TIME)

%    This function was generated by the Symbolic Math Toolbox version 6.3.
%    21-Nov-2017 10:30:42
[ng,nc] = size(udg);
nch = 9;
nd = 2;
nl1 = nl(:,1);
nl2 = nl(:,2);
one = ones(ng,1);
param1 = param{1};
param2 = param{2};
param3 = param{3};
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
t3 = 1.0./uh1.^2;
t4 = uh6.^2;
t5 = t4.*(1.0./2.0);
t6 = uh7.^2;
t7 = t6.*(1.0./2.0);
t8 = uh8.^2;
t9 = t8.*(1.0./2.0);
t10 = 1.0./uh1;
t11 = uh3.^2;
t12 = param1-1.0;
t13 = t2.*t3.*(1.0./2.0);
t14 = t3.*t11.*(1.0./2.0);
t15 = uh4.^2;
t16 = t3.*t15.*(1.0./2.0);
t17 = t13+t14+t16;
t18 = t17.*uh1;
t19 = t5+t7+t9+t18-uh5;
t20 = uh6.*uh7;
t21 = t20-t10.*uh2.*uh3;
t22 = t10.*uh2.*uh7;
t23 = t22-t10.*uh3.*uh6;
t24 = param3.^2;
t25 = one.*param2;
t26 = 1.0./uh1.^3;
t27 = t2.*t26;
t28 = t11.*t26;
t29 = t15.*t26;
t30 = t27+t28+t29;
t31 = t13+t14+t16-t30.*uh1;
t32 = t12.*t31;
t33 = param1.*t2.*uh2;
t34 = param1.*t11.*uh2;
t35 = param1.*t15.*uh2;
t36 = uh1.*uh3.*uh6.*uh7.*2.0;
t37 = uh1.*uh4.*uh6.*uh8.*2.0;
t38 = param1.*t4.*uh1.*uh2;
t39 = param1.*t6.*uh1.*uh2;
t40 = param1.*t8.*uh1.*uh2;
t41 = t33+t34+t35+t36+t37+t38+t39+t40-t2.*uh2-t11.*uh2-t15.*uh2-t6.*uh1.*uh2.*2.0-t8.*uh1.*uh2.*2.0-param1.*uh1.*uh2.*uh5.*2.0;
t42 = param1.*t11.*uh3;
t43 = param1.*t2.*uh3;
t44 = param1.*t15.*uh3;
t45 = uh1.*uh2.*uh6.*uh7.*2.0;
t46 = uh1.*uh4.*uh7.*uh8.*2.0;
t47 = param1.*t4.*uh1.*uh3;
t48 = param1.*t6.*uh1.*uh3;
t49 = param1.*t8.*uh1.*uh3;
t50 = t42+t43+t44+t45+t46+t47+t48+t49-t2.*uh3-t11.*uh3-t15.*uh3-t4.*uh1.*uh3.*2.0-t8.*uh1.*uh3.*2.0-param1.*uh1.*uh3.*uh5.*2.0;
fh = [nl1.*uh2+nl2.*uh3+param2.*(u1-uh1);-nl2.*t21+nl1.*(-t5+t7+t9+t2.*t10-t12.*t19)+param2.*(u2-uh2);-nl1.*t21+nl2.*(t5-t7+t9+t10.*t11-t12.*t19)+param2.*(u3-uh3);-nl1.*(uh6.*uh8-t10.*uh2.*uh4)-nl2.*(uh7.*uh8-t10.*uh3.*uh4)+param2.*(u4-uh4);param2.*(u5-uh5)-nl1.*t3.*t41.*(1.0./2.0)-nl2.*t3.*t50.*(1.0./2.0);-nl2.*t23+nl1.*uh9+param2.*(u6-uh6);nl1.*t23+nl2.*uh9+param2.*(u7-uh7);param2.*(u8-uh8)+nl1.*(t10.*uh2.*uh8-t10.*uh4.*uh6)+nl2.*(t10.*uh3.*uh8-t10.*uh4.*uh7);param2.*(u9-uh9)+nl1.*t24.*uh6+nl2.*t24.*uh7];
if nargout > 1
    fh_udg = [t25;zero;zero;zero;zero;zero;zero;zero;zero;zero;t25;zero;zero;zero;zero;zero;zero;zero;zero;zero;t25;zero;zero;zero;zero;zero;zero;zero;zero;zero;t25;zero;zero;zero;zero;zero;zero;zero;zero;zero;t25;zero;zero;zero;zero;zero;zero;zero;zero;zero;t25;zero;zero;zero;zero;zero;zero;zero;zero;zero;t25;zero;zero;zero;zero;zero;zero;zero;zero;zero;t25;zero;zero;zero;zero;zero;zero;zero;zero;zero;t25];
end
if nargout > 2
    t51 = t3.*uh2.*uh7;
    t52 = t51-t3.*uh3.*uh6;
    t53 = param1.*uh2.*uh3.*2.0;
    t54 = uh1.*uh6.*uh7.*2.0;
    t55 = t53+t54-uh2.*uh3.*2.0;
    t56 = t8.*uh1.*2.0;
    t57 = param1.*uh1.*uh5.*2.0;
    t58 = nl1.*t10.*uh2;
    t59 = nl2.*t10.*uh3;
    t60 = t12.*uh6;
    t61 = nl1.*one.*t10.*uh4;
    t62 = uh1.*uh4.*uh8.*2.0;
    t63 = nl2.*one.*t10.*uh4;
    t64 = uh8-t12.*uh8;
    t65 = nl1.*uh6;
    t66 = nl2.*uh7;
    t67 = -param2+t58+t59;
    t68 = one.*t67;
    t69 = nl1.*one;
    t70 = nl2.*one;
    fh_uh = [-t25;-one.*(nl1.*(t32+t2.*t3)+nl2.*t3.*uh2.*uh3);-one.*(nl2.*(t32+t3.*t11)+nl1.*t3.*uh2.*uh3);-one.*(nl1.*t3.*uh2.*uh4+nl2.*t3.*uh3.*uh4);-one.*(nl1.*t3.*(t6.*uh2.*-2.0-t8.*uh2.*2.0+param1.*t4.*uh2+param1.*t6.*uh2+param1.*t8.*uh2-param1.*uh2.*uh5.*2.0+uh3.*uh6.*uh7.*2.0+uh4.*uh6.*uh8.*2.0).*(1.0./2.0)+nl2.*t3.*(t4.*uh3.*-2.0-t8.*uh3.*2.0+param1.*t4.*uh3+param1.*t6.*uh3+param1.*t8.*uh3-param1.*uh3.*uh5.*2.0+uh2.*uh6.*uh7.*2.0+uh4.*uh7.*uh8.*2.0).*(1.0./2.0)-nl1.*t26.*t41-nl2.*t26.*t50);nl2.*one.*t52;-nl1.*one.*t52;-one.*(nl1.*(t3.*uh2.*uh8-t3.*uh4.*uh6)+nl2.*(t3.*uh3.*uh8-t3.*uh4.*uh7));zero;t69;one.*(-param2+t59+nl1.*(t10.*uh2.*2.0-t10.*t12.*uh2));one.*(nl1.*t10.*uh3-nl2.*t10.*t12.*uh2);t61;one.*(nl1.*t3.*(t2.*3.0+t11+t15+t56+t57-param1.*t2.*3.0-param1.*t11-param1.*t15+t6.*uh1.*2.0-param1.*t4.*uh1-param1.*t6.*uh1-param1.*t8.*uh1).*(1.0./2.0)-nl2.*t3.*t55.*(1.0./2.0));-nl2.*one.*t10.*uh7;nl1.*one.*t10.*uh7;nl1.*one.*t10.*uh8;zero;t70;one.*(nl2.*t10.*uh2-nl1.*t10.*t12.*uh3);one.*(-param2+t58+nl2.*(t10.*uh3.*2.0-t10.*t12.*uh3));t63;one.*(nl2.*t3.*(t2+t11.*3.0+t15+t56+t57-param1.*t2-param1.*t11.*3.0-param1.*t15+t4.*uh1.*2.0-param1.*t4.*uh1-param1.*t6.*uh1-param1.*t8.*uh1).*(1.0./2.0)-nl1.*t3.*t55.*(1.0./2.0));nl2.*one.*t10.*uh6;-nl1.*one.*t10.*uh6;nl2.*one.*t10.*uh8;zero;zero;-nl1.*one.*t10.*t12.*uh4;-nl2.*one.*t10.*t12.*uh4;t68;-one.*(nl1.*t3.*(uh2.*uh4.*-2.0+param1.*uh2.*uh4.*2.0+uh1.*uh6.*uh8.*2.0).*(1.0./2.0)+nl2.*t3.*(uh3.*uh4.*-2.0+param1.*uh3.*uh4.*2.0+uh1.*uh7.*uh8.*2.0).*(1.0./2.0));zero;zero;-one.*(nl1.*t10.*uh6+nl2.*t10.*uh7);zero;zero;nl1.*one.*t12;nl2.*one.*t12;zero;one.*(-param2+nl1.*param1.*t10.*uh2+nl2.*param1.*t10.*uh3);zero;zero;zero;zero;zero;-one.*(t66+nl1.*(t60+uh6));-one.*(nl1.*uh7+nl2.*(t60-uh6));-nl1.*one.*uh8;-one.*(nl1.*t3.*(t62+uh1.*uh3.*uh7.*2.0+param1.*uh1.*uh2.*uh6.*2.0).*(1.0./2.0)+nl2.*t3.*(uh1.*uh2.*uh7.*2.0-uh1.*uh3.*uh6.*4.0+param1.*uh1.*uh3.*uh6.*2.0).*(1.0./2.0));-one.*(param2-t59);-nl1.*one.*t10.*uh3;-t61;nl1.*one.*t24;zero;-one.*(nl2.*uh6-nl1.*(uh7-t12.*uh7));-one.*(t65+nl2.*(uh7+t12.*uh7));-nl2.*one.*uh8;-one.*(nl2.*t3.*(t62+uh1.*uh2.*uh6.*2.0+param1.*uh1.*uh3.*uh7.*2.0).*(1.0./2.0)+nl1.*t3.*(uh1.*uh2.*uh7.*-4.0+uh1.*uh3.*uh6.*2.0+param1.*uh1.*uh2.*uh7.*2.0).*(1.0./2.0));-nl2.*one.*t10.*uh2;-one.*(param2-t58);-t63;nl2.*one.*t24;zero;nl1.*one.*t64;nl2.*one.*t64;-one.*(t65+t66);-one.*(nl1.*t3.*(uh1.*uh2.*uh8.*-4.0+uh1.*uh4.*uh6.*2.0+param1.*uh1.*uh2.*uh8.*2.0).*(1.0./2.0)+nl2.*t3.*(uh1.*uh3.*uh8.*-4.0+uh1.*uh4.*uh7.*2.0+param1.*uh1.*uh3.*uh8.*2.0).*(1.0./2.0));zero;zero;t68;zero;zero;zero;zero;zero;zero;t69;t70;zero;-t25];
end
fh = reshape(fh,ng,nch);
fh_udg = reshape(fh_udg,ng,nch,nc);
fh_uh = reshape(fh_uh,ng,nch,nch);