function [fh,fh_udg,fh_uh] = fhat(nl,p,udg,uh,param,time)
%FHAT flux function
%   [fh,fhu,fhq,fhm] = fhat(nl,p,u,q,m,param)
%
%      NL(N,ND)              Normal N points
%      P(N,ND)               Coordinates for N points
%      U(N,NC)               Unknown vector for N points with NC components
%      Q(N,NC,ND)            Flux vector for N points with NC components in the
%                            coordinate directions
%      M(N,NC)               Hybrid unkowns for N points with NC components
%      PARAM                 Parameter list
%      FH(N,NC):              Volume flux at N points
%      FHU(N,NC,NC):         Jacobian of the flux flux vector w.r.t. U
%      FHQ(N,NC,NC,ND):      Jacobian of the flux flux vector w.r.t. Q
%      FHM(N,NC,NC):         Jacobian of the flux flux vector w.r.t. Lambda

[ng,nch] = size(uh);

u = udg(:,1:nch);

[f,f_uhdg] = flux(p,uh,param);

% [An,Anm] = getan(nl,uh,param,2);
% 
% fh = permute(mapContractK(f,nl,2,3,1,2,[],1) + mapContractK(An,u-uh,2,3,1,2,[],1),[2 1]);
% 
% fh_udg = An;
% 
% fh_uh = permute(mapContractK(f_uhdg,nl,[2 4],3,1,2,[],1)+mapContractK(Anm,u-uh,[2 4],3,1,2,[],1),[3 1 2])-An;

tau = param{3};
fh = permute(mapContractK(f,nl,2,3,1,2,[],1),[2 1]) + tau*(u-uh);
fh_udg = zeros(ng,nch,nch);
fh_uh = permute(mapContractK(f_uhdg,nl,[2 4],3,1,2,[],1),[3 1 2]);
for i=1:nch
    fh_udg(:,i,i) = tau;
    fh_uh(:,i,i) = fh_uh(:,i,i) - tau;
end

% [ng,nc] = size(udg);
% nch = 4;
% nd = 2;
% nl1 = nl(:,1);
% nl2 = nl(:,2);
% one = ones(ng,1);
% param1 = param{1};
% param3 = param{3};
% param4 = param{4};
% u1 = udg(:,1);
% u2 = udg(:,2);
% u3 = udg(:,3);
% u4 = udg(:,4);
% uh1 = uh(:,1);
% uh2 = uh(:,2);
% uh3 = uh(:,3);
% uh4 = uh(:,4);
% zero = zeros(ng,1);
% t2 = 1.0./param4;
% t3 = uh2.^2;
% t4 = uh3.^2;
% t5 = t2.*uh1.*2.0e1;
% t6 = exp(t5);
% t7 = t6+1.0;
% t8 = log(t7);
% t9 = 1.0./t8;
% t10 = nl1.*uh2;
% t11 = nl2.*uh3;
% t12 = t10+t11;
% t13 = one.*param3;
% t14 = 1.0./param4.^2;
% t15 = 1.0./t8.^2;
% t16 = nl1.*t3.*3.0e1;
% t17 = nl1.*t4.*1.0e1;
% t18 = nl2.*uh2.*uh3.*2.0e1;
% t19 = t16+t17+t18-nl1.*param1.*t3.*1.0e1-nl1.*param1.*t4.*1.0e1;
% t20 = 1.0./t7;
% t21 = nl2.*t3.*1.0e1;
% t22 = nl2.*t4.*3.0e1;
% t23 = nl1.*uh2.*uh3.*2.0e1;
% t24 = t21+t22+t23-nl2.*param1.*t3.*1.0e1-nl2.*param1.*t4.*1.0e1;
% t25 = t3+t4;
% t26 = param1-1.0;
% fh = [t10+t11+param3.*(u1-uh1);-t2.*(nl1.*param4.*uh4-param3.*param4.*u2+param3.*param4.*uh2-nl1.*param1.*param4.*uh4)+t2.*t9.*t19;-t2.*(nl2.*param4.*uh4-param3.*param4.*u3+param3.*param4.*uh3-nl2.*param1.*param4.*uh4)+t2.*t9.*t24;param3.*(u4-uh4)+param1.*t2.*t9.*t12.*uh4.*2.0e1-t12.*t14.*t15.*t25.*t26.*2.0e2];
% if nargout > 1
%     fh_udg = [t13;zero;zero;zero;zero;t13;zero;zero;zero;zero;t13;zero;zero;zero;zero;t13];
% end
% if nargout > 2
%     t27 = nl1.*uh3.*2.0e1;
%     t28 = nl2.*uh2.*2.0e1;
%     fh_uh = [-t13;one.*t6.*t14.*t15.*t19.*t20.*-2.0e1;one.*t6.*t14.*t15.*t20.*t24.*-2.0e1;one.*(1.0./param4.^3.*t6.*1.0./t8.^3.*t12.*t20.*t25.*t26.*8.0e3-param1.*t6.*t12.*t14.*t15.*t20.*uh4.*4.0e2);nl1.*one;-one.*(param3-t2.*t9.*(nl1.*uh2.*6.0e1+nl2.*uh3.*2.0e1-nl1.*param1.*uh2.*2.0e1));one.*t2.*t9.*(t27+t28-nl2.*param1.*uh2.*2.0e1);-one.*(nl1.*param1.*t2.*t9.*uh4.*-2.0e1+nl1.*t14.*t15.*t25.*t26.*2.0e2+t12.*t14.*t15.*t26.*uh2.*4.0e2);nl2.*one;one.*t2.*t9.*(t27+t28-nl1.*param1.*uh3.*2.0e1);-one.*(param3-t2.*t9.*(nl1.*uh2.*2.0e1+nl2.*uh3.*6.0e1-nl2.*param1.*uh3.*2.0e1));-one.*(nl2.*param1.*t2.*t9.*uh4.*-2.0e1+nl2.*t14.*t15.*t25.*t26.*2.0e2+t12.*t14.*t15.*t26.*uh3.*4.0e2);zero;-one.*t2.*(nl1.*param4-nl1.*param1.*param4);-one.*t2.*(nl2.*param4-nl2.*param1.*param4);-one.*(param3-param1.*t2.*t9.*t12.*2.0e1)];
% end
% fh = reshape(fh,ng,nch);
% fh_udg = reshape(fh_udg,ng,nch,nc);
% fh_uh = reshape(fh_uh,ng,nch,nch);
