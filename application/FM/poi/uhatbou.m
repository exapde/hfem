function uh = uhatbou(ib,ui,nl,p,udg,param,time)        
% fhg = fbou(ib, uinf, nlg1, pg1, udg1, uhg, app.arg, app.time);
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
%      FHQ(N,NC,NC):         Jacobian of the flux flux vector w.r.t. Q


tau   = param{end};
switch ib
    case 1  % Dirichlet        
        uh = ui;
    case 2   % Neumman
        qn = udg(:,2).*nl(:,1) + udg(:,3).*nl(:,2) - ui;        
        uh = qn + tau*udg(:,1);                        
    otherwise
        error('unknown boundary type');
end

