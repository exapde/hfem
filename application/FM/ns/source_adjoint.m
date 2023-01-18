function [sr,sr_udg] = source_adjoint(p,udg,param,time)

[ng,nc] = size(udg);

a  = 0.0;
sr = a.*udg(:,1);
sr_udg = zeros(ng,nc); 
sr_udg(:,1) = a;




