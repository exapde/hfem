function [sr,sr_udg] = source_adjoint(p,udg,param,time)

[ng,nc] = size(udg);

%a = (1/(2*pi*0.0012))*exp(-0.5*(p(:,1).^2/(0.02^2) + (p(:,2)-0.25).^2/(0.06^2)));
a  = 0.0;
sr = a.*udg(:,1);
sr_udg = zeros(ng,nc); 
sr_udg(:,1) = a;




