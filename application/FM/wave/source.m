function [sr,sr_udg] = source(pg,udg,param,time)

[ng,nc] = size(udg);
nd = size(pg,2);
nch = nc/(nd+1);

% freq = param{3};
% fc = 400*sin(freq*2*pi*time).*exp(-pi*freq*time^2/8);

sr = zeros(ng,nch);
%sr = fc*exp(-100*((pg(:,1)+3).^2+(pg(:,2)-3).^2));

sr_udg = zeros(ng,nch,nc); 
