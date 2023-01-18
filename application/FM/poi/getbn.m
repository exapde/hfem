function Bn = getbn(ib, nl, udg, param)

tau   = param{end};
ng = size(udg,1);
nch = 1;
Bn = zeros(ng,nch,nch);
Bn(:,1,1) = tau;


