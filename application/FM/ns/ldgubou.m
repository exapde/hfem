function uh = ldgubou(ib,ui,nl,p,udg,param,time)
% UH(:,:,i) = ubou(ib,uinf,nl,p,u,param,time);

%[ng,nc] = size(udg);
nch = size(ui,2);

switch ib
    case 1  % far field
        uh = udg(:,1:nch);
    case 2  % adiabatic wall
        uh = udg(:,1:nch);
        uh(:,2:1:(nch-1)) = 0;
    otherwise
        error('unknown boundary type');
end




