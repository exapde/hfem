function uh = ldgubou(ib,ui,nl,p,udg,param,time)
% UH(:,:,i) = ubou(ib,uinf,nl,p,u,param,time);

%[ng,nc] = size(udg);
nch = 1;

switch ib
    case 1  % Dirichlet
        uh(:,1) = ui;            
    case 2  % neumman
        uh = udg(:,nch);
    case 3
        x = p(:,1);
        y = p(:,2);
        ui = sin(x)*sin(y);        
        uh = ui;
    otherwise
        error('unknown boundary type');
end




