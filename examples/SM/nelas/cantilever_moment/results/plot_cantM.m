% Width default
set(0,'defaultlinelinewidth',2)

% Chargement resultats HDG (res.WB)
%load('cantM_k2_16_t10_DL.mat');
nloads = size(res.WA1,1);

% load increment vector :
if isfield(res,'Loa')
    X = res.Loa;
else
    X = linspace(1./nloads,1,nloads);
end

% Resultats de reference
Xpi = 2.*pi*X;
refUA =-12./Xpi .* sin(Xpi) + 12.;
refWA = 12./Xpi .* (1. -cos(Xpi));

%WA = res.WA1(:,3);
%norm(refWA-WA)
WA = 0.5*(res.WA1(:,3)+res.WA2(:,3))-0.05;
%norm(refWA-WA)
UA = 12.- 0.5*(res.WA1(:,1)+res.WA2(:,1));
%UA = 10. - res.WA2(:,1);

% Plot results
figure(1); clf; plot(refUA,X,'-', UA,X,'-', refWA,X,'-', WA,X,'-'); ylim([0 1.]); xlim([0 15.]); grid on;
title('Deflections'); xlabel('Tip Deflections'); ylabel('Normalized Bending Moment');
legend('Reference-X','HDG Solution-X','Reference-Z','HDG Solution-Z','Location','southeast');

