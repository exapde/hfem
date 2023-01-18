% Width default
set(0,'defaultlinelinewidth',1.5)

% Resultats de reference
Xref = zeros(1,23); refWA=zeros(1,32);
Xref(1) = 0.0877; refWA(1) = 0.693;  Xref(19) = 0.3245; refWA(19)= 16.590;
Xref(2) = 0.1980; refWA(2) = 1.638;  Xref(20) = 0.2717; refWA(20)= 17.094;
Xref(3) = 0.3473; refWA(3) = 3.087;  Xref(21) = 0.2272; refWA(21)= 17.657;
Xref(4) = 0.4686; refWA(4) = 4.477;  Xref(22) = 0.1940; refWA(22)= 18.299;
Xref(5) = 0.5647; refWA(5) = 5.802;  Xref(23) = 0.1750; refWA(23)= 19.028;
Xref(6) = 0.6381; refWA(6) = 7.057;  Xref(24) = 0.1729; refWA(24)= 19.852;
Xref(7) = 0.6908; refWA(7) = 8.237;  Xref(25) = 0.1905; refWA(25)= 20.771;
Xref(8) = 0.7246; refWA(8) = 9.339;  Xref(26) = 0.2303; refWA(26)= 21.780;
Xref(9) = 0.7412; refWA(9) = 10.358; Xref(27) = 0.2950; refWA(27)= 22.875;
Xref(10)= 0.7421; refWA(10)= 11.293; Xref(28) = 0.3871; refWA(28)= 24.049;
Xref(11)= 0.7286; refWA(11)= 12.141; Xref(29) = 0.4443; refWA(29)= 24.663;
Xref(12)= 0.7023; refWA(12)= 12.903; Xref(30) = 0.5093; refWA(30)= 25.293;
Xref(13)= 0.6649; refWA(13)= 13.583; Xref(31) = 0.5826; refWA(31)= 25.940;
Xref(14)= 0.6182; refWA(14)= 14.188; Xref(32) = 0.6644; refWA(32)= 26.601;
Xref(15)= 0.5643; refWA(15)= 14.728; Xref(33) = 0.7551; refWA(33)= 27.276;
Xref(16)= 0.5055; refWA(16)= 15.217; Xref(34) = 0.8549; refWA(34)= 27.964;
Xref(17)= 0.4442; refWA(17)= 15.676; Xref(35) = 0.9643; refWA(35)= 28.663;
Xref(18)= 0.3830; refWA(18)= 16.125; Xref(36) = 1.0835; refWA(36)= 29.374;


% Chargement resultats HDG (res.WB)
% load('k2_8_t1.mat');

nloads = size(res.WA1,1);

% load increment vector :
if isfield(res,'Loa')
    X = res.Loa;
else
    X = linspace(1./nloads,1,nloads);
end

R  = 2540. ;
WA = R - 0.5*(res.WB1(:,2) + res.WB2(:,2));

% Plot results
figure(1); clf; plot(refWA,Xref,'-', WA, X,'s');
ylim([0. 1.1]); xlim([0.,30.]); grid on;
%title('Comparaisons of radial deflections');
xlabel('Radial displacements of points A'); ylabel('Normalized pulling force P/Pmax');
legend('V_A 16x16 Shells Elem', 'V_A 8x12 HDG-k2 Elem', ...
      'Location','northwest');


