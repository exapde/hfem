% Width default
set(0,'defaultlinelinewidth',1.5)

% Resultats de reference

Xref(1) = 0.0517 ; refWA(1) = 1.846;  Xref(19) = -0.1001; refWA(19)= 14.520;
Xref(2) = 0.1182 ; refWA(2) = 5.271;  Xref(20) = -0.1142; refWA(20)= 14.451;
Xref(3) = 0.1583 ; refWA(3) = 8.257;  Xref(21) = -0.1247; refWA(21)= 14.862;
Xref(4) = 0.1837 ; refWA(4) = 10.799; Xref(22) = -0.1288; refWA(22)= 15.778;
Xref(5) = 0.1914 ; refWA(5) = 11.904; Xref(23) = -0.1271; refWA(23)= 16.961;
Xref(6) = 0.1953 ; refWA(6) = 12.892; Xref(24) = -0.1196; refWA(24)= 18.320;
Xref(7) = 0.1950 ; refWA(7) = 13.752; Xref(25) = -0.1055; refWA(25)= 19.817;
Xref(8) = 0.1901 ; refWA(8) = 14.472; Xref(26) = -0.0825; refWA(26)= 21.420;
Xref(9) = 0.1806 ; refWA(9) = 15.050; Xref(27) = -0.0484; refWA(27)= 23.100;
Xref(10)= 0.1671 ; refWA(10)= 15.501; Xref(28) = -0.0006; refWA(28)= 24.824;
Xref(11)= 0.1323 ; refWA(11)= 16.145; Xref(29) = 0.0626 ; refWA(29)= 26.565;
Xref(12)= 0.0923 ; refWA(12)= 16.602; Xref(30) = 0.1427 ; refWA(30)= 28.302;
Xref(13)= 0.0504 ; refWA(13)= 16.915; Xref(31) = 0.2403 ; refWA(31)= 30.023;
Xref(14)= 0.0083 ; refWA(14)= 17.008; Xref(32) = 0.3559 ; refWA(32)= 31.720;
Xref(15)= -0.0312; refWA(15)= 16.697; Xref(33) = 0.4898 ; refWA(33)= 33.388;
Xref(16)= -0.0622; refWA(16)= 15.780; Xref(34) = 0.6417 ; refWA(34)= 35.024;
Xref(17)= -0.0739; refWA(17)= 15.206; Xref(35) = 0.8114 ; refWA(35)= 36.626;
Xref(18)= -0.0861; refWA(18)= 14.767; Xref(36) = 1.0313	; refWA(36)= 38.450;


% Chargement resultats HDG (res.WB)
load('k2_8_t1_Dl10_AL_h635.mat');

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
ylim([-0.2 1.1]); xlim([0.,40.]); grid on;
%title('Comparaisons of radial deflections');
xlabel('Radial displacements of points A'); ylabel('Normalized pulling force P/Pmax');
legend('V_A 16x16 Shells Elem', 'V_A 8x12 HDG-k2 Elem', ...
      'Location','northwest');


