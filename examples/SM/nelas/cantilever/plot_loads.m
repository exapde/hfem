% Resultats de reference
refWA=zeros(1,12); refWB = zeros(1,12);
refWA(1) = 2.455;  refWB(1) = 3.370;
refWA(2) = 4.277;  refWB(2) = 5.876;
refWA(3) = 5.649;  refWB(3) = 7.736;
refWA(4) = 6.725;  refWB(4) = 9.160;
refWA(5) = 7.602;  refWB(5) = 10.288;
refWA(6) = 8.340;  refWB(6) = 11.213;
refWA(7) = 8.974;  refWB(7) = 11.992;
refWA(8) = 9.529;  refWB(8) = 12.661;
refWA(9) = 10.023; refWB(9) = 13.247;
refWA(10)= 10.468; refWB(10)= 13.768;
refWA(11)= 10.876; refWB(11)= 14.240;
refWA(12)= 11.257; refWB(12)= 14.674;
refWA(13)= 11.620; refWB(13)= 15.081;
refWA(14)= 11.970; refWB(14)= 15.469;
refWA(15)= 12.310; refWB(15)= 15.842;
refWA(16)= 12.642; refWB(16)= 16.202;
refWA(17)= 12.966; refWB(17)= 16.550;
refWA(18)= 13.282; refWB(18)= 16.886;
refWA(19)= 13.590; refWB(19)= 17.212;
refWA(20)= 13.891; refWB(20)= 17.528;


% Chargement resultats HDG (res.WB)
load('ann_k1_10loads.mat');

Xref = linspace(0.05,1,20);
X    = linspace(0.1, 1,10);

% Plot results
figure(1); clf; plot(refWA,Xref,'-', res.WA2(:,3),X,'s', refWB,Xref,'-', res.WB2(:,3),X,'s'); ylim([0 1.]); grid on;
title('Comparaisons of vertical deflections'); xlabel('Vertical deflections Z_B'); ylabel('Line Force');
legend('Reference','HDG Solution','Location','northwest');

