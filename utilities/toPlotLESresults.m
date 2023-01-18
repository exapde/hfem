
close all
figure(1)
plot(x2d{1}(:,1),Cf2d_avg,'.')
hold on
plot(linspace(-0.1,1.1,100),zeros(100,1),'-.')
xlim([-0.1,1.1]);
ylim([-0.02,0.02]);
xlabel('$x/c$','Interpreter','latex','FontSize',16);
ylabel('$C_f$','Interpreter','latex','FontSize',16);
title('Time-average skin-friction coefficient','Interpreter','latex','FontSize',16);
set(gcf,'color','w');

close all
figure(1)
plot(x2d{1}(:,1),Cp2d_avg,'.')
xlim([-0.1,1.1]);
xlabel('$x/c$','Interpreter','latex','FontSize',16);
ylabel('$C_p$','Interpreter','latex','FontSize',16);
title('Time-average pressure coefficient','Interpreter','latex','FontSize',16);
set(gcf,'color','w');

close all

figure(1)
subplot(2,2,1)
plot(sFieldLowerSide,deltaS_s2tAvg(iLowerSide),'.')
hold on
plot(sFieldLowerSide,thetaS_s2tAvg(iLowerSide),'.')
xlim([-0.1,1.1]);
ylim([0,0.015]);%ylim([0,max(delta_s2tAvg)*1.15]);
xlabel('$x/c$','Interpreter','latex','FontSize',16)
ylabel('Thickness','Interpreter','latex','FontSize',16)
title('Boundary layer thickness - Pressure side','Interpreter','latex','FontSize',16)
legend1 = legend('Displacement thickness','Momentum thickness');
set(legend1,'FontSize',14,'Interpreter','latex');
set(gcf,'color','w');

subplot(2,2,2)
plot(sFieldLowerSide,HS_s2tAvg(iLowerSide),'.')
xlim([-0.1,1.1]);
ylim([0,12]);%ylim([0,max(H_s2tAvg)*1.15]);
xlabel('$x/c$','Interpreter','latex','FontSize',16)
ylabel('$H_1$','Interpreter','latex','FontSize',16)
title('Shape parameter - Pressure side','Interpreter','latex','FontSize',16)
set(gcf,'color','w');

subplot(2,2,3)
plot(sFieldLowerSide, NS1(iLowerSide),'.')
hold on
plot(sFieldLowerSide, N_envelopeMethod_lowerSide,'.');
xlim([-0.1,1.1]);
ylim([0,14]);
% ylim([0,max(max(N_envelopeMethod_lowerSide(:)),max(N_envelopeMethod_upperSide(:)))*1.15]);
xlabel('$x/c$','Interpreter','latex','FontSize',16)
ylabel('$N_1$','Interpreter','latex','FontSize',16)
title('Amplification factor of streamwise instabilities - Pressure side','Interpreter','latex','FontSize',16)
legend1 = legend('LES','$e^N$ envelope method');
set(legend1,'FontSize',14,'Interpreter','latex');
set(gcf,'color','w');

subplot(2,2,4)
plot(sFieldLowerSide,NS2(iLowerSide),'.')
xlim([-0.1,1.1]);
ylim([0,14]);
xlabel('$x/c$','Interpreter','latex','FontSize',16)
ylabel('$N_2$','Interpreter','latex','FontSize',16)
title('Amplification factor of cross-flow instabilities - Pressure side','Interpreter','latex','FontSize',16)
legend1 = legend('LES');
set(legend1,'FontSize',14,'Interpreter','latex');
set(gcf,'color','w');


figure(2)
subplot(2,2,1)
plot(sFieldUpperSide,deltaS_s2tAvg(iUpperSide),'-')
hold on
plot(sFieldUpperSide,thetaS_s2tAvg(iUpperSide),'-')
xlim([-0.1,1.1]);
ylim([0,0.015]);%ylim([0,max(delta_s2tAvg)*1.15]);
xlabel('$x/c$','Interpreter','latex','FontSize',16)
ylabel('Thickness','Interpreter','latex','FontSize',16)
title('Boundary layer thickness - Suction side','Interpreter','latex','FontSize',16)
legend1 = legend('Displacement thickness','Momentum thickness');
set(legend1,'FontSize',14,'Interpreter','latex');
set(gcf,'color','w');

subplot(2,2,2)
plot(sFieldUpperSide,HS_s2tAvg(iUpperSide),'-')
xlim([-0.1,1.1]);
ylim([0,12]);%ylim([0,max(H_s2tAvg)*1.15]);
xlabel('$x/c$','Interpreter','latex','FontSize',16)
ylabel('$H_1$','Interpreter','latex','FontSize',16)
title('Shape parameter - Suction side','Interpreter','latex','FontSize',16)
set(gcf,'color','w');

subplot(2,2,3)
plot(sFieldUpperSide, NS1(iUpperSide),'-')
hold on
plot(sFieldUpperSide, N_envelopeMethod_upperSide,'-');
xlim([-0.1,1.1]);
ylim([0,14]);
% ylim([0,max(max(N_envelopeMethod_lowerSide(:)),max(N_envelopeMethod_upperSide(:)))*1.15]);
xlabel('$x/c$','Interpreter','latex','FontSize',16)
ylabel('$N_1$','Interpreter','latex','FontSize',16)
title('Amplification factor of streamwise instabilities - Suction side','Interpreter','latex','FontSize',16)
legend1 = legend('LES results','Predicted by $e^N$ envelope method');
set(legend1,'FontSize',14,'Interpreter','latex');
set(gcf,'color','w');

subplot(2,2,4)
plot(sFieldUpperSide,NS2(iUpperSide),'-')
xlim([-0.1,1.1]);
ylim([0,14]);
xlabel('$x/c$','Interpreter','latex','FontSize',16)
ylabel('$N_2$','Interpreter','latex','FontSize',16)
title('Amplification factor of cross-flow instabilities - Suction side','Interpreter','latex','FontSize',16)
legend1 = legend('LES results');
set(legend1,'FontSize',14,'Interpreter','latex');
set(gcf,'color','w');



% close all
% 
% figure(1)
% subplot(1,2,1)
% plot(sFieldLowerSide,deltaS_s2tAvg(iLowerSide),'-')
% hold on
% plot(sFieldUpperSide,deltaS_s2tAvg(iUpperSide),'-')
% xlim([-0.1,1.1]);
% ylim([0,0.015]);%ylim([0,max(delta_s2tAvg)*1.15]);
% xlabel('$x/c$','Interpreter','latex','FontSize',16)
% ylabel('$\delta_1$','Interpreter','latex','FontSize',16)
% title('Displacement thickness','Interpreter','latex','FontSize',16)
% legend1 = legend('Pressure side','Suction side');
% set(legend1,'FontSize',14,'Interpreter','latex');
% set(gcf,'color','w');
% 
% subplot(1,2,2)
% plot(sFieldLowerSide,thetaS_s2tAvg(iLowerSide),'-')
% hold on
% plot(sFieldUpperSide,thetaS_s2tAvg(iUpperSide),'-')
% xlim([-0.1,1.1]);
% ylim([0,0.008]);%ylim([0,max(theta_s2tAvg)*1.15]);
% xlabel('$x/c$','Interpreter','latex','FontSize',16)
% ylabel('$\theta_1$','Interpreter','latex','FontSize',16)
% title('Momentum thickness','Interpreter','latex','FontSize',16)
% legend1 = legend('Pressure side','Suction side');
% set(legend1,'FontSize',14,'Interpreter','latex');
% set(gcf,'color','w');
% 
% figure(2)
% plot(sFieldLowerSide,HS_s2tAvg(iLowerSide),'-')
% hold on
% plot(sFieldUpperSide,HS_s2tAvg(iUpperSide),'-')
% xlim([-0.1,1.1]);
% ylim([0,12]);%ylim([0,max(H_s2tAvg)*1.15]);
% xlabel('$x/c$','Interpreter','latex','FontSize',16)
% ylabel('$H_1$','Interpreter','latex','FontSize',16)
% title('Shape parameter','Interpreter','latex','FontSize',16)
% legend1 = legend('Pressure side','Suction side');
% set(legend1,'FontSize',14,'Interpreter','latex');
% set(gcf,'color','w');
% 
% figure(3)
% subplot(1,2,1)
% plot(sFieldLowerSide, NS1(iLowerSide),'-')
% hold on
% plot(sFieldUpperSide, NS1(iUpperSide),'-')
% hold on
% plot(sFieldLowerSide, N_envelopeMethod_lowerSide,'-');
% hold on
% plot(sFieldUpperSide, N_envelopeMethod_upperSide,'-');
% xlim([-0.1,1.1]);
% ylim([0,14]);
% % ylim([0,max(max(N_envelopeMethod_lowerSide(:)),max(N_envelopeMethod_upperSide(:)))*1.15]);
% xlabel('$x/c$','Interpreter','latex','FontSize',16)
% ylabel('$N_1$','Interpreter','latex','FontSize',16)
% title('Amplification factor of streamwise instabilities','Interpreter','latex','FontSize',16)
% legend1 = legend('LES - Pressure side','LES - Suction side','Envelope method - Pressure side','Envelope method - Suction side');
% set(legend1,'FontSize',14,'Interpreter','latex');
% set(gcf,'color','w');
% 
% subplot(1,2,2)
% plot(sFieldLowerSide,NS2(iLowerSide),'-')
% hold on
% plot(sFieldUpperSide,NS2(iUpperSide),'-')
% xlim([-0.1,1.1]);
% ylim([0,14]);
% xlabel('$x/c$','Interpreter','latex','FontSize',16)
% ylabel('$N_2$','Interpreter','latex','FontSize',16)
% title('Amplification factor of cross-flow instabilities','Interpreter','latex','FontSize',16)
% legend1 = legend('LES - Pressure side','LES - Suction side');
% set(legend1,'FontSize',14,'Interpreter','latex');
% set(gcf,'color','w');


figure(4)
subplot(2,3,1)
plot(u_pert_InNormal_s2tRMS(iUnique25PctLower,in25PctLower,1)/ueSmag_s2tAvg(iUnique25PctLower), nDist_InNormal(iUnique25PctLower,in25PctLower)/max(nDist_InNormal(iUnique25PctLower,in25PctLower)),'-');
xlabel('$u^{''}_{1,RMS} / \bar{u}_{1,e}$','Interpreter','latex','FontSize',16);
ylabel('$n / n_e$','Interpreter','latex','FontSize',16);
title('Pressure side ($x/c=0.25$)','Interpreter','latex','FontSize',16);
set(gcf,'color','w');

subplot(2,3,2)
plot(u_pert_InNormal_s2tRMS(iUnique50PctLower,in50PctLower,1)/ueSmag_s2tAvg(iUnique50PctLower), nDist_InNormal(iUnique50PctLower,in50PctLower)/max(nDist_InNormal(iUnique50PctLower,in50PctLower)),'-');
xlabel('$u^{''}_{1,RMS} / \bar{u}_{1,e}$','Interpreter','latex','FontSize',16);
ylabel('$n / n_e$','Interpreter','latex','FontSize',16);
title('Pressure side ($x/c=0.50$)','Interpreter','latex','FontSize',16);
set(gcf,'color','w');

subplot(2,3,3)
plot(u_pert_InNormal_s2tRMS(iUnique75PctLower,in75PctLower,1)/ueSmag_s2tAvg(iUnique75PctLower), nDist_InNormal(iUnique75PctLower,in75PctLower)/max(nDist_InNormal(iUnique75PctLower,in75PctLower)),'-');
xlabel('$u^{''}_{1,RMS} / \bar{u}_{1,e}$','Interpreter','latex','FontSize',16);
ylabel('$n / n_e$','Interpreter','latex','FontSize',16);
title('Pressure side ($x/c=0.75$)','Interpreter','latex','FontSize',16);
set(gcf,'color','w');

subplot(2,3,4)
plot(u_pert_InNormal_s2tRMS(iUnique25PctUpper,in25PctUpper,1)/ueSmag_s2tAvg(iUnique25PctUpper), nDist_InNormal(iUnique25PctUpper,in25PctUpper)/max(nDist_InNormal(iUnique25PctUpper,in25PctUpper)),'-');
xlabel('$u^{''}_{1,RMS} / \bar{u}_{1,e}$','Interpreter','latex','FontSize',16);
ylabel('$n / n_e$','Interpreter','latex','FontSize',16);
title('Suction side ($x/c=0.25$)','Interpreter','latex','FontSize',16);
set(gcf,'color','w');

subplot(2,3,5)
plot(u_pert_InNormal_s2tRMS(iUnique50PctUpper,in50PctUpper,1)/ueSmag_s2tAvg(iUnique50PctUpper), nDist_InNormal(iUnique50PctUpper,in50PctUpper)/max(nDist_InNormal(iUnique50PctUpper,in50PctUpper)),'-');
xlabel('$u^{''}_{1,RMS} / \bar{u}_{1,e}$','Interpreter','latex','FontSize',16);
ylabel('$n / n_e$','Interpreter','latex','FontSize',16);
title('Suction side ($x/c=0.50$)','Interpreter','latex','FontSize',16);
set(gcf,'color','w');

subplot(2,3,6)
plot(u_pert_InNormal_s2tRMS(iUnique75PctUpper,in75PctUpper,1)/ueSmag_s2tAvg(iUnique75PctUpper), nDist_InNormal(iUnique75PctUpper,in75PctUpper)/max(nDist_InNormal(iUnique75PctUpper,in75PctUpper)),'-');
xlabel('$u^{''}_{1,RMS} / \bar{u}_{1,e}$','Interpreter','latex','FontSize',16);
ylabel('$n / n_e$','Interpreter','latex','FontSize',16);
title('Suction side ($x/c=0.75$)','Interpreter','latex','FontSize',16);
set(gcf,'color','w');


figure(5)
subplot(1,2,1)
plot(u1_InNormal_s2tAvg(iUnique25PctLower,in25PctLower)/u1_InNormal_s2tAvg(iUnique25PctLower,in25PctLower(end)), nDist_InNormal(iUnique25PctLower,in25PctLower)/max(nDist_InNormal(iUnique25PctLower,in25PctLower)),'-');
hold on
plot(u1_InNormal_s2tAvg(iUnique50PctLower,in50PctLower)/u1_InNormal_s2tAvg(iUnique50PctLower,in50PctLower(end)), nDist_InNormal(iUnique50PctLower,in50PctLower)/max(nDist_InNormal(iUnique50PctLower,in50PctLower)),'-');
hold on
plot(u1_InNormal_s2tAvg(iUnique75PctLower,in75PctLower)/u1_InNormal_s2tAvg(iUnique75PctLower,in75PctLower(end)), nDist_InNormal(iUnique75PctLower,in75PctLower)/max(nDist_InNormal(iUnique75PctLower,in75PctLower)),'-');
xlim([-0.1,1.1]); ylim([0,1]);
xlabel('$\bar{u}_1 / \bar{u}_{1,e}$','Interpreter','latex','FontSize',16);
ylabel('$n / n_e$','Interpreter','latex','FontSize',16);
title('Average streamwise velocity profiles - Pressure side','Interpreter','latex','FontSize',16);
legend1 = legend(['$x/c$ = 0.25 (H = ',num2str(HS_s2tAvg(iLowerSide(find(sFieldLowerSide<0.25,1,'first')))),')'], ...
    ['$x/c$ = 0.50 (H = ',num2str(HS_s2tAvg(iLowerSide(find(sFieldLowerSide<0.50,1,'first')))),')'], ...
    ['$x/c$ = 0.75 (H = ',num2str(HS_s2tAvg(iLowerSide(find(sFieldLowerSide<0.75,1,'first')))),')']);
set(legend1,'FontSize',14,'Interpreter','latex');
set(gcf,'color','w');

subplot(1,2,2)
plot(u1_InNormal_s2tAvg(iUnique25PctUpper,in25PctUpper)/u1_InNormal_s2tAvg(iUnique25PctUpper,in25PctUpper(end)), nDist_InNormal(iUnique25PctUpper,in25PctUpper)/max(nDist_InNormal(iUnique25PctUpper,in25PctUpper)),'-');
hold on
plot(u1_InNormal_s2tAvg(iUnique50PctUpper,in50PctUpper)/u1_InNormal_s2tAvg(iUnique50PctUpper,in50PctUpper(end)), nDist_InNormal(iUnique50PctUpper,in50PctUpper)/max(nDist_InNormal(iUnique50PctUpper,in50PctUpper)),'-');
hold on
plot(u1_InNormal_s2tAvg(iUnique75PctUpper,in75PctUpper)/u1_InNormal_s2tAvg(iUnique75PctUpper,in75PctUpper(end)), nDist_InNormal(iUnique75PctUpper,in75PctUpper)/max(nDist_InNormal(iUnique75PctUpper,in75PctUpper)),'-');
xlim([-0.1,1.1]); ylim([0,1]);
xlabel('$\bar{u}_1 / \bar{u}_{1,e}$','Interpreter','latex','FontSize',16);
ylabel('$n / n_e$','Interpreter','latex','FontSize',16);
title('Average streamwise velocity profiles - Suction side','Interpreter','latex','FontSize',16);
legend1 = legend(['$x/c$ = 0.25 (H = ',num2str(HS_s2tAvg(iUpperSide(find(sFieldUpperSide>0.25,1,'first')))),')'], ...
    ['$x/c$ = 0.50 (H = ',num2str(HS_s2tAvg(iUpperSide(find(sFieldUpperSide>0.50,1,'first')))),')'], ...
    ['$x/c$ = 0.75 (H = ',num2str(HS_s2tAvg(iUpperSide(find(sFieldUpperSide>0.75,1,'first')))),')']);
set(legend1,'FontSize',14,'Interpreter','latex');
set(gcf,'color','w');


% figure(5)
% subplot(2,2,1)
% plot(u1_InNormal_s2tAvg(iUnique25PctLower,in25PctLower)/u1_InNormal_s2tAvg(iUnique25PctLower,in25PctLower(end)), nDist_InNormal(iUnique25PctLower,in25PctLower)/max(nDist_InNormal(iUnique25PctLower,in25PctLower)),'-');
% hold on
% plot(u1_InNormal_s2tAvg(iUnique25PctUpper,in25PctUpper)/u1_InNormal_s2tAvg(iUnique25PctUpper,in25PctUpper(end)), nDist_InNormal(iUnique25PctUpper,in25PctUpper)/max(nDist_InNormal(iUnique25PctUpper,in25PctUpper)),'-');
% xlim([-0.1,1.1]); ylim([0,1]);
% xlabel('$\bar{u}_1 / \bar{u}_{1,e}$','Interpreter','latex','FontSize',16);
% ylabel('$n / n_e$','Interpreter','latex','FontSize',16);
% title('Average streamwise velocity profile at $x/c=0.25$','Interpreter','latex','FontSize',16);
% legend1 = legend(['Pressure side (H = ',num2str(HS_s2tAvg(iLowerSide(find(sFieldLowerSide<0.25,1,'first')))),')'], ['Suction side (H = ',num2str(HS_s2tAvg(iUpperSide(find(sFieldUpperSide>0.25,1,'first')))),')']);
% set(legend1,'FontSize',14,'Interpreter','latex');
% set(gcf,'color','w');
% 
% subplot(2,2,2)
% plot(u1_InNormal_s2tAvg(iUnique50PctLower,in50PctLower)/u1_InNormal_s2tAvg(iUnique50PctLower,in50PctLower(end)), nDist_InNormal(iUnique50PctLower,in50PctLower)/max(nDist_InNormal(iUnique50PctLower,in50PctLower)),'-');
% hold on
% plot(u1_InNormal_s2tAvg(iUnique50PctUpper,in50PctUpper)/u1_InNormal_s2tAvg(iUnique50PctUpper,in50PctUpper(end)), nDist_InNormal(iUnique50PctUpper,in50PctUpper)/max(nDist_InNormal(iUnique50PctUpper,in50PctUpper)),'-');
% xlim([-0.1,1.1]); ylim([0,1]);
% xlabel('$\bar{u}_1 / \bar{u}_{1,e}$','Interpreter','latex','FontSize',16);
% ylabel('$n / n_e$','Interpreter','latex','FontSize',16);
% title('Average streamwise velocity profile at $x/c=0.50$','Interpreter','latex','FontSize',16);
% legend1 = legend(['Pressure side (H = ',num2str(HS_s2tAvg(iLowerSide(find(sFieldLowerSide<0.50,1,'first')))),')'], ['Suction side (H = ',num2str(HS_s2tAvg(iUpperSide(find(sFieldUpperSide>0.50,1,'first')))),')']);
% set(legend1,'FontSize',14,'Interpreter','latex');
% set(gcf,'color','w');
% 
% subplot(2,2,3)
% plot(u1_InNormal_s2tAvg(iUnique75PctLower,in75PctLower)/u1_InNormal_s2tAvg(iUnique75PctLower,in75PctLower(end)), nDist_InNormal(iUnique75PctLower,in75PctLower)/max(nDist_InNormal(iUnique75PctLower,in75PctLower)),'-');
% hold on
% plot(u1_InNormal_s2tAvg(iUnique75PctUpper,in75PctUpper)/u1_InNormal_s2tAvg(iUnique75PctUpper,in75PctUpper(end)), nDist_InNormal(iUnique75PctUpper,in75PctUpper)/max(nDist_InNormal(iUnique75PctUpper,in75PctUpper)),'-');
% xlim([-0.1,1.1]); ylim([0,1]);
% xlabel('$\bar{u}_1 / \bar{u}_{1,e}$','Interpreter','latex','FontSize',16);
% ylabel('$n / n_e$','Interpreter','latex','FontSize',16);
% title('Average streamwise velocity profile at $x/c=0.75$','Interpreter','latex','FontSize',16);
% legend1 = legend(['Pressure side (H = ',num2str(HS_s2tAvg(iLowerSide(find(sFieldLowerSide<0.75,1,'first')))),')'], ['Suction side (H = ',num2str(HS_s2tAvg(iUpperSide(find(sFieldUpperSide>0.75,1,'first')))),')']);
% set(legend1,'FontSize',14,'Interpreter','latex');
% set(gcf,'color','w');


% figure(6)
% subplot(2,2,1)
% plot(pPert_InNormal_s2tRMS(iUnique25PctLower,in25PctLower) / pe_s2tAvg(iUnique25PctLower), nDist_InNormal(iUnique25PctLower,in25PctLower)/nDist_InNormal(iUnique25PctLower,in25PctLower(end)),'-');
% hold on
% plot(pPert_InNormal_s2tRMS(iUnique25PctUpper,in25PctUpper) / pe_s2tAvg(iUnique25PctUpper), nDist_InNormal(iUnique25PctUpper,in25PctUpper)/nDist_InNormal(iUnique25PctUpper,in25PctUpper(end)),'-');
% xlabel('$p^{''}_{RMS} / \bar{p}_e$','Interpreter','latex','FontSize',16);
% ylabel('$n / n_e$','Interpreter','latex','FontSize',16);
% title('Fluctuating pressure profile at $x/c=0.25$','Interpreter','latex','FontSize',16);
% legend1 = legend('Suction side', 'Pressure side');
% set(legend1,'FontSize',14,'Interpreter','latex');
% set(gcf,'color','w');
% subplot(2,2,2)
% plot(pPert_InNormal_s2tRMS(iUnique50PctLower,in50PctLower) / pe_s2tAvg(iUnique50PctLower), nDist_InNormal(iUnique50PctLower,in50PctLower)/nDist_InNormal(iUnique50PctLower,in50PctLower(end)),'-');
% hold on
% plot(pPert_InNormal_s2tRMS(iUnique50PctUpper,in50PctUpper) / pe_s2tAvg(iUnique50PctUpper), nDist_InNormal(iUnique50PctUpper,in50PctUpper)/nDist_InNormal(iUnique50PctUpper,in50PctUpper(end)),'-');
% xlabel('$p^{''}_{RMS} / \bar{p}_e$','Interpreter','latex','FontSize',16);
% ylabel('$n / n_e$','Interpreter','latex','FontSize',16);
% title('Fluctuating pressure profile at $x/c=0.50$','Interpreter','latex','FontSize',16);
% legend1 = legend('Suction side', 'Pressure side');
% set(legend1,'FontSize',14,'Interpreter','latex');
% set(gcf,'color','w');
% subplot(2,2,3)
% plot(pPert_InNormal_s2tRMS(iUnique75PctLower,in75PctLower) / pe_s2tAvg(iUnique75PctLower), nDist_InNormal(iUnique75PctLower,in75PctLower)/nDist_InNormal(iUnique75PctLower,in75PctLower(end)),'-');
% hold on
% plot(pPert_InNormal_s2tRMS(iUnique75PctUpper,in75PctUpper) / pe_s2tAvg(iUnique75PctUpper), nDist_InNormal(iUnique75PctUpper,in75PctUpper)/nDist_InNormal(iUnique75PctUpper,in75PctUpper(end)),'-');
% xlabel('$p^{''}_{RMS} / \bar{p}_e$','Interpreter','latex','FontSize',16);
% ylabel('$n / n_e$','Interpreter','latex','FontSize',16);
% title('Fluctuating pressure profile at $x/c=0.75$','Interpreter','latex','FontSize',16);
% legend1 = legend('Suction side', 'Pressure side');
% set(legend1,'FontSize',14,'Interpreter','latex');
% set(gcf,'color','w');
