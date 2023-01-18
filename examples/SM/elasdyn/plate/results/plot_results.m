
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        PLOTTING CONVERGENCE FOR k=1,2,3            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Width default
set(0,'defaultlinelinewidth',1.5)

% Setting slopes triangles
triX = [1, 2, 2, 1]; 
tri1 = [1, 1./2, 1, 1]; tri2 = [1, 1./4, 1, 1];
tri3 = [1, 1./8, 1, 1]; tri4 = [1, 1./16,1, 1];
TXTx = [1.4 2.1]; TXTy = [1.5 0.4];

% Loading data
load('/home/terrana/Codes/DIGASO/problem/SM/elasdyn/plate/results/C01_k1_t2.mat')
h1 = L2errors.h; U1 = L2errors.Uerror(:,1); Q1 = L2errors.Qerror(:,1);
load('/home/terrana/Codes/DIGASO/problem/SM/elasdyn/plate/results/C01_k2_t2.mat')
h2 = L2errors.h; U2 = L2errors.Uerror(:,1); Q2 = L2errors.Qerror(:,1);
load('/home/terrana/Codes/DIGASO/problem/SM/elasdyn/plate/results/C01_k3_t2.mat')
h3 = L2errors.h; U3 = L2errors.Uerror(:,1); Q3 = L2errors.Qerror(:,1);

% Convergence plot - Velocities
figure(2); clf; loglog(1./h1,U1,'-s',1./h2,U2,'-s',1./h3,U3,'-s'); grid on; xlim([1.,40.]);
title('Convergence for error on velocities v'); xlabel('1/h'); ylabel('||v-v_h||')
% Positioning slope triangles
hold on; Cx = 12. ; Cy = 1.e-3; 
loglog(Cx*triX,Cy*tri2,'k','LineWidth',1); text(Cx*TXTx,Cy*TXTy,{'1','2'});
hold on; Cx = 8. ; Cy = 3.5e-5; 
loglog(Cx*triX,Cy*tri3,'k','LineWidth',1); text(Cx*TXTx,Cy*TXTy,{'1','3'});
hold on; Cx = 7. ; Cy = 2.5e-6; 
loglog(Cx*triX,Cy*tri4,'k','LineWidth',1); text(Cx*TXTx,Cy*TXTy,{'1','4'});
legend('k=1','k=2','k=3','Location','southwest');

% Convergence plot - Deformation Gradient
figure(3); clf; loglog(1./h1,Q1,'-s',1./h2,Q2,'-s',1./h3,Q3,'-s'); grid on; xlim([1.,40.]);
title('Convergence for error on F'); xlabel('1/h'); ylabel('||F-F_h||')
% Positioning slope triangles
hold on; Cx = 10. ; Cy = 1.e-2; 
loglog(Cx*triX,Cy*tri2,'k','LineWidth',1); text(Cx*TXTx,Cy*TXTy,{'1','2'});
hold on; Cx = 8. ; Cy = 1.5e-4; 
loglog(Cx*triX,Cy*tri3,'k','LineWidth',1); text(Cx*TXTx,Cy*TXTy,{'1','3'});
hold on; Cx = 6. ; Cy = 1.5e-5; 
loglog(Cx*triX,Cy*tri4,'k','LineWidth',1); text(Cx*TXTx,Cy*TXTy,{'1','4'});
legend('k=1','k=2','k=3','Location','southwest');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        PLOTTING FOR A GIVEN k, AND VARYING h       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setting slopes triangles (Lower)
triX = [1, 2, 1, 1]; 
tri1 = [1, 1, 2, 1]; tri2 = [1, 1, 4, 1];
tri3 = [1, 1, 8, 1]; tri4 = [1, 1,16, 1];
TXTx = [1.4 0.9]; TXTy = [0.7 2];

% Loading data
load('/home/terrana/Codes/DIGASO/problem/SM/elasdyn/plate/results/C01_k3_t2_h01.mat')
h1 = L2errors.h; U1 = L2errors.Uerror(:,1); Q1 = L2errors.Qerror(:,1);
load('/home/terrana/Codes/DIGASO/problem/SM/elasdyn/plate/results/C01_k3_t2_h001.mat')
h2 = L2errors.h; U2 = L2errors.Uerror(:,1); Q2 = L2errors.Qerror(:,1);
load('/home/terrana/Codes/DIGASO/problem/SM/elasdyn/plate/results/C01_k3_t2_h0001.mat')
h3 = L2errors.h; U3 = L2errors.Uerror(:,1); Q3 = L2errors.Qerror(:,1);
load('/home/terrana/Codes/DIGASO/problem/SM/elasdyn/plate/results/C01_k3_t2_h00001.mat')
h4 = L2errors.h; U4 = L2errors.Uerror(:,1); Q4 = L2errors.Qerror(:,1);
load('/home/terrana/Codes/DIGASO/problem/SM/elasdyn/plate/results/C01_k3_t2_h000001.mat')
h5 = L2errors.h; U5 = L2errors.Uerror(:,1); Q5 = L2errors.Qerror(:,1);
load('/home/terrana/Codes/DIGASO/problem/SM/elasdyn/plate/results/C01_k3_t2_h0000001.mat')
h6 = L2errors.h; U6 = L2errors.Uerror(:,1); Q6 = L2errors.Qerror(:,1);


% Convergence plot - Velocities
figure(2); clf; loglog(1./h1,U1,'-s',1./h2,U2,'-s',1./h3,U3,'-s',...
    1./h4,U4,'-s',1./h5,U5,'-s',1./h6,U6,'-s'); grid on; xlim([1.,25.]);
title('Convergence of v for different thicknesses (k=1)'); xlabel('1/h'); ylabel('||v-v_h||')
% Positioning slope triangles
%hold on; Cx = 6. ; Cy = 6.e-6; % (for k=1)
%hold on; Cx = 6. ; Cy = 6.e-8; % (for k=2)
%hold on; Cx = 2.5 ; Cy = 3.e-8; % (for k=3)
hold on; Cx = 5.3 ; Cy = 3.e-6; % (for k=3)
loglog(Cx*triX,Cy*tri2,'k','LineWidth',1); text(Cx*TXTx,Cy*TXTy,{'1','2'});
legend('t=1.e-1','t=1.e-2','t=1.e-3','t=1.e-4','t=1.e-5','t=1.e-6','Location','southwest');

% Convergence plot - Deformation Gradient
figure(3); clf; loglog(1./h1,Q1,'-s',1./h2,Q2,'-s',1./h3,Q3,'-s',...
    1./h4,Q4,'-s',1./h5,Q5,'-s',1./h6,Q6,'-s'); grid on; xlim([1.,25.]);
title('Convergence of F for different thicknesses (k=1)'); xlabel('1/h'); ylabel('||F-F_h||')
% Positioning slope triangles
%hold on; Cx = 4. ; Cy = 8.e-5; % (for k=1)
%hold on; Cx = 4. ; Cy = 8.e-7; % (for k=2)
hold on; Cx = 3. ; Cy = 1.e-6; % (for k=3)
%hold on; Cx = 8. ; Cy = 3.e-2; % (for k=3)
loglog(Cx*triX,Cy*tri2,'k','LineWidth',1); text(Cx*TXTx,Cy*TXTy,{'1','2'});
legend('t=1.e-1','t=1.e-2','t=1.e-3','t=1.e-4','t=1.e-5','t=1.e-6','Location','southwest');


hold on; Cx = 5.5 ; Cy = 2.5e-3; 
loglog(Cx*triX,Cy*tri3,'k','LineWidth',1); text(Cx*TXTx,Cy*TXTy,{'1','3'});
