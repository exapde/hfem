% %% All Explicit dt = 0.001
% load('/home/lkolkman/Dropbox (MIT)/DIGASO/err_p4_nref1_K05_dt0001.mat')
load('C:\Users\Lauren\Dropbox (MIT)\DIGASO\err_p4_nref1_K05_dt0001.mat')
% load('Sq_Ref_Solution_Simplified.mat')
% 
% Calculate Error in Comparison to the Implicit Mesh
% errexp1 = calerror_IM_reference(UDGexpAll,meshEX,master,UDGimpAll_REF,meshIM_REF,time);
errexp1 = erruexp;
errimp1 = erruimp;

%% All Explicit dt = 0.001
% load('/home/lkolkman/Dropbox (MIT)/DIGASO/err_p4_nref1_K05_dt0004.mat')
load('C:\Users\Lauren\Dropbox (MIT)\DIGASO\err_p4_nref1_K05_dt0004.mat')
% load('Sq_Ref_Solution_Simplified.mat')

% Calculate Error in Comparison to the Implicit Mesh
% errexp1 = calerror_IM_reference(UDGexpAll,meshEX,master,UDGimpAll_REF,meshIM_REF,time);
errexp2 = erruexp;
errimp2 = erruimp;



%% dt = 0.0017
load('/home/lkolkman/Dropbox (MIT)/DIGASO/err_p4_nref1_K05_dt0005.mat')
% load('C:\Users\Lauren\Dropbox (MIT)\DIGASO\err_p4_nref1_K05_dt0005.mat')

% Calculate Error in Comparison to the Implicit Mesh
% errexp2 = calerror_IM_reference(UDGexpAll,meshEX,master,UDGimpAll_REF,meshIM_REF,time);
% errimp2 = calerror_IM_reference(UDGimpAll,meshIM,master,UDGimpAll_REF,meshIM_REF,time);
errexp3 = erruexp;
errimp3 = erruimp;


%% dt = 0.003
% load('/home/lkolkman/Dropbox (MIT)/DIGASO/err_p4_nref1_K05_dt0007.mat')
load('C:\Users\Lauren\Dropbox (MIT)\DIGASO\err_p4_nref1_K05_dt0007.mat')

% Calculate Error in Comparison to the Implicit Mesh
% errexp3 = calerror_IM_reference(UDGexpAll,meshEX,master,UDGimpAll_REF,meshIM_REF,time);
% errimp3 = calerror_IM_reference(UDGimpAll,meshIM,master,UDGimpAll_REF,meshIM_REF,time);
errexp4 = erruexp;
errimp4 = erruimp;


%% dt = 0.0035
% load('/home/lkolkman/Dropbox (MIT)/DIGASO/err_p4_nref1_K05_dt0010.mat')
load('C:\Users\Lauren\Dropbox (MIT)\DIGASO\err_p4_nref1_K05_dt0010.mat')

% Calculate Error in Comparison to the Implicit Mesh
% errexp4 = calerror_IM_reference(UDGexpAll,meshEX,master,UDGimpAll_REF,meshIM_REF,time);
% errimp4 = calerror_IM_reference(UDGimpAll,meshIM,master,UDGimpAll_REF,meshIM_REF,time);
errimp5 = erruimp;


%%  dt = 0.004
% load('/home/lkolkman/Dropbox (MIT)/DIGASO/err_p4_nref1_K05_dt0020.mat')
load('C:\Users\Lauren\Dropbox (MIT)\DIGASO\err_p4_nref1_K05_dt0020.mat')

% Calculate Error in Comparison to the Implicit Mesh
% errimp5 = calerror_IM_reference(UDGimpAll,meshIM,master,UDGimpAll_REF,meshIM_REF,time);

errimp6 = erruimp;


%% Computed Error
% err1 = sqrt(errexp1.^2 + errimp1.^2);
err2 = sqrt(errexp2.^2 + errimp2.^2);
err3 = sqrt(errexp3.^2 + errimp3.^2);
err4 = sqrt(errexp4.^2 + errimp4.^2);
err5 = errimp5;
err6 = errimp6;


All_Error = [err2, err3, err4, err5 err6];

V = All_Error(1,:);
Q = sqrt(All_Error(2,:).^2 + All_Error(3,:).^2);
U = All_Error(4,:);


%% Plot

dt = [0.0001; 0.0004; 0.0005; 0.0007; 0.0010; 0.0020];
dt2 = [0.0004; 0.0005; 0.0007; 0.0010; 0.0020];

dt = dt2;
%%%%%%%%%%%%%%%%%% v Plot %%%%%%%%%%%%%%%%%%
figure()
loglog(dt2,All_Error(1,:),'-rs');

title('v convergence');
xlabel('log(dt)');
ylabel('log(Error)');
axis tight

Vm = zeros(8,1);

Vm(1) = (log10(V(1)) - log10(V(2)))/(log10(dt(1)) - log10(dt(2)));
Vm(2) = (log10(V(2)) - log10(V(3)))/(log10(dt(2)) - log10(dt(3)));
Vm(3) = (log10(V(3)) - log10(V(4)))/(log10(dt(3)) - log10(dt(4)));
Vm(4) = (log10(V(4)) - log10(V(5)))/(log10(dt(4)) - log10(dt(5)));
% Vm(5) = (log10(V(5)) - log10(V(6)))/(log10(dt(5)) - log10(dt(6)));
% Vm(6) = (log10(V(6)) - log10(V(7)))/(log10(dt(6)) - log10(dt(7)));
% Vm(7) = (log10(V(7)) - log10(V(8)))/(log10(dt(7)) - log10(dt(8)));
% Vm(8) = (log10(V(8)) - log10(V(9)))/(log10(dt(8)) - log10(dt(9)));

%%%%%%%%%%%%%%%%%% u Plot %%%%%%%%%%%%%%%%%%
figure()
loglog(dt2,All_Error(4,:),'-s');

title('u convergence');
xlabel('log(dt)');
ylabel('log(Error)');
axis tight

Um = zeros(8,1);

Um(1) = (log10(U(1)) - log10(U(2)))/(log10(dt(1)) - log10(dt(2)));
Um(2) = (log10(U(2)) - log10(U(3)))/(log10(dt(2)) - log10(dt(3)));
Um(3) = (log10(U(3)) - log10(U(4)))/(log10(dt(3)) - log10(dt(4)));
Um(4) = (log10(U(4)) - log10(U(5)))/(log10(dt(4)) - log10(dt(5)));
% Um(5) = (log10(U(5)) - log10(U(6)))/(log10(dt(5)) - log10(dt(6)));
% Um(6) = (log10(U(6)) - log10(U(7)))/(log10(dt(6)) - log10(dt(7)));
% Um(7) = (log10(U(7)) - log10(U(8)))/(log10(dt(7)) - log10(dt(8)));
% Um(8) = (log10(U(8)) - log10(U(9)))/(log10(dt(8)) - log10(dt(9)));

%%%%%%%%%%%%%%%%%% u Plot %%%%%%%%%%%%%%%%%%
figure()
loglog(dt2,Q,'-s','Color',[148/255,0,211/255]);

title('q convergence');
xlabel('log(dt)');
ylabel('log(Error)');
axis tight

Qm = zeros(8,1);

Qm(1) = (log10(Q(1)) - log10(Q(2)))/(log10(dt(1)) - log10(dt(2)));
Qm(2) = (log10(Q(2)) - log10(Q(3)))/(log10(dt(2)) - log10(dt(3)));
Qm(3) = (log10(Q(3)) - log10(Q(4)))/(log10(dt(3)) - log10(dt(4)));
Qm(4) = (log10(Q(4)) - log10(Q(5)))/(log10(dt(4)) - log10(dt(5)));
% Qm(5) = (log10(Q(5)) - log10(Q(6)))/(log10(dt(5)) - log10(dt(6)));
% Qm(6) = (log10(Q(6)) - log10(Q(7)))/(log10(dt(6)) - log10(dt(7)));
% Qm(7) = (log10(Q(7)) - log10(Q(8)))/(log10(dt(7)) - log10(dt(8)));
% Qm(8) = (log10(Q(8)) - log10(Q(9)))/(log10(dt(8)) - log10(dt(9)));

