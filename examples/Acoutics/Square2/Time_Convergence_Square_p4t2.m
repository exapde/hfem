%% All Explicit dt = 0.001
load('/home/lkolkman/Dropbox (MIT)/research_2/Error_Plots/Time_Convergence/p4_nref2/dt_0.0005_t2.mat')
% load('C:\Users\Lauren\Dropbox (MIT)\research_2\Error_Plots\Time_Convergence\p4_nref2\dt_0.0005.mat')
load('Sq_Ref_Solution_Simplified.mat')

% Calculate Error in Comparison to the Implicit Mesh
% errexp1 = calerror_IM_reference(UDGexpAll,meshEX,master,UDGimpAll_REF,meshIM_REF,time);
errexp1 = erruexp;

%% All Explicit dt = 0.001
load('/home/lkolkman/Dropbox (MIT)/research_2/Error_Plots/Time_Convergence/p4_nref2/dt_0.001_t2.mat')
% load('C:\Users\Lauren\Dropbox (MIT)\research_2\Error_Plots\Time_Convergence\p4_nref2\dt_0.001.mat')
load('Sq_Ref_Solution_Simplified.mat')

% Calculate Error in Comparison to the Implicit Mesh
% errexp1 = calerror_IM_reference(UDGexpAll,meshEX,master,UDGimpAll_REF,meshIM_REF,time);
errexp2 = erruexp;



%% dt = 0.0017
load('/home/lkolkman/Dropbox (MIT)/research_2/Error_Plots/Time_Convergence/p4_nref2/dt_0.0015_t2.mat')
% load('C:\Users\Lauren\Dropbox (MIT)\research_2\Error_Plots\Time_Convergence\p4_nref2\dt_0.0015.mat')

% Calculate Error in Comparison to the Implicit Mesh
% errexp2 = calerror_IM_reference(UDGexpAll,meshEX,master,UDGimpAll_REF,meshIM_REF,time);
% errimp2 = calerror_IM_reference(UDGimpAll,meshIM,master,UDGimpAll_REF,meshIM_REF,time);
errexp3 = erruexp;
errimp3 = erruimp;


%% dt = 0.003
load('/home/lkolkman/Dropbox (MIT)/research_2/Error_Plots/Time_Convergence/p4_nref2/dt_0.002.mat')
% load('C:\Users\Lauren\Dropbox (MIT)\research_2\Error_Plots\Time_Convergence\p4_nref2\dt_0.002.mat')

% Calculate Error in Comparison to the Implicit Mesh
% errexp3 = calerror_IM_reference(UDGexpAll,meshEX,master,UDGimpAll_REF,meshIM_REF,time);
% errimp3 = calerror_IM_reference(UDGimpAll,meshIM,master,UDGimpAll_REF,meshIM_REF,time);
errexp4 = erruexp;
errimp4 = erruimp;


%% dt = 0.0035
load('/home/lkolkman/Dropbox (MIT)/research_2/Error_Plots/Time_Convergence/p4_nref2/dt_0.0025_t2.mat')
% load('C:\Users\Lauren\Dropbox (MIT)\research_2\Error_Plots\Time_Convergence\p4_nref2\dt_0.0025.mat')

% Calculate Error in Comparison to the Implicit Mesh
% errexp4 = calerror_IM_reference(UDGexpAll,meshEX,master,UDGimpAll_REF,meshIM_REF,time);
% errimp4 = calerror_IM_reference(UDGimpAll,meshIM,master,UDGimpAll_REF,meshIM_REF,time);
errexp5 = erruexp;
errimp5 = erruimp;


%%  dt = 0.004
load('/home/lkolkman/Dropbox (MIT)/research_2/Error_Plots/Time_Convergence/p4_nref2/dt_0.003_t2.mat')
% load('C:\Users\Lauren\Dropbox (MIT)\research_2\Error_Plots\Time_Convergence\p4_nref2\dt_0.003.mat')

% Calculate Error in Comparison to the Implicit Mesh
% errimp5 = calerror_IM_reference(UDGimpAll,meshIM,master,UDGimpAll_REF,meshIM_REF,time);

errexp6 = erruexp;
errimp6 = erruimp;

%% 
load('/home/lkolkman/Dropbox (MIT)/research_2/Error_Plots/Time_Convergence/p4_nref2/dt_0.0035_t2.mat')
% load('C:\Users\Lauren\Dropbox (MIT)\research_2\Error_Plots\Time_Convergence\p4_nref2\dt_0.0035.mat')

% Calculate Error in Comparison to the Implicit Mesh
% errimp5 = calerror_IM_reference(UDGimpAll,meshIM,master,UDGimpAll_REF,meshIM_REF,time);
errexp7 = erruexp;
errimp7 = erruimp;

%% All Implicit dt = 0.016
load('/home/lkolkman/Dropbox (MIT)/research_2/Error_Plots/Time_Convergence/p4_nref2/dt_0.004_t2.mat')
% load('C:\Users\Lauren\Dropbox (MIT)\research_2\Error_Plots\Time_Convergence\p4_nref2\dt_0.004.mat')

% Calculate Error in Comparison to the Implicit Mesh
% errimp5 = calerror_IM_reference(UDGimpAll,meshIM,master,UDGimpAll_REF,meshIM_REF,time);
errimp8 = erruimp;

%% All Implicit dt = 0.016
load('/home/lkolkman/Dropbox (MIT)/research_2/Error_Plots/Time_Convergence/p4_nref2/dt_0.008_t2.mat')
% load('C:\Users\Lauren\Dropbox (MIT)\research_2\Error_Plots\Time_Convergence\p4_nref2\dt_0.008.mat')

% Calculate Error in Comparison to the Implicit Mesh
% errimp5 = calerror_IM_reference(UDGimpAll,meshIM,master,UDGimpAll_REF,meshIM_REF,time);
errimp9 = erruimp;

%% Computed Error
err1 = errexp1;
err2 = errexp2;
err3 = sqrt(errexp3.^2 + errimp3.^2);
err4 = sqrt(errexp4.^2 + errimp4.^2);
err5 = sqrt(errexp5.^2 + errimp5.^2);
err6 = sqrt(errexp6.^2 + errimp6.^2);
err7 = sqrt(errexp7.^2 + errimp7.^2);
err8 = errimp8;
err9 = errimp9;

All_Error = [err1, err2, err3, err4, err5 err6 err7 err8 err9];

V = All_Error(1,:);
Q = sqrt(All_Error(2,:).^2 + All_Error(3,:).^2);
U = All_Error(4,:);


%% Plot

dt = [0.0005; 0.001; 0.0015; 0.002; 0.0025; 0.003; 0.0035; 0.004; 0.008];

%%%%%%%%%%%%%%%%%% v Plot %%%%%%%%%%%%%%%%%%
figure()
loglog(dt,All_Error(1,:),'-s');

title('v convergence');
xlabel('log(dt)');
ylabel('log(Error)');
axis tight

Vm = zeros(8,1);

Vm(1) = (log10(V(1)) - log10(V(2)))/(log10(dt(1)) - log10(dt(2)));
Vm(2) = (log10(V(2)) - log10(V(3)))/(log10(dt(2)) - log10(dt(3)));
Vm(3) = (log10(V(3)) - log10(V(4)))/(log10(dt(3)) - log10(dt(4)));
Vm(4) = (log10(V(4)) - log10(V(5)))/(log10(dt(4)) - log10(dt(5)));
Vm(5) = (log10(V(5)) - log10(V(6)))/(log10(dt(5)) - log10(dt(6)));
Vm(6) = (log10(V(6)) - log10(V(7)))/(log10(dt(6)) - log10(dt(7)));
Vm(7) = (log10(V(7)) - log10(V(8)))/(log10(dt(7)) - log10(dt(8)));
Vm(8) = (log10(V(8)) - log10(V(9)))/(log10(dt(8)) - log10(dt(9)));

%%%%%%%%%%%%%%%%%% u Plot %%%%%%%%%%%%%%%%%%
figure()
loglog(dt,All_Error(4,:),'-s');

title('u convergence');
xlabel('log(dt)');
ylabel('log(Error)');
axis tight

Um = zeros(8,1);

Um(1) = (log10(U(1)) - log10(U(2)))/(log10(dt(1)) - log10(dt(2)));
Um(2) = (log10(U(2)) - log10(U(3)))/(log10(dt(2)) - log10(dt(3)));
Um(3) = (log10(U(3)) - log10(U(4)))/(log10(dt(3)) - log10(dt(4)));
Um(4) = (log10(U(4)) - log10(U(5)))/(log10(dt(4)) - log10(dt(5)));
Um(5) = (log10(U(5)) - log10(U(6)))/(log10(dt(5)) - log10(dt(6)));
Um(6) = (log10(U(6)) - log10(U(7)))/(log10(dt(6)) - log10(dt(7)));
Um(7) = (log10(U(7)) - log10(U(8)))/(log10(dt(7)) - log10(dt(8)));
Um(8) = (log10(U(8)) - log10(U(9)))/(log10(dt(8)) - log10(dt(9)));

%%%%%%%%%%%%%%%%%% u Plot %%%%%%%%%%%%%%%%%%
figure()
loglog(dt,Q,'-s');

title('q convergence');
xlabel('log(dt)');
ylabel('log(Error)');
axis tight

Qm = zeros(8,1);

Qm(1) = (log10(Q(1)) - log10(Q(2)))/(log10(dt(1)) - log10(dt(2)));
Qm(2) = (log10(Q(2)) - log10(Q(3)))/(log10(dt(2)) - log10(dt(3)));
Qm(3) = (log10(Q(3)) - log10(Q(4)))/(log10(dt(3)) - log10(dt(4)));
Qm(4) = (log10(Q(4)) - log10(Q(5)))/(log10(dt(4)) - log10(dt(5)));
Qm(5) = (log10(Q(5)) - log10(Q(6)))/(log10(dt(5)) - log10(dt(6)));
Qm(6) = (log10(Q(6)) - log10(Q(7)))/(log10(dt(6)) - log10(dt(7)));
Qm(7) = (log10(Q(7)) - log10(Q(8)))/(log10(dt(7)) - log10(dt(8)));
Qm(8) = (log10(Q(8)) - log10(Q(9)))/(log10(dt(8)) - log10(dt(9)));

