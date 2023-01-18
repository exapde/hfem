%% All Explicit dt = 0.0004
load('/home/lkolkman/Dropbox (MIT)/research_2/Error_Plots/Lmesh/Time/dt_0.0004.mat')
% load('C:\Users\Lauren\Dropbox (MIT)\research_2\Error_Plots\Lmesh\Time\dt_0.0004.mat')
load('L_ref_Solution1.1_Simplified.mat')

UDGimpAll_REF = [UDGimp_REF, PDGimp_REF];

% Calculate Error in Comparison to the Implicit Mesh
UDGexpAll = [UDGexp, PDGexp];
UDGimpAll = [UDGimp, PDGimp];
errexp1 = calerror_IM_reference(UDGexpAll,meshEX,master,UDGimpAll_REF,meshIM_REF,time);


%% dt = 0.0009
load('/home/lkolkman/Dropbox (MIT)/research_2/Error_Plots/Lmesh/Time/dt_0.0009.mat')
% load('C:\Users\Lauren\Dropbox (MIT)\research_2\Error_Plots\Lmesh\Time\dt_0.0009.mat')

% Calculate Error in Comparison to the Implicit Mesh
UDGexpAll = [UDGexp, PDGexp];
UDGimpAll = [UDGimp, PDGimp];
errexp2 = calerror_IM_reference(UDGexpAll,meshEX,master,UDGimpAll_REF,meshIM_REF,time);
errimp2 = calerror_IM_reference(UDGimpAll,meshIM,master,UDGimpAll_REF,meshIM_REF,time);



%% dt = 0.0013
load('/home/lkolkman/Dropbox (MIT)/research_2/Error_Plots/Lmesh/Time/dt_0.0013_run2.mat')
% load('C:\Users\Lauren\Dropbox (MIT)\research_2\Error_Plots\Lmesh\Time\dt_0.0013_run2.mat')

% Calculate Error in Comparison to the Implicit Mesh
UDGexpAll = [UDGexp, PDGexp];
UDGimpAll = [UDGimp, PDGimp];
errexp3 = calerror_IM_reference(UDGexpAll,meshEX,master,UDGimpAll_REF,meshIM_REF,time);
errimp3 = calerror_IM_reference(UDGimpAll,meshIM,master,UDGimpAll_REF,meshIM_REF,time);


%% dt = 0.003
load('/home/lkolkman/Dropbox (MIT)/research_2/Error_Plots/Lmesh/Time/dt_0.003.mat')
% load('C:\Users\Lauren\Dropbox (MIT)\research_2\Error_Plots\Lmesh\Time\dt_0.003.mat')

% Calculate Error in Comparison to the Implicit Mesh
UDGexpAll = [UDGexp, PDGexp];
UDGimpAll = [UDGimp, PDGimp];
errexp4 = calerror_IM_reference(UDGexpAll,meshEX,master,UDGimpAll_REF,meshIM_REF,time);
errimp4 = calerror_IM_reference(UDGimpAll,meshIM,master,UDGimpAll_REF,meshIM_REF,time);


%% All Implicit dt = 0.004
load('/home/lkolkman/Dropbox (MIT)/research_2/Error_Plots/Lmesh/Time/dt_0.004_test.mat')
% load('C:\Users\Lauren\Dropbox (MIT)\research_2\Error_Plots\Lmesh\Time\dt_0.004_test.mat')

% Calculate Error in Comparison to the Implicit Mesh
UDGexpAll = [UDGexp, PDGexp];
UDGimpAll = [UDGimp, PDGimp];
errimp5 = calerror_IM_reference(UDGimpAll,meshIM,master,UDGimpAll_REF,meshIM_REF,time);


%% Computed Error
err1 = errexp1;
err2 = sqrt(errexp2.^2 + errimp2.^2);
err3 = sqrt(errexp3.^2 + errimp3.^2);
err4 = sqrt(errexp4.^2 + errimp4.^2);
err5 = errimp5;


All_Error = [err1, err2, err3, err4, err5];

V = All_Error(1,:);
Q = sqrt(All_Error(2,:).^2 + All_Error(3,:).^2);
U = All_Error(4,:);


%% Plot

dt = [0.0004; 0.0009; 0.0013; 0.003; 0.004];

%%%%%%%%%%%%%%%%%% v Plot %%%%%%%%%%%%%%%%%%
figure()
loglog(dt,All_Error(1,:),'-rs');

title('v convergence');
xlabel('log(dt)');
ylabel('log(Error)');
axis tight

Vm = zeros(4,1);

Vm(1) = (log10(V(1)) - log10(V(2)))/(log10(dt(1)) - log10(dt(2)));
Vm(2) = (log10(V(2)) - log10(V(3)))/(log10(dt(2)) - log10(dt(3)));
Vm(3) = (log10(V(3)) - log10(V(4)))/(log10(dt(3)) - log10(dt(4)));
Vm(4) = (log10(V(4)) - log10(V(5)))/(log10(dt(4)) - log10(dt(5)));

%%%%%%%%%%%%%%%%%% u Plot %%%%%%%%%%%%%%%%%%
figure()
loglog(dt,All_Error(4,:),'-s');

title('u convergence');
xlabel('log(dt)');
ylabel('log(Error)');
axis tight

Um = zeros(4,1);

Um(1) = (log10(U(1)) - log10(U(2)))/(log10(dt(1)) - log10(dt(2)));
Um(2) = (log10(U(2)) - log10(U(3)))/(log10(dt(2)) - log10(dt(3)));
Um(3) = (log10(U(3)) - log10(U(4)))/(log10(dt(3)) - log10(dt(4)));
Um(4) = (log10(U(4)) - log10(U(5)))/(log10(dt(4)) - log10(dt(5)));

%%%%%%%%%%%%%%%%%% q Plot %%%%%%%%%%%%%%%%%%
figure()
loglog(dt,Q,'-ms');

title('q convergence');
xlabel('log(dt)');
ylabel('log(Error)');
axis tight

Qm = zeros(4,1);

Qm(1) = (log10(Q(1)) - log10(Q(2)))/(log10(dt(1)) - log10(dt(2)));
Qm(2) = (log10(Q(2)) - log10(Q(3)))/(log10(dt(2)) - log10(dt(3)));
Qm(3) = (log10(Q(3)) - log10(Q(4)))/(log10(dt(3)) - log10(dt(4)));
Qm(4) = (log10(Q(4)) - log10(Q(5)))/(log10(dt(4)) - log10(dt(5)));


