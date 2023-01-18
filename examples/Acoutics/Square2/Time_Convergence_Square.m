%% All Explicit dt = 0.0005
load('/home/lkolkman/Dropbox (MIT)/research_2/Error_Plots/Time_Convergence/p3nref8_dt0.0005.mat')
% load('C:\Users\Lauren\Dropbox (MIT)\research_2\Error_Plots\Time_Convergence\p3nref8_dt0.0005.mat')
load('Sq_Ref_Solution_Simplified.mat')

% Calculate Error in Comparison to the Implicit Mesh
% errexp1 = calerror_IM_reference(UDGexpAll,meshEX,master,UDGimpAll_REF,meshIM_REF,time);
errexp1 = erruexp;

%% dt = 0.0009
load('/home/lkolkman/Dropbox (MIT)/research_2/Error_Plots/Time_Convergence/p3nref8_dt0.0009.mat')
% load('C:\Users\Lauren\Dropbox (MIT)\research_2\Error_Plots\Time_Convergence\p3nref8_dt0.0009.mat')

% Calculate Error in Comparison to the Implicit Mesh
% errexp2 = calerror_IM_reference(UDGexpAll,meshEX,master,UDGimpAll_REF,meshIM_REF,time);
% errimp2 = calerror_IM_reference(UDGimpAll,meshIM,master,UDGimpAll_REF,meshIM_REF,time);
errexp2 = erruexp;
errimp2 = erruimp;


%% dt = 0.0013
load('/home/lkolkman/Dropbox (MIT)/research_2/Error_Plots/Time_Convergence/p3nref8_dt0.0013.mat')
% load('C:\Users\Lauren\Dropbox (MIT)\research_2\Error_Plots\Time_Convergence\p3nref8_dt0.0013.mat')

% Calculate Error in Comparison to the Implicit Mesh
% errexp3 = calerror_IM_reference(UDGexpAll,meshEX,master,UDGimpAll_REF,meshIM_REF,time);
% errimp3 = calerror_IM_reference(UDGimpAll,meshIM,master,UDGimpAll_REF,meshIM_REF,time);
errexp3 = erruexp;
errimp3 = erruimp;


%% dt = 0.0014
load('/home/lkolkman/Dropbox (MIT)/research_2/Error_Plots/Time_Convergence/p3nref8_dt0.0014.mat')
% load('C:\Users\Lauren\Dropbox (MIT)\research_2\Error_Plots\Time_Convergence\p3nref8_dt0.0014.mat')

% Calculate Error in Comparison to the Implicit Mesh
% errexp4 = calerror_IM_reference(UDGexpAll,meshEX,master,UDGimpAll_REF,meshIM_REF,time);
% errimp4 = calerror_IM_reference(UDGimpAll,meshIM,master,UDGimpAll_REF,meshIM_REF,time);
errexp4 = erruexp;
errimp4 = erruimp;


%% All Implicit dt = 0.0016
load('/home/lkolkman/Dropbox (MIT)/research_2/Error_Plots/Time_Convergence/p3nref8_dt0.0016_run2.mat')
% load('C:\Users\Lauren\Dropbox (MIT)\research_2\Error_Plots\Time_Convergence\p3nref8_dt0.0016.mat')

% Calculate Error in Comparison to the Implicit Mesh
% errimp5 = calerror_IM_reference(UDGimpAll,meshIM,master,UDGimpAll_REF,meshIM_REF,time);
errimp5 = erruimp;


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

dt = [0.0005; 0.0009; 0.0013; 0.0014; 0.0016];

%%%%%%%%%%%%%%%%%% v Plot %%%%%%%%%%%%%%%%%%
figure()
loglog(dt,All_Error(1,:),'-s');

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

%%%%%%%%%%%%%%%%%% u Plot %%%%%%%%%%%%%%%%%%
figure()
loglog(dt,Q,'-s');

title('q convergence');
xlabel('log(dt)');
ylabel('log(Error)');
axis tight

Qm = zeros(4,1);

Qm(1) = (log10(Q(1)) - log10(Q(2)))/(log10(dt(1)) - log10(dt(2)));
Qm(2) = (log10(Q(2)) - log10(Q(3)))/(log10(dt(2)) - log10(dt(3)));
Qm(3) = (log10(Q(3)) - log10(Q(4)))/(log10(dt(3)) - log10(dt(4)));
Qm(4) = (log10(Q(4)) - log10(Q(5)))/(log10(dt(4)) - log10(dt(5)));
