%% P = 1 Errors

% v and q
E1_im = zeros(3,5);
E1_ex = zeros(3,5);

E1_ex(:,1) = [0.0159657422380087;0.0119027806584292;0.0119027806584292]; %nref = 1
E1_ex(:,2) = [0.00230060365836086;0.00285434280677394;0.00285434280677393]; %nref = 2
E1_ex(:,3) = [0.000380565372435871;0.000685692968094465;0.000685692968094475]; %nref = 4
E1_ex(:,4) = [7.92594237896768e-05;0.000167680834333829;0.000167680834333867]; %nref = 8
E1_ex(:,5) = [1.884580365455571e-05;4.161118468376498e-05;4.161118468373459e-05]; %nref = 16

E1_im(:,1) = [0.000925810574525559;0.00110657420286225;0.00110657420286219]; %nref = 1
E1_im(:,2) = [9.53949044269054e-05;0.000262607114452550;0.000262607114452549]; %nref = 2
E1_im(:,3) = [1.47712122728205e-05;5.97005342172845e-05;5.97005342172732e-05]; %nref = 4
E1_im(:,4) = [5.47975519410531e-06;1.33744852787009e-05;1.33744852787045e-05]; %nref = 8
E1_im(:,5) = [1.187977758299471e-06;3.098136151677355e-06;3.098136151672767e-06]; %nref = 16


% vstar

Vstar1_im = zeros(1,5);
Vstar1_ex = zeros(1,5);

Vstar1_ex(:,1) = 0.016514483450718; %nref = 1
Vstar1_ex(:,2) = 0.002176774221594; %nref = 2
Vstar1_ex(:,3) = 2.749013154256933e-04; %nref = 4
Vstar1_ex(:,4) = 3.476420675100419e-05; %nref = 8
Vstar1_ex(:,5) = 4.515771155457245e-06; %nref = 16

Vstar1_im(:,1) = 0.001014087840031; %nref = 1
Vstar1_im(:,2) = 8.990515379765029e-05; %nref = 2
Vstar1_im(:,3) = 5.968534312706003e-06; %nref = 4
Vstar1_im(:,4) = 4.218775624711878e-06; %nref = 8
Vstar1_im(:,5) = 8.296576226000192e-07; %nref = 16

% U

U1_im = zeros(1,5);
U1_ex = zeros(1,5);

U1_ex(:,1) = [0.00199960830608309]; %nref = 1
U1_ex(:,2) = [0.000299864333661615]; %nref = 2
U1_ex(:,3) = [4.58408795953520e-05]; %nref = 4
U1_ex(:,4) = [8.13605518908165e-06]; %nref = 8
U1_ex(:,5) = [1.754762510959323e-06]; %nref = 16

U1_im(:,1) = [5.53826930473025e-05]; %nref = 1
U1_im(:,2) = [8.62492171439193e-06]; %nref = 2
U1_im(:,3) = [1.98715706094116e-06]; %nref = 4
U1_im(:,4) = [4.34401688619851e-07]; %nref = 8
U1_im(:,5) = [9.770737449925677e-08]; %nref = 16

% Ustar

Ustar1_im = zeros(1,4);
Ustar1_ex = zeros(1,4);

Ustar1_ex(:,1) = 0.002086125785247; %nref = 1
Ustar1_ex(:,2) = 2.992146959705180e-04; %nref = 2
Ustar1_ex(:,3) = 4.025711256999723e-05; %nref = 4
Ustar1_ex(:,4) = 5.192132558236530e-06; %nref = 8
Ustar1_ex(:,5) = 6.579133570594379e-07; %nref = 16

Ustar1_im(:,1) = 3.847190631128837e-05; %nref = 1
Ustar1_im(:,2) = 4.765698100929722e-06; %nref = 2
Ustar1_im(:,3) = 1.146336529286778e-06; %nref = 4
Ustar1_im(:,4) = 2.008721461787488e-07; %nref = 8
Ustar1_im(:,5) = 2.770800296632217e-08; %nref = 16

%% P = 2 Errors

% v and q
E2_im = zeros(3,4);
E2_ex = zeros(3,4);


E2_ex(:,1) = [0.000377407453428187;0.00128517331338548;0.00128517331338536]; %nref = 1
E2_ex(:,2) = [4.60270524970059e-05;0.000138061537472420;0.000138061537472370]; %nref = 2
E2_ex(:,3) = [6.20516571749317e-06;1.61966700411837e-05;1.61966700412782e-05]; %nref = 4
E2_ex(:,4) = [8.07876148800009e-07;1.96759977584070e-06;1.96759977681994e-06]; %nref = 8

E2_im(:,1) = [7.05885472292468e-05;6.72097428432049e-05;6.72097428432408e-05]; %nref = 1
E2_im(:,2) = [2.26119950724244e-06;7.50348567773919e-06;7.50348567768770e-06]; %nref = 2
E2_im(:,3) = [2.53941263356812e-07;7.37246547022843e-07;7.37246547072749e-07]; %nref = 4
E2_im(:,4) = [2.71033135219596e-08;7.93789264303954e-08;7.93789264620394e-08]; %nref = 8


% vstar

Vstar2_im = zeros(1,4);
Vstar2_ex = zeros(1,4);

Vstar2_ex(:,1) = 0.000270023174821213; %nref = 1
Vstar2_ex(:,2) = 1.418612167951484e-05; %nref = 2
Vstar2_ex(:,3) = 8.475741238675893e-07; %nref = 4
Vstar2_ex(:,4) = 5.619165969168567e-08; %nref = 8

Vstar2_im(:,1) = 7.00906294798797e-05; %nref = 1
Vstar2_im(:,2) = 1.613164189409937e-06; %nref = 2
Vstar2_im(:,3) = 1.463602095513888e-07; %nref = 4
Vstar2_im(:,4) = 7.138252959514678e-09; %nref = 8

% U

U2_im = zeros(1,4);
U2_ex = zeros(1,4);

U2_ex(:,1) = 5.48259425802488e-05; %nref = 1
U2_ex(:,2) = 4.58914775494636e-06; %nref = 2
U2_ex(:,3) = [5.62441575176994e-07]; %nref = 4
U2_ex(:,4) = [7.28109159184251e-08]; %nref = 8

U2_im(:,1) = 3.12132225876817e-06; %nref = 1
U2_im(:,2) = 1.97539347284733e-07; %nref = 2
U2_im(:,3) = [1.75897178220378e-08]; %nref = 4
U2_im(:,4) = [2.04772217477330e-09]; %nref = 8

% Ustar

Ustar2_im = zeros(1,4);
Ustar2_ex = zeros(1,4);

Ustar2_ex(:,1) = 6.65138233288823e-05; %nref = 1
Ustar2_ex(:,2) = 3.52779279237385e-06; %nref = 2
Ustar2_ex(:,3) = 2.019284675901002e-07; %nref = 4
Ustar2_ex(:,4) = 1.221602763996971e-08; %nref = 8

Ustar2_im(:,1) = 2.97949734792521e-06; %nref = 1
Ustar2_im(:,2) = 1.516876757937891e-07; %nref = 2
Ustar2_im(:,3) = 7.059796046760352e-09; %nref = 4
Ustar2_im(:,4) = 2.930597464710975e-10; %nref = 8

%% P = 3 Errors

% v and q
E3_im = zeros(3,4);
E3_ex = zeros(3,4);

E3_ex(:,1) = [4.94707156989462e-05;7.12473271026033e-05;7.12473271025985e-05]; %nref = 1
E3_ex(:,2) = [3.09593621215007e-06;4.83255395664021e-06;4.83255395661897e-06]; %nref = 2
E3_ex(:,3) = [1.92258286580976e-07;3.05465097657568e-07;3.05465098065865e-07]; %nref = 4
% E3_ex(:,4) = [1.18312755074593e-08;1.91644365156702e-08;]; %nref = 8
E3_ex(:,4) = [1.18312755074593e-08;1.91644365156702e-08;1.91644372575802e-08]; %nref = 8

E3_im(:,1) = [4.71827858805950e-06;3.75588958570770e-06;3.75588958573254e-06]; %nref = 1
E3_im(:,2) = [1.76174416028252e-07;3.49228481930775e-07;3.49228481942034e-07]; %nref = 2
E3_im(:,3) = [1.36503162592177e-08;1.76798601024587e-08;1.76798600474657e-08]; %nref = 4
% E3_im(:,4) = [5.17534731982276e-10;6.03286062595280e-10;]; %nref = 8
E3_im(:,4) = [5.17534731982276e-10;6.03286062595280e-10;6.03285502696562e-10]; %nref = 8

% vstar

Vstar3_im = zeros(1,4);
Vstar3_ex = zeros(1,4);

Vstar3_ex(:,1) = 1.454059508886479e-05; %nref = 1
Vstar3_ex(:,2) = 6.729956798024829e-07; %nref = 2
Vstar3_ex(:,3) = 4.043577143762468e-08; %nref = 4
% Vstar3_ex(:,4) = ; %nref = 8
Vstar3_ex(:,4) = 1.895871524719472e-09; %nref = 8

Vstar3_im(:,1) = 4.653922042205659e-06; %nref = 1
Vstar3_im(:,2) = 1.680904204671066e-07; %nref = 2
Vstar3_im(:,3) = 1.328725684404459e-08; %nref = 4
% Vstar3_im(:,4) = ; %nref = 8
Vstar3_im(:,4) = 7.857515003975513e-10; %nref = 8

% U

U3_im = zeros(1,4);
U3_ex = zeros(1,4);

U3_ex(:,1) = [6.09192697560230e-06]; %nref = 1
U3_ex(:,2) = [4.17544861502881e-07]; %nref = 2
U3_ex(:,3) = [2.64905776423895e-08]; %nref = 4
% U3_ex(:,4) = []; %nref = 8
U3_ex(:,4) = [1.68917269781842e-09]; %nref = 8

U3_im(:,1) = [2.34923477912767e-07]; %nref = 1
U3_im(:,2) = [1.52039867995823e-08]; %nref = 2
U3_im(:,3) = [5.85208025180835e-10]; %nref = 4
% U3_im(:,4) = []; %nref = 8
U3_im(:,4) = [3.41470171859376e-11]; %nref = 8

% Ustar

Ustar3_im = zeros(1,4);
Ustar3_ex = zeros(1,4);

Ustar3_ex(:,1) = 3.125657004180897e-06; %nref = 1
Ustar3_ex(:,2) = 9.406670142465745e-08; %nref = 2
Ustar3_ex(:,3) = 2.970477229632685e-09; %nref = 4
% Ustar3_ex(:,4) = ; %nref = 8
Ustar3_ex(:,4) = 3.265993512802381e-10; %nref = 8

Ustar3_im(:,1) = 1.930983689492355e-07; %nref = 1
Ustar3_im(:,2) = 1.280988687525945e-08; %nref = 2
Ustar3_im(:,3) = 2.904869015919024e-10; %nref = 4
% Ustar3_im(:,4) = ; %nref = 8
Ustar3_im(:,4) = 1.293317765758507e-11; %nref = 8


%% P = 4 Errors

% v and q
E4_im = zeros(3,4);
E4_ex = zeros(3,4);

E4_ex(:,1) = [4.35756820642283e-06;6.05784426073026e-06;6.05784426073273e-06]; %nref = 1
E4_ex(:,2) = [1.34977156160827e-07;1.63145393217267e-07;1.63145393313830e-07]; %nref = 2

E4_im(:,1) = [1.83247947252935e-07;2.15557316385206e-07;2.15557316456171e-07]; %nref = 1
E4_im(:,2) = [8.86229671025610e-09;5.44917748598753e-09;5.44917757749544e-09]; %nref = 2


% vstar

Vstar4_im = zeros(1,4);
Vstar4_ex = zeros(1,4);

Vstar4_ex(:,1) = 1.431645329939087e-06; %nref = 1
Vstar4_ex(:,2) = 2.445523260931923e-08; %nref = 2

Vstar4_im(:,1) = 1.822980033021460e-07; %nref = 1
Vstar4_im(:,2) = 8.912599420965596e-09; %nref = 2

% U

U4_im = zeros(1,4);
U4_ex = zeros(1,4);

U4_ex(:,1) = [6.96632009709185e-07]; %nref = 1
U4_ex(:,2) = [2.08506661487927e-08]; %nref = 2
% U4_ex(:,3) = ; %nref = 4
% U4_ex(:,4) = ; %nref = 8

U4_im(:,1) = [1.12283398927964e-08]; %nref = 1
U4_im(:,2) = [1.56851521322319e-10]; %nref = 2
% U4_im(:,3) = ; %nref = 4
% U4_im(:,4) = ; %nref = 8

% Ustar

Ustar4_im = zeros(1,4);
Ustar4_ex = zeros(1,4);

Ustar4_ex(:,1) = 1.992993121021680e-07; %nref = 1
Ustar4_ex(:,2) = 2.537613746448653e-09; %nref = 2
% Ustar4_ex(:,3) = ; %nref = 4
% Ustar4_ex(:,4) = ; %nref = 8

Ustar4_im(:,1) = 1.092558094167440e-08; %nref = 1
Ustar4_im(:,2) = 1.297793143260934e-10; %nref = 2
% Ustar4_im(:,3) = ; %nref = 4
% Ustar4_im(:,4) = ; %nref = 8


%% Compute Total Errors

% Compute Total Q Error
Qexp1 = sqrt(E1_ex(2,:).^2 + E1_ex(3,:).^2);
Qexp2 = sqrt(E2_ex(2,:).^2 + E2_ex(3,:).^2);
Qexp3 = sqrt(E3_ex(2,:).^2 + E3_ex(3,:).^2);
Qexp4 = sqrt(E4_ex(2,:).^2 + E4_ex(3,:).^2);


Qim1 = sqrt(E1_im(2,:).^2 + E1_im(3,:).^2);
Qim2 = sqrt(E2_im(2,:).^2 + E2_im(3,:).^2);
Qim3 = sqrt(E3_im(2,:).^2 + E3_im(3,:).^2);
Qim4 = sqrt(E4_im(2,:).^2 + E4_im(3,:).^2);


Q1 = sqrt(Qexp1.^2 + Qim1.^2); % p = 1
Q2 = sqrt(Qexp2.^2 + Qim2.^2); % p = 2
Q3 = sqrt(Qexp3.^2 + Qim3.^2); % p = 3
Q4 = sqrt(Qexp4.^2 + Qim4.^2); % p = 4


% Compute Total V Error
V1 = sqrt(E1_ex(1,:).^2 + E1_im(1,:).^2);
V2 = sqrt(E2_ex(1,:).^2 + E2_im(1,:).^2);
V3 = sqrt(E3_ex(1,:).^2 + E3_im(1,:).^2);
V4 = sqrt(E4_ex(1,:).^2 + E4_im(1,:).^2);


% Compute Total Vstar Error
Vstar1 = sqrt(Vstar1_ex(1,:).^2 + Vstar1_im(1,:).^2);
Vstar2 = sqrt(Vstar2_ex(1,:).^2 + Vstar2_im(1,:).^2);
Vstar3 = sqrt(Vstar3_ex(1,:).^2 + Vstar3_im(1,:).^2);
Vstar4 = sqrt(Vstar4_ex(1,:).^2 + Vstar4_im(1,:).^2);


% Compute Total V Error
U1 = sqrt(U1_ex(1,:).^2 + U1_im(1,:).^2);
U2 = sqrt(U2_ex(1,:).^2 + U2_im(1,:).^2);
U3 = sqrt(U3_ex(1,:).^2 + U3_im(1,:).^2);
U4 = sqrt(U4_ex(1,:).^2 + U4_im(1,:).^2);


% Compute Total Vstar Error
Ustar1 = sqrt(Ustar1_ex(1,:).^2 + Ustar1_im(1,:).^2);
Ustar2 = sqrt(Ustar2_ex(1,:).^2 + Ustar2_im(1,:).^2);
Ustar3 = sqrt(Ustar3_ex(1,:).^2 + Ustar3_im(1,:).^2);
Ustar4 = sqrt(Ustar4_ex(1,:).^2 + Ustar4_im(1,:).^2);

%% Plots

% h lengths
hmin_est3 = [0.07,0.035,0.0175];
hmin_est = [0.07,0.035,0.0175, 0.0088];
hmin_est2 = [0.07,0.035,0.0175, 0.0088, 0.0044];

%%%%%%%%%%%%%%%%%% U Plot %%%%%%%%%%%%%%%%%%
figure()
loglog(hmin_est2,U1,'-s');
hold on
loglog(hmin_est,U2,'-s');
loglog(hmin_est,U3,'-s');
loglog(hmin_est,U4,'-s');

title('u convergence');
xlabel('log(h)');
ylabel('log(Error)');
legend('p=1','p=2','p=3','p=4');
axis tight


%%%%%%%%%%%%%%%%%% V Plot %%%%%%%%%%%%%%%%%%
figure()
loglog(hmin_est2,V1,'-s');
hold on
loglog(hmin_est,V2,'-s');
loglog(hmin_est,V3,'-s');
loglog(hmin_est,V4,'-s');

title('v convergence');
xlabel('log(h)');
ylabel('log(Error)');
legend('p=1','p=2','p=3','p=4');
axis tight


%%%%%%%%%%%%%%%%%% Q Plot %%%%%%%%%%%%%%%%%%
figure()
loglog(hmin_est2,Q1,'-s');
hold on
loglog(hmin_est,Q2,'-s');
loglog(hmin_est,Q3,'-s');
loglog(hmin_est,Q4,'-s');

title('q convergence');
xlabel('log(h)');
ylabel('log(Error)');
legend('p=1','p=2','p=3','p=4');
axis tight

%%%%%%%%%%%%%%%%%% Ustar Plot %%%%%%%%%%%%%%%%%%
figure()
loglog(hmin_est2,Ustar1,'-s');
hold on
loglog(hmin_est,Ustar2,'-s');
loglog(hmin_est,Ustar3,'-s');
loglog(hmin_est,Ustar4,'-s');

title('u* convergence');
xlabel('log(h)');
ylabel('log(Error)');
legend('p=1','p=2','p=3','p=4');
axis tight

%%%%%%%%%%%%%%%%%% Vstar Plot %%%%%%%%%%%%%%%%%%
figure()
loglog(hmin_est2,Vstar1,'-s');
hold on
loglog(hmin_est,Vstar2,'-s');
loglog(hmin_est,Vstar3,'-s');
loglog(hmin_est,Vstar4,'-s');

title('v* convergence');
xlabel('log(h)');
ylabel('log(Error)');
legend('p=1','p=2','p=3','p=4');
axis tight

%%%%%%%%%%%%%%%%%% Vstar-V Plots %%%%%%%%%%%%%%%%%%

% p = 1
figure()
loglog(hmin_est2,V1,'-s');
hold on
loglog(hmin_est2,Vstar1,'-s');

title('v convergence, p = 1');
xlabel('log(h)');
ylabel('log(Error)');
legend('v','vstar');

% p = 2
figure()
loglog(hmin_est,V2,'-s');
hold on
loglog(hmin_est,Vstar2,'-s');

title('v convergence, p = 2');
xlabel('log(h)');
ylabel('log(Error)');
legend('v','vstar');

% p = 3
figure()
loglog(hmin_est,V3,'-s');
hold on
loglog(hmin_est,Vstar3,'-s');

title('v convergence, p = 3');
xlabel('log(h)');
ylabel('log(Error)');
legend('v','vstar');

% p = 4
figure()
loglog(hmin_est,V4,'-s');
hold on
loglog(hmin_est,Vstar4,'-s');

title('v convergence, p = 4');
xlabel('log(h)');
ylabel('log(Error)');
legend('v','vstar');


%% Slope Calculations


%%%%%%%%%%%%%%%%%% u slopes %%%%%%%%%%%%%%%%%%
U1m = zeros(4,1);
U2m = zeros(4,1);
U3m = zeros(4,1);
U4m = zeros(4,1);

% p = 1
U1m(1) = (log10(U1(1)) - log10(U1(2)))/(log10(hmin_est(1)) - log10(hmin_est(2)));
U1m(2) = (log10(U1(2)) - log10(U1(3)))/(log10(hmin_est(2)) - log10(hmin_est(3)));
U1m(3) = (log10(U1(3)) - log10(U1(4)))/(log10(hmin_est(3)) - log10(hmin_est(4)));
U1m(4) = (log10(U1(4)) - log10(U1(5)))/(log10(hmin_est2(4)) - log10(hmin_est2(5)));
% p = 2
U2m(1) = (log10(U2(1)) - log10(U2(2)))/(log10(hmin_est(1)) - log10(hmin_est(2)));
U2m(2) = (log10(U2(2)) - log10(U2(3)))/(log10(hmin_est(2)) - log10(hmin_est(3)));
U2m(3) = (log10(U2(3)) - log10(U2(4)))/(log10(hmin_est(3)) - log10(hmin_est(4)));
% p = 3
U3m(1) = (log10(U3(1)) - log10(U3(2)))/(log10(hmin_est(1)) - log10(hmin_est(2)));
U3m(2) = (log10(U3(2)) - log10(U3(3)))/(log10(hmin_est(2)) - log10(hmin_est(3)));
U3m(3) = (log10(U3(3)) - log10(U3(4)))/(log10(hmin_est(3)) - log10(hmin_est(4)));
% p = 4
U4m(1) = (log10(U4(1)) - log10(U4(2)))/(log10(hmin_est(1)) - log10(hmin_est(2)));

%%%%%%%%%%%%%%%%%% v slopes %%%%%%%%%%%%%%%%%%
V1m = zeros(4,1);
V2m = zeros(4,1);
V3m = zeros(4,1);
V4m = zeros(4,1);

% p = 1
V1m(1) = (log10(V1(1)) - log10(V1(2)))/(log10(hmin_est(1)) - log10(hmin_est(2)));
V1m(2) = (log10(V1(2)) - log10(V1(3)))/(log10(hmin_est(2)) - log10(hmin_est(3)));
V1m(3) = (log10(V1(3)) - log10(V1(4)))/(log10(hmin_est(3)) - log10(hmin_est(4)));
V1m(4) = (log10(V1(4)) - log10(V1(5)))/(log10(hmin_est2(4)) - log10(hmin_est2(5)));
% p = 2
V2m(1) = (log10(V2(1)) - log10(V2(2)))/(log10(hmin_est(1)) - log10(hmin_est(2)));
V2m(2) = (log10(V2(2)) - log10(V2(3)))/(log10(hmin_est(2)) - log10(hmin_est(3)));
V2m(3) = (log10(V2(3)) - log10(V2(4)))/(log10(hmin_est(3)) - log10(hmin_est(4)));
% p = 3
V3m(1) = (log10(V3(1)) - log10(V3(2)))/(log10(hmin_est(1)) - log10(hmin_est(2)));
V3m(2) = (log10(V3(2)) - log10(V3(3)))/(log10(hmin_est(2)) - log10(hmin_est(3)));
V3m(3) = (log10(V3(3)) - log10(V3(4)))/(log10(hmin_est(3)) - log10(hmin_est(4)));
% p = 4
V4m(1) = (log10(V4(1)) - log10(V4(2)))/(log10(hmin_est(1)) - log10(hmin_est(2)));


%%%%%%%%%%%%%%%%%% Q slopes %%%%%%%%%%%%%%%%%%
Q1m = zeros(4,1);
Q2m = zeros(4,1);
Q3m = zeros(4,1);
Q4m = zeros(4,1);

% p = 1
Q1m(1) = (log10(Q1(1)) - log10(Q1(2)))/(log10(hmin_est(1)) - log10(hmin_est(2)));
Q1m(2) = (log10(Q1(2)) - log10(Q1(3)))/(log10(hmin_est(2)) - log10(hmin_est(3)));
Q1m(3) = (log10(Q1(3)) - log10(Q1(4)))/(log10(hmin_est(3)) - log10(hmin_est(4)));
Q1m(4) = (log10(Q1(4)) - log10(Q1(5)))/(log10(hmin_est2(4)) - log10(hmin_est2(5)));
% p = 2
Q2m(1) = (log10(Q2(1)) - log10(Q2(2)))/(log10(hmin_est(1)) - log10(hmin_est(2)));
Q2m(2) = (log10(Q2(2)) - log10(Q2(3)))/(log10(hmin_est(2)) - log10(hmin_est(3)));
Q2m(3) = (log10(Q2(3)) - log10(Q2(4)))/(log10(hmin_est(3)) - log10(hmin_est(4)));
% p = 3
Q3m(1) = (log10(Q3(1)) - log10(Q3(2)))/(log10(hmin_est(1)) - log10(hmin_est(2)));
Q3m(2) = (log10(Q3(2)) - log10(Q3(3)))/(log10(hmin_est(2)) - log10(hmin_est(3)));
Q3m(3) = (log10(Q3(3)) - log10(Q3(4)))/(log10(hmin_est(3)) - log10(hmin_est(4)));
% p = 4
Q4m(1) = (log10(Q4(1)) - log10(Q4(2)))/(log10(hmin_est(1)) - log10(hmin_est(2)));


%%%%%%%%%%%%%%%%%% Ustar slopes %%%%%%%%%%%%%%%%%%
Ustar1m = zeros(4,1);
Ustar2m = zeros(4,1);
Ustar3m = zeros(4,1);
Ustar4m = zeros(4,1);

% p = 1
Ustar1m(1) = (log10(Ustar1(1)) - log10(Ustar1(2)))/(log10(hmin_est(1)) - log10(hmin_est(2)));
Ustar1m(2) = (log10(Ustar1(2)) - log10(Ustar1(3)))/(log10(hmin_est(2)) - log10(hmin_est(3)));
Ustar1m(3) = (log10(Ustar1(3)) - log10(Ustar1(4)))/(log10(hmin_est(3)) - log10(hmin_est(4)));
Ustar1m(4) = (log10(Ustar1(4)) - log10(Ustar1(5)))/(log10(hmin_est2(4)) - log10(hmin_est2(5)));
% p = 2
Ustar2m(1) = (log10(Ustar2(1)) - log10(Ustar2(2)))/(log10(hmin_est(1)) - log10(hmin_est(2)));
Ustar2m(2) = (log10(Ustar2(2)) - log10(Ustar2(3)))/(log10(hmin_est(2)) - log10(hmin_est(3)));
Ustar2m(3) = (log10(Ustar2(3)) - log10(Ustar2(4)))/(log10(hmin_est(3)) - log10(hmin_est(4)));
% p = 3
Ustar3m(1) = (log10(Ustar3(1)) - log10(Ustar3(2)))/(log10(hmin_est(1)) - log10(hmin_est(2)));
Ustar3m(2) = (log10(Ustar3(2)) - log10(Ustar3(3)))/(log10(hmin_est(2)) - log10(hmin_est(3)));
Ustar3m(3) = (log10(Ustar3(3)) - log10(Ustar3(4)))/(log10(hmin_est(3)) - log10(hmin_est(4)));
% p = 4
Ustar4m(1) = (log10(Ustar4(1)) - log10(Ustar4(2)))/(log10(hmin_est(1)) - log10(hmin_est(2)));

%%%%%%%%%%%%%%%%%% Vstar slopes %%%%%%%%%%%%%%%%%%
Vstar1m = zeros(4,1);
Vstar2m = zeros(4,1);
Vstar3m = zeros(4,1);
Vstar4m = zeros(4,1);

% p = 1
Vstar1m(1) = (log10(Vstar1(1)) - log10(Vstar1(2)))/(log10(hmin_est(1)) - log10(hmin_est(2)));
Vstar1m(2) = (log10(Vstar1(2)) - log10(Vstar1(3)))/(log10(hmin_est(2)) - log10(hmin_est(3)));
Vstar1m(3) = (log10(Vstar1(3)) - log10(Vstar1(4)))/(log10(hmin_est(3)) - log10(hmin_est(4)));
Vstar1m(4) = (log10(Vstar1(4)) - log10(Vstar1(5)))/(log10(hmin_est2(4)) - log10(hmin_est2(5)));
% p = 2
Vstar2m(1) = (log10(Vstar2(1)) - log10(Vstar2(2)))/(log10(hmin_est(1)) - log10(hmin_est(2)));
Vstar2m(2) = (log10(Vstar2(2)) - log10(Vstar2(3)))/(log10(hmin_est(2)) - log10(hmin_est(3)));
Vstar2m(3) = (log10(Vstar2(3)) - log10(Vstar2(4)))/(log10(hmin_est(3)) - log10(hmin_est(4)));
% p = 3
Vstar3m(1) = (log10(Vstar3(1)) - log10(Vstar3(2)))/(log10(hmin_est(1)) - log10(hmin_est(2)));
Vstar3m(2) = (log10(Vstar3(2)) - log10(Vstar3(3)))/(log10(hmin_est(2)) - log10(hmin_est(3)));
Vstar3m(3) = (log10(Vstar3(3)) - log10(Vstar3(4)))/(log10(hmin_est(3)) - log10(hmin_est(4)));
% p = 4
Vstar4m(1) = (log10(Vstar4(1)) - log10(Vstar4(2)))/(log10(hmin_est(1)) - log10(hmin_est(2)));


%%%%%%%%%%%%%%%%%% All Values %%%%%%%%%%%%%%%%%%
All_Vm = [V1m V2m V3m V4m];
All_Qm = [Q1m Q2m Q3m Q4m];
All_Um = [U1m U2m U3m U4m];
All_Vstarm = [Vstar1m Vstar2m Vstar3m Vstar4m];
All_Ustarm = [Ustar1m Ustar2m Ustar3m Ustar4m];