load('sole040delta25uinf0.mat')

% beta = app.arg{end-1}(1);
% Etmax = 400*1e3;
% delta_peek = 0.25/100;
% windlevel = 0;
% iterwire;

load('sole0470delta50uinf200.mat');
beta = app.arg{10}(1);
windvel = 200;
delta_peek = 0.5/100;
Einf = [480:10:600]*1e3;
for jj = 1:length(Einf)                    
    Etmax = Einf(jj);
    beta = 1.2*beta;
    iterwire;            
    fn = ['sole0' num2str(Etmax/1e3) 'delta' num2str(delta_peek*1e4) 'uinf' num2str(windvel) '.mat'];
    save(fn,'UDG','UH','app','Etmax','delta_peek','windvel');    
end       


load('sole0200delta50uinf25.mat');
beta = app.arg{10}(1);
Uinf = [150 100 75 50 25 10 5];
delta_peek = 0.5/100;
Einf = 200;
for jj = 7:length(Uinf)                    
    windvel = Uinf(jj);
    beta = beta/1.75;
    iterwire;            
    fn = ['sole0' num2str(Etmax/1e3) 'delta' num2str(delta_peek*1e4) 'uinf' num2str(windvel) '.mat'];
    save(fn,'UDG','UH','app','Etmax','delta_peek','windvel');    
end       


load('sole0600delta50uinf200.mat');
beta = app.arg{10}(1);
Uinf = [150 100 75 50 25 10 5];
delta_peek = 0.5/100;
Einf = 600;
for jj = 1:length(Uinf)                    
    windvel = Uinf(jj);
    beta = beta/1.25;
    iterwire;            
    fn = ['sole0' num2str(Etmax/1e3) 'delta' num2str(delta_peek*1e4) 'uinf' num2str(windvel) '.mat'];
    save(fn,'UDG','UH','app','Etmax','delta_peek','windvel');    
end       

load('sole0400delta50uinf5.mat');
beta = app.arg{10}(1);
Uinf = [4 3 2 1 0];
delta_peek = 0.5/100;
Einf = 400;
for jj = 1:length(Uinf)                    
    windvel = Uinf(jj);
    beta = beta/1.05;
    iterwire;            
    fn = ['sole0' num2str(Etmax/1e3) 'delta' num2str(delta_peek*1e4) 'uinf' num2str(windvel) '.mat'];
    save(fn,'UDG','UH','app','Etmax','delta_peek','windvel');    
end       

load('sole0300delta50uinf5.mat');
beta = app.arg{10}(1);
Uinf = [4 3 2 1 0];
delta_peek = 0.5/100;
Einf = 300;
for jj = 1:length(Uinf)                    
    windvel = Uinf(jj);
    beta = beta/1.05;
    iterwire;            
    fn = ['sole0' num2str(Etmax/1e3) 'delta' num2str(delta_peek*1e4) 'uinf' num2str(windvel) '.mat'];
    save(fn,'UDG','UH','app','Etmax','delta_peek','windvel');    
end       

load('sole0200delta50uinf5.mat');
beta = app.arg{10}(1);
Uinf = [4 3 2 1 0];
delta_peek = 0.5/100;
Einf = 200;
for jj = 1:length(Uinf)                    
    windvel = Uinf(jj);
    beta = beta/1.05;
    iterwire;            
    fn = ['sole0' num2str(Etmax/1e3) 'delta' num2str(delta_peek*1e4) 'uinf' num2str(windvel) '.mat'];
    save(fn,'UDG','UH','app','Etmax','delta_peek','windvel');    
end       

load('sole0100delta50uinf5.mat');
beta = app.arg{10}(1);
Uinf = [4 3 2 1 0];
delta_peek = 0.5/100;
Einf = 100;
for jj = 1:length(Uinf)                    
    windvel = Uinf(jj);
    beta = beta/1.05;
    iterwire;            
    fn = ['sole0' num2str(Etmax/1e3) 'delta' num2str(delta_peek*1e4) 'uinf' num2str(windvel) '.mat'];
    save(fn,'UDG','UH','app','Etmax','delta_peek','windvel');    
end       

load('sole0100delta50uinf0.mat');
beta = app.arg{10}(1);
windvel = 0;
delta_peek = 0.5/100;
Einf = [150:100:650 700:50:800]*1e3;
for jj = 1:length(Einf)                    
    Etmax = Einf(jj);
    beta = beta*1.25;
    [Etmax delta_peek windvel]
    iterwire;            
    fn = ['sole0' num2str(Etmax/1e3) 'delta' num2str(delta_peek*1e4) 'uinf' num2str(windvel) '.mat'];
    save(fn,'UDG','UH','app','Etmax','delta_peek','windvel');    
end       

load('sole1390delta50uinf0.mat');
beta = app.arg{10}(1);
windvel = 0;
delta_peek = 0.5/100;
Einf = [400:25:800]*1e3;
for jj = 1:length(Einf)                    
    Etmax = Einf(jj);
    beta = beta*1.015;
    [Etmax delta_peek windvel]
    iterwire;            
    fn = ['sole1' num2str(Etmax/1e3) 'delta' num2str(delta_peek*1e4) 'uinf' num2str(windvel) '.mat'];
    save(fn,'UDG','UH','app','Etmax','delta_peek','windvel');    
end       


load('sole0100delta50uinf50.mat');
beta = app.arg{10}(1);
Uinf = [150 100 75 50 25 10 5];
delta_peek = 0.5/100;
Einf = 100;
for jj = 5:length(Uinf)                    
    windvel = Uinf(jj);
    beta = beta/1.25;
    iterwire;            
    fn = ['sole0' num2str(Etmax/1e3) 'delta' num2str(delta_peek*1e4) 'uinf' num2str(windvel) '.mat'];
    save(fn,'UDG','UH','app','Etmax','delta_peek','windvel');    
end       

load('sole0200delta50uinf50.mat');
beta = app.arg{10}(1);
Uinf = [150 100 75 50 25 10 5];
delta_peek = 0.5/100;
Einf = 200;
for jj = 5:length(Uinf)                    
    windvel = Uinf(jj);
    beta = beta/1.25;
    iterwire;            
    fn = ['sole0' num2str(Etmax/1e3) 'delta' num2str(delta_peek*1e4) 'uinf' num2str(windvel) '.mat'];
    save(fn,'UDG','UH','app','Etmax','delta_peek','windvel');    
end       

Einf = [100:100:600]*1e3;
Uinf = [200 150 100 75 50 25 10 5];
delta_peek = 0.5/100;
I1 = zeros(length(Uinf),length(Einf));
for ii = 1:length(Einf)
    Etmax = Einf(ii);
    for jj = 1:length(Uinf)   
        windvel = Uinf(jj);
        fn = ['sole0' num2str(Etmax/1e3) 'delta' num2str(delta_peek*1e4) 'uinf' num2str(windvel) '.mat'];
        load(fn);
        %streamer;
        I1(jj,ii) = current(master,mesh,UDG,-5,Etmax);
    end
end

figure(5);clf;
plot(Uinf,1e6*I1(:,1),'b-',Uinf,1e6*I1(:,2),'r-',...
     Uinf,1e6*I1(:,3),'g-',Uinf,1e6*I1(:,4),'k-',...
     Uinf,1e6*I1(:,5),'m-',Uinf,1e6*I1(:,6),'c-',...
     'MarkerSize',8,'LineWidth',1.5);
xlabel('u_\infty (m/s)','FontSize',20); 
ylabel('Corona current per unit length (\muA/m)','FontSize',20);
legend({'$E_\infty = 100 \mbox{ KV/m}$','$E_\infty = 200 \mbox{ KV/m}$',...
        '$E_\infty = 300 \mbox{ KV/m}$','$E_\infty = 400 \mbox{ KV/m}$',...
        '$E_\infty = 500 \mbox{ KV/m}$','$E_\infty = 600 \mbox{ KV/m}$'},...
        'interpreter','latex','location','NorthWest','FontSize',18);
set(gca,'FontSize',18); box on;
axis tight;
%axis([0 200 0 200]);
set(gca, 'LooseInset', get(gca, 'TightInset'));
fn = ['wire_current_vs_velocity3.png'];
print('-dpng',fn);

figure(6);clf;
plot(Einf/1e3,1e6*I1(1,:),'b-',Einf/1e3,1e6*I1(3,:),'r-',...
     Einf/1e3,1e6*I1(5,:),'g-',Einf/1e3,1e6*I1(6,:),'k-',...
     Einf/1e3,1e6*I1(7,:),'m-',Einf/1e3,1e6*I1(8,:),'c-',...
     'MarkerSize',8,'LineWidth',1.5);
xlabel('E_\infty (KV/m)','FontSize',20); 
ylabel('Corona current per unit length (\muA/m)','FontSize',20);
legend({'$u_\infty = 200 \mbox{ m/s}$','$u_\infty = 100 \mbox{ m/s}$',...
        '$u_\infty = 50 \mbox{ m/s}$','$u_\infty = 25 \mbox{ m/s}$',...
        '$u_\infty = 10 \mbox{ m/s}$','$u_\infty = 5 \mbox{ m/s}$'},...
        'interpreter','latex','location','NorthWest','FontSize',18);
set(gca,'FontSize',18); box on;
axis tight;
%axis([0 200 0 200]);
set(gca, 'LooseInset', get(gca, 'TightInset'));
fn = ['wire_current_vs_externalfield1.png'];
print('-dpng',fn);

Einf = [30 40 50 60:20:400 410:10:600]*1e3;
windvel = 200;
delta_peek = 0.5/100;
I2 = zeros(length(Einf),1);
for ii = 1:length(Einf)
    Etmax = Einf(ii);             
    fn = ['sole0' num2str(Etmax/1e3) 'delta' num2str(delta_peek*1e4) 'uinf' num2str(windvel) '.mat'];
    load(fn);
    I2(ii) = current(master,mesh,UDG,-5,Etmax);    
end

figure(7);clf;
plot(Einf/1e3,1e6*I2(:),'b-','MarkerSize',8,'LineWidth',1.5);
xlabel('E_\infty (KV/m)','FontSize',20); 
ylabel('Corona current per unit length (\muA/m)','FontSize',20);
legend({'$u_\infty = 200 \mbox{ m/s}$'},...
        'interpreter','latex','location','NorthWest','FontSize',18);
set(gca,'FontSize',18); box on;
axis tight;
axis([30 600 0 2000]);
set(gca, 'LooseInset', get(gca, 'TightInset'));
fn = ['wire_current_vs_externalfield2.png'];
print('-dpng',fn);


