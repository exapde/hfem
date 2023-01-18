syms u1 u2 u3 u4
syms uh1
syms uinf1
syms x1 x2 x3 
syms nl1 nl2 nl3
syms time
syms param1 param2

appname = 'poisson';

param = [param1 param2];
if nd==2    
    udg = [u1 u2 u3];
    uh = uh1;
    uinf = uinf1;
    pg = [x1 x2];  
    nl = [nl1 nl2];    
else
    udg = [u1 u2 u3 u4];
    uh = uh1;
    uinf = uinf1;
    pg = [x1 x2 x3];  
    nl = [nl1 nl2 nl3];    
end

ncu = length(uh);
nc = length(udg);
nch = ncu;

kappa  = param(1);
tau  = param(2);
if nd==2                                               
    u    = udg(1);
    qx   = udg(2);
    qy   = udg(3);    
    f = [kappa*qx, kappa*qy];        
    fhat = simplify(f(1)*nl(1) + f(2)*nl(2) + tau*(u-uh));
else                                         
    u    = udg(1);
    qx   = udg(2);
    qy   = udg(3);    
    qz   = udg(4);    
    f = [kappa*qx, kappa*qy, kappa*qz];        
    fhat = simplify(f(1)*nl(1) + f(2)*nl(2) + f(3)*nl(3) + tau*(u-uh));
end

filename1 = ['flux_' appname num2str(nd) 'd' '.c'];
filename2 = ['source_' appname num2str(nd) 'd'  '.c'];
filename3 = ['fhat_' appname num2str(nd) 'd' '.c'];

genccode; % generate source codes

filename1 = ['fluxonly_' appname num2str(nd) 'd' '.c'];
filename2 = ['sourceonly_' appname num2str(nd) 'd'  '.c'];
filename3 = ['fhatonly_' appname num2str(nd) 'd' '.c'];

genccode_withoutjac; % generate source codes

% generate an application file
gid = fopen('fluxes.c','w');
fid = fopen('../../fluxes_template.c','r');    
tline = fgetl(fid); 
while ischar(tline)        
    str = tline;
    str = strrep(str, 'appname', appname);    
    fprintf(gid, '%s\n', str);    
    tline = fgetl(fid);            
end
fclose(fid);
fclose(gid);




