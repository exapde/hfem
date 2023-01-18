function mesh = mkmesh_naca0012foil(xf, yf, porder, nelem, d)

shift = 0.001;                   

dxu = xf(2)-xf(1);
dyu = yf(2)-yf(1);
dxl = xf(end-1)-xf(end);
dyl = yf(end-1)-yf(end);
tau = acos((dxu*dxl+dyu*dyl)/sqrt(dxu*dxu+dyu*dyu)/sqrt(dxl*dxl+dyl*dyl));
n   = 2-tau/pi;

chord = max(xf)-min(xf);
xmo = 0.5*(max(xf)+min(xf));
ymo = yf(1);
xn = 2.0*(1.0+shift)*(xf-xmo)/chord-shift;
yn = 2.0*(1.0+shift)*(yf-ymo)/chord;

% Near circle
[x,y] = trefftz_inv( xn, yn, n, 1/n, 0,1);

% Center
xm = 0.5*(max(x)+min(x));
ym = 0.5*(max(y)+min(y));
fc = (max(x)-min(x)+max(y)-min(y))/4;
x = (x-xm)/fc;
y = (y-ym)/fc;

% Now a circle
[A,B] = GetTG( 200, x, y);

t = linspace(0,2*pi,nelem*porder+1);
xc = 0.5*cos(t); 
yc = 0.5*sin(t); 

% Now a near circle
[xg, yg] = TG( 2.0*xc, 2.0*yc, A, B);

% Finally the real thing
[xgf, ygf] = trefftz( xg*fc+xm, yg*fc+ym, n, 1/n, 0);

xgf = (xgf+shift)*chord/(2.0*(1.0+shift)) + xmo;
ygf = ygf*chord/(2.0*(1.0+shift)) + ymo;

pcg = [xgf(:) ygf(:)]; 
%pcg = pcg(end:-1:1,:);
tcg = zeros(nelem,porder+1);
for i = 1:porder+1
    tcg(:,i) = (i:porder:(porder*(nelem-1)+i))';
end

nodetype=0;
elemtype=1;
plocal = masternodes(porder,1,elemtype,nodetype);
shapft = mkshape(porder,plocal,plocal,elemtype);
dshapft = shapft(:,:,2)';

pdg = zeros(porder+1,2,nelem);
ndg = zeros(porder+1,2,nelem);
for i = 1:nelem
    pn = pcg(tcg(i,:),:);
    dpg = dshapft*pn;        
    jac = sqrt(dpg(:,1).^2+dpg(:,2).^2);
    nl   = [dpg(:,2),-dpg(:,1)];
    nl   = bsxfun(@rdivide, nl, jac);
    pdg(:,:,i) = pn;
    ndg(:,:,i) = nl;
    if i>1
        n1 = ndg(end,:,i-1);
        n2 = ndg(1,:,i);
        nn = 0.5*(n1+n2);
        ndg(end,:,i-1) = nn;
        ndg(1,:,i) = nn;
    end
end


% ne = porder*nelem;
% n1 = reshape(ndg(porder+1,:,nelem),[1 2]); 
% ncg = ndg(1:porder,:,:);
% ncg = reshape(permute(ncg,[1 3 2]),[ne 2]);
% ncg = [ncg; n1];
% 
% dp = d*plocal;
% X = zeros(ne+1,porder+1);
% Y = X;
% for i = 1:(ne+1)
%     for j=1:(porder+1)
%         X(i,j) = pcg(i,1) + ncg(i,1)*dp(j);
%         Y(i,j) = pcg(i,2) + ncg(i,2)*dp(j);
%     end
% %     if (i==1) || (i==ne+1)
% %         for j=1:(porder+1)
% %             X(i,j) = pcg(i,1) + 6*dp(j);
% %             Y(i,j) = 0;
% %         end
% %     end    
% end
% 
% [p,t,dgnodes]=cart2dg(elemtype,porder,X,Y);
% 
% bndexpr = {'all(p(:,2)<1e-6)','true'};     
% mesh = mkmesh(p,t,porder,bndexpr,elemtype,nodetype);
% mesh.dgnodes = dgnodes;
% 
% figure(1); clf; plot(pcg(:,1),pcg(:,2),'-');
% figure(2); clf; meshplot(mesh,1);

p1 = porder+1;
dp = d*plocal;
dgnodes = zeros((porder+1),(porder+1),2,nelem);
peg = zeros(porder+1,2,nelem);
for i = 1:nelem    
    pt = pdg(:,:,i);      
    nl = ndg(:,:,i);
    peg(:,:,i) = pt + dp(porder+1)*nl;
    for j=1:p1
        dgnodes(:,j,:,i) = reshape(pt + dp(j)*nl,[(porder+1) 1 2]);           
    end
end
dgnodes = reshape(dgnodes,[(porder+1)*(porder+1),2,nelem]);

x = reshape(pdg(porder+1,:,nelem),[1 2]); 
pcg = pdg(1:porder,:,:);
pcg = reshape(permute(pcg,[1 3 2]),[porder*nelem 2]);
pcg = [pcg; x];

x = reshape(peg(porder+1,:,nelem),[1 2]); 
peg = peg(1:porder,:,:);
peg = reshape(permute(peg,[1 3 2]),[porder*nelem 2]);
peg = [peg; x];

% ind = 1:porder:(porder*nelem+1);
% p = [pcg(ind,:); peg(ind,:)];
% m = length(ind);
% n = 2;
% t = [1 2 m+2 m+1];
% t=kron(t,ones(n-1,1))+kron(ones(size(t)),(0:n-2)'*m);
% t=kron(t,ones(m-1,1))+kron(ones(size(t)),(0:m-2)');
% 
% %figure(1); clf; simpplot(p,t); 
% bndexpr = {'all(p(:,2)<1e-6)','true'};     
% mesh = mkmesh(p,t,porder,bndexpr,elemtype,nodetype);
% mesh.dgnodes = dgnodes;

ind = 1:porder:(porder*nelem+1);
p = [pcg(ind,:); peg(ind,:); [max(xf)+4*d 0]; [max(xf)+4*d d/2]; [max(xf)+4*d -d/2]];
m = length(ind);
n = 2;
t = [1 2 m+2 m+1];
t=kron(t,ones(n-1,1))+kron(ones(size(t)),(0:n-2)'*m);
t=kron(t,ones(m-1,1))+kron(ones(size(t)),(0:m-2)');
t=[t; [2*m+1 1 m+1 2*m+2]; [1 2*m+1 2*m+3 2*m]];

bndexpr = {'all(p(:,2)<1e-6)','true'};     
mesh = mkmesh(p,t,porder,bndexpr,elemtype,nodetype);
mesh.dgnodes(:,:,1:(end-2)) = dgnodes;

%figure(1); clf; simpplot(p,t); 

% xgf = dgnodes(:,1,:);
% ygf = dgnodes(:,2,:);
% figure(1); clf; plot(xf,yf,'-',xgf(:),ygf(:),'o'); axis equal;
% figure(2); clf; plot(pcg(:,1),pcg(:,2),'*',pdg(:,1),pdg(:,2),'o'); axis equal;
% figure(3); clf; plot(p(:,1),p(:,2),'*'); axis equal;



function [x,y] = trefftz( x1, y1, n, cx, cy)
z1 = complex( x1, y1);
cc = complex( cx, cy);
A = ((z1-cc)./(z1+cc)).^n;
z = ((1+A)./(1-A))*n*cc;
x = real(z);
y = imag(z);

function [x,y] = trefftz_inv( x1, y1, n, cx, cy, track)
z1 = complex( x1, y1);
cc = complex( cx, cy);
A = ((z1-n*cc)./(z1+n*cc));
if track
   R = abs(A); 
   T = angle(A);
   for k = 2:size(T,1)
       d(1) = T(k) + 2*pi;
       d(2) = T(k);
       d(3) = T(k) - 2*pi;
       [minv,j] = min(abs(d-T(k-1)));
       T(k) = d(j);
   end
   B = ((R).^(1/n)).*exp(i*T/n);
else
   B = A.^(1/n);
end
z = ((1+B)./(1-B))*cc;
x = real(z);
y = imag(z);

function [x,y] = TG( xc, yc, A, B)
N = size(A,2)-1;
zc = complex( xc, yc);
e = complex(zeros(size(zc)),zeros(size(zc)));
for j=1:N+1
    e = e + (A(j)+i*B(j)).*zc.^(1-j);
end
ix = isnan(e);
e(ix) = 0;
z = zc.*exp(e);
x = real(z);
y = imag(z);


function [A,B] = GetTG( N, x, y)
lr = log(sqrt(x.^2 + y.^2));
th = atan2(y,x);
for k=2:size(th,1)
    if (th(k) < th(k-1)) 
        th(k) = th(k) + 2*pi;
    end
end
thi = [th(1:end-1) - 2*pi; th(1:end); th(2:end)+2*pi];
lri = lr([1:end-1,1:end,2:end]);
A = ones(1,N+1);
B = zeros(1,N+1);
Y = complex(zeros(1,2*N), zeros(1,2*N));
tt = 0:pi/N:2*pi;
tt = tt(1:end-1);
Anew = 0*A;
Bnew = B; 
while norm(A-Anew) + norm(B-Bnew) > 1.e-15,
    A = Anew;
    B = Bnew;
    B(1) = th(1) - sum(B(2:N+1));
    B(N+1) = 0;
    Y(1) = 2*N*B(1);
    Y(2:N) = N*(B(2:N)+ i*A(2:N));
    Y(N+1) = 2*N*B(N+1);
    Y(N+2:2*N) = conj(Y(N:-1:2));
    zt = tt + ifft(Y);
    rr = spline(thi, lri, zt);
    Y = fft(rr);
    Anew(1) = real(Y(1))/(2*N);
    Anew(2:N) = real(Y(2:N))/N;
    Anew(N+1) = real(Y(N+1))/(2*N);
    Bnew(1) = B(1);
    Bnew(2:N) = -imag(Y(2:N))/N;
    Bnew(N+1) = 0;
end




