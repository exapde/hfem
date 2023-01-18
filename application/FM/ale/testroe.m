
syms r ru rv rE gam nx ny v1 v2

gam1 = gam-1;
r1   = 1./r;
u    = ru.*r1;
v    = rv.*r1;
E    = rE.*r1;
q    = 0.5*(u.*u+v.*v);
p    = gam1*(rE-r.*q);
h    = E+p.*r1;
c2   = gam* p.*r1;
c    = sqrt(c2);
                                        
f1 = [ru, ru.*u+p, rv.*u,   ru.*h].';
f2 = [rv, ru.*v,   rv.*v+p, rv.*h].';

g1 = [r*v1 ru*v1 rv*v1 rE*v1].';
g2 = [r*v2 ru*v2 rv*v2 rE*v2].';
f1 = f1 - g1;
f2 = f2 - g2;

f1 = simplify(f1);
f2 = simplify(f2);
df1= simplify([diff(f1,'r') diff(f1,'ru') diff(f1,'rv') diff(f1,'rE')]);
df2= simplify([diff(f2,'r') diff(f2,'ru') diff(f2,'rv') diff(f2,'rE')]);
dfn= simplify(df1*nx + df2*ny);
lam= eig(dfn);
lam= simplify(lam);

lam1 = (nx*u + ny*v) - (nx*v1 + ny*v2);
lam2 = (nx*u + ny*v) - (nx*v1 + ny*v2);
lam3 = (nx*u + ny*v) - (nx*v1 + ny*v2) - (gam*(gam-1)*(rE - 0.5*(u^2+v^2))/r)^(1/2);
lam4 = (nx*u + ny*v) - (nx*v1 + ny*v2) + (gam*(gam-1)*(rE - 0.5*(u^2+v^2))/r)^(1/2);
 