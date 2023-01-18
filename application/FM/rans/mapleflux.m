with(linalg):
with(codegen):
with(CodeGeneration):

FileTools[Remove]("saflux2d.c"):;

U:=vector(15):
r:=U[1]:   rx:=U[6]:   ry:=U[11]:
ru:=U[2]: rux:=U[7]:  ruy:=U[12]:
rv:=U[3]: rvx:=U[8]:  rvy:=U[13]:
rE:=U[4]: rEx:=U[9]:  rEy:=U[14]:
rN:=U[5]: rNx:=U[10]: rNy:=U[15]:


# Define Flux

p:=(gam-1)*(rE-1/2*(ru^2+rv^2)/r):
Fix:=evalm([ru,ru^2/r+p,ru*rv/r,ru*(rE+p)/r,rN*ru/r]):
Fiy:=evalm([rv,ru*rv/r,rv^2/r+p,rv*(rE+p)/r,rN*rv/r]):

u:=ru/r:
v:=rv/r:
E:=rE/r:
N:=rN/r:

ux:=(rux-rx*u)/r:
vx:=(rvx-rx*v)/r:
Ex:=(rEx-rx*E)/r:
Nx:=(rNx-rx*N)/r:

uy:=(ruy-ry*u)/r:
vy:=(rvy-ry*v)/r:
Ey:=(rEy-ry*E)/r:
Ny:=(rNy-ry*N)/r:

txx:=2/3*(2*ux-vy):
txy:=vx+uy:
tyy:=2/3*(2*vy-ux):
ex:=Ex-u*ux-v*vx:
ey:=Ey-u*uy-v*vy:

chi:=rN*Re:
fv1:=chi^3/(chi^3+cv1^3):
muN:=rN*fv1:
muT:=muN+1/Re:
 fs:=muT/sigm:

Fvx:=evalm([0,txx*muT,txy*muT,(txx*u+txy*v+gam/Pr*ex)*muT,fs*Nx]):
Fvy:=evalm([0,txy*muT,tyy*muT,(txy*u+tyy*v+gam/Pr*ey)*muT,fs*Ny]):

Jix:=jacobian(Fix,U):
Jiy:=jacobian(Fiy,U):
Jvx:=jacobian(Fvx,U):
Jvy:=jacobian(Fvy,U):

# Generate C code 

Jix:=convert(evalm(transpose(Jix)),vector):
Jiy:=convert(evalm(transpose(Jiy)),vector):
Jvx:=convert(evalm(transpose(Jvx)),vector):
Jvy:=convert(evalm(transpose(Jvy)),vector):

Fix:=makevoid(makeproc(Fix,[U,Fix])):
Fiy:=makevoid(makeproc(Fiy,Fiy)):
Fvx:=makevoid(makeproc(Fvx,Fvx)):
Fvy:=makevoid(makeproc(Fvy,Fvy)):
Jix:=makevoid(makeproc(Jix,Jix)):
Jiy:=makevoid(makeproc(Jiy,Jiy)):
Jvx:=makevoid(makeproc(Jvx,Jvx)):
Jvy:=makevoid(makeproc(Jvy,Jvy)):

saflux2d:=makevoid(joinprocs([Fix,Fiy,Fvx,Fvy,Jix,Jiy,Jvx,Jvy])):

saflux2d:=declare(U::'array'(1..15,numeric),saflux2d):
saflux2d:=declare(Fix::'array'(1..5,numeric),saflux2d):
saflux2d:=declare(Fiy::'array'(1..5,numeric),saflux2d):
saflux2d:=declare(Fvx::'array'(1..5,numeric),saflux2d):
saflux2d:=declare(Fvy::'array'(1..5,numeric),saflux2d):
saflux2d:=declare(Jix::'array'(1..75,numeric),saflux2d):
saflux2d:=declare(Jiy::'array'(1..75,numeric),saflux2d):
saflux2d:=declare(Jvx::'array'(1..75,numeric),saflux2d):
saflux2d:=declare(Jvy::'array'(1..75,numeric),saflux2d):

C(saflux2d,optimize,deducetypes=false,output="saflux2d.c"):


