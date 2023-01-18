with(linalg):
with(codegen):
with(CodeGeneration):

FileTools[Remove]("sasource2d.c"):;

U:=vector(15):
r:=U[1]:   rx:=U[6]:   ry:=U[11]:
ru:=U[2]: rux:=U[7]:  ruy:=U[12]:
rv:=U[3]: rvx:=U[8]:  rvy:=U[13]:
rE:=U[4]: rEx:=U[9]:  rEy:=U[14]:
rN:=U[5]: rNx:=U[10]: rNy:=U[15]:

#param:=vector(2):
#unprotect(Re):
#Re:=param[1]:
#Pr:=param[2]:

#ga:=7/5;
#cv1:=7.1;
#cb1:=0.1355;
#cb2:=0.622;
#sigma:=2/3;
#kappa:=0.41;
#cw1:=3.2391;
#cw2:=0.3;
#cw3:=2;

# Define Source

u:=ru/r:
v:=rv/r:
N:=rN/r:
Nx:=(rNx - rx*N)/r:
Ny:=(rNy - ry*N)/r:
vx:=(rvx - rx*v)/r:
uy:=(ruy - ry*u)/r:

%Ds:=(cb2/sigm)*r*(Nx+Ny)*(Nx+Ny):
Ds:=(cb2/sigm)*(rNx*rNx+rNy*rNy)/r:

Om:=((uy-vx)*(uy-vx) + 1.0e-20)^(0.5):
chi:=rN*Re:
fv1:=chi^3/(chi^3+cv1^3):
fv2:=1.0-chi/(1.0+chi*fv1):
Nk:=N/(kappa^2*d^2):
Sh:=Om + Nk*fv2:
Sp:=cb1*Sh*rN:

rd:=Nk/Sh:
gd:=rd + cw2*(rd^(6.0) - rd):
fw:=gd*((1+cw3^6)/(gd^6 + cw3^6))^(1.0/6.0):
Sd:=cw1*fw*(rN/d)^2/r:

Sm:=Ds+Sp-Sd:


#Omega:=((uy-vx)*(uy-vx) + 1e-20)^(1/2): 
#chi:=rN*Re:
#fv1:=chi^3/(chi^3+cv1^3):
#fv2:=1.0-chi/(1.0+chi*fv1):
#ft1:=r*kappa^2*dd^2:
#Omega0:=Omega+rN*fv2/ft1:
#Omega1:=(Omega0^6 + (0.05*Omega)^6)^(1/6):
#rb:=rN/(Omega1*ft1):
#ra:=rb*cw4/((rb^6 + cw4^6)^(1/6)):
#g:=ra+cw2*(ra^6-ra):
#ftt:=(1+cw36)/(g^6+cw36):
#fw:=g*(ftt^(1/6)):
#Sa:=cb1*rN*Omega1:
#Sb:=(cb2/sigma)*(rNx*rNx+rNy*rNy)/r:
#Sc:=cw1*fw*(rN/dd)^2/r:
#Sm:=Sa+Sb-Sc:


S:=evalm([Sm]):
JS:=jacobian(S,U):

# Generate C code 

S:=makevoid(makeproc(S,[U,S])):

JS:=convert(evalm(transpose(JS)),vector):
JS:=makevoid(makeproc(JS,JS)):

sasource2d:=makevoid(joinprocs([S,JS])):

sasource2d:=declare(U::'array'(1..15,numeric),sasource2d):
sasource2d:=declare(S::'array'(1..1,numeric),sasource2d):
sasource2d:=declare(JS::'array'(1..15,numeric),sasource2d):

C(sasource2d,optimize,deducetypes=false,output="sasource2d.c"):


