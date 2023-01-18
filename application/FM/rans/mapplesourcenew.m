with(linalg):
with(codegen):
with(CodeGeneration):

FileTools[Remove]("sasource2dnew.c"):;

U:=vector(15):
r:=U[1]:   rx:=U[6]:   ry:=U[11]:
ru:=U[2]: rux:=U[7]:  ruy:=U[12]:
rv:=U[3]: rvx:=U[8]:  rvy:=U[13]:
rE:=U[4]: rEx:=U[9]:  rEy:=U[14]:
rN:=U[5]: rNx:=U[10]: rNy:=U[15]:


# Define Source

u:=ru/r:
v:=rv/r:
N:=rN/r:
Nx:=(rNx - rx*N)/r:
Ny:=(rNy - ry*N)/r:
vx:=(rvx - rx*v)/r:
uy:=(ruy - ry*u)/r:

chi:=rN*Re:
shi:=log(1 + exp(b*chi))/b:
phi:=exp(b*chi)/(1+exp(b*chi)):
#chix:=rNx*Re:
#chiy:=rNy*Re:
#shix:=phi*chix:
#shiy:=phi*chiy:

#Ds:=(cb2/sigm)*r*(Nx+Ny)*(Nx+Ny):
Ds:=(cb2/sigm)*(rNx*rNx+rNy*rNy)/r:
#Ds:=(cb2/(sigm*Re*Re))*(shix*shix+shiy*shiy)/r:
#Ds:=(phi*phi*cb2/sigm)*(rNx*rNx+rNy*rNy)/r:

Om:=((uy-vx)*(uy-vx) + 1.0e-20)^(0.5):
fv1:=shi^3/(shi^3+cv1^3):
fv2:=1.0-shi/(1.0+shi*fv1):
Nk:=(shi/(r*Re))/(kappa^2*d^2):
Sh:=Om + Nk*fv2:
Sp:=cb1*Sh*(shi/Re):

rd:=Nk/Sh:
gd:=rd + cw2*(rd^(6.0) - rd):
fw:=gd*((1+cw3^6)/(gd^6 + cw3^6))^(1.0/6.0):
Sd:=cw1*fw*(rN/d)^2/r:

Sm:=Ds+Sp-Sd:

S:=evalm([Sm]):
JS:=jacobian(S,U):

# Generate C code 

S:=makevoid(makeproc(S,[U,S])):

JS:=convert(evalm(transpose(JS)),vector):
JS:=makevoid(makeproc(JS,JS)):

sasource2dnew:=makevoid(joinprocs([S,JS])):

sasource2dnew:=declare(U::'array'(1..15,numeric),sasource2dnew):
sasource2dnew:=declare(S::'array'(1..1,numeric),sasource2dnew):
sasource2dnew:=declare(JS::'array'(1..15,numeric),sasource2dnew):

C(sasource2dnew,optimize,deducetypes=false,output="sasource2dnew.c"):


