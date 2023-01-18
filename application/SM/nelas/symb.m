% Plane Strain Neo-Hookean Material Models

syms F11 F12 F21 F22 pres mu kappa lambda

F=[F11 F12; F21 F22];

JF=det(F);
IF=inv(F);
TF=F.';
VF=IF.';
CF=F.'*F;

wcase=1;

% Strain energy function and the volumetric component of the first
% Piola-Kirchoff Tensor
switch (wcase)
    case 1
        W  = (mu/2)*(trace(CF)+1-3) - mu*log(JF) + (kappa/2)*(JF-1)^2;
        Wi = (mu/2)*(trace(CF)+1-3) - mu*log(JF);
        Wv = (kappa/2)*(JF-1)^2;
        Pv = pres*JF*VF;
        g  = JF-1; 
    case 2
        W  = (mu/2)*(trace(CF)+1-3) - mu*log(JF) + (kappa/2)*(log(JF))^2;
        Wi = (mu/2)*(trace(CF)+1-3) - mu*log(JF);
        Wv = (kappa/2)*(log(JF))^2;
        Pv = pres*VF;
        g  = log(JF);
    case 3
        W  = (mu/2)*(JF^(-2/3)*(trace(CF)+1)-3) + (kappa/2)*(JF-1)^2;
        Wi = (mu/2)*(JF^(-2/3)*(trace(CF)+1)-3);
        Wv = (kappa/2)*(JF-1)^2;
        Pv = pres*JF*VF;
        g  = JF-1; 
    case 4
        W  = (mu/2)*(JF^(-2/3)*(trace(CF)+1)-3) + (kappa/2)*(log(JF))^2;
        Wi = (mu/2)*(JF^(-2/3)*(trace(CF)+1)-3);
        Wv = (kappa/2)*(log(JF))^2;
        Pv = pres*VF;
        g  = log(JF);
    otherwise
        error('Not a valid model');
end

% Isochoric component of the first Piola-Kirchoff Tensor
Pi(1,1)=diff(Wi,'F11');
Pi(1,2)=diff(Wi,'F12');
Pi(2,1)=diff(Wi,'F21');
Pi(2,2)=diff(Wi,'F22');

% the first Piola-Kirchoff Tensor
P = Pi + Pv;

% Derivativies of the first Piola-Kirchoff Tensor with respect to the
% deformation gradient
D=P.'; D=D(:); 
dDdF=[diff(D(1),'F11');  diff(D(2),'F11'); diff(D(1),'F12'); diff(D(2),'F12');...
      diff(D(1),'F21');  diff(D(2),'F21'); diff(D(1),'F22'); diff(D(2),'F22');...
      diff(D(3),'F11');  diff(D(4),'F11'); diff(D(3),'F12'); diff(D(4),'F12');...
      diff(D(3),'F21');  diff(D(4),'F21'); diff(D(3),'F22'); diff(D(4),'F22')];
dPdF=simplify(dDdF);  

% Derivativies of the first Piola-Kirchoff Tensor with respect to the
% pressure
dPdp=diff(P,'pres');
dPdp=simplify(dPdp);  

dgdF = simplify([diff(g,'F11') diff(g,'F12') diff(g,'F21') diff(g,'F22')]);




