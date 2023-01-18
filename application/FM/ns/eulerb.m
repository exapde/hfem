function fn = eulerb(up,np,ib,ui,p,param,time)
%EULERB Calculate the boundary flux for the Euler equations.
%   FN=EULERB(UP,NP,IB,UI,P,PARAM,TIME)
%
%      UP(np,3):   np plus states
%      NP(np,2):   np normal plus vectors 
%      IB:         Boundary type
%                  - IB: 1 Far-Field (Radiation)
%                  - IB: 2 Solid Wall(Reflection)
%                  - IB: 3 Non Homogenous Far-Filed (Incoming Wave)
%      UI(3):      Infinity state associated with IB
%      P(np,2):    np x,y coordinates
%      PARAM{1}:   Cell array containing either the wave speed c=PARAM{1}
%      TIME:       Time
%      FN(np,3):   np normal fluxes (f plus)  
%                          
% - Written by: J. Peraire
% 
if ib == 1                 % Far field
    um = repmat(ui,size(up,1),1);
elseif  ib == 2            % Reflect
    un = up(:,2).*np(:,1)+up(:,3).*np(:,2);
    um = [up(:,1),up(:,2)-2*un.*np(:,1),up(:,3)-2*un.*np(:,2),up(:,4)];
end 

fn = euleri_roe(up,um,np,p,param,time);
