function [IMEXVhat] = getVhatIMEX(mesh,UH)
%GETVHATIMEX gets the hybrid variable Vhat from the implicit solution on
%the IMEX boundary. This hybrid variable will later be used by the Explicit
%resolution as a Dirichlet boundary condition.
%This code works only for the scalar wave equation in its present state.
%
% SYNTAX:  [IMEXQHAT] = GETVHATIMEX(MASTER,MESH,APP,UDG)
%
% INPUTS:
%    MESH      - Mesh structure (of Implicit domain)
%    UDG       - Vector of Hybrid unknowns from Implicit domain
%
% OUTPUTS:
%    IMEXVhat  - Hybrid variable on the IMEX boundary. 
%
% SEE ALSO: GETQHATIMEX
% 
% Author(s): Lauren Kolkman, Sebastien Terrana
% April 2018

% Get some dimensions
nch = size(UH,1);
[npf,~] = size(mesh.perm);
% IMEX faces
fimex = mesh.fimex;
nfimex = length(fimex);
% Initialize field IMEXVhat
IMEXVhat = zeros(nch,npf,nfimex);

% Filling IMEXVhat with hybrid variables on the IMEX faces
for i=1:nfimex
    pos = (fimex(i)-1)*npf+1;
    IMEXVhat(:,:,i) = UH(:,pos:pos+npf-1);
end

end