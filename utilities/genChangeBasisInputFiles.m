
% This script computes and writes to a binary file the change of basis
% matrix from orthogonal to nodal basis (and viceversa) for different
% polynomial orders and element types

% We assume that the orthogonal basis are ordered in non-decreasing
% polynomial order.

porderMin = 0;
porderMax = 10;
nodetype = 0;
endian = 'native'; % Options: 'native', 'ieee-le' , 'ieee-be'

if nodetype ~= 0; error('There may be underlying assumptions in C++ code on nodetype == 0.'); end

for porder = porderMin:porderMax
    
    % Line:
    [plocal,~,~,~,~,~,~] = masternodes(porder,1,0,nodetype);
    orth2nodal = tensorproduct(plocal,porder);
    nodal2orth = inv(orth2nodal);
    fileID = fopen(['basisChange_line_p', num2str(porder),'.bin'],'w');
    fwrite(fileID,orth2nodal,'double',endian);
    fwrite(fileID,nodal2orth,'double',endian);
    fclose(fileID);
    
    % Tris:
    [plocal,~,~,~,~,~,~] = masternodes(porder,2,0,nodetype);
    orth2nodal = koornwinder(plocal,porder);
    nodal2orth = inv(orth2nodal);
    fileID = fopen(['basisChange_tri_p', num2str(porder),'.bin'],'w');
    fwrite(fileID,orth2nodal,'double',endian);
    fwrite(fileID,nodal2orth,'double',endian);
    fclose(fileID);
    
    % Quads:
    [plocal,~,~,~,~,~,~] = masternodes(porder,2,1,nodetype);
    orth2nodal = tensorproduct(plocal,porder);
    nodal2orth = inv(orth2nodal);
    fileID = fopen(['basisChange_quad_p', num2str(porder),'.bin'],'w');
    fwrite(fileID,orth2nodal,'double',endian);
    fwrite(fileID,nodal2orth,'double',endian);
    fclose(fileID);
    
    % Tets:
    [plocal,~,~,~,~,~,~] = masternodes(porder,3,0,nodetype);
    orth2nodal = koornwinder(plocal,porder);
    nodal2orth = inv(orth2nodal);
    fileID = fopen(['basisChange_tet_p', num2str(porder),'.bin'],'w');
    fwrite(fileID,orth2nodal,'double',endian);
    fwrite(fileID,nodal2orth,'double',endian);
    fclose(fileID);
    
    % Hexas:
    [plocal,~,~,~,~,~,~] = masternodes(porder,3,1,nodetype);
    orth2nodal = tensorproduct(plocal,porder);
    nodal2orth = inv(orth2nodal);
    fileID = fopen(['basisChange_hex_p', num2str(porder),'.bin'],'w');
    fwrite(fileID,orth2nodal,'double',endian);
    fwrite(fileID,nodal2orth,'double',endian);
    fclose(fileID);
end
