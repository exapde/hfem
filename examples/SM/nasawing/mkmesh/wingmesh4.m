
% Get Datas from Nastran file
run read_nastran2.m

tol = 1e-6;
p = round(GRID/1e-8)*1e-8;

% Increase Thickness for visualization purposes
% thickness = 2*thickness;

% Nunber of Quads-Shells
nquad = size(CQUAD,1);
nv    = size(CQUAD,2);
% Number of nodes
np = size(p,1);

% Create connectivity edge-->quads and edge-->points
[edge2p,edge2quad,quad2edge] = mkQ2e(CQUAD);

% Break Shells at intersections to create new shells
[quad2sh,shell2th,e2b,ib] = break_shells(edge2quad,quad2edge,quad2sh,shell2th);

% Check Breaking of shells
% writeINPshells('VisuWingShBreak', GRID, CQUAD, quad2sh, thickness);

% Create topological corners (beam connection)
corner = mkCorners(CQUAD,edge2p,e2b,np);

% Number of shells, beams, corners
nsh = size(shell2th,1);
nbe = length(ib);
nco = length(corner);

% Create volumes plates from shells
hex2sh = []; sh2quad = {};
pp = cell(nsh,1);
tp = cell(nsh,1);
for i = 1:nsh
    % Shell to Quads connectivity
    sh2quad{i} = find(quad2sh==i);
    if (isempty(sh2quad{i})); continue; end
    % Extrude each shell to 3D plates
    [pp{i}, tp{i}] = mkplate(p, CQUAD(sh2quad{i},:), shell2th(i), 0);
end


% Create beams at shell intersections
for i = 1:nbe
    [pb{i}, tb{i}, ~, pp] = mkbeam_new(p, CQUAD, pp, sh2quad, shell2th, ib{i}, tol);
end

% Create corners cells at beams intersections
% nco = 4; % A SUPPRIMER !!!!!
for i = 1:nco
    [pc{i}, tc{i}, pp, pb] = mkCornerCell(p, pp, tp, pb, tb, ib, corner{i}, shell2th);
end

ptot = []; ttot = []; ntot = 0;
% Add nodes of each plates
for i = 1:nsh
    ptot   = [ptot; pp{i}];
    ttot   = [ttot; tp{i}+ntot];
    hex2sh = [hex2sh; i*ones(size(tp{i},1),1)];
    ntot = ntot+length(pp{i});
end

% Add nodes of each beams
for i = 1:nbe
    ptot   = [ptot; pb{i}];
    ttot   = [ttot; tb{i}+ntot];
    hex2sh = [hex2sh; (nsh+i)*ones(size(tb{i},1),1)];
    ntot = ntot+length(pb{i});
end

% Add nodes of each corners
for i = 1:nco
    ptot   = [ptot; pc{i}];
    ttot   = [ttot; tc{i}+ntot];
    hex2sh = [hex2sh; nsh+nbe+1];
    ntot = ntot+length(pc{i});
end

% Suppression of dupplicated nodes
[ptot,ttot]=fixmesh2(ptot,ttot);

% Output Visu Paraview
writeINPshells('VisuWingShExtruCorners', ptot, ttot, hex2sh, zeros(size(ttot,1),1));


return;

