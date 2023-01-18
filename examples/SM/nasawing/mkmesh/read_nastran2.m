close all
clear all

fid = fopen('CRM_FEM.bdf','r');
BDF = textscan(fid, '%s','delimiter',sprintf('\n'));
BDF = BDF{1};
fclose all;

% 'a' the list of lines where the string 'GRID ' appears
a = strfind(BDF,'GRID '); a = find(not(cellfun('isempty', a)));
GRID = char(BDF(a));
GRID = [str2num(GRID(:,25:32)),str2num(GRID(:,33:40)),str2num(GRID(:,41:end))];
X = GRID(:,1); Y = GRID(:,2); Z = GRID(:,3);

% 'a' the list of lines where the string 'CQUAD ' appears
a = strfind(BDF,'CQUAD'); a = find(not(cellfun('isempty', a)));
CQUAD = char(BDF(a));
CQUAD = [str2num(CQUAD(:,9:16)),str2num(CQUAD(:,17:24)),str2num(CQUAD(:,25:32)),str2num(CQUAD(:,33:40)),str2num(CQUAD(:,41:48)),str2num(CQUAD(:,49:56))];

% Reordering of quads (usefull if for some reason the quads are not counted
% in a monotonous way)
IdQ = CQUAD(:,1);
if (max(abs(IdQ-[1:size(CQUAD,1)]'))~=0)
    CQUAD(IdQ,:) = CQUAD;
    fprintf('Quads have been reordered')
end

% Quad to Shell connectivity
quad2sh = CQUAD(:,2);

% Store only Quad Connectivities in CQUAD
CQUAD = CQUAD(:,3:end);

% 'a' the list of lines where the string 'PSHELL ' appears
a = strfind(BDF,'PSHELL'); a = find(not(cellfun('isempty', a)));
thickness = char(BDF(a));
% Thickness for shells and quads
shell2th  = str2num(thickness(:,25:32));
thickness = shell2th(quad2sh);


%%%%%%  --- INP PARAVIEW VISUALIZATION ---  %%%%%
% writeINPshells('VisuWing', GRID, CQUAD, quad2sh, thickness);

return;

%%%%%%  --- MATLAB VISUALIZATION PART  ---  %%%%%

% Find The four groups of shells (for visualization purpose)
a = strfind(BDF,'upper skin'); a = find(not(cellfun('isempty', a))); US = char(BDF(a)); US = str2num(US(:,20:end)); 
a = strfind(BDF,'lower skin'); a = find(not(cellfun('isempty', a))); LS = char(BDF(a)); LS = str2num(LS(:,20:end));
a = strfind(BDF,'spar'); a = find(not(cellfun('isempty', a))); SPAR = char(BDF(a)); SPAR = str2num(SPAR(:,14:end));
a = strfind(BDF,'rib'); a = find(not(cellfun('isempty', a))); RIB = char(BDF(a)); RIB = str2num(RIB(:,14:end));
ID = [US*0+0;LS*0+1;SPAR*0+2;RIB*0+3];
ID = ID(quad2sh);

% Full Shell model
figure(1); clf;
patch(X(CQUAD(ID==0,:))',Y(CQUAD(ID==0,:))',Z(CQUAD(ID==0,:))'+4,thickness(ID==0)','edgecolor','k')
patch(X(CQUAD(ID==1,:))',Y(CQUAD(ID==1,:))',Z(CQUAD(ID==1,:))'-4,thickness(ID==1)','edgecolor','k')
patch(X(CQUAD(ID==2,:))',Y(CQUAD(ID==2,:))',Z(CQUAD(ID==2,:))',thickness(ID==2)','edgecolor','k')
patch(X(CQUAD(ID==3,:))',Y(CQUAD(ID==3,:))',Z(CQUAD(ID==3,:))',thickness(ID==3)','edgecolor','k')
axis equal
axis tight
colormap(jet);
colorbar
view([75 27])

% Spars only Shell model
figure(2); clf;
patch(X(CQUAD(ID==2,:))',Y(CQUAD(ID==2,:))',Z(CQUAD(ID==2,:))',thickness(ID==2)','edgecolor','k')
axis equal
axis tight
colormap(jet);
colorbar
view([75 27])

% Ribs only Shell model
figure(3); clf;
patch(X(CQUAD(ID==3,:))',Y(CQUAD(ID==3,:))',Z(CQUAD(ID==3,:))',thickness(ID==3)','edgecolor','k')
axis equal
axis tight
colormap(jet);
colorbar
view([75 27])

clear('ID','X','Y','Z','a');
