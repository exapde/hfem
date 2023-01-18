fprintf('  --> Initializing DIGASO')
d0=fileparts([pwd,filesep]);
addpath([d0,'/mesh/distmesh']);
addpath([d0,'/utilities']);
addpath([d0,'/mesh/mkmesh']);
addpath([d0,'/mesh/cmesh']);
addpath([d0,'/mesh/airfoilTools']);
addpath([d0,'/mesh/airfoilTools/geometries']);
addpath([d0,'/mesh']);
addpath([d0,'/kernel']);
addpath([d0,'/master']);
addpath([d0,'/postprocessing']);
addpath([d0,'/preprocessing']);
fprintf(' \n Done. \n\n');
clear d0
