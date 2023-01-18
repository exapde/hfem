function mesh = mkmesh_annular(n,nr,nt,R,r,h,porder,elemtype)
%MKMESH_ANNULAR Creates 3D mesh data structure for a slit annular plate.
%   MESH=MKMESH_ANNULAR(M,N,PORDER)
%
%      MESH:      Mesh structure
%      N:         Number of points in the angular direction 
%      NR:        Number of points in the radial direction
%      NT:        Number of points in the thickness direction
%      R:         Exterior radius
%      r:         Interior radius
%      h:         Thickness of the annular plate
%      PORDER:    Polynomial Order of Approximation (default=1)
%      ELEMTYPE:  Flag determining element type
%                 Flag = 0 tri/tet elements (default)
%                 Flag = 1 quad/hex elements
%   See also: CUBEMESH, MKMESH

if nargin<7, porder = 1; end
if nargin<8, elemtype=0; end

% Undeformed unit square mesh
[p,t] = cubemesh(nr,n,nt,elemtype);

% Definition of limit conditions
bndexpr = {'all(p(:,1)<1e-6)','all(p(:,1)>max(p0(:,1))-1e-6)', ...
           'all(p(:,2)<1e-6)','all(p(:,2)>max(p0(:,2))-1e-6)', ...
           'all(p(:,3)<1e-6)','all(p(:,3)>max(p0(:,3))-1e-6)'};  

nodetype = 0;
mesh = mkmesh(p,t,porder,bndexpr,elemtype,nodetype);


RR = r + (R-r)*mesh.p(:,1);
t = mesh.p(:,2);
mesh.p(:,1) = RR.*cos(2*pi*t);
mesh.p(:,2) = RR.*sin(2*pi*t);
mesh.p(:,3) = h*mesh.p(:,3);

RR = r + (R-r)*mesh.dgnodes(:,1,:);
t = mesh.dgnodes(:,2,:);
mesh.dgnodes(:,1,:) = RR.*cos(2*pi*t);
mesh.dgnodes(:,2,:) = RR.*sin(2*pi*t);
mesh.dgnodes(:,3,:) = h*mesh.dgnodes(:,3,:);



% xzmap = [[0,1,1.2,1+0.5]' [0,0,1,1]'];
% 
% pxz = mesh.p(:,[1 3]);
% pxz(:,1) = (pxz(:,1)-pmin(1))/(pmax(1)-pmin(1));
% pxz(:,2) = (pxz(:,2)-pmin(3))/(pmax(3)-pmin(3));
% pxz = mapp(pxz,xzmap);
% pxz(:,1) = pmin(1) + (pmax(1)-pmin(1))*pxz(:,1);
% pxz(:,2) = pmin(3) + (pmax(3)-pmin(3))*pxz(:,2);
% mesh.p(:,1) = pxz(:,1);
% mesh.p(:,3) = pxz(:,2);
% 
% x = mesh.dgnodes(:,1,:);
% z = mesh.dgnodes(:,3,:);
% pxz = [x(:) z(:)];
% pxz(:,1) = (pxz(:,1)-pmin(1))/(pmax(1)-pmin(1));
% pxz(:,2) = (pxz(:,2)-pmin(3))/(pmax(3)-pmin(3));
% pxz = mapp(pxz,xzmap);
% pxz(:,1) = pmin(1) + (pmax(1)-pmin(1))*pxz(:,1);
% pxz(:,2) = pmin(3) + (pmax(3)-pmin(3))*pxz(:,2);
% mesh.dgnodes(:,1,:) = reshape(pxz(:,1), size(mesh.dgnodes(:,1,:)));
% mesh.dgnodes(:,3,:) = reshape(pxz(:,2), size(mesh.dgnodes(:,1,:)));








