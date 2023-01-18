function [p,t] = mkcgpt(dgnodes)

% remove duplicate nodes in mesh.p1
[ns,dim,nt]=size(dgnodes);
A = reshape(permute(dgnodes,[1,3,2]),[ns*nt,dim]);
snap = 1e-8;
A = round(A/snap)*snap;
B = flip(A,1);
[~,I] = unique(B,'rows'); 
B = B(sort(I),:);
B = flip(B,1);
[~,b] = ismember(A,B,'rows');

% CG mesh
p = B;
t = reshape(b,[ns nt])';



