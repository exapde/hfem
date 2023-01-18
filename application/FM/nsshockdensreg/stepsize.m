function s = stepsize(udg,uh,param)

gam = param{1};

r    = udg(:,1,:);
ru   = udg(:,2,:);
rv   = udg(:,3,:);
rE   = udg(:,4,:);

r1   = 1./r;
u    = ru.*r1;
v    = rv.*r1;    
E    = rE.*r1;
q    = 0.5*(u.*u+v.*v);
p    = (gam-1)*(rE-r.*q);    
h    = E+p.*r1;    

s = 1;
if (min(r(:))<=0)
    fprintf('Negative density in UDG:   %e \n', min(r(:)));
    s = 0;
end
if (min(p(:))<=0)
    fprintf('Negative pressure in UDG:  %e \n', min(p(:)));
    s = 0;
end
if (min(E(:))<=0)
    fprintf('Negative energy in UDG:  %e \n', min(E(:)));
    s = 0;
end
if (min(h(:))<=0)
    fprintf('Negative enthalpy in UDG:  %e \n', min(h(:)));
    s = 0;
end

% if (min(r(:))<=0) || (min(rE(:))<=0) || (min(E(:))<=0) || (min(p(:))<=0) || (min(h(:))<=0)
%     warning('Negative states in UDG');
%     s = 0;
%     [min(r(:)) min(rE(:)) min(E(:)) min(p(:)) min(h(:))]
%     ind = find(p(:)==min(p(:)));
%     [rE(ind) r(ind) q(ind) u(ind) v(ind) ru(ind) rv(ind)]    
% end
    
r    = uh(1,:);
ru   = uh(2,:);
rv   = uh(3,:);
rE   = uh(4,:);

r1   = 1./r;
u    = ru.*r1;
v    = rv.*r1;    
E    = rE.*r1;
q    = 0.5*(u.*u+v.*v);
p    = (gam-1)*(rE-r.*q);    
h    = E+p.*r1;    

s = 1;
if (min(r(:))<=0)
    fprintf('Negative density in UH:   %e \n', min(r(:)));
    s = 0;
end
if (min(p(:))<=0)
    fprintf('Negative pressure in UH:  %e \n', min(p(:)));
    s = 0;
end
if (min(E(:))<=0)
    fprintf('Negative energy in UH:  %e \n', min(E(:)));
    s = 0;
end
if (min(h(:))<=0)
    fprintf('Negative enthalpy in UH:  %e \n', min(h(:)));
    s = 0;
end

% [min(r(:)) min(rE(:)) min(E(:)) min(p(:)) min(h(:))]
% if (min(r(:))<=0) || (min(rE(:))<=0) || (min(E(:))<=0) || (min(p(:))<=0) || (min(h(:))<=0)
%     warning('Negative states in UH');
%     s = 0;
% end
% 
