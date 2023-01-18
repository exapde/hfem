function [] = mapEqual(A,B)
tol = 1e-6;

if(numel(A) ~= numel(B))
    size(A)
    size(B)
    error('numel(A) ~= numel(B) at mapEqual.m');
elseif (ndims(A) ~= ndims(B))
    size(A)
    size(B)
    error('ndims(A) ~= ndims(B) at mapEqual.m');
elseif norm(size(A)-size(B)) > 1e-6
    size(A)
    size(B)     
    error('norm(size(A)-size(B)) > 1e-6 at mapEqual.m');
end

amax = max(abs(A(:)));
maxDiff = max( abs( A(:) - B(:) ) );
if(maxDiff > tol*amax & maxDiff > tol)
    size(A)
    size(B)
    fprintf('WARNING: A and B are different (within a tolerance) %e %e',maxDiff,amax);
    pause
end
