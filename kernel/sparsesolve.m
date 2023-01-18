function DUH = sparsesolve(AE, FE, f, f2f, iterative, hybrid)
% SPARSESOLVE  Solve the linear system with eliminating dofs on the domain boundary
%
%   DUH = SPARSESOLVE(AE, FE, f, f2f, iterative, hybrid)
%
%      AE:        Matrix in Matlab-sparse format
%      FE:        Vector
%      f:         face-to-element connectivities
%      f2f:       face-to-face connectivities
%      iterative: flag to determine solver used
%      hybrid:    flag     
%
%      DUH:      Output vector

if strcmp(hybrid,'edg') || strcmp(hybrid,'hedg') 
    DUH = AE\FE;
    return;
end

nf  = size(f,1);            % Number of faces
ndf = size(AE,2);           % Totaf number of DOF on all faces
ncf = ndf/nf;               % Number of DOF per face
fext= (find(f(:,end) < 0)); % External faces
fint= (find(f(:,end) >-1)); % Internal faces
nex = length(fext);         % Number of Ext faces
nin = length(fint);         % Number of Int faces

fext = ncf*(fext-1) * ones(1,ncf) + ones(nex,1) * (1:ncf);
fext = reshape(fext',[ncf*nex,1]);
fint = ncf*(fint-1) * ones(1,ncf) + ones(nin,1) * (1:ncf);
fint = reshape(fint',[ncf*nin,1]);

% eliminate dofs on the boundary
[D,G] = mkDG(AE(:,:), FE(fext), f, f2f, fext, fint);
% [D,G] = mkDG(AE(1:n,:), FE(1:n), f, f2f, m);
% D = AE(1:n,1:n)\AE(1:n,n+1:end);
% G = AE(1:n,1:n)\FE(1:n);
% [max(abs(D(:)-D1(:))) max(abs(G(:)-G1(:)))]
% pause

% obtain the system for the interior dofs
FE = (FE(fint)-AE(fint,fext)*G);
AE = (AE(fint,fint)-AE(fint,fext)*D);

if iterative == 0
    % solve for dofs on the interior faces
    DUHI = AE\FE;               
else                 
    % compute block-Jacobi preconditioner    
    ME  = zeros(ncf,ncf,nf-m);    
    for i = 1:(nf-m)      
        r = (i-1)*ncf+1:i*ncf;    
        ME(:,:,i) = inv(full(AE(r,r)));  
        FE(r) = ME(:,:,i)*FE(r); 
    end    
    
    % solve for dofs on the interior faces using gmres
    DUHI = sparsegmres(AE, ME, FE);           
end

% compute dofs on the boundary
DUHB = G - D*DUHI;    
    
% get the output vector
DUH = sparse(ndf,1);
DUH(fint) = DUHI;
DUH(fext) = DUHB;
%DUH  = [DUHB; DUHI];                           


function [D,G] = mkDG(AE, F, f, f2f, fext, fint)

nf  = size(f,1);
ndf = size(AE,2);
ncf = ndf/nf;
nef = (1:ncf)';
facext = (find(f(:,end) < 0));
m   = length(facext);

ind = [];
for i = 1:m                      % loop through each boundary faces
    fe = facext(i);                % current external face
    fi = f2f(:,fe);              % neighboring faces to face i 
    bf = (fi > 0);               % Checking nonzero faces
    fi = fi(bf);
    bf = (f(fi,end) < 0);        % Checking neighboring faces on the boundary  
    fj = fi(bf);
    if any(ismember(fj,ind))==0  % do if faces in fj are not in the list ind
        ind = [ind; fj];         % list of all boundary faces considered so far
        ff = sort(fj);           % sort faces in fj
        dd = bsxfun(@plus,(ff-1)'*ncf,nef);
        AE(dd(:),dd(:)) = inv(AE(dd(:),dd(:))); 
    end
end
G = AE(fext,fext)*F(:);
D = AE(fext,fext)*AE(fext,fint);


function [x, flags, iter, rev] = sparsegmres(A, M, b, x, restart, tol, maxit)

% check the number of input arguments
if nargin < 2
  error('Not enough input arguments.');
end

% get the dimension from b
N = size(b(:),1);

% Default parameters
if nargin < 4 || isempty(x),       x = zeros(N,1);  end;
if nargin < 5 || isempty(restart), restart = 20;    end;
if nargin < 6 || isempty(tol),     tol = 1e-7;      end;
if nargin < 7 || isempty(maxit),   maxit = min(100*restart,N);  end;

% initialization
nrmb   = norm(b); 
flags  = inf; 
iter   = 0; 
cycle  = 0;

% allocate memory
hh = zeros(restart+1,restart);
v  = zeros(N,restart+1);
e1 = zeros(restart+1,1); e1(1)=1;
rev = zeros(restart,1);

n1 = size(M,1);
n2 = size(M,3);

while (1) 
    % perform matrix-vector multiplication              
    d = reshape(A*x,[n1 n2]);
%     for i=1:n2
%         d(:,i) = M(:,:,i)*d(:,i);
%     end
    d = arraymul(M,d);
    
    % compute the residual vector      
    r = b(:) - d(:);        
    beta = norm(r);
    v(:,1) = r/beta;

    res  = beta;
    iter = iter+1;
    rev(iter) = res;    
    for j = 1:restart

        % set flag=0 if convergence
        if res/nrmb <= tol
            flags = 0;
            break;
        end        
        
        % set flag=1 if reaching maximum iteration 
        if iter>=maxit
            flags = 1;
            break;
        end                
                 
        % perform matrix-vector multiplication
        d = reshape(A*v(:,j),[n1 n2]);
        d = arraymul(M,d);
%         for i=1:n2
%             d(:,i) = M(:,:,i)*d(:,i);
%         end
    
        % Arnoldi process (i.e., GS orthogonalization on the Krylov supspace)        
        v(:,j+1) = d(:);
        for i = 1:j        
            hh(i,j) = v(:,i)'*v(:,j+1);
            v(:,j+1) = v(:,j+1) - hh(i,j)*v(:,i);        
        end
        hh(j+1,j) = norm(v(:,j+1));     
        if (hh(j+1,j) ~= 0.0)          
            v(:,j+1) = v(:,j+1)/hh(j+1,j);   
        else            
            break;
        end

        % solve the reduced system 
        y = hh(1:j+1,1:j)\(beta*e1(1:j+1));
        
        % compute the residual norm
        res = norm(beta*e1(1:j+1)-hh(1:j+1,1:j)*y);
        
        iter = iter + 1; 
        rev(iter) = res;
    end 
      
    % compute the solution    
    x(:) = x(:) + v(:,1:length(y))*y;    
    
    cycle = cycle + 1;

    % stop if converging or reaching the maximum iteration
    if flags < inf, 
        fprintf('gmres(%d) converges at outer iteration %d to a solution with relative residual %g\n', [restart cycle res/nrmb]);       
        break; 
    end;     
end
rev = rev/nrmb;

