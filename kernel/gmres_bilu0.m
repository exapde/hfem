function [x, flags, iter, rev] = gmres_bilu0(mesh, M, A, b, x, restart, tol, maxit)

% check the number of input arguments
if nargin < 4
  error('Not enough input arguments.');
end

% get the dimension from b
N = size(b(:),1);

% Default parameters
if nargin < 5 || isempty(x),       x = zeros(N,1);  end;
if nargin < 6 || isempty(restart), restart = 20;    end;
if nargin < 7 || isempty(tol),     tol = 1e-6;      end;
if nargin < 8 || isempty(maxit),   maxit = N;       end;

b = applyBILU0(mesh, M, b);

% initialization
nrmb   = norm(b(:)); 
flags  = inf; 
iter   = 0; 
cycle  = 0;

% allocate memory
hh = zeros(restart+1,restart);
v  = zeros(N,restart+1);
e1 = zeros(restart+1,1); e1(1)=1;
rev = zeros(restart,1);

while (1) 
    % perform matrix-vector multiplication        
    r = matvec(mesh, A,x);
    
    % Apply preconditioner
    r = applyBILU0(mesh, M, r);
    
    r = b - r;
    
    % Normalize the residual
    beta = norm(r(:));
    %if beta 
    
    v(:,1) = r(:)/beta;    
    
    res  = beta;
%     iter = iter+1;
%     rev(iter) = res;    
    for j = 1:restart            
        
        % perform matrix-vector multiplication        
        r = matvec(mesh,A,reshape(v(:,j),size(x)));
        
        % Apply preconditioner
        r = applyBILU0(mesh, M, r);        
                        
        % Arnoldi process (i.e., GS orthogonalization on the Krylov supspace)        
        v(:,j+1) = r(:);
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
        disp(res/nrmb)
        
        iter = iter + 1; 
        rev(iter) = res;
        
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
        
    end 
      
    % compute the solution    
    x(:) = x(:) + v(:,1:length(y))*y;    
    
    cycle = cycle + 1;

    % stop if converging or reaching the maximum iteration
    if flags < inf,      
        break; 
    end;     
end
fprintf('gmres(%d) converges at iteration %d to a solution with relative residual %g\n', [restart iter res/nrmb]);       
rev = rev/nrmb;
