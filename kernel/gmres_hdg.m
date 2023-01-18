function [x, flags, iter, rev] = gmres_hdg(A, b, f2f, x, restart, tol, maxit)

% check the number of input arguments
if nargin < 2
  error('Not enough input arguments.');
end

% get the dimension from b
N = size(b(:),1);

% Default parameters
if nargin < 4 || isempty(x),       x = zeros(N,1);  end;
if nargin < 5 || isempty(restart), restart = 20;    end;
if nargin < 6 || isempty(tol),     tol = 1e-6;      end;
if nargin < 7 || isempty(maxit),   maxit = N;       end;

% initialization
nrmb   = norm(b); 
flags  = inf; 
iter   = 0; 
cycle  = 0;

sz = size(A);
nb = sz(1);
nf = size(f2f,2);
ind = (f2f(:) == 0);
f2f(ind) = nf+1;
f2f(1,:) = [];
nn = size(f2f,1);

% allocate memory
hh = zeros(restart+1,restart);
v  = zeros(N,restart+1);
e1 = zeros(restart+1,1); e1(1)=1;
rev = zeros(restart,1);
d = zeros(nb,nf);

% A: [(npf*nch)x(npf*nch)*(2*nfe-1) nf]
% f: [(npf*nch)*(2*nfe-1) nf]

% b, x, A, v
% temporary f 

x = reshape(x,[nb nf]);
while (1) 
    % perform matrix-vector multiplication              
    s = cat(2,reshape(x,[nb nf]),zeros(nb,1));
    f = reshape(s(:,f2f(:)),[nb*nn nf]);
    d = arraymul(A,f);    
    d = x + d; % d = M*A*x
    
    % compute the residual vector      
    r = b(:) - d(:); % r = M*b - M*Ax = M*(b-A*x)       
    beta = norm(r);
    v(:,1) = r/beta;

    % v(:,1) = b - A*x;
    
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
        t = reshape(v(:,j),[nb nf]);
        s = cat(2,t,zeros(nb,1));        
        f = reshape(s(:,f2f(:)),[nb*nn nf]);
        d = arraymul(A,f);    
        d = t + d;
        
        % Arnoldi process (i.e., GS orthogonalization on the Krylov supspace)        
        v(:,j+1) = d(:);
%         for i = 1:j        
%             hh(i,j) = v(:,i)'*v(:,j+1);
%             v(:,j+1) = v(:,j+1) - hh(i,j)*v(:,i);        
%         end
        hh(1:j,j) = v(:,1:j)'*v(:,j+1);
        v(:,j+1) = v(:,j+1) - v(:,1:j)*hh(1:j,j);        
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

% function [x, flags, iter, rev] = gmres(A, b, f2f, x, restart, tol, maxit)
% 
% % check the number of input arguments
% if nargin < 2
%   error('Not enough input arguments.');
% end
% 
% % get the dimension from b
% N = size(b(:),1);
% 
% % Default parameters
% if nargin < 4 || isempty(x),       x = zeros(N,1);  end;
% if nargin < 5 || isempty(restart), restart = 20;    end;
% if nargin < 6 || isempty(tol),     tol = 1e-4;      end;
% if nargin < 7 || isempty(maxit),   maxit = N;       end;
% 
% [x, flags, iter, rev] = gmres1(A, b, f2f, x, restart, tol, maxit);
% %[x, flags, iter, rev] = gmres2(A, b, f2f, x, restart, tol, maxit);
% 
% function [x, flags, iter, rev] = gmres1(A, b, f2f, x, restart, tol, maxit)
% % 
% % % check the number of input arguments
% % if nargin < 2
% %   error('Not enough input arguments.');
% % end
% % 
% % % get the dimension from b
% % N = size(b(:),1);
% % 
% % % Default parameters
% % if nargin < 4 || isempty(x),       x = zeros(N,1);  end;
% % if nargin < 5 || isempty(restart), restart = 20;    end;
% % if nargin < 6 || isempty(tol),     tol = 1e-6;      end;
% % if nargin < 7 || isempty(maxit),   maxit = N;       end;
% 
% N = size(b(:),1);
% 
% % initialization
% nrmb   = norm(b); 
% flags  = inf; 
% iter   = 0; 
% cycle  = 0;
% 
% sz = size(A);
% nb = sz(1);
% nf = size(f2f,2);
% ind = (f2f(:) == 0);
% f2f(ind) = nf+1;
% f2f(1,:) = [];
% nn = size(f2f,1);
% 
% % allocate memory
% hh = zeros(restart+1,restart);
% cc = hh;
% v  = zeros(N,restart+1);
% e1 = zeros(restart+1,1); e1(1)=1;
% rev = zeros(restart,1);
% d = zeros(nb,nf);
% 
% np = max(1,matlabpool('size'));
% nm = floor(nf/np);
% ns = nm*ones(np,1);
% for j=1:nf-nm*np
%     ns(j) = ns(j)+1;
% end
% fpart = zeros(nf,1);
% for i=2:np
%     ind = sum(ns(1:i-1))+1:sum(ns(1:i));
%     fpart(ind) = i-1;
% end
% for i=1:np
%     inp{i} = find(fpart==i-1);    
% end
% for i=np:-1:1
%     An{i} = A(:,:,inp{i});
%     %A(:,:,inp{i}) = [];
% end
% clear A;
% 
% x = reshape(x,[nb nf]);
% while (1) 
%     % perform matrix-vector multiplication              
%     %s = cat(2,reshape(x,[nb nf]),zeros(nb,1));
%     s = [x zeros(nb,1)];
%     f = reshape(s(:,f2f(:)),[nb*nn nf]);    
% %     parfor i=1:nf
% %         d(:,i) = A(:,:,i)*f(:,i);        
% %     end    
% %     d = d+x;
%     for i = 1:np
%         fn{i} = f(:,inp{i});
%     end
%     parfor i=1:np
%         %dn{i} =  arraymul(An{i},fn{i});
%         for k=1:size(fn{i},2)
%             dn{i}(:,k) = An{i}(:,:,k)*fn{i}(:,k);
%         end
%     end        
%     for i = 1:np
%         d(:,inp{i}) = x(:,inp{i}) + dn{i};
%     end    
%     
%     % compute the residual vector      
%     r = b(:) - d(:);        
%     beta = norm(r);
%     v(:,1) = r/beta;
% 
%     res  = beta;
%     iter = iter+1;
%     rev(iter) = res;   
%     for j = 1:restart
% 
%         % set flag=0 if convergence
%         if res/nrmb <= tol
%             flags = 0;
%             break;
%         end        
%         
%         % set flag=1 if reaching maximum iteration 
%         if iter>=maxit
%             flags = 1;
%             break;
%         end                
%                  
%         % perform matrix-vector multiplication
%         t = reshape(v(:,j),[nb nf]);
%         %s = cat(2,t,zeros(nb,1));        
%         s = [t zeros(nb,1)];
%         f = reshape(s(:,f2f(:)),[nb*nn nf]);
% %         parfor i=1:nf
% %             d(:,i) = A(:,:,i)*f(:,i);                
% %         end    
% %         d = t + d;
%         for i = 1:np
%             fn{i} = f(:,inp{i});
%         end
%         parfor i=1:np
%             %dn{i} =  arraymul(An{i},fn{i});
%             for k=1:size(fn{i},2)
%                 dn{i}(:,k) = An{i}(:,:,k)*fn{i}(:,k);
%             end
%         end    
%         for i = 1:np
%             d(:,inp{i}) = t(:,inp{i}) + dn{i};
%         end     
%         
%         % Arnoldi process (i.e., GS orthogonalization on the Krylov supspace)        
%         v(:,j+1) = d(:);
% %         for i = 1:j        
% %             hh(i,j) = v(:,i)'*v(:,j+1);
% %             v(:,j+1) = v(:,j+1) - hh(i,j)*v(:,i);        
% %         end
%         hh(1:j,j) = v(:,1:j)'*v(:,j+1);
%         v(:,j+1) = v(:,j+1) - v(:,1:j)*hh(1:j,j);
% %         cc(1:j,j) = v(:,1:j)'*v(:,j+1);
% %         v(:,j+1) = v(:,j+1) - v(:,1:j)*cc(1:j,j);
% %         hh(1:j,j) = hh(1:j,j) + cc(1:j,j);
%         hh(j+1,j) = norm(v(:,j+1));             
%         if (hh(j+1,j) ~= 0.0)          
%             v(:,j+1) = v(:,j+1)/hh(j+1,j);   
%         else            
%             break;
%         end
% 
%         % solve the reduced system 
%         y = hh(1:j+1,1:j)\(beta*e1(1:j+1));
%         
%         % compute the residual norm
%         res = norm(beta*e1(1:j+1)-hh(1:j+1,1:j)*y);
%         
%         iter = iter + 1; 
%         rev(iter) = res;
%     end 
%       
%     % compute the solution    
%     x(:) = x(:) + v(:,1:length(y))*y;    
%     
%     cycle = cycle + 1;
% 
%     % stop if converging or reaching the maximum iteration
%     if flags < inf, 
%         fprintf('gmres(%d) converges at outer iteration %d to a solution with relative residual %g\n', [restart cycle res/nrmb]);       
%         break; 
%     end;     
% end
% rev = rev/nrmb;
% 
% 
% function [x, flags, iter, rev] = gmres2(A, b, f2f, x, restart, tol, maxit)
% 
% % initialization
% nrmb   = norm(b); 
% flags  = inf; 
% iter   = 0; 
% cycle  = 0;
% 
% sz = size(A);
% nb = sz(1);
% nf = size(f2f,2);
% ind = (f2f(:) == 0);
% f2f(ind) = nf+1;
% f2f(1,:) = [];
% nn = size(f2f,1);
% 
% % allocate memory
% hh = zeros(restart+1,restart);
% %v  = zeros(N,restart+1);
% e1 = zeros(restart+1,1); e1(1)=1;
% rev = zeros(restart,1);
% d = zeros(nb,nf);
% 
% np = max(1,matlabpool('size'));
% nm = floor(nf/np);
% ns = nm*ones(np,1);
% for j=1:nf-nm*np
%     ns(j) = ns(j)+1;
% end
% fpart = zeros(nf,1);
% for i=2:np
%     ind = sum(ns(1:i-1))+1:sum(ns(1:i));
%     fpart(ind) = i-1;
% end
% for i=1:np
%     inp{i} = find(fpart==i-1);    
%     v{i} = zeros(length(inp{i})*nb,restart+1);
% end
% for i=np:-1:1
%     An{i} = A(:,:,inp{i});
%     A(:,:,inp{i}) = [];
% end
% 
% x = reshape(x,[nb nf]);
% b = reshape(b,[nb nf]);
% for i = 1:np
%     bn{i} = b(:,inp{i});    
% end
% 
% while (1) 
%     % perform matrix-vector multiplication              
%     s = [x zeros(nb,1)];
%     f = reshape(s(:,f2f(:)),[nb*nn nf]);    
%     for i = 1:np
%         fn{i} = f(:,inp{i});
%         xn{i} = x(:,inp{i});
%     end
%     beta = 0;
%     parfor i=1:np
%         rn{i} = bn{i} - xn{i} - arraymul(An{i},fn{i});        
%         beta = beta + sum(rn{i}(:).^2);
%     end        
%     beta = sqrt(beta);
%     for i=1:np
%         v{i}(:,1) = rn{i}(:)/beta;
%     end
%     
%     res  = beta;
%     iter = iter+1;
%     rev(iter) = res;    
%     for j = 1:restart
% 
%         % set flag=0 if convergence
%         if res/nrmb <= tol
%             flags = 0;
%             break;
%         end        
%         
%         % set flag=1 if reaching maximum iteration 
%         if iter>=maxit
%             flags = 1;
%             break;
%         end                
%                  
%         % perform matrix-vector multiplication
%         s = zeros(nb, nf+1);        
%         for i=1:np
%             s(:,inp{i}) = reshape(v{i}(:,j),nb,[]);
%         end
%         f = reshape(s(:,f2f(:)),[nb*nn nf]);
%         for i = 1:np
%             fn{i} = f(:,inp{i});
%             tn{i} = s(:,inp{i});
%         end
%         nrm = 0;
%         tm = zeros(j,1);
%         parfor i=1:np
%             v{i}(:,j+1) =  reshape(tn{i} + arraymul(An{i},fn{i}),[],1);
%             tm = tm + v{i}(:,1:j)'*v{i}(:,j+1);            
%         end            
%         parfor i=1:np
%             v{i}(:,j+1) = v{i}(:,j+1) - v{i}(:,1:j)*tm;
%             nrm = nrm + sum(v{i}(:,j+1).^2);
%         end
%         hh(1:j,j) = tm;
%         hh(j+1,j) = sqrt(nrm);
%         
%         if (hh(j+1,j) ~= 0.0)          
%             for i=1:np
%                 v{i}(:,j+1) = v{i}(:,j+1)/hh(j+1,j);   
%             end
%         else            
%             break;
%         end
%         
%         % solve the reduced system 
%         y = hh(1:j+1,1:j)\(beta*e1(1:j+1));
%         
%         % compute the residual norm
%         res = norm(beta*e1(1:j+1)-hh(1:j+1,1:j)*y);
%         
%         iter = iter + 1; 
%         rev(iter) = res;
%     end 
%       
%     % compute the solution    
%     parfor i=1:np
%         xn{i}(:) = xn{i}(:) + v{i}(:,1:length(y))*y;    
%     end
%     for i=1:np
%         x(:,inp{i}) = xn{i};
%     end
%     
%     cycle = cycle + 1;
% 
%     % stop if converging or reaching the maximum iteration
%     if flags < inf, 
%         fprintf('gmres(%d) converges at outer iteration %d to a solution with relative residual %g\n', [restart cycle res/nrmb]);       
%         break; 
%     end;     
% end
% rev = rev/nrmb;
