function [K,F] = mkDenseBlockSystem(AE, FE, fe1, fe2, elcon)

if nargin == 5
    [K,F] = mkBlockSystem1(AE, FE, fe1, fe2, elcon);    
else
    [K,F] = mkBlockSystem2(AE, FE, fe1, fe2);    
end

function [K,F] = mkBlockSystem1(AE, FE, f, t2f, elcon)
% mkBlockMatrix returns 4-dimensional array for the Jacobian matrix 
%  
%   AE    :  5-dimensional array for the elemental matrices
%   FE    :  3-dimensional array for the elemental vectors
%   f     :  Face to element connectivity
%   t2f   :  Element to face connectivity
%
%   K     :  4-dimensional array for the Jacobian matrix
%   F     :  2-dimensional array for the RHS vector
%   f2f   :  Face to face connectivity

nch = size(AE,1);   % number of components of UH
ndf = size(AE,2);   % number of points per face times number of faces per element
ne  = size(AE,5);   % number of elements
nf  = size(f,1);    % number of faces
nfe = size(t2f,2);  % number of faces per element
npf = ndf/nfe;      % number of points per face
ncf = nch*npf;      % number of components of UH times number of points per face
nbf = 2*nfe-1;      % number of neighboring faces

FE  = reshape(FE,[nch npf nfe ne]);
AE  = reshape(AE,[nch npf nfe nch npf nfe ne]);
elcon = reshape(elcon,[npf nfe ne]);


K   = zeros(ncf,ncf,nbf-1,nf);
F   = zeros(ncf,nf);
for i = 1:nf
    fi = f(i,end-1:end); % obtain two elements sharing the same face i      
    if fi(2)>0           % face i is an interior face
        kf = t2f(fi,:);         % obtain neighboring faces 
        i1 = find(kf(1,:)==i);  % obtain the index of face i in the 1st element
        i2 = find(kf(2,:)==i);  % obtain the index of face i in the 2nd element            
                        
        j1 = elcon(:,i1,fi(1)) - (i-1)*npf;        
        j2 = elcon(:,i2,fi(2)) - (i-1)*npf;        
        
        % the first block        
        k = 1;
        M = reshape(AE(:,j1,i1,:,j1,i1,fi(1)) + AE(:,j2,i2,:,j2,i2,fi(2)), [ncf ncf]);        
        F(:,i)   = reshape(FE(:,j1,i1,fi(1)) + FE(:,j2,i2,fi(2)), [ncf 1]);                
        
        for is=1:nfe % loop over each faces of the 1st element
            if is ~= i1  
                k = k + 1;                 
                j = kf(1,is);
                j3 = elcon(:,is,fi(1)) - (j-1)*npf;                
                K(:,:,k-1,i) = reshape(AE(:,j1,i1,:,j3,is,fi(1)), [ncf ncf]);                                                                    
            end
        end
        
        for is=1:nfe % loop over faces of the 2nd element
            if is ~= i2                                                
                k = k + 1;                 
                j = kf(2,is);
                j4 = elcon(:,is,fi(2)) - (j-1)*npf;                
                K(:,:,k-1,i) = reshape(AE(:,j2,i2,:,j4,is,fi(2)), [ncf ncf]);                                                    
            end
        end        
    else % face i is a boundary face
        kf = t2f(fi(1),:);      % obtain neighboring faces 
        i1 = find(kf(1,:)==i);  % obtain the index of face i in the 1st element         
        j1 = elcon(:,i1,fi(1)) - (i-1)*npf;                
        
        % the first block        
        k = 1;
        M = reshape(AE(:,j1,i1,:,j1,i1,fi(1)), [ncf ncf]);        
        F(:,i)   = reshape(FE(:,j1,i1,fi(1)), [ncf 1]);                
        
        for is=1:nfe % loop over each faces of the 1st element
            if is ~= i1  
                k = k + 1;          
                j = kf(1,is);
                j3 = elcon(:,is,fi(1)) - (j-1)*npf;
                K(:,:,k-1,i) = reshape(AE(:,j1,i1,:,j3,is,fi(1)), [ncf ncf]);                                                    
            end
        end        
    end    
    
    Mi = inv(M);
    F(:,i) = Mi*F(:,i);      
    for j=1:nbf-1
       K(:,:,j,i) = Mi*K(:,:,j,i);
    end    
end
K = reshape(K,[ncf,ncf*(nbf-1),nf]);


function [K,F] = mkBlockSystem2(AE, FE, fe1, fe2)
% mkBlockMatrix returns 4-dimensional array for the Jacobian matrix 
%  
%   AE    :  5-dimensional array for the elemental matrices
%   FE    :  3-dimensional array for the elemental vectors
%   fe1   :  Face to element connectivity
%   fe2   :  Face to element connectivity
%
%   K     :  4-dimensional array for the Jacobian matrix
%   F     :  2-dimensional array for the RHS vector

nch = size(AE,1);   % number of components of UH
ndf = size(AE,2);   % number of points per face times number of faces per element
ne  = size(AE,5);   % number of elements
nf  = size(fe1,1);  % number of faces
nfe = size(fe1,2)-1;% number of faces per element
npf = ndf/nfe;      % number of points per face
ncf = nch*npf;      % number of components of UH times number of points per face
nbf = 2*nfe-1;      % number of neighboring faces

FE  = reshape(FE,[nch npf nfe ne]);
AE  = reshape(AE,[nch npf nfe nch npf nfe ne]);

K   = zeros(ncf,ncf,nbf,nf);
F   = zeros(ncf,nf);

for i = 1:nf
    fi = [fe1(i,1) fe2(i,1)]; % obtain two elements sharing the same face i      
    if fi(2)>0           % face i is an interior face
        i1 = fe1(i,2);
        i2 = fe2(i,2);                
        
        % the first block        
        k = 1;
        K(:,:,k,i) = reshape(AE(:,:,i1,:,:,i1,fi(1)) + AE(:,:,i2,:,:,i2,fi(2)), [ncf ncf]);        
        F(:,i)   = reshape(FE(:,:,i1,fi(1)) + FE(:,:,i2,fi(2)), [ncf 1]);                               
        
        for is=2:nfe             
            k = k + 1;                                 
            K(:,:,k,i) = reshape(AE(:,:,i1,:,:,fe1(i,is+1),fi(1)), [ncf ncf]);                                                                                
        end
        
        for is=2:nfe % loop over remaining faces of the 2nd element             
            k = k + 1;                                 
            K(:,:,k,i) = reshape(AE(:,:,i2,:,:,fe2(i,is+1),fi(2)), [ncf ncf]);                                                                    
        end        
    else % face i is a boundary face
        i1 = fe1(i,2);
        
        % the first block        
        k = 1;
        K(:,:,k,i) = reshape(AE(:,:,i1,:,:,i1,fi(1)), [ncf ncf]);        
        F(:,i)   = reshape(FE(:,:,i1,fi(1)), [ncf 1]);                        
        
        for is=2:nfe % loop over each faces of the 1st element            
            k = k + 1;                          
            K(:,:,k,i) = reshape(AE(:,:,i1,:,:,fe1(i,is+1),fi(1)), [ncf ncf]);                                                                    
        end        
    end    
end


