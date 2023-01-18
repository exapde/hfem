function [f,t2f] = mkt2f(t)
%MKT2F Compute Face Connectivity and Triangle to Face Connetivity.
%   [F,T2F]=MKT2F(T)
%
%      T:         Triangle indices (NT,3)
%      F:         Face connectivity (NF,4) (for boundary edges F(:,4)=0)
%      T2F:       Triangle to Face Connectivity (NT,3)
%
%   See also: MKT2T.

[nt,dim1] = size(t);
dim = dim1-1;

switch dim1-1
    case 2
        edg = [2,3; 3,1; 1,2];
    case 3
        edg = [2,3,4; 1,4,3; 1,2,4; 1,3,2];
    otherwise
        error('Only can handle dim=2 or dim=3');
end

t2t = mkt2t(t);
nb = sum(sum(t2t <= 0));
f = zeros((dim1*nt+nb)/2,dim1+1);
t2f = zeros(size(t));
jf = 0;
for i=1:nt
    for j=1:dim1
        if t2t(i,j) > i || t2t(i,j) <=0
            ie = t2t(i,j);
            jf = jf + 1;
            
            f(jf,1:dim) = t(i,edg(j,:));
            f(jf,dim1) = i;
            f(jf,dim1+1) = ie;
            t2f(i,j) = jf;
        
            if ie > 0
                t2f(ie,find(t(ie,:)==(sum(t(ie,:))-sum(f(jf,1:dim))))) = jf;
            end
  
        end
    end
end
clear t2t;

% Reorder faces - First interior then boundary

nf = size(f,1);
[a,mp] = sort(f(:,dim1+1) == 0);
f = f(mp,:);

a = zeros(nf,1);
a(mp) = (1:nf)';
t2f = reshape(t2f,dim1*nt,1);
ii = find(t2f > 0);
t2f(ii) = a(t2f(ii));
ii = find(t2f < 0);
t2f(ii) = -a(-t2f(ii));
t2f = reshape(t2f,nt,dim1);





