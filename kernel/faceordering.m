function P = faceordering(A,f2f)

f2f(:,1) = [];

np = size(A,1);
nb = size(A,3);
nf = size(A,4);
C  = zeros(nf,nb);

for i = 1:nf
    for j = 1:nb
        C(i,j) = norm(reshape(A(:,:,j,i),[np*np,1]));
    end
end

w = zeros(nf,1);
for k = 1:nf % for each face k
    D = zeros(nb,nb);
    for i = 1:nb
        for j = 1:nb
            if i ~= j         
                ia = f2f(k,i); % face ia  
                ib = f2f(k,j); % face ib                  
                if ia>0 && ib>0
                    ja = f2f(ia,:)==k;                     
                    D(i,j) = C(ia,ja)*C(k,j);
                end                
            end
        end
    end
    w(k) = norm(D(:));
end
P=0;

% P = zeros(nf,1);
% for i = 1:nf
%     [~,P(i)] = min(w);
%     w(P(i)) = inf;    
%     for k = 1:nb % update the weight of neighbors of P(i)
%         fb = f2f(P(i),k); % face fb is a neighbor of P(i)
%         if fb==0
%             break;
%         end
%         if ismember(fb,P(1:i))==0 % face fb is not pivoted yet
%             D = zeros(nb,nb);
%             for m = 1:nb
%                 for n = 1:nb                   
%                     ia = f2f(fb,m);                         
%                     ib = f2f(fb,n);
%                     if (m ~= n) && ia>0 && ib>0 && ismember(ia,P(1:i))==0 && ismember(ib,P(1:i))==0
%                         ja = f2f(ia,:)==fb;
%                         D(i,j) = C(ia,ja)*C(fb,n);
%                     end                    
%                 end
%             end
%         end
%     end
% end
% 
