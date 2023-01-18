function plotpartition(p,t,fedg,elem2cpu,ent2cpu,extintelem,extintent,extintelempts,hybrid)

bcol = [1 1 0;... % yellow 
        1 0 1;... % magneta
        0 1 1;... % cyan
        1 0 0;... % red
        0 1 0;... % green
        0 0 1;... % blue
        1,0.4,0.6;...
        0.4,0.6,1;...
       ];
   
np = max(elem2cpu)+1;

% plot nonoverlapping subdomains
figure(1); clf;
hold on;        
for i=0:np-1
    ind = elem2cpu==i;
    ti = t(:,ind);
    simpplot(p,ti',[],bcol(i+1,:));                       
end
% for it=1:size(t,2)
%     pmid=mean(p(t(:,it),:),1);
%     txtpars={'fontname','times','fontsize',20,'horizontala','center','BackgroundColor',[1,1,1]};
%     text(pmid(1),pmid(2),num2str(it),txtpars{:});
% end
if strcmp(hybrid,'hdg')
    for it=1:size(fedg,1)
        pmid=mean(p(fedg(it,1:2),:),1);
        i = ent2cpu(it)+1;
        txtpars={'fontname','times','fontsize',20,'horizontala','center','BackgroundColor',bcol(7-i,:)};
        text(pmid(1),pmid(2),num2str(it),txtpars{:});
    end
elseif strcmp(hybrid,'edg')        
    for it=1:size(fedg,1)
        pmid=fedg(it,:);
        i = ent2cpu(it)+1;
        txtpars={'fontname','times','fontsize',20,'horizontala','center','BackgroundColor',bcol(7-i,:)};
        text(pmid(1),pmid(2),num2str(it),txtpars{:});
    end
end

hold off;
axis equal;      
axis tight;
axis on;  


% plot overlapping subdomains
for n=1:np               
    figure(n+1); clf;    
    hold on;        
    simpplot(p,t',[],'none');   
    
    tn = t(:,extintelem{n}(1:sum(extintelempts(n,1:2)))+1);
    simpplot(p,tn',[],bcol(n,:));                           
        
    tn = t(:,extintelem{n}(sum(extintelempts(n,1:2))+1:end)+1);
    simpplot(p,tn',[],[0 0 0]);                           
    
    if strcmp(hybrid,'hdg')
        for it=1:size(fedg,1)
            pmid=mean(p(fedg(it,1:2),:),1);
            if ismember(it,extintent{n}+1)
                txtpars={'fontname','times','fontsize',20,'horizontala','center','BackgroundColor',bcol(7-n,:)};
            else
                txtpars={'fontname','times','fontsize',20,'horizontala','center','BackgroundColor',[1 1 1]};   
            end
            text(pmid(1),pmid(2),num2str(it),txtpars{:});
        end
    elseif strcmp(hybrid,'edg')                
        for it=1:size(fedg,1)
            pmid=fedg(it,:);
            if ismember(it,extintent{n}+1)
                txtpars={'fontname','times','fontsize',20,'horizontala','center','BackgroundColor',bcol(7-n,:)};                
            else
                txtpars={'fontname','times','fontsize',20,'horizontala','center','BackgroundColor',[1 1 1]};                
            end
            text(pmid(1),pmid(2),num2str(it),txtpars{:});
        end
    end

    hold off;
    axis equal;      
    axis tight;
    axis on;  
    
%     intent = find(ent2cpu==n-1);
%     extintent{n}'+1
%     intent'
%     pause    
end

