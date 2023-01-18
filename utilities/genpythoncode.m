
% generate matlab code for the flux function
if exist('f','var') && exist('udg','var')
    
    f_udg  = jacobian(f,udg);
    f = reshape(f,[ncu,nd]);
    f_udg = reshape(f_udg,[ncu,nd,nc]);


    % generate temporary C file
    ccode(f(1:end,:),'file','tmp1.c')
    ccode(f_udg(1:end,:,:),'file','tmp2.c')

    delete(filename1);
    gid = fopen(filename1,'wt');
    
    str = 'from numpy import *';
    fprintf(gid, '%s\n\n', str);                  
    str = 'def flux(pg,udg,param,time):';
    fprintf(gid, '%s\n\n', str);         
    
    str = 'ng = pg.shape[0]';
    fprintf(gid, '\t%s\n', str);         
    str = 'nd = pg.shape[1]';
    fprintf(gid, '\t%s\n', str);    
    str = 'nc = udg.shape[1]';
    fprintf(gid, '\t%s\n', str);         
    str = ['ncu = ' num2str(ncu)];
    fprintf(gid, '\t%s\n', str);             
    str = 'f = zeros((ng,ncu,nd))';
    fprintf(gid, '\t%s\n', str);             
    str = 'f_udg = zeros((ng,ncu,nd,nc))';
    fprintf(gid, '\t%s\n', str);             
    fprintf(gid, '\n');      

    for i=1:length(pg)
        str = ['x' num2str(i) ' = pg[:,' num2str(i-1) ']'];
        fprintf(gid, '\t%s\n', str);                  
    end
    for i=1:length(udg)
        str = ['u' num2str(i) ' = udg[:,' num2str(i-1) ']'];
        fprintf(gid, '\t%s\n', str);                  
    end
    for i=1:length(param)
        str = ['param' num2str(i) ' = param[' num2str(i-1) ']'];
        fprintf(gid, '\t%s\n', str);                  
    end    
    fprintf(gid, '\n');      
    
    fid = fopen('tmp1.c','r');    
    tline = fgetl(fid); 
    i=1;       
    while ischar(tline)        
        %str = tline;
        str = strrep(tline, ';', '');
        str = strrep(str, 't0 = ', '');
        str = strrep(str, 'A0[', 'f[:,');
        str = strrep(str, '][', ',');           
        fprintf(gid, '\t%s\n', str);                  
        tline = fgetl(fid);        
        i=i+1;        
    end
    fclose(fid);
    fprintf(gid, '\n');  
        
    fid = fopen('tmp2.c','r');    
    tline = fgetl(fid); 
    i=1;       
    while ischar(tline)        
        %str = tline;
        str = strrep(tline, ';', '');
        str = strrep(str, 't0 = ', '');
        str = strrep(str, 'A0[', 'f_udg[:,');
        str = strrep(str, '][', ',');             
        fprintf(gid, '\t%s\n', str);                  
        tline = fgetl(fid);        
        i=i+1;        
    end
    fclose(fid);
    fprintf(gid, '\n');  
        
    str = 'return f, f_udg';
    fprintf(gid, '\t%s', str);                  
    fprintf(gid, '\n'); 
    
    fclose(gid);
    
%     delete('tmp1.c');
%     delete('tmp2.c');
end


% generate matlab code for the source function
if exist('s','var') && exist('udg','var')
    s_udg  = jacobian(s,udg);
    s = reshape(s,[ncu,1]);
    s_udg = reshape(s_udg,[ncu,nc]);
    
        % generate temporary C file
    ccode(s,'file','tmp1.c')
    ccode(s_udg(1:end,:),'file','tmp2.c')

    delete(filename2);
    gid = fopen(filename2,'wt');
    
    str = 'from numpy import *';
    fprintf(gid, '%s\n\n', str);                  
    str = 'def source(pg,udg,param,time):';
    fprintf(gid, '%s\n\n', str);         
    
    str = 'ng = pg.shape[0]';
    fprintf(gid, '\t%s\n', str);         
    str = 'nd = pg.shape[1]';
    fprintf(gid, '\t%s\n', str);    
    str = 'nc = udg.shape[1]';
    fprintf(gid, '\t%s\n', str);         
    str = ['ncu = ' num2str(ncu)];
    fprintf(gid, '\t%s\n', str);             
    str = 's = zeros((ng,ncu))';
    fprintf(gid, '\t%s\n', str);             
    str = 's_udg = zeros((ng,ncu,nc))';
    fprintf(gid, '\t%s\n', str);             
    fprintf(gid, '\n');      

    for i=1:length(pg)
        str = ['x' num2str(i) ' = pg[:,' num2str(i-1) ']'];
        fprintf(gid, '\t%s\n', str);                  
    end
    for i=1:length(udg)
        str = ['u' num2str(i) ' = udg[:,' num2str(i-1) ']'];
        fprintf(gid, '\t%s\n', str);                  
    end
    for i=1:length(param)
        str = ['param' num2str(i) ' = param[' num2str(i-1) ']'];
        fprintf(gid, '\t%s\n', str);                  
    end    
    fprintf(gid, '\n');      
    
    fid = fopen('tmp1.c','r');    
    tline = fgetl(fid); 
    i=1;       
    while ischar(tline)        
        %str = tline;
        str = strrep(tline, ';', '');
        str = strrep(str, 't0 = ', '');
        str = strrep(str, 'A0[', 's[:,');        
        fprintf(gid, '\t%s\n', str);                  
        tline = fgetl(fid);        
        i=i+1;        
    end
    fclose(fid);
    fprintf(gid, '\n');  
        
    fid = fopen('tmp2.c','r');    
    tline = fgetl(fid); 
    i=1;       
    while ischar(tline)        
        %str = tline;
        str = strrep(tline, ';', '');
        str = strrep(str, 't0 = ', '');
        str = strrep(str, 'A0[', 's_udg[:,');
        str = strrep(str, '][', ',');             
        fprintf(gid, '\t%s\n', str);                  
        tline = fgetl(fid);        
        i=i+1;        
    end
    fclose(fid);
    fprintf(gid, '\n');  
        
    str = 'return s, s_udg';
    fprintf(gid, '\t%s', str);                  
    fprintf(gid, '\n'); 
    
    fclose(gid);
    
    delete('tmp1.c');
    delete('tmp2.c');
elseif exist('filename2','var')
    gid = fopen(filename2,'wt');
    
    str = 'from numpy import *';
    fprintf(gid, '%s\n\n', str);                  
    str = 'def source(pg,udg,param,time):';
    fprintf(gid, '%s\n\n', str);         
    
    str = 'ng = pg.shape[0]';
    fprintf(gid, '\t%s\n', str);         
    str = 'nd = pg.shape[1]';
    fprintf(gid, '\t%s\n', str);    
    str = 'nc = udg.shape[1]';
    fprintf(gid, '\t%s\n', str);         
    str = ['ncu = ' num2str(ncu)];
    fprintf(gid, '\t%s\n', str);             
    str = 's = zeros((ng,ncu))';
    fprintf(gid, '\t%s\n', str);             
    str = 's_udg = zeros((ng,ncu,nc))';
    fprintf(gid, '\t%s\n', str);             
    fprintf(gid, '\n');      

    str = 'return s, s_udg';
    fprintf(gid, '\t%s', str);                  
    fprintf(gid, '\n'); 
        
    fclose(gid);    
end

% generate matlab code for the fhat function
if exist('fh','var') && exist('udg','var') && exist('uh','var')
    
    fhat_udg   = jacobian(fhat,udg);
    fhat_uh = jacobian(fhat,uh);    
    fhat = reshape(fhat,[ncu,1]);
    fhat_udg = reshape(fhat_udg,[ncu,nc]);
    fhat_uh = reshape(fhat_uh,[ncu,ncu]);

    % generate temporary C file
    ccode(fhat(1:end),'file','tmp1.c')
    ccode(fhat_udg(1:end,:),'file','tmp2.c')
    ccode(fhat_uh(1:end,:),'file','tmp3.c')
    
    delete(filename3);
    gid = fopen(filename3,'wt');
    
    str = 'from numpy import *';
    fprintf(gid, '%s\n\n', str);                  
    str = 'def fhat(nl,pg,udg,uh,param,time):';
    fprintf(gid, '%s\n\n', str);         
    
    str = 'ng = pg.shape[0]';
    fprintf(gid, '\t%s\n', str);         
    str = 'nd = pg.shape[1]';
    fprintf(gid, '\t%s\n', str);    
    str = 'nc = udg.shape[1]';
    fprintf(gid, '\t%s\n', str);         
    str = ['ncu = ' num2str(ncu)];
    fprintf(gid, '\t%s\n', str);             
    str = 'fh = zeros((ng,ncu))';
    fprintf(gid, '\t%s\n', str);             
    str = 'fh_udg = zeros((ng,ncu,nc))';
    fprintf(gid, '\t%s\n', str);             
    str = 'fh_uh = zeros((ng,ncu,ncu))';
    fprintf(gid, '\t%s\n', str);             
    fprintf(gid, '\n');      

    for i=1:length(pg)
        str = ['x' num2str(i) ' = pg[:,' num2str(i-1) ']'];
        fprintf(gid, '\t%s\n', str);                  
    end
    for i=1:length(nl)
        str = ['nl' num2str(i) ' = nl[:,' num2str(i-1) ']'];
        fprintf(gid, '\t%s\n', str);                  
    end
    for i=1:length(udg)
        str = ['u' num2str(i) ' = udg[:,' num2str(i-1) ']'];
        fprintf(gid, '\t%s\n', str);                  
    end    
    for i=1:length(uh)
        str = ['uh' num2str(i) ' = uh[:,' num2str(i-1) ']'];
        fprintf(gid, '\t%s\n', str);                  
    end
    for i=1:length(param)
        str = ['param' num2str(i) ' = param[' num2str(i-1) ']'];
        fprintf(gid, '\t%s\n', str);                  
    end    
    fprintf(gid, '\n');      
    
    fid = fopen('tmp1.c','r');    
    tline = fgetl(fid); 
    i=1;       
    while ischar(tline)        
        %str = tline;
        str = strrep(tline, ';', '');
        str = strrep(str, 't0 = ', '');
        str = strrep(str, 'A0[', 'fh[:,');
        str = strrep(str, '][', ',');           
        fprintf(gid, '\t%s\n', str);                  
        tline = fgetl(fid);        
        i=i+1;        
    end
    fclose(fid);
    fprintf(gid, '\n');  
        
    fid = fopen('tmp2.c','r');    
    tline = fgetl(fid); 
    i=1;       
    while ischar(tline)        
        %str = tline;
        str = strrep(tline, ';', '');
        str = strrep(str, 't0 = ', '');
        str = strrep(str, 'A0[', 'fh_udg[:,');
        str = strrep(str, '][', ',');             
        fprintf(gid, '\t%s\n', str);                  
        tline = fgetl(fid);        
        i=i+1;        
    end
    fclose(fid);
    fprintf(gid, '\n');  
        
    fid = fopen('tmp3.c','r');    
    tline = fgetl(fid); 
    i=1;       
    while ischar(tline)        
        %str = tline;
        str = strrep(tline, ';', '');
        str = strrep(str, 't0 = ', '');
        str = strrep(str, 'A0[', 'fh_uh[:,');
        str = strrep(str, '][', ',');             
        fprintf(gid, '\t%s\n', str);                  
        tline = fgetl(fid);        
        i=i+1;        
    end
    fclose(fid);
    fprintf(gid, '\n');  
    
    str = 'return fh, fh_udg, fh_uh';
    fprintf(gid, '\t%s', str);                  
    fprintf(gid, '\n'); 
    
    fclose(gid);       
    
%     delete('tmp1.c');
%     delete('tmp2.c');
%     delete('tmp3.c');
end

% generate matlab code for the fbou function
if exist('fb','var') && exist('udg','var') && exist('uh','var')
    
    delete(filename4);
    gid = fopen(filename4,'wt');
    
    str = 'from numpy import *';
    fprintf(gid, '%s\n\n', str);                  
    str = 'def fbou(ib,uinf,nl,pg,udg,uh,param,time):';
    fprintf(gid, '%s\n\n', str);         
    
    str = 'ng = pg.shape[0]';
    fprintf(gid, '\t%s\n', str);         
    str = 'nd = pg.shape[1]';
    fprintf(gid, '\t%s\n', str);    
    str = 'nc = udg.shape[1]';
    fprintf(gid, '\t%s\n', str);         
    str = ['ncu = ' num2str(ncu)];
    fprintf(gid, '\t%s\n', str);             
    str = 'fh = zeros((ng,ncu))';
    fprintf(gid, '\t%s\n', str);             
    str = 'fh_udg = zeros((ng,ncu,nc))';
    fprintf(gid, '\t%s\n', str);             
    str = 'fh_uh = zeros((ng,ncu,ncu))';
    fprintf(gid, '\t%s\n', str);             
    fprintf(gid, '\n');      

    for i=1:length(pg)
        str = ['x' num2str(i) ' = pg[:,' num2str(i-1) ']'];
        fprintf(gid, '\t%s\n', str);                  
    end
    for i=1:length(nl)
        str = ['nl' num2str(i) ' = nl[:,' num2str(i-1) ']'];
        fprintf(gid, '\t%s\n', str);                  
    end
    for i=1:length(uinf)
        str = ['uinf' num2str(i) ' = uinf[:,' num2str(i-1) ']'];
        fprintf(gid, '\t%s\n', str);                  
    end        
    for i=1:length(udg)
        str = ['u' num2str(i) ' = udg[:,' num2str(i-1) ']'];
        fprintf(gid, '\t%s\n', str);                  
    end    
    for i=1:length(uh)
        str = ['uh' num2str(i) ' = uh[:,' num2str(i-1) ']'];
        fprintf(gid, '\t%s\n', str);                  
    end
    for i=1:length(param)
        str = ['param' num2str(i) ' = param[' num2str(i-1) ']'];
        fprintf(gid, '\t%s\n', str);                  
    end    
    fprintf(gid, '\n');              
    
    for k=1:length(fb)
        k
        fh = fb{k};
        fh_udg = jacobian(fh,udg);
        fh_uh = jacobian(fh,uh);
       
        ccode(fh(1:end),'file','tmp1.c')
        ccode(fh_udg(1:end,:),'file','tmp2.c')
        ccode(fh_uh(1:end,:),'file','tmp3.c')  
        
        if k==1
            str = ['if ib ==' num2str(k) ':'];
        else
            str = ['elif ib ==' num2str(k) ':'];
        end                
        fprintf(gid, '\t%s', str);                  
        fprintf(gid, '\n');  
        
        fid = fopen('tmp1.c','r');    
        tline = fgetl(fid); 
        i=1;       
        while ischar(tline)        
            %str = tline;
            str = strrep(tline, ';', '');
            str = strrep(str, 't0 = ', '');
            str = strrep(str, 'A0[', 'fh[:,');
            str = strrep(str, '][', ',');           
            fprintf(gid, '\t\t%s\n', str);                  
            tline = fgetl(fid);        
            i=i+1;        
        end
        fclose(fid);
        fprintf(gid, '\n');  

        fid = fopen('tmp2.c','r');    
        tline = fgetl(fid); 
        i=1;       
        while ischar(tline)        
            %str = tline;
            str = strrep(tline, ';', '');
            str = strrep(str, 't0 = ', '');
            str = strrep(str, 'A0[', 'fh_udg[:,');
            str = strrep(str, '][', ',');             
            fprintf(gid, '\t\t%s\n', str);                  
            tline = fgetl(fid);        
            i=i+1;        
        end
        fclose(fid);
        fprintf(gid, '\n');  

        fid = fopen('tmp3.c','r');    
        tline = fgetl(fid); 
        i=1;       
        while ischar(tline)        
            %str = tline;
            str = strrep(tline, ';', '');
            str = strrep(str, 't0 = ', '');
            str = strrep(str, 'A0[', 'fh_uh[:,');
            str = strrep(str, '][', ',');             
            fprintf(gid, '\t\t%s\n', str);                  
            tline = fgetl(fid);        
            i=i+1;        
        end
        fclose(fid);
        fprintf(gid, '\n');          
    end
    
    str = 'return fh, fh_udg, fh_uh';
    fprintf(gid, '\t%s', str);                  
    fprintf(gid, '\n'); 
    
    fclose(gid);           
end
