
% generate c code for the flux function
if exist('f','var') && exist('udg','var')
        
    f = reshape(f,[ncu,nd]);    

    % generate temporary C file
    ccode(f(:),'file','tmp1.c')    

    delete(filename1);
    gid = fopen(filename1,'wt');
            
    if nd==2
        str = ['void fluxonly_' appname '2d(double *f, double *pg, double *udg, double *param, double time, int ng, int nc, int ncu, int nd, int ncd)'];
    elseif nd==3
        str = ['void fluxonly_' appname '3d(double *f, double *pg, double *udg, double *param, double time, int ng, int nc, int ncu, int nd, int ncd)'];
    end
    fprintf(gid, '%s\n', str); 
    str = '{';
    fprintf(gid, '%s\n', str);         
    
    for i=1:length(param)
        str = ['double param' num2str(i) ' = param[' num2str(i-1) '];'];
        fprintf(gid, '\t%s\n', str);                  
    end        

%     str = 'double ';
%     for i=1:length(pg)-1
%         str = [str 'x' num2str(i) ','];        
%     end
%     str = [str 'x' num2str(length(pg)) ';'];        
%     fprintf(gid, '\t%s\n', str);         
%     
%     str = 'double ';
%     for i=1:length(udg)-1
%         str = [str 'u' num2str(i) ','];        
%     end
%     str = [str 'u' num2str(length(udg)) ';'];        
%     fprintf(gid, '\t%s\n', str);         
        
    str = 'for (int i = 0; i <ng; i++) {';
    fprintf(gid, '\n\t%s\n', str);         
        
    for i=1:length(pg)
        str = ['double x' num2str(i) ' = pg[' num2str(i-1) '*ng+i];'];
        fprintf(gid, '\t\t%s\n', str);                  
    end
    for i=1:length(udg)
        str = ['double u' num2str(i) ' = udg[' num2str(i-1) '*ng+i];'];
        fprintf(gid, '\t\t%s\n', str);                  
    end
    fprintf(gid, '\n');      
       
    fid = fopen('tmp1.c','r');    
    tline = fgetl(fid); 
    i=1; a1 = 0;       
    while ischar(tline)        
        str = tline;
        
        i1 = strfind(str,'[');        
        i2 = strfind(str,']');        
        if isempty(i1)==0    
            a2 = str2num(str((i1+1):(i2-1)));                        
            for j = a1:(a2-1)                
                strj = ['f[' num2str(j) '*ng+i] = 0.0;'];
                fprintf(gid, '\t\t%s\n', strj);                  
            end
            a1 = a2+1;              
        end
                
        str = strrep(str, '  ', '');
        str = strrep(str, 'A0[', 'f[');
        str = strrep(str, '][0]', '*ng+i]');                          
        if isempty(i1)==1
            str = ['double ' str];
        end
        
        fprintf(gid, '\t\t%s\n', str);    
        tline = fgetl(fid);        
        i=i+1;   
    end
    fclose(fid);
    fprintf(gid, '\n');  
        
    str = '}';
    fprintf(gid, '\t%s\n', str);             
        
    str = '}';
    fprintf(gid, '%s\n', str);             
    fprintf(gid, '\n'); 
    
    fclose(gid);            
    
    delete('tmp1.c');    
end

% generate c code for the source function
if exist('s','var') && exist('udg','var')    
    s = reshape(s,[ncu,1]);    
    
    % generate temporary C file
    ccode(s,'file','tmp1.c')    

    delete(filename2);
    gid = fopen(filename2,'wt');
    
    if nd==2
        str = ['void sourceonly_' appname '2d(double *s, double *pg, double *udg, double *param, double time, int ng, int nc, int ncu, int nd, int ncd)'];
    elseif nd==3
        str = ['void sourceonly_' appname '3d(double *s, double *pg, double *udg, double *param, double time, int ng, int nc, int ncu, int nd, int ncd)'];
    end
    fprintf(gid, '%s\n', str); 
    str = '{';
    fprintf(gid, '%s\n', str);         
    
    for i=1:length(param)
        str = ['double param' num2str(i) ' = param[' num2str(i-1) '];'];
        fprintf(gid, '\t%s\n', str);                  
    end        

    str = 'for (int i = 0; i <ng; i++) {';
    fprintf(gid, '\n\t%s\n', str);         
        
    for i=1:length(pg)
        str = ['double x' num2str(i) ' = pg[' num2str(i-1) '*ng+i];'];
        fprintf(gid, '\t\t%s\n', str);                  
    end
    for i=1:length(udg)
        str = ['double u' num2str(i) ' = udg[' num2str(i-1) '*ng+i];'];
        fprintf(gid, '\t\t%s\n', str);                  
    end
    fprintf(gid, '\n');      
        
    fid = fopen('tmp1.c','r');    
    tline = fgetl(fid);         
    i=1; a1 = 0;       
    while ischar(tline)        
        str = tline;
        
        i1 = strfind(str,'[');        
        i2 = strfind(str,']');        
        if isempty(i1)==0    
            a2 = str2num(str((i1+1):(i2-1)));                        
            for j = a1:(a2-1)                
                strj = ['s[' num2str(j) '*ng+i] = 0.0;'];
                fprintf(gid, '\t\t%s\n', strj);                  
            end
            a1 = a2+1;              
        end
                
        str = strrep(str, '  ', '');
        str = strrep(str, 'A0[', 's[');
        str = strrep(str, '][0]', '*ng+i]');                          
        if isempty(i1)==1
            str = ['double ' str];
        end
        
        fprintf(gid, '\t\t%s\n', str);    
        tline = fgetl(fid);        
        i=i+1;   
    end
    fclose(fid);
    fprintf(gid, '\n');  
        
    str = '}';
    fprintf(gid, '\t%s\n', str);             
        
    str = '}';
    fprintf(gid, '%s\n', str);             
    fprintf(gid, '\n');     
    
    fclose(gid);
    
    delete('tmp1.c');    
elseif exist('filename2','var')
    gid = fopen(filename2,'wt');
    
    if nd==2
        str = ['void sourceonly_' appname '2d(double *s, double *pg, double *udg, double *param, double time, int ng, int nc, int ncu, int nd, int ncd)'];
    elseif nd==3
        str = ['void sourceonly_' appname '3d(double *s, double *pg, double *udg, double *param, double time, int ng, int nc, int ncu, int nd, int ncd)'];
    end
    
    fprintf(gid, '%s\n', str); 
    str = '{';
    fprintf(gid, '%s\n', str);         
    
    str = 'int i;';
    fprintf(gid, '\t%s\n', str);         
    str = 'for (i = 0; i <ng*ncu; i++)';
    fprintf(gid, '\t%s\n', str);                 
    str = 's[i] = 0.0;';
    fprintf(gid, '\t\t%s\n', str);             
    
    str = '}';
    fprintf(gid, '%s\n', str);             
    fprintf(gid, '\n');         
        
    fclose(gid);    
end


% generate c code for the fhat function
if exist('fhat','var') && exist('udg','var') && exist('uh','var')
    
    fhat = reshape(fhat,[ncu,1]);
    % generate temporary C file
    ccode(fhat(:),'file','tmp1.c')
    
    delete(filename3);
    gid = fopen(filename3,'wt');
    
    if nd==2
        str = ['void fhatonly_' appname '2d(double *fh, double *pg, double *udg, double *uh, double *nl, double *param, double time, int ng, int nc, int ncu, int nd, int ncd)'];
    elseif nd==3
        str = ['void fhatonly_' appname '3d(double *fh, double *pg, double *udg, double *uh, double *nl, double *param, double time, int ng, int nc, int ncu, int nd, int ncd)'];
    end
    fprintf(gid, '%s\n', str); 
    str = '{';
    fprintf(gid, '%s\n', str);         
    
    for i=1:length(param)
        str = ['double param' num2str(i) ' = param[' num2str(i-1) '];'];
        fprintf(gid, '\t%s\n', str);                  
    end        

    str = 'for (int i = 0; i <ng; i++) {';
    fprintf(gid, '\n\t%s\n', str);         
    
    for i=1:length(pg)
        str = ['double x' num2str(i) ' = pg[' num2str(i-1) '*ng+i];'];
        fprintf(gid, '\t\t%s\n', str);                  
    end
    for i=1:length(udg)
        str = ['double u' num2str(i) ' = udg[' num2str(i-1) '*ng+i];'];
        fprintf(gid, '\t\t%s\n', str);                  
    end
    for i=1:length(uh)
        str = ['double uh' num2str(i) ' = uh[' num2str(i-1) '*ng+i];'];
        fprintf(gid, '\t\t%s\n', str);                  
    end
    for i=1:length(nl)
        str = ['double nl' num2str(i) ' = nl[' num2str(i-1) '*ng+i];'];
        fprintf(gid, '\t\t%s\n', str);                  
    end
    fprintf(gid, '\n');      
        
    fid = fopen('tmp1.c','r');    
    tline = fgetl(fid); 
    i=1; a1 = 0;       
    while ischar(tline)        
        str = tline;
        
        i1 = strfind(str,'[');        
        i2 = strfind(str,']');        
        if isempty(i1)==0    
            a2 = str2num(str((i1+1):(i2-1)));                        
            for j = a1:(a2-1)                
                strj = ['fh[' num2str(j) '*ng+i] = 0.0;'];
                fprintf(gid, '\t\t%s\n', strj);                  
            end
            a1 = a2+1;              
        end
                
        str = strrep(str, '  ', '');
        str = strrep(str, 'A0[', 'fh[');
        str = strrep(str, '][0]', '*ng+i]');                          
        if isempty(i1)==1
            str = ['double ' str];
        end
        
        fprintf(gid, '\t\t%s\n', str);    
        tline = fgetl(fid);        
        i=i+1;   
    end        
    fclose(fid);
    fprintf(gid, '\n');  
        
    str = '}';
    fprintf(gid, '\t%s\n', str);             
    
    str = '}';
    fprintf(gid, '%s\n', str);             
    fprintf(gid, '\n');         
    
    fclose(gid);           
    delete('tmp1.c');
end

return;

% generate C code for the fbou function
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
