
icase = 6;
nproc = app.nfile;
if nproc == 1
    filename = appser.fileout;
else
    filename = apppar.fileout;

end
gam = app.arg{1};
[UDG,UH] = getsolfrombinaryfile(filename,nproc,master.npv,app.nc,master.npf,app.nch,app.hybrid);


switch (icase)
    case 1 % scalar case
        time = 7;
        % --- Evaluate L2 errors in time ---
        UDGnum = zeros(size(UDG));
        UDGnum(:,1,:) = 2 + sin(mesh.dgnodes(:,1,:) + mesh.dgnodes(:,2,:) - ...
            (UDG(:,2,:)+UDG(:,3,:))*time./UDG(:,1,:));
        UDGnum(:,2,:) = UDG(:,2,:)./UDG(:,1,:);  
        UDGnum(:,3,:) = UDG(:,3,:)./UDG(:,1,:);
        UDGnum(:,4,:) = UDG(:,4,:)./UDG(:,1,:);
        UDGnum(:,6,:) = UDG(:,6,:);
        UDGnum(:,7,:) = UDG(:,7,:);
        UDGnum(:,8,:) = UDG(:,8,:);
        vv = UDGnum(:,2,:).*UDGnum(:,2,:) + UDGnum(:,3,:).*UDGnum(:,3,:) + UDGnum(:,4,:).*UDGnum(:,4,:);
        bb = UDG(:,6,:).*UDG(:,6,:) + UDG(:,7,:).*UDG(:,7,:) + UDG(:,8,:).*UDG(:,8,:);
        UDGnum(:,5,:) = (gam-1)*(UDG(:,5,:) - UDGnum(:,1,:).*vv/2 -bb/2);
      %  UDGnum(:,9,:) = UDG(:,9,:);
        err = calerror(UDGnum,mesh,master,@exactsol,time)
        
    case 2 
        clf; scaplot(mesh,UDG(:,1,:),[],2,0); 
        axis on; axis square; axis tight; colorbar; colormap jet;
        drawnow
        
    case 3
        figure(1)
        clf; scaplot(mesh,UDG(:,1,:),[],2,0); 
        axis on; axis square; axis tight; colorbar; colormap jet;
        title('density');
        drawnow
    
        figure(2)
        clf; scaplot(mesh,UDG(:,4,:),[],2,0); 
        axis on; axis square; axis tight; colorbar; colormap jet;
        title('energy');
        drawnow
    
        p1 = size(UDG(:,1,:),1);
        p3 = size(UDG(:,1,:),3);
        for ie = 1:p3
            for ip = 1:p1
                c = sqrt(gam*5/UDG(ip,1,ie));
                ux = UDG(ip,2,ie); 
                uy = UDG(ip,3,ie); 
                normu = norm([ux uy]);
                machnum(ip,1,ie) = normu/c;
                bx = UDG(ip,5,ie); 
                by = UDG(ip,6,ie); 
                normb = norm([bx by])^2;
                magnepress(ip,1,ie) = normb/2;
            end
        end
        figure(3)
        clf; scaplot(mesh,machnum,[],2,0); 
        axis on; axis square; axis tight; colorbar; colormap jet;
        title('Mach number');
        xlim([0.3 0.7]);
        ylim([0.3 0.7]);
        drawnow
    
        figure(4)
        clf; scaplot(mesh,magnepress,[],2,0); 
        axis on; axis square; axis tight; colorbar; colormap jet;
        title('Magnetic pressure');
        drawnow
        
    case 4
        err = calerror(UDG,mesh,master,@exactsolAlfven);
        disp('Yes!');
        figure(1)
        clf; scaplot(mesh,UDG(:,1,:),[],2,0); 
        axis on; axis tight; colorbar; colormap jet;
        title('rho');
        drawnow
        
        figure(2)
        clf; scaplot(mesh,UDG(:,2,:),[],2,0); 
        axis on; axis tight; colorbar; colormap jet;
        title('rho*ux');
        drawnow 
        
        figure(3)
        clf; scaplot(mesh,UDG(:,3,:),[],2,0); 
        axis on; axis tight; colorbar; colormap jet;
        title('rho*uy');
        drawnow        
        
        figure(4)
        clf; scaplot(mesh,UDG(:,4,:),[],2,0); 
        axis on; axis tight; colorbar; colormap jet;
        title('rho*uz');
        drawnow      
        
        figure(5)
        clf; scaplot(mesh,UDG(:,5,:),[],2,0); 
        axis on; axis tight; colorbar; colormap jet;
        title('E');
        drawnow  
        
        figure(6)
        clf; scaplot(mesh,UDG(:,6,:),[],2,0); 
        axis on; axis tight; colorbar; colormap jet;
        title('Bx');
        drawnow     
        
        figure(7)
        clf; scaplot(mesh,UDG(:,7,:),[],2,0); 
        axis on; axis tight; colorbar; colormap jet;
        title('By');
        drawnow  
        
        figure(8)
        clf; scaplot(mesh,UDG(:,8,:),[],2,0); 
        axis on; axis tight; colorbar; colormap jet;
        title('Bz');
        drawnow
        
        figure(9)
        clf; scaplot(mesh,UDG(:,9,:),[],2,0); 
        axis on; axis tight; colorbar; colormap jet;
        title('psi');
        drawnow    
    case 5
        err = calerror(UDG,mesh,master,@exactsolvortex);
        disp('Yes!');
%         figure(1)
%         clf; scaplot(mesh,UDG(:,1,:),[],2,0); 
%         axis on; axis tight; colorbar; colormap jet;
%         title('rho');
%         drawnow
%         
%         figure(2)
%         clf; scaplot(mesh,UDG(:,2,:),[],2,0); 
%         axis on; axis tight; colorbar; colormap jet;
%         title('rho*ux');
%         drawnow 
%         
%         figure(3)
%         clf; scaplot(mesh,UDG(:,3,:),[],2,0); 
%         axis on; axis tight; colorbar; colormap jet;
%         title('rho*uy');
%         drawnow        
%         
%         figure(4)
%         clf; scaplot(mesh,UDG(:,4,:),[],2,0); 
%         axis on; axis tight; colorbar; colormap jet;
%         title('rho*uz');
%         drawnow      
%         
%         figure(5)
%         clf; scaplot(mesh,UDG(:,5,:),[],2,0); 
%         axis on; axis tight; colorbar; colormap jet;
%         title('E');
%         drawnow  
        
        figure(6)
        clf; scaplot(mesh,UDG(:,6,:),[-0.1592,0.1592],2,0); 
        axis on; axis tight; colorbar; colormap jet;
        title('Bx');
        drawnow     
        
%         figure(7)
%         clf; scaplot(mesh,UDG(:,7,:),[],2,0); 
%         axis on; axis tight; colorbar; colormap jet;
%         title('By');
%         drawnow  
%         
%         figure(8)
%         clf; scaplot(mesh,UDG(:,8,:),[],2,0); 
%         axis on; axis tight; colorbar; colormap jet;
%         title('Bz');
%         drawnow
%         
%         figure(9)
%         clf; scaplot(mesh,UDG(:,9,:),[],2,0); 
%         axis on; axis tight; colorbar; colormap jet;
%         title('psi');
%         drawnow  
    case 6
        disp('Yes!');
        %f = figure('visible','off');
        scaplot(mesh,UDG(:,1,:),[0.9,2.1],2,0); 
        axis on; axis tight; colormap jet;
        %print -dpng KH120.png 
        %close(f)
        %title('rho');
        drawnow
        
%         figure(2)
%         clf; scaplot(mesh,UDG(:,2,:),[],2,0); 
%         axis on; axis tight; colorbar; colormap jet;
%         title('rho*ux');
%         drawnow 
%         
%         figure(3)
%         clf; scaplot(mesh,UDG(:,3,:),[],2,0); 
%         axis on; axis tight; colorbar; colormap jet;
%         title('rho*uy');
%         drawnow        
% %         
%         figure(4)
%         clf; scaplot(mesh,UDG(:,4,:),[],2,0); 
%         axis on; axis tight; colorbar; colormap jet;
%         title('rho*uz');
%         drawnow      
%         
%         figure(5)
%         clf; scaplot(mesh,UDG(:,5,:),[],2,0); 
%         axis on; axis tight; colorbar; colormap jet;
%         title('E');
%         drawnow  
%         
%         figure(6)
%         clf; scaplot(mesh,UDG(:,6,:),[],2,0); 
%         axis on; axis tight; colorbar; colormap jet;
%         title('Bx');
%         drawnow     
%         
%         figure(7)
%         clf; scaplot(mesh,UDG(:,7,:),[],2,0); 
%         axis on; axis tight; colorbar; colormap jet;
%         title('By');
%         drawnow  
%         
%         figure(8)
%         clf; scaplot(mesh,UDG(:,8,:),[],2,0); 
%         axis on; axis tight; colorbar; colormap jet;
%         title('Bz');
%         drawnow
%         
%         figure(9)
%         clf; scaplot(mesh,UDG(:,9,:),[],2,0); 
%         axis on; axis tight; colorbar; colormap jet;
%         title('psi');
%         drawnow  
%         figure(10)
%         clf; scaplot(mesh,sqrt(UDG(:,6,:).^2+UDG(:,7,:).^2)./UDG(:,8,:),[],2,0); 
%         axis on; axis tight; colorbar; colormap jet;
%         title('KH');
%         drawnow         
    otherwise 

end