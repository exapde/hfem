
for k=2:4
    fn = ['hedgp' num2str(k)];
    load(fn);
    
%     figure(1); clf;
%     set(axes, 'FontSize', 16, 'LineWidth', 1.0, 'TickLength', [0.015 0.015]);
%     scaplot(mesh,eulereval(UDG(:,1:4,:),'p',1.4),[],2,0); axis off;
%     axis([-0.5 1.5 -0.75 0.75]);
%     colorbar('FontSize',14);
%     fn = ['nacapreshedgp' num2str(porder)];
%     print('-depsc',fn);
% 
%     figure(2); clf;
%     set(axes, 'FontSize', 16, 'LineWidth', 1.0, 'TickLength', [0.015 0.015]);
%     scaplot(mesh,eulereval(UDG(:,1:4,:),'M',1.4),[],2,0); axis off;
%     axis([-0.5 1.5 -0.75 0.75]);
%     colorbar('FontSize',14);   
%     fn = ['nacamachhedgp' num2str(porder)];
%     print('-depsc',fn);
% 
%     figure(3); clf;
%     set(axes, 'FontSize', 16, 'LineWidth', 1.0, 'TickLength', [0.015 0.015]);
%     scaplot(mesh,UDG(:,5,:),[],2,0); axis off;
%     axis([-0.5 2.5 -0.75 0.75]);
%     colorbar('FontSize',14);
%     fn = ['nacaeddyhedgp' num2str(porder)];
%     print('-depsc',fn);
%     
    if k==2
        [Cp2,Cf2,x2]=getsurfacedata(master,mesh,app,UDG,UH,2);
    elseif k==3
        [Cp3,Cf3,x3]=getsurfacedata(master,mesh,app,UDG,UH,2);
    elseif k==4
        [Cp4,Cf4,x4]=getsurfacedata(master,mesh,app,UDG,UH,2);
    end
end

[LCp,UCp,Cpxfoil,Cfxfoil]=naca0012_exp_data;
figure(4); clf;
set(axes, 'FontSize', 18, 'LineWidth', 1.0, 'TickLength', [0.015 0.015], 'GridLineStyle', '-', ...
    'XTick', [0 0.2 0.4 0.6 0.8 1],'YTick', [-0.5 -0.25 0 0.25 0.5 0.75 1], 'YDir', 'Reverse');
hold on;
plot(UCp(:,1),UCp(:,3), 'ok', 'LineWidth', 1.5);
plot(x2(1:end-1), -Cp2(1:end-1), '-b','LineWidth', 1.5);
plot(x3(1:end-1), -Cp3(1:end-1), '--r','LineWidth', 1.5);
plot(x4(1:end-1), -Cp4(1:end-1), '-.g','LineWidth', 1.5);
xlabel('$x/c$','FontSize', 24, 'Interpreter', 'latex');
ylabel('$c_p$','FontSize', 24, 'Interpreter', 'latex');
legend({'Experiment', '$p=2$', '$p=3$', '$p=4$'},'Location','SE','Interpreter','latex');
axis([0 1 -0.5 1]);
box on; 
grid on;
hold off;
fn = ['nacacoefhedgp' num2str(porder)];
print('-depsc',fn);




