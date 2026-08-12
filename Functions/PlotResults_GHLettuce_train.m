UserColor    = 1/255*[0,0,0];

% % convert outside Co2 and humidity to desired units
% d(:,2)    = co2dens2ppm(d(:,3), d(:,2))*1e-3;
% d(:,4)    = vaporDens2rh(d(:,3), d(:,4));

if ops.nx==4
   q = 0; 
else
   q = ops.nu;
end

figure(1);clf;
for ll=1:ops.nx-q  %statenya diloop
    subplot(3,ops.nx-q,ll)
    stairs(Data(:,ll),'linewidth',1.5);grid;hold on
    axis tight
    if ll==1
        ylabel('$x_1$ (g/m$^2$)','fontsize',12,'interpreter','latex')
    elseif ll==2
        ylabel('$x_2$ (kg/m$^3$)','fontsize',12,'interpreter','latex')
        ylim([0 1.1*max(Data(:,ll))])
    elseif ll==3
        ylabel('$x_3$ ($^\circ$C)','fontsize',12,'interpreter','latex')
        ylim([0 1.1*max(Data(:,ll))])
    elseif ll==4
        ylabel('$x_4$ (kg/m$^3$)','fontsize',12,'interpreter','latex')
        ylim([0 1.1*max(Data(:,ll))])
    end
end

for ll=1:ops.nd
    subplot(3,ops.nx-q,ops.nx-q+ll)
    stairs(ops.t/c,d(:,ll),'linewidth',1.5);grid
    axis tight
    if ll==1
        ylabel('$d_1$ (W/m$^2$)','fontsize',13,'interpreter','latex')
    elseif ll==2
        ylabel('$d_2$ (kg/m$^3$)','fontsize',13,'interpreter','latex')
        ylim([0 1.1*max(d(:,ll))])
    elseif ll==3
        ylabel('$d_3$ ($^\circ$C)','fontsize',13,'interpreter','latex')
        ylim([0 1.1*max(d(:,ll))])
    elseif ll==4
        ylabel('$d_4$ (kg/m$^3$)','fontsize',13,'interpreter','latex')
         ylim([0 1.1*max(d(:,ll))])
        xlabel('Time (days)','fontsize',13,'interpreter','latex')
    end
end
for ll=1:ops.nu
    subplot(3,ops.nx-q,ops.nx-q+ops.nd+ll)
    stairs(ops.t/c,u(:,ll),'linewidth',1.5);grid;hold on
    axis tight
    if ll==1
        ylabel('$u_1$ (mg/m$^2$/s)','fontsize',13,'interpreter','latex')
        xlabel('Time (days)','fontsize',13,'interpreter','latex')
    elseif ll==2
        ylabel('$u_2$ (mm/s)','fontsize',13,'interpreter','latex')
        xlabel('Time (days)','fontsize',13,'interpreter','latex')
    elseif ll==3
        ylabel('$u_3$ (W/m$^2$)','fontsize',13,'interpreter','latex')
        xlabel('Time (days)','fontsize',13,'interpreter','latex')
    end
end

