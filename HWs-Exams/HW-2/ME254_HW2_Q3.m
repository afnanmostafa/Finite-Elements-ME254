%% Afnan Mostafa
%% ME 441 FEM HW 2 Q 3
%% 10/12/2023

%% ~~~ clear space ~~~ %%
clear
clc
close all
rng('shuffle')

%% ~~~  ~~~ %%
[X, Y] = meshgrid(0.5:0.005:1,0.5:0.005:1);

%% exact solution
Z = surf(X,Y,sin(X.^2+Y.^2));

%% 4 nodal values
Z1 = sin(0.5^2+0.5^2);
Z2 = sin(1^2+0.5^2);
Z3 = sin(1^2+1^2);
Z4 = sin(0.5^2+1^2);

%% 4 shape functions
N1 = 4.*X.*Y-4.*X-4.*Y+4;
N2 = -4.*X.*Y+4.*X+2.*Y-2;
N3 = 4.*X.*Y-2.*X-2.*Y+1;
N4 = -4.*X.*Y+2.*X+4.*Y-2;

%% Field variable equation
ZZ = Z1*N1 + Z2*N2 + Z3*N3 + Z4*N4;

%% co-efficients of a0, a1, a2, a3 (from the submitted HW paper)
fitZ = -0.9694 + 1.9584*X + 1.9584*Y - 2.0384*X.*Y;

%% plot exact
hold on
surf(X,Y,fitZ);

%% Plot shape functions (use plot_shapeFunc=1 for plotting shape fucntions)
plot_shapeFunc=0;

if plot_shapeFunc==1
    figure
    
    xlim([0.5,1])
    ylim([0.5,1])
    s1=subplot(2,2,1);
    surf(X,Y,N1)
    box on
    set(gca,'FontName','Garamond','FontSize',14,'FontWeight','bold',...
        'LineWidth',1.5,'XMinorTick','off',...
        'YMinorTick','off','GridAlpha',0.07,...
        'GridLineStyle','--','LineWidth',1);
    
    ylabel('Y',...
        'FontName','Garamond','FontSize',16)
    xlabel('X',...
        'FontName','Garamond','FontSize',16);
    zlabel('N_1',...
        'FontName','Garamond','FontSize',16);
    
    s2=subplot(2,2,2);
    surf(X,Y,N2)
    box on
    set(gca,'FontName','Garamond','FontSize',14,'FontWeight','bold',...
        'LineWidth',1.5,'XMinorTick','off',...
        'YMinorTick','off','GridAlpha',0.07,...
        'GridLineStyle','--','LineWidth',1);
    
    ylabel('Y',...
        'FontName','Garamond','FontSize',16)
    xlabel('X',...
        'FontName','Garamond','FontSize',16);
    zlabel('N_2',...
        'FontName','Garamond','FontSize',16);
    
    s3=subplot(2,2,3);
    surf(X,Y,N3)
    box on
    set(gca,'FontName','Garamond','FontSize',14,'FontWeight','bold',...
        'LineWidth',1.5,'XMinorTick','off',...
        'YMinorTick','off','GridAlpha',0.07,...
        'GridLineStyle','--','LineWidth',1);
    
    ylabel('Y',...
        'FontName','Garamond','FontSize',16)
    xlabel('X',...
        'FontName','Garamond','FontSize',16);
    zlabel('N_3',...
        'FontName','Garamond','FontSize',16);
    
    s4=subplot(2,2,4);
    surf(X,Y,N4)
    
    box on
    set(gca,'FontName','Garamond','FontSize',14,'FontWeight','bold',...
        'LineWidth',1.5,'XMinorTick','off',...
        'YMinorTick','off','GridAlpha',0.07,...
        'GridLineStyle','--','LineWidth',1);
    
    ylabel('Y',...
        'FontName','Garamond','FontSize',16)
    xlabel('X',...
        'FontName','Garamond','FontSize',16);
    zlabel('N_4',...
        'FontName','Garamond','FontSize',16);
    
    hold on
    box on
end