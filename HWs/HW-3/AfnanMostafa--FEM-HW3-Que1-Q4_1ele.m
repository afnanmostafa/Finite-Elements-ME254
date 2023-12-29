%% FEM, HW 3, Problem 1
%% Afnan Mostafa
%% 11/28/2023

%% %%%%%%%%%%%%%%%%%%%% clearing space %%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
close all
rng('shuffle')

%% %%%%%%%%%%%%%%%%%%%% material properties %%%%%%%%%%%%%%%%%%%%%

E_0 = 1e6;
nu = 0.5;
t = 0.001;

%% %%%%%%%%%%%%%%%%%%%% symbolic math %%%%%%%%%%%%%%%%%%%%%%%%%%%

syms n e x y u2 u3 v3 u4 v4

%%           |
%            V
%  4 ___________________3
%  |                    |
%  |                    |
%  |                    |
%  |                    |
%  |         #1         |                
%  |                    |
%  |                    |
%  1____________________2
%  /\                  /\
%  ____________________oo__

%% %%%%%%%%%%%%%%%%%%%% coordinate: bottom left %%%%%%%%%%%%%%%%%

xyMat = [0 0; 1 0; 1 1; 0 1];

%% %%%%%%%%%%%%%%%%%%%% coordinate: middle %%%%%%%%%%%%%%%%%%%%%%

% xyMat = [-0.5 -0.5; 0.5 -0.5; 0.5 0.5; -0.5 0.5];

%% %%%%%%%%%%%%%%%%%%%% Jacobian Matrix %%%%%%%%%%%%%%%%%%%%%%%%%

J = 1/4*[-(1-n) (1-n) (1+n) -(1+n);
    -(1-e) -(1+e) (1+e) (1-e)]*xyMat;
invJ = inv(J);
Zer = zeros(size(invJ));

%% %%%%%%%%%%%%%%%%%%%% Alpha Matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%

alphaMat = [1 0 0 0; 0 0 0 1; 0 1 1 0];

%% %%%%%%%%%%%%%%%%%%%% Beta Matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

betaMat = [[invJ] [Zer]; [Zer] [invJ]];

%% %%%%%%%%%%%%%%%%%%%% Gamma Matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%

gammaMat = (1/4)*[-1+n 0 1-n 0 1+n 0 -1-n 0;
    -1+e 0 -1-e 0 1+e 0 1-e 0;
    0 -1+n 0 1-n 0 1+n 0 -1-n;
    0 -1+e 0 -1-e 0 1+e 0 1-e];

%% %%%%%%%%%%%%%%%%%%%% Strain-Displacement Matrix %%%%%%%%%%%%%%

B = alphaMat*betaMat*gammaMat;
B_t = transpose(B);

%% %%%%%%%%%%%%%%%%%%%% Constitutive Matrix: plane stress %%%%%%%

E = (E_0/(1-(nu)^2))*[1 nu 0; nu 1 0; 0 0 (1-nu)/2];

%% %%%%%%%%%%%%%%%%%%%% Stiffness Matrix %%%%%%%%%%%%%%%%%%%%%%%%

integrand = eval(B_t*E*B*t*det(J));

%% %%%%%%%%%%%%%%%%%%%% Gauss Quadrature (2nd order) %%%%%%%%%%%%

gaussPoints = [-1/sqrt(3) -1/sqrt(3);
    1/sqrt(3) -1/sqrt(3)
    1/sqrt(3) 1/sqrt(3)
    -1/sqrt(3) 1/sqrt(3)];

for i=1:length(integrand)
    for j=1:length(integrand)
        func = integrand(i,j);
        c=0;
        for k=1:4
            e = gaussPoints(k,1);
            n = gaussPoints(k,2);
            gq(k) = eval(subs(func));
        end
        integral(i,j) = gq(1)+gq(2)+gq(3)+gq(4);
    end
end

%% %%%%%%%%%%%%%%%%%%%% disp, force matrices %%%%%%%%%%%%%%%%%%%

disp = [0;0;u2;0;u3;v3;u4;v4];
force = [0;0.5;0;0.5;0;-0.5;0;-0.5];

%% %%%%%%%%%%%%%%%%%%%% apply BCs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp_BC = [u2;u3;v3;u4;v4];
force_BC = [0;0;-0.5;0;-0.5];
new_integral = integral;

%% %%%%%%%%%%%% remove K singularity %%%%%%%%%%%%%%%%%%%%%%%%%%%

new_integral(:,1) = [];
new_integral(:,1) = [];
new_integral(:,2) = [];
new_integral(1,:) = [];
new_integral(1,:) = [];
new_integral(2,:) = [];

%% %%%%%%%%%%%%%%%%%%%% solve Kd = F %%%%%%%%%%%%%%%%%%%%%%%%%%%

nod_disp = inv(new_integral)*force_BC;

allDisp = [
    0 0; nod_disp(1) 0; nod_disp(2) nod_disp(3); nod_disp(4) nod_disp(5)];

%% %%%%%%%%%%%%%%%%%%%% print out nodal displ %%%%%%%%%%%%%%%%%%

[id_v] = find(ismember(abs(allDisp(:,2)), max(abs(allDisp(:,2)))));
[id_u] = find(ismember(abs(allDisp(:,1)), max(abs(allDisp(:,1)))));

sprintf('Max Vertical displacement occurs at node %d: %0.4f units', id_v, allDisp(id_v,2))
sprintf('Max Horizontal displacement occurs at node %d: %0.4f units', id_u, allDisp(id_u,1))

%%           |
%            V
%  4 ___________________3
%  |                    |
%  |                    |
%  |                    |
%  |        #1          |
%  |                    |                
%  |                    |
%  |                    |
%  1____________________2
%  /\                  /\
%  ____________________oo__

%% %%%%%%%%%%% plot original and deformed systems %%%%%%%%%%%%

plotDeform=1;
if plotDeform
    hold on
    % plots the original system
    l1 = line([0,1],[0,0],'LineWidth',3,'Color','k');
    line([1,1],[0,1],'LineWidth',3,'Color','k');
    line([0,0],[0,1],'LineWidth',3,'Color','k');
    line([1,0],[1,1],'LineWidth',3,'Color','k');
    x = [0; 1; 1; 0; 0];
    y = [0; 0; 1; 1; 0];
    disp_u = [allDisp(1,1); allDisp(2,1); allDisp(3,1); allDisp(4,1); allDisp(1,1)];

    disp_v = [allDisp(1,2); allDisp(2,2); allDisp(3,2); allDisp(4,2); allDisp(1,2)];

    defX = x + disp_u;
    defY = y + disp_v;
    box on
    % plots the deformed system
    p1 = plot(defX, defY, 'r-', 'LineWidth', 1.5);
    set(gca,'FontName','Garamond','FontSize',18,'FontWeight','bold',...
        'LineWidth',2,'XMinorTick','off',...
        'YMinorTick','off','GridAlpha',0.07,...
        'GridLineStyle','--','LineWidth',2);
    title('Deformation Plot in Real Space (global coordinate system)');
    xlabel('X');
    ylabel('Y');
    legend([l1 p1],{'Original System', 'Deformed System'},'Location','southeast',...
        'Color',[0.941176470588235 0.941176470588235 0.941176470588235]);
    % set(gcf,'units','points','position',[100,100,1024,700])
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%