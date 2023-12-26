% ================================================
% ME 441: FEM, Final Project
% Afnan Mostafa
% 12/09/2023
% =================================================
%
%===================================================
% this script needs 7 functions: 
% 1. IntegrandStiffMatQ4.m,
% 2. JacobianMatQ4.m (called within IntegrandStiffMatQ4), 
% 3. GaussQuadQ4.m, 
% 4. globalizeStiffMat.m, 
% 5. construct_colMat.m, 
% 6. nodalStress.m, and 
% 7. show_displacements.m
%====================================================

%% %%%%%%%%%%%%%%%%%%%% clearing space %%%%%%%%%%%%%%%%%%%%%

clear
clc
close all
rng('shuffle')
tic

%% %%%%%%%%%%%%%%%%%%%% read data files %%%%%%%%%%%%%%%%%%%%

load coord.txt      %% x-y coordinate
load elemconn.txt   %% element node numbers
load bc_code.txt    %% boundary conditions
load loads.txt      %% force
load matprop.txt    %% material properties

%% %%%%%%%%%%%%%%%%%% determine work size %%%%%%%%%%%%%%%%%%

NumNodes = size(coord,1);
NumElements = size(elemconn,1);
DOF = size(loads,2);    %% u and v at each node

%% %%%%%%%%%% elastic modulus and Poisson's ratio %%%%%%%%%%

nu = matprop(2);
E_0 = matprop(3);

%% %%%%%%%%%%%%%%%%%%%% symbolic math %%%%%%%%%%%%%%%%%%%%%%

syms n e x y 
syms u v [1 2*NumNodes]

%% %%%%%%%%%%%%%%%% plane stress or strain %%%%%%%%%%%%%%%%%

% Thickness does not vary, so t = t. For plane strain, t=1
planechoice = input('Enter 1 for plane stress, 2 for plane strain: ');
if planechoice == 1
    t = matprop(1);
    isPlaneStrain = 0;
    isPlaneStress = 1;
elseif planechoice == 2
    t = 1;
    isPlaneStrain = 1;
    isPlaneStress = 0;
else % any other entry
    disp('invalid entry, selecting "Plane Stress" by default\n')
    t = matprop(1);
    isPlaneStrain = 0;
    isPlaneStress = 1;
end

%% %%%%%%%%%%%%% full or reduced integration %%%%%%%%%%%%%

intchoice = input(['Enter 1 for reduced integration,' ...
    ' 2 for full integration: ']);
order = intchoice;

%% %%%%%%%%%%%%%%%%%%%% xy matrix %%%%%%%%%%%%%%%%%%%%%%%

% storing coords
xy = cell(1,NumElements);
for n = 1:NumElements
    xy{1,n} = coord(elemconn(n,:),:);
end

%% %%%%%%%%%%%%%% cell creation for storage %%%%%%%%%%%%%%

% cell declaration for stiffness matrices and integrands
kmat_all = cell(1,NumElements);
intgrd_all = cell(1,NumElements);
B_all = cell(1,NumElements);
stiffMatSet = cell(NumElements,1);

%% %% function callback (1,2,3): get elemental stiff mat. %%

for p=1:NumElements
    [integrand,B,B_t,E] = IntegrandStiffMatQ4(xy{1,p},t,E_0,nu,...
        isPlaneStrain,isPlaneStress);
    intgrd_all{1,p} = integrand;
    B_all{1,p} = B;
    [stiffMat] = GaussQuadQ4(order,intgrd_all{1,p});
    kmat_all{1,p} = stiffMat;
end

%% %%% function callback (4): to get global stiff mat. %%%%%

for u=1:NumElements
    [mat_cp] = globalizeStiffMat(kmat_all{1,u},elemconn(u,:),...
        size(elemconn,2),NumNodes,DOF);
    stiffMatSet{u,1} = mat_cp;
end
globalMat = zeros(size(mat_cp));
for v=1:NumElements
    globalMat = globalMat + stiffMatSet{v,1};
end
globalMat_backup = globalMat;

%% %%%%%%% function callback (5): displacements %%%%%%%%%%%

[allU, allV] = construct_colMat(bc_code,1,0); %fixed=1, free=0

%% %%%%%%%%%%%%%% unified disp. matrices %%%%%%%%%%%%%%%%%%

unified_disp = sym(zeros(size(bc_code,1)*size(bc_code,2),1));
unified_disp(1:2:end) = allU;
unified_disp(2:2:end) = allV;
recover_disp = unified_disp;

%% %%%%%%%%%%%%%%%%% unified forces %%%%%%%%%%%%%%%%%%%%%%

unified_F = sym(zeros(size(loads,1)*size(loads,2),1));
unified_F(1:2:end) = loads(:,1);
unified_F(2:2:end) = loads(:,2);

%% %%%%%%%%%%%%%%%%%%%% apply BCs %%%%%%%%%%%%%%%%%%%%%%%%%

remove_DOF = find(unified_disp == 0);
unified_F(remove_DOF, :) = [];
unified_disp(remove_DOF, :) = [];
globalMat(remove_DOF, :) = [];
globalMat(:, remove_DOF) = [];

%% %%%%%%%%%%%%%%%%%% d = inv(K)*f %%%%%%%%%%%%%%%%%%%%%%%%

nod_disp = globalMat\unified_F;

%% %%%%%%% recover the displacements for all nodes %%%%%%%%

format long
solved_disp = setdiff(1:size(recover_disp,1),remove_DOF);
recover_disp(solved_disp) = vpa((nod_disp));

%% %%%%%%%%%%%%% check for reaction forces %%%%%%%%%%%%%%%%

vert_edge_nodes = find(coord(:,1) == 0);
horiz_edge_nodes = find(coord(:,2) == 0);
forces = globalMat_backup*recover_disp;
forcesX = forces(1:2:end); 
forcesY = forces(2:2:end);
reactionForceBase = ... 
    double(subs(sym(sqrt((sum(forcesX(vert_edge_nodes)))^2 ...
+ (sum(forcesY(horiz_edge_nodes)))^2))));

%% %%%%%%%%%%%% ux, uy, ut for plotting %%%%%%%%%%%%%%%%%%%

ux = double(recover_disp(1:2:end));
uy = double(recover_disp(2:2:end));
ut = sqrt((ux).^2 + (uy).^2);

%% %%%%% function callback (7): displacement contours %%%%%

show_displacements(elemconn,coord,ux,ux,uy,NumNodes,...
    'x Displacement Contours');

show_displacements(elemconn,coord,uy,ux,uy,NumNodes,...
    'y Displacement Contours');

show_displacements(elemconn,coord,ut,ux,uy,NumNodes,...
    'Total Displacement Contours');

%% %%% function callback (6): generate stress matrices %%%%

[sx,sy,sxy,sz,VM] = nodalStress(B_all,elemconn,recover_disp,E,...
    NumNodes,planechoice,nu);

%% %%%%%% function callback (7):  stress contours %%%%%%%%%%

show_displacements(elemconn,coord,sx,ux,uy,NumNodes,...
    'Normal Stress Contours along x');
show_displacements(elemconn,coord,sy,ux,uy,NumNodes,...
    'Normal Stress Contours along y');
show_displacements(elemconn,coord,sxy,ux,uy,NumNodes,...
    'Shear Stress Contours along xy');
show_displacements(elemconn,coord,VM,ux,uy,NumNodes,...
    'von Mises Stress');
show_displacements(elemconn,coord,sz,ux,uy,NumNodes,...
    'Normal Stress Contours along z');

disp('Done')
toc

% =========================================================
% ~~~~~~~~~~~~~~~~~~~ END OF SCRIPT ~~~~~~~~~~~~~~~~~~~~~~
% =========================================================
