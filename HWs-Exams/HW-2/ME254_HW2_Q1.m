clear clc
clear all
close all

rng('shuffle')

%% given parameters
theta12 = 150;
theta14 = 45;
theta13 = 270;
E = 1e+9;
neu = 0.3;
A = .0025; %% 5cm x 5cm
L12 = 5/(sind(30));
L14 = 5/(sind(45));

%% Bar 1-2
x = theta12;
rot_matrix = [(cosd(x))^2 sind(x)*cosd(x) -(cosd(x))^2 -sind(x)*cosd(x);
    sind(x)*cosd(x) (sind(x))^2 -sind(x)*cosd(x) -(sind(x))^2;
    -(cosd(x))^2 -sind(x)*cosd(x) (cosd(x))^2 sind(x)*cosd(x);
    -sind(x)*cosd(x) -(sind(x))^2 sind(x)*cosd(x) (sind(x))^2];

K12 = ((E*A)/L12)*(rot_matrix);

%% Bar 1-4
x = theta14;
rot_matrix2 = [(cosd(x))^2 sind(x)*cosd(x) -(cosd(x))^2 -sind(x)*cosd(x);
    sind(x)*cosd(x) (sind(x))^2 -sind(x)*cosd(x) -(sind(x))^2;
    -(cosd(x))^2 -sind(x)*cosd(x) (cosd(x))^2 sind(x)*cosd(x);
    -sind(x)*cosd(x) -(sind(x))^2 sind(x)*cosd(x) (sind(x))^2];

K14 = ((E*A)/L14)*(rot_matrix2);

%% Beam 1-3
x = theta13;
rot_matrix3 = [cosd(x) sind(x) 0 0 0 0;
    -sind(x) cosd(x) 0 0 0 0;
    0 0 1 0 0 0;
    0 0 0 cosd(x) sind(x) 0;
    0 0 0 -sind(x) cosd(x) 0;
    0 0 0 0 0 1];

L = 5;
I = (0.05*0.05^3)/12; %% I = (b*h^3)/12

K13 = [0 0 0 0 0 0;
    0 12 6*L 0 -12 6*L;
    0 6*L  4*L^2 0 -6*L 2*L^2;
    0 0 0 0 0 0;
    0 -12 -6*L 0 12 -6*L;
    0 6*L 2*L^2 0 -6*L 4*L^2].*((E*I)/(L^3));

new_K13 = transpose(rot_matrix3)*K13*rot_matrix3;

%% reshaping matrix
%% Bar 12
gl_K12 = K12;
gl_K12(end+1,:) = 0;
gl_K12(end+1,:) = 0;
gl_K12(:,end+1) = 0;
gl_K12(:,end+1) = 0;

%% Bar 14
gl_K14 = K14;
gl_K14(end+2,:) = 0;
gl_K14(:,end+2) = 0;
gl_K14(5,:) = gl_K14(3,:);
gl_K14(6,:) = gl_K14(4,:);
gl_K14(3,:) = 0;
gl_K14(4,:) = 0;
gl_K14(:,5) = gl_K14(:,3);
gl_K14(:,6) = gl_K14(:,4);
gl_K14(:,3) = 0;
gl_K14(:,4) = 0;

%% Beam 13 globalization to 10x10
gl_K13 = new_K13;
gl_K13(end+4,:) = 0;
gl_K13(:,end+4) = 0;
gl_K13(:,8) = gl_K13(:,6);
gl_K13(:,7) = gl_K13(:,5);
gl_K13(:,6) = gl_K13(:,4);
gl_K13(:,4) = 0;
gl_K13(:,5) = 0;
gl_K13(8,:) = gl_K13(6,:);
gl_K13(7,:) = gl_K13(5,:);
gl_K13(6,:) = gl_K13(4,:);
gl_K13(4,:) = 0;
gl_K13(5,:) = 0;

%% BAR: rotation matrix globalization to 10x10
%% Bar 12 globalization
gl_K12_1 = gl_K12;
gl_K12_1(end+4,:) = 0;
gl_K12_1(:,end+4) = 0;
gl_K12_1(:,5) = gl_K12_1(:,4);
gl_K12_1(:,4) = gl_K12_1(:,3);
gl_K12_1(:,3) = 0;
gl_K12_1(5,:) = gl_K12_1(4,:);
gl_K12_1(4,:) = gl_K12_1(3,:);
gl_K12_1(3,:) = 0;

%% Bar 14 globalization
gl_K14_1 = gl_K14;
gl_K14_1(end+4,:) = 0;
gl_K14_1(:,end+4) = 0;
gl_K14_1(:,9)=gl_K14_1(:,5);
gl_K14_1(:,10)=gl_K14_1(:,6);
gl_K14_1(:,5) = gl_K14_1(:,4);
gl_K14_1(:,4) = gl_K14_1(:,3);
gl_K14_1(:,3) = 0;
gl_K14_1(:,6) = 0;
gl_K14_1(9,:)=gl_K14_1(5,:);
gl_K14_1(10,:)=gl_K14_1(6,:);
gl_K14_1(5,:)=0;
gl_K14_1(6,:)=0;

%% FINAL GLOBAL MATRIX (10x10)
K_global = gl_K12_1 + gl_K14_1 + gl_K13;

%% checking for singularity
% DOFs = inv(K_global);

%% applying boundary conditions

% K_global_BC = [K_global(1,1) K_global(1,2) K_global(1,3) K_global(1,4) K_global(1,9);
%     K_global(2,1) K_global(2,2) K_global(2,3) K_global(2,4) K_global(2,9);
%     K_global(3,1) K_global(3,2) K_global(3,3) K_global(3,4) K_global(3,9);
%     K_global(4,1) K_global(4,2) K_global(4,3) K_global(4,4) K_global(4,9);
%     K_global(9,1) K_global(9,2) K_global(9,3) K_global(9,4) K_global(9,9)];

K_global_BC = [K_global(1,1) K_global(1,2) K_global(1,3);
    K_global(2,1) K_global(2,2) K_global(2,3);
    K_global(3,1) K_global(3,2) K_global(3,3)];

DOFs = inv(K_global_BC);

% syms M

Force = [3000; -5000; 0];

answer = DOFs*Force;