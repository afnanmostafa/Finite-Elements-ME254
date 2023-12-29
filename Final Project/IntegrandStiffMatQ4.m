function [integrand,B,B_t,E] = IntegrandStiffMatQ4(xyMat,t,E_0,nu,isPlaneStrain,isPlaneStress)
%% written by Afnan Mostafa as part of ME 441 at UR
%IntegrandStiffMatQ4 evaluates the "integrand" inside the integral of
%  stiffness matrix of Q4 elements in isoparametric space
%
%  takes the coords matrix, thickness, Poisson's ratio, Young's Mod, and
%  either plane strain or stress condition and then calls JacobianMatQ4
%  matrix to compute strain-displacement matrix (B) and integrand
%
%  JacobianMatQ4 takes the beta matrix (4x4) and Jacobian (2x2) and then
%  does: transpose(B)*E*B*thickness*determinant(J)
%  Please Note: B ~= (not equals) betaMat
%  B = strain-displacement matrix, betaMat = matrix of derivatives of N
%
%  inputs: [x y] matrix, thickness (t), Poisson's ratio (nu), Young's Mod
%         (E_0), either plane strain or stress condition (just use 1 or 0)
%         ex1: IntegrandStiffMatQ4(xyMat,0.001,1e6,0.5,0,1)
%         ex2: IntegrandStiffMatQ4(xyMat,0.001,1e6,0.5,1,0)
%         ex3: IntegrandStiffMatQ4(xyMat,0.001,1e6,0.5,0,0) (plane stress
%         by default)
%         ex4: IntegrandStiffMatQ4(xyMat,0.001,1e6,0.5,1,1) (ERROR)
%
%  outputs: integrand, B, transpose of B, constitutive matrix

%% %%%%%%%%%%%% sanity check for plane conds.%%%%%%%%%%%%%

if sum([isPlaneStrain, isPlaneStress]) == 2
    error("Can't use both plane strain and plane stress, use any one")
elseif sum([isPlaneStrain, isPlaneStress]) == 1
    % do nothing
else
   warning('Choosing plane stress condition by default')
end

%% %%%%%%%%%% sanity check for for symbolic e,n %%%%%%%%%%

% if sum([strcmp(class(n),'sym'), strcmp(class(e),'sym')]) == 2
%     % do nothing
% elseif sum([strcmp(class(n),'sym'), strcmp(class(e),'sym')]) < 2
    syms n e
% end

%% %%%%%%%%%%%%% call JacobianMatQ4 function %%%%%%%%%%%%

[J,~,betaMat] = JacobianMatQ4(xyMat);

%% %%%%%%%%%%%%%%%%%%%% Alpha Matrix %%%%%%%%%%%%%%%%%%%%

alphaMat = [1 0 0 0; 0 0 0 1; 0 1 1 0];

%% %%%%%%%%%%%%%%%%%%%% Beta Matrix %%%%%%%%%%%%%%%%%%%%%

%betaMat = [[invJ] [Zer]; [Zer] [invJ]];  %% no need to redefine

%% %%%%%%%%%%%%%%%%%%%% Gamma Matrix %%%%%%%%%%%%%%%%%%%%%

gammaMat = (1/4)*[-1+n 0 1-n 0 1+n 0 -1-n 0;
    -1+e 0 -1-e 0 1+e 0 1-e 0;
    0 -1+n 0 1-n 0 1+n 0 -1-n;
    0 -1+e 0 -1-e 0 1+e 0 1-e];

%% %%%%%%%%%%%%%%%%% Strain-Displacement Matrix %%%%%%%%%%%%

B = alphaMat*betaMat*gammaMat;
B_t = transpose(B);

%% %%%%%%%%%%%% Constitutive Matrix: plane stress %%%%%%%%%%%

if isPlaneStress == 1 && isPlaneStrain == 0
    E = (E_0/(1-(nu)^2))*[1 nu 0; nu 1 0; 0 0 (1-nu)/2];
elseif isPlaneStrain == 1 && isPlaneStress == 0
    E = (E_0/((1+nu)*(1-2*nu)))*[1-nu nu 0; nu 1-nu 0; 0 0 0.5-nu];
else
    disp('Choosing plane stress condition by default')
    E = (E_0/(1-(nu)^2))*[1 nu 0; nu 1 0; 0 0 (1-nu)/2];
end

%% %%%%%%%%%%%%%%%%% Stiffness Matrix %%%%%%%%%%%%%%%%%%%%%%%

integrand = eval(B_t*E*B*t*det(J));

end

% =========================================================
% ~~~~~~~~~~~~~~~~~~ END OF FUNCTION ~~~~~~~~~~~~~~~~~~~~~~
% =========================================================