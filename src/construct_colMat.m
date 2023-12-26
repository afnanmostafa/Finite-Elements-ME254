function [allU, allV] = construct_colMat(bc_code,fixedDOFBool,freeDOFBool)
% written by Afnan Mostafa as part of ME 441 Profject at UR
%construct_colMat generates column matrix for displacement and load
%
%   inputs:  bc_code: text file with two columns (u v or fx fy),
%            fixedDOFBool: what boolean is set for fixed DOF (1 or 0) 
%            freeDOFBool: what boolean is set for free DOF (1 or 0) 
%   outputs: freshly-baked column matrices (u and v or fx or fy)
%
% it assumes the first and second columns to be u and v (or fx and fy)
% respectively, please modify the function for otherwise
%

%% %%%%% generate displacement  %%%%%%%%%

% convert to symbolic
allU = sym(bc_code(:,1));   %% modify here
allV = sym(bc_code(:,2));   %% modify here

% get indices of free dof
free_dofU = find(allU == freeDOFBool);
free_dofV = find(allV == freeDOFBool);

% replace free DOFs (initial value = 0) with symbolic variables
allU(free_dofU) = arrayfun(@(x) sym(['u' num2str(x)]), free_dofU);
allV(free_dofV) = arrayfun(@(x) sym(['v' num2str(x)]), free_dofV);

% replace fixed DOFs (initial value = 1) with symbolic variable 0
allU(allU == fixedDOFBool) = sym('0');
allV(allV == fixedDOFBool) = sym('0');
end

% =========================================================
% ~~~~~~~~~~~~~~~~~~ END OF FUNCTION ~~~~~~~~~~~~~~~~~~~~~~
% =========================================================