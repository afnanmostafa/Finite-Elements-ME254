function [mat3] = globalizeStiffMat(mat,posNodes,eleNodes,totNodes,DoF)
%% written by Afnan Mostafa as part of ME 441 at UR
%globalizeStiffMat = globalizes element matrix if given nodes (CCW)
%
% inputs:
%   mat         = matrix for globalization, 
%   posNodes    = nodal positions (CCW) from bottom left, 
%   eleNodes    = how many nodes in an element, 
%   totNodes    = total nodes in the entire system, and
%   DoF         = 2 for u, v (per node)
%
% output: global stiffness matrix
%
%% %%%%%%%%%%%%%% define size %%%%%%%%%%%%%

mat_cp = mat;
diffDOF = totNodes*DoF - length(mat);
mat_cp(end+diffDOF,:) = 0;
mat_cp(:,end+diffDOF) = 0;
mat2 = zeros(size(mat_cp));

%% %%%%%%%%%%%%%% rearrange columns %%%%%%%%%%%%%

posDisp = [(posNodes.*2)-1; posNodes.*2];
p = 0;
for m = 1:eleNodes
    for n = 1:DoF
        mat2(:,posDisp(n,m)) = mat_cp(:,n+(2*p));
    end
    p = p+1;
end

%% %%%%%%%%%%%%%%% rearrange rows %%%%%%%%%%%%%%%

mat3 = zeros(size(mat2));
p = 0;
for m = 1:eleNodes
    for n = 1:DoF
        mat3(posDisp(n,m),:) = mat2(n+(2*p),:);
    end
    p = p+1;
end
end

% =========================================================
% ~~~~~~~~~~~~~~~~~~ END OF FUNCTION ~~~~~~~~~~~~~~~~~~~~~~
% =========================================================