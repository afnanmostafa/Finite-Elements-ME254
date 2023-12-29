function [J,invJ,betaMat] = JacobianMatQ4(xyMat)
%% written by Afnan Mostafa as part of ME 441 at UR
%JacobianMat evaluates the Jacobian matrix for Q4 in isoparametric space
%
%   takes the [x y] matrix (4x2) and multiplies it with a prefactor and
%   another matrix (2x4) that consists of the derivatives of shape
%   functions w.r.t. eta and n. Also, it calculates the inverse of Jacobian
%   and then assembles the Beta matrix needed for stiffness matrix of Q4
%   elements in isoparametric space.
%
% input: [x y] matrix
% outputs: Jacobian matrix, inverse Jacobian matrix, and Beta Matrix

%%

syms n e

%% %%%%%%%%%%%% main body function %%%%%%%%%%%
J = (1/4)*[-(1-n) (1-n) (1+n) -(1+n);
    -(1-e) -(1+e) (1+e) (1-e)]*xyMat;

invJ = inv(J);

Zer = zeros(size(invJ));
betaMat = [invJ Zer; Zer invJ];

end

% =========================================================
% ~~~~~~~~~~~~~~~~~~ END OF FUNCTION ~~~~~~~~~~~~~~~~~~~~~~
% =========================================================