function [stiffMat] = GaussQuadQ4(order,integrand)
%% written by Afnan Mostafa as part of ME 441 at UR
%GaussQuadQ4 evaluates the gauss integral for Q4 in isoparametric space to
%obtain the stiffness matrix of a Q4 element.
%
%   takes the GQ order (either 1 or 2 or 3) and integrand obtained from the
%   IntegrandStiffMatQ4.m file (function file) and then performs the GQ
%   integration to obtain stiffness matrix (integral)
%
% input: GQ order and integrand matrix in terms of e,n (xi, eta)
% output: Stiffness Matrix

%% %%%%%%%%%%% sanity check for for symbolic e,n %%%%%%%%%%%%

% if sum([strcmp(class(n),'sym'), strcmp(class(e),'sym')]) == 2
% 
% elseif sum([strcmp(class(n),'sym'), strcmp(class(e),'sym')]) < 2
    syms n e
% end

%% %%%%%%% gauss points and weights for 2d integration %%%%%%%

switch (order)
    case 1
        gaussPt = 0; wt = 4;
    case 2
        gaussPt = [
            -1/sqrt(3) -1/sqrt(3);
            1/sqrt(3) -1/sqrt(3);
            1/sqrt(3) 1/sqrt(3);
            -1/sqrt(3) 1/sqrt(3)];
        wt = [
            1 1;
            1 1;
            1 1;
            1 1];
    case 3
        gaussPt = [
            0 0;
            0 -sqrt(3/5);
            0 sqrt(3/5)
            sqrt(3/5) 0;
            -sqrt(3/5) 0;
            sqrt(3/5) sqrt(3/5);
            sqrt(3/5) -sqrt(3/5);
            -sqrt(3/5) sqrt(3/5);
            -sqrt(3/5) -sqrt(3/5)];
        wt = [
            8/9 8/9;
            8/9 5/9;
            8/9 5/9;
            5/9 8/9;
            5/9 8/9;
            5/9 5/9;
            5/9 5/9;
            5/9 5/9;
            5/9 5/9];
        % case 4
    otherwise
        error("Can't do more than 3rd order GQ")
end

%% %%%%%%%%%%%%% apply gauss quadrature %%%%%%%%%%%%%%%%

gq = [];
integral = zeros(length(integrand),length(integrand));
for i=1:length(integrand)
    for j=1:length(integrand)
        func = integrand(i,j);
        if order == 1
            gq = eval(subs(func,[e,n],[0,0]))*wt;
        else
            for k=1:length(gaussPt)
                % e = gaussPt(k,1);
                % n = gaussPt(k,2);
                gq(k) = eval(subs(func,[e,n],[gaussPt(k,1),...
                    gaussPt(k,2)]))*wt(k,1)*wt(k,2);
            end
        end
        integral(i,j) = sum([gq]);
    end
end
stiffMat = integral;
end

% =========================================================
% ~~~~~~~~~~~~~~~~~~ END OF FUNCTION ~~~~~~~~~~~~~~~~~~~~~~
% =========================================================