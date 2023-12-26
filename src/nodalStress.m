function [sx,sy,sxy,sz,VM] = nodalStress(B_all,elemconn,recover_disp,E,NumNodes,planechoice,nu)
% written by Afnan Mostafa as part of ME 441 Project at UR
%nodalStress generates column matrix for nodal stresses
%
%   inputs:  B_all: strain-displacement matrix for all elements,
%            elemconn: element connection (nodal connectors), 
%            recover_disp: full displacment matrix (u, v),
%            E: constitutive matrix,
%            nu: Poisson's ratio,
%            NumNodes: No. of nodes, and
%            planechoice: plane stress/strain
%   outputs: freshly-baked stress column matrices (sx, sy, sxy, sz, VM)
%
% sx=normal stress in x; sy=normal stress in y; sxy=shear stress in xy,
% sz=normal stress in z; VM=von Mises stress
%
%%

syms e n

%% %%%% get elemental displacment matrix  %%%%%%%

element_size = size(B_all,2);
elem_disp = cell(1,element_size);
for g = 1:element_size
    for gg = 1:size(elemconn,2)
        elem_disp{1,g}((2*gg)-1,1) = recover_disp(2*elemconn(g,gg)-1);  
        elem_disp{1,g}((2*gg),1) = recover_disp(2*elemconn(g,gg));
    end
end

%% %%% get stresses at Gauss Point (centroid for reduced Q4) %%%

stress_GP = cell(1,element_size);
for b = 1:element_size
    stress_GP{1,b} = eval(subs((E*B_all{1,b}*elem_disp{1,b}),[e,n],[0,0]));
end
stress_GPmat = cell2mat(stress_GP);

%% %%%%% distribute centroidal stress to the nodes %%%%%%%%

%%%%%% way 1
% bot_edge = find(coord(:,2) == 0);
% left_edge = find(coord(:,1) == 0);
% curved_edge = find(coord(:,1) > (0.7*max(coord(:,1))) |...
%     coord(:,2) > (0.7*max(coord(:,2))));
% 
% nodes_edge = [bot_edge; left_edge; curved_edge];
% nodes_edge = unique(nodes_edge);
% 

%%%%%%% way 2
[~, ~, uidx] = unique(elemconn(:));
counts = accumarray(uidx, 1);

%% %%%%%%%%%% generate stress matrices %%%%%%%%%%%%%%%

sx = zeros(NumNodes,1);
sy = zeros(NumNodes,1); 
sxy = zeros(NumNodes,1);
sz = zeros(NumNodes,1);

for yy = 1:element_size
    for zz = 1:size(elemconn,2)
        sx(elemconn(yy,zz),1) = sx(elemconn(yy,zz),1) + ...
            (stress_GPmat(1,yy)/(counts(elemconn(yy,zz),1)));
        
        sy(elemconn(yy,zz),1) = sy(elemconn(yy,zz),1) + ...
            (stress_GPmat(2,yy)/(counts(elemconn(yy,zz),1)));
        
        sxy(elemconn(yy,zz),1) = sxy(elemconn(yy,zz),1) + ...
            (stress_GPmat(3,yy)/(counts(elemconn(yy,zz),1)));
    end
end

if planechoice == 1
    sigma1 = ((sx+sy)./2) + sqrt((((sx-sy).^2)./4)+sxy.^2);
    sigma2 = ((sx+sy)./2) - sqrt((((sx-sy).^2)./4)+sxy.^2);
    tau_max = (sx-sy)./2;
    VM = sqrt((sigma1).^2 - (sigma1.*sigma2) + (sigma2).^2);
    sz = zeros(NumNodes,1);
elseif planechoice == 2
    sz = nu.*(sx+sy);
    VM = sqrt((((sx-sy).^2) + ((sy-sz).^2) + ((sx-sz).^2) + 6*((sxy).^2))./2);
end

end

% =========================================================
% ~~~~~~~~~~~~~~~~~~ END OF FUNCTION ~~~~~~~~~~~~~~~~~~~~~~
% =========================================================