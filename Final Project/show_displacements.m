function show_displacements(elemconn,coord,ut,ux,uy,NumNodes,customTitle)

figure;

%% %%%%%%%%%%%%%%%%%%%%% show mesh %%%%%%%%%%%%%%%%%%%%%

mesh=zeros(NumNodes,1);
subplot(2,1,1)
p=trisurf ( elemconn, coord(:,1), coord(:,2),mesh');
title('FEA Mesh');
axis equal tight;
view(2);
set(p,'edgecolor','k');
set(p,'LineWidth',1);
set(p,'Marker','d');
set(gca, 'FontName','Garamond','FontSize',12, 'FontWeight','bold')
% ylim([0,1])
% xlim([0,1])
xlabel('X (m)')
ylabel('Y (m)')
elmcon = reshape(elemconn, size(elemconn,1)*size(elemconn,2),1);
unq_elmcon = unique(elmcon, 'stable');

%% %%%%%%%%%%%%%%%% plot node numbers %%%%%%%%%%%%%%%%%

for h=1:NumNodes
        text(coord(unq_elmcon(h,1),1)+0.01,coord(unq_elmcon(h,1),2)+0.03,num2str(unq_elmcon(h,1)), 'Color','r')
end

%% %%%%%%%%%%%%%%%% Show displacements %%%%%%%%%%%%%%%%

subplot(2,1,2)
scale=1.;
displ=zeros(NumNodes,2);
displ(:,1)=(ux)*scale;
displ(:,2)=(uy)*scale;
distorted=coord+displ;
%   p=trisurf ( elemconn, distorted(:,1), distorted(:,2)+ut(:,2), ut' );
p=trisurf (elemconn, distorted(:,1), distorted(:,2), ut' );
axis equal tight;
view(2);
shading interp;
colorbar;
set(p,'edgecolor','k');
set(p,'LineWidth',1);
set(gcf, 'color', 'white');
set(gca, 'FontName','Garamond','FontSize',16,'FontWeight','bold')
title(customTitle, 'FontSize', 13)
% ylim([0,1])
% xlim([0,1])
xlabel('X (m)')
ylabel('Y (m)')
set(gcf,'units','points','position',[100,80,1000,600])
end

% =========================================================
% ~~~~~~~~~~~~~~~~~~ END OF FUNCTION ~~~~~~~~~~~~~~~~~~~~~~
% =========================================================