% FINAL PLOTS 
clc 
clear
close all

%% DEFINE GENERATE THE MESH 

mesh = 'almond'; %options : 'messi' , 'almond' , 'sphere' 
refine = 0;

switch mesh
    case 'messi'
        messi = 'cabeza lowpoly messi.stl';
        load 99_samples_Rough_mesh_messi.mat
        %load('100_samples_Rough_mesh_messi.mat', 'C')
        TR = stlread(messi);
        P = TR.Points; P = P/400;
        tri = TR.ConnectivityList;  
    case 'almond'
        load 101_samples_Rough_mesh_almond.mat
        almond = 'nasaAlmond.STL';
        TR = stlread(almond);
        P = TR.Points; P = P / 1000;            
        tri = TR.ConnectivityList;  
    case 'sphere'
        load('101_samples_Rough_mesh_sphere.mat')
        [ P, tri ] = generateSphereMesh( 4, 'oct' ); P = P';
        TR = triangulation(tri,P);
    otherwise 
        fprintf('Wrong mesh name!!!!!!!!!! \n')
end 

if refine == 1
    [P, tri] = barycentric_refine(P, tri); % BARYCENTRIC REFINE 
    [P, tri] = barycentric_refine(P, tri); % BARYCENTRIC REFINE 
    TR = triangulation(tri, P);   % new refined triangulation
end 


figure()
trisurf( tri, P(:,1), P(:,2), P(:,3), 'EdgeColor', 'k' );
axis equal 




%% FACTOR 0.1

P_rough = PS(:,:,1,1);
Covariance = C(:,:,1);
alpha = correlation_length(1)^2/2;
nbins = 25;
i0 = 500;                  % reference node
C_ref = Covariance(i0,:);           % covariance from node i0 to all nodes
% Geodesic-like distance using mesh edges
TR2 = triangulation(tri,P);
E = edges(TR2);
edge_w = vecnorm(P(E(:,1),:) - P(E(:,2),:),2,2);
Gmesh = graph(E(:,1),E(:,2),edge_w,size(P,1));
dist_geo = distances(Gmesh,i0).';
dist = dist_geo;
%bin
C_num = C_ref(:)/C_ref(i0);
edges_bin = linspace(0,max(dist),nbins+1);
[~,~,binID] = histcounts(dist,edges_bin);
valid = binID > 0;
dist_bin = accumarray(binID(valid),dist(valid),[],@mean);
C_bin    = accumarray(binID(valid),C_num(valid),[],@mean);
% Smooth theoretical curve
dist_smooth = linspace(0,max(dist_bin),10000);
Cth_smooth = exp(-dist_smooth.^2/(8*alpha*1));

figure('Color','w','Units','inches','Position',[1 1 10 7])

t = tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
Pall = [P; P_rough];
xl = [min(Pall(:,1)) max(Pall(:,1))];
yl = [min(Pall(:,2)) max(Pall(:,2))];
zl = [min(Pall(:,3)) max(Pall(:,3))];

% ---- Correlation / field plot
field = abs(Covariance(i0,:)./max(abs(Covariance(i0,:))));

% ---- Original mesh
ax1 = nexttile(1);
trisurf(tri,P(:,1),P(:,2),P(:,3), ...
    'EdgeColor','none','LineWidth',0.4,'Facecolor','interp');
title({'\textbf{Nominal}','\textbf{surface}'}, ...
    'Interpreter','latex','FontSize',25)
%title('Nominal Surface','Interpreter','latex','FontSize',20)
axis equal tight off
xlim(xl); ylim(yl); zlim(zl)
view(35,25)
light
colormap(ax1,"white")

%==============================
if strcmp(mesh,'almond')
    zmin = min(P(:,2));
    zmax = max(P(:,2));
    height_mesh = zmax - zmin;
x1 = 0.09;
x2 = 0.28;
y  = 0.70;
cap = 0.015;   % cap half-height
% Main horizontal line
annotation('line',[x1 x2],[y y],...
    'Color','k','LineWidth',2.5)
% Left cap
annotation('line',[x1 x1],[y-cap y+cap],...
    'Color','k','LineWidth',2.5)
% Right cap
annotation('line',[x2 x2],[y-cap y+cap],...
    'Color','k','LineWidth',2.5)
% Centered text
annotation('textbox',...
    [(x1+x2)/2-0.04, y-0.05, 0.08, 0.04],...
    'String',sprintf('%.2f m',height_mesh),...
    'Color','k',...
    'EdgeColor','none',...
    'HorizontalAlignment','center',...
    'FontSize',15,...
    'Interpreter','latex')
else 
    zmin = min(P(:,3));
    zmax = max(P(:,3));
    height_mesh = zmax - zmin;
    annotation('line',[0.105 0.105],[0.55 0.8],...
        'Color','k','LineWidth',2.5)
    annotation('line',[0.097 0.113],[0.8 0.8],...
        'Color','k','LineWidth',2.5)
    annotation('line',[0.097 0.113],[0.55 0.55],...
        'Color','k','LineWidth',2.5)
    annotation('textbox',[0.02 0.63 0.08 0.05],...
        'String',sprintf('%.2f m',height_mesh),...
        'Color','k',...
        'EdgeColor','none',...
        'FontSize',15,'Interpreter','latex')
end 
%===============================

% ---- Rough mesh
ax2 = nexttile(2);
trisurf(tri,P_rough(:,1),P_rough(:,2),P_rough(:,3), ...
    'EdgeColor','none','LineWidth',0.4,'Facecolor','interp');
title({'\textbf{Rough surface}','\textbf{realization}'}, ...
    'Interpreter','latex','FontSize',25)
%title('Rough Surface Realization','Interpreter','latex','FontSize',20)
axis equal tight off
xlim(xl); ylim(yl); zlim(zl)
view(35,25)
light
colormap(ax2,"white")

% ---- Covariance field
ax3 = nexttile(3);
trisurf(tri,P(:,1),P(:,2),P(:,3),field, ...
    'EdgeColor','none','LineWidth',0.4,'Facecolor','interp');
title({'\textbf{Normalized}','\textbf{Covariance}'}, ...
    'Interpreter','latex','FontSize',25)
axis equal tight off
xlim(xl); ylim(yl); zlim(zl)
view(35,25)
colormap(ax3,"hot")
clim([0 1])
cb = colorbar;
cb.Ticks = [0 0.5 1];
cb.TickLabels = {'$0$','$0.5$','$1$'};

cb.TickLabelInterpreter = 'latex';
cb.FontSize = 12;
cb.Color = 'k';
cb.LineWidth = 1.5;

% ---- Large bottom covariance validation plot
ax4 = nexttile(4,[1 3]);
plot(dist_smooth,Cth_smooth,'k-','LineWidth',2); hold on
plot(dist_bin,abs(C_bin),'ro','MarkerFaceColor','r')
title(sprintf('Correlation length:  $\\ell = %.5f$ m', correlation_length(1)), ...
    'Interpreter','latex', ...
    'FontSize',25)
xlabel('\textbf{Geodesic distance (meters)}','Interpreter','latex','FontSize',20)
ylabel('\textbf{Normalized covariance}','Interpreter','latex','FontSize',20)
legend('Theory','Monte Carlo','Interpreter','latex','FontSize',20)
grid on

%% 0.1 (2)
figure('Color','w','Units','inches','Position',[1 1 10 4])

% Main covariance plot
ax4 = axes;
plot(dist_smooth,Cth_smooth,'k-','LineWidth',2); hold on
plot(dist_bin,abs(C_bin),'ro','MarkerFaceColor','r')

title(sprintf('Correlation length:  $\\ell = %.5f$ m', correlation_length(1)), ...
    'Interpreter','latex','FontSize',40)

xlabel('\textbf{Geodesic distance (meters)}', ...
    'Interpreter','latex','FontSize',40)

ylabel('\textbf{Normalized covariance}', ...
    'Interpreter','latex','FontSize',30)

legend('Theory','Monte Carlo', ...
    'Interpreter','latex','FontSize',20, ...
    'Location','northeast')

grid on
ylim([0 1])
xlim([0 max(dist_smooth)])

% Save main axis position
mainPos = ax4.Position;

% Inset positions inside the main plot
w = 0.16;
h = 0.28;
y0 = mainPos(2) + 0.1*mainPos(4);

xA = mainPos(1) + 0.2*mainPos(3);
xB = mainPos(1) + 0.5*mainPos(3);
xC = mainPos(1) + 0.8*mainPos(3);

% Common limits
Pall = [P; P_rough];
xl = [min(Pall(:,1)) max(Pall(:,1))];
yl = [min(Pall(:,2)) max(Pall(:,2))];
zl = [min(Pall(:,3)) max(Pall(:,3))];

% ---- Nominal surface inset
ax1 = axes('Position',[xA y0 w h]);
trisurf(tri,P(:,1),P(:,2),P(:,3), ...
    'EdgeColor','none','FaceColor',[0.7 0.7 0.7]);
axis equal tight off
xlim(xl); ylim(yl); zlim(zl)
view(35,25)
light
title({'\textbf{Nominal}','\textbf{surface}'}, ...
    'Interpreter','latex','FontSize',25)

% ---- Rough surface inset
ax2 = axes('Position',[xB y0 w h]);
trisurf(tri,P_rough(:,1),P_rough(:,2),P_rough(:,3), ...
    'EdgeColor','none','FaceColor',[0.7 0.7 0.7]);
axis equal tight off
xlim(xl); ylim(yl); zlim(zl)
view(35,25)
light
title({'\textbf{Rough surface}','\textbf{realization}'}, ...
    'Interpreter','latex','FontSize',25)

% ---- Covariance field inset
ax3 = axes('Position',[xC y0 w h]);
trisurf(tri,P(:,1),P(:,2),P(:,3),field, ...
    'EdgeColor','none','FaceColor','interp');
axis equal tight off
xlim(xl); ylim(yl); zlim(zl)
view(35,25)
colormap(ax3,'hot')
clim([0 1])
title({'\textbf{Normalized}','\textbf{Covariance}'}, ...
    'Interpreter','latex','FontSize',25)

cb = colorbar(ax3);
cb.Ticks = [0 0.5 1];
cb.TickLabels = {'$0$','$0.5$','$1$'};
cb.TickLabelInterpreter = 'latex';
cb.FontSize = 9;
%% FACTOR 1
P_rough = PS(:,:,1,2);
Covariance = C(:,:,2);
alpha = correlation_length(2)^2/2;

C_ref = Covariance(i0,:);           % covariance from node i0 to all nodes
% Geodesic-like distance using mesh edges
TR2 = triangulation(tri,P);
E = edges(TR2);
edge_w = vecnorm(P(E(:,1),:) - P(E(:,2),:),2,2);
Gmesh = graph(E(:,1),E(:,2),edge_w,size(P,1));
dist_geo = distances(Gmesh,i0).';
dist = dist_geo;
%bin
C_num = C_ref(:)/C_ref(i0);
edges_bin = linspace(0,max(dist),nbins+1);
[~,~,binID] = histcounts(dist,edges_bin);
valid = binID > 0;
dist_bin = accumarray(binID(valid),dist(valid),[],@mean);
C_bin    = accumarray(binID(valid),C_num(valid),[],@mean);
% Smooth theoretical curve
dist_smooth = linspace(0,max(dist_bin),10000);
Cth_smooth = exp(-dist_smooth.^2/(8*alpha*1));

figure('Color','w','Units','inches','Position',[1 1 10 7])

t = tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
Pall = [P; P_rough];
xl = [min(Pall(:,1)) max(Pall(:,1))];
yl = [min(Pall(:,2)) max(Pall(:,2))];
zl = [min(Pall(:,3)) max(Pall(:,3))];

field = abs(Covariance(i0,:)./max(abs(Covariance(i0,:))));

% ---- Original mesh
ax1 = nexttile(1);
trisurf(tri,P(:,1),P(:,2),P(:,3), ...
    'EdgeColor','none','LineWidth',0.4,'Facecolor','interp');
title({'\textbf{Nominal}','\textbf{surface}'}, ...
    'Interpreter','latex','FontSize',25)
%title('Nominal Surface','Interpreter','latex','FontSize',20)
axis equal tight off
xlim(xl); ylim(yl); zlim(zl)
view(35,25)
light
colormap(ax1,"white")

%==============================
%==============================
if strcmp(mesh,'almond')
    zmin = min(P(:,1));
    zmax = max(P(:,1));
    height_mesh = zmax - zmin;
x1 = 0.09;
x2 = 0.28;
y  = 0.70;
cap = 0.015;   % cap half-height
% Main horizontal line
annotation('line',[x1 x2],[y y],...
    'Color','k','LineWidth',2.5)
% Left cap
annotation('line',[x1 x1],[y-cap y+cap],...
    'Color','k','LineWidth',2.5)
% Right cap
annotation('line',[x2 x2],[y-cap y+cap],...
    'Color','k','LineWidth',2.5)
% Centered text
annotation('textbox',...
    [(x1+x2)/2-0.04, y-0.05, 0.08, 0.04],...
    'String',sprintf('%.2f m',height_mesh),...
    'Color','k',...
    'EdgeColor','none',...
    'HorizontalAlignment','center',...
    'FontSize',15,...
    'Interpreter','latex')
else 
    zmin = min(P(:,3));
    zmax = max(P(:,3));
    height_mesh = zmax - zmin;
    annotation('line',[0.105 0.105],[0.55 0.8],...
        'Color','k','LineWidth',2.5)
    annotation('line',[0.097 0.113],[0.8 0.8],...
        'Color','k','LineWidth',2.5)
    annotation('line',[0.097 0.113],[0.55 0.55],...
        'Color','k','LineWidth',2.5)
    annotation('textbox',[0.02 0.63 0.08 0.05],...
        'String',sprintf('%.2f m',height_mesh),...
        'Color','k',...
        'EdgeColor','none',...
        'FontSize',15,'Interpreter','latex')
end 
%===============================
%===============================


% ---- Rough mesh
ax2 = nexttile(2);
trisurf(tri,P_rough(:,1),P_rough(:,2),P_rough(:,3), ...
    'EdgeColor','none','LineWidth',0.4,'Facecolor','interp');
title({'\textbf{Rough surface}','\textbf{realization}'}, ...
    'Interpreter','latex','FontSize',25)
%title('Rough Surface Realization','Interpreter','latex','FontSize',20)
axis equal tight off
xlim(xl); ylim(yl); zlim(zl)
view(35,25)
light
colormap(ax2,"white")

% ---- Covariance field
ax3 = nexttile(3);
trisurf(tri,P(:,1),P(:,2),P(:,3),field, ...
    'EdgeColor','none','LineWidth',0.4,'Facecolor','interp');
title({'\textbf{Normalized}','\textbf{Covariance}'}, ...
    'Interpreter','latex','FontSize',25)
axis equal tight off
xlim(xl); ylim(yl); zlim(zl)
view(35,25)
colormap(ax3,"hot")
clim([0 1])
cb = colorbar;
cb.Ticks = [0 0.5 1];
cb.TickLabels = {'$0$','$0.5$','$1$'};

cb.TickLabelInterpreter = 'latex';
cb.FontSize = 12;
cb.Color = 'k';
cb.LineWidth = 1.5;

% ---- Large bottom covariance validation plot
ax4 = nexttile(4,[1 3]);
plot(dist_smooth,Cth_smooth,'k-','LineWidth',2); hold on
plot(dist_bin,abs(C_bin),'ro','MarkerFaceColor','r')
title(sprintf('Correlation length:  $\\ell = %.5f$ m', correlation_length(2)), ...
    'Interpreter','latex', ...
    'FontSize',25)
xlabel('\textbf{Geodesic distance (meters)}','Interpreter','latex','FontSize',20)
ylabel('\textbf{Normalized covariance}','Interpreter','latex','FontSize',20)
legend('Theory','Monte Carlo','Interpreter','latex','FontSize',20)
grid on

%% 1 (2)

figure('Color','w','Units','inches','Position',[1 1 10 4])

% Main covariance plot
ax4 = axes;
plot(dist_smooth,Cth_smooth,'k-','LineWidth',2); hold on
plot(dist_bin,abs(C_bin),'ro','MarkerFaceColor','r')

title(sprintf('Correlation length:  $\\ell = %.5f$ m', correlation_length(1)), ...
    'Interpreter','latex','FontSize',40)

xlabel('\textbf{Geodesic distance (meters)}', ...
    'Interpreter','latex','FontSize',30)

ylabel('\textbf{Normalized covariance}', ...
    'Interpreter','latex','FontSize',30)

legend('Theory','Monte Carlo', ...
    'Interpreter','latex','FontSize',20, ...
    'Location','northeast')

grid on
ylim([0 1])
xlim([0 max(dist_smooth)])

% Save main axis position
mainPos = ax4.Position;

% Inset positions inside the main plot
w = 0.16;
h = 0.28;
y0 = mainPos(2) + 0.1*mainPos(4);

xA = mainPos(1) + 0.2*mainPos(3);
xB = mainPos(1) + 0.5*mainPos(3);
xC = mainPos(1) + 0.8*mainPos(3);

% Common limits
Pall = [P; P_rough];
xl = [min(Pall(:,1)) max(Pall(:,1))];
yl = [min(Pall(:,2)) max(Pall(:,2))];
zl = [min(Pall(:,3)) max(Pall(:,3))];

% ---- Nominal surface inset
ax1 = axes('Position',[xA y0 w h]);
trisurf(tri,P(:,1),P(:,2),P(:,3), ...
    'EdgeColor','none','FaceColor',[0.7 0.7 0.7]);
axis equal tight off
xlim(xl); ylim(yl); zlim(zl)
view(35,25)
light
title({'\textbf{Nominal}','\textbf{surface}'}, ...
    'Interpreter','latex','FontSize',25)

% ---- Rough surface inset
ax2 = axes('Position',[xB y0 w h]);
trisurf(tri,P_rough(:,1),P_rough(:,2),P_rough(:,3), ...
    'EdgeColor','none','FaceColor',[0.7 0.7 0.7]);
axis equal tight off
xlim(xl); ylim(yl); zlim(zl)
view(35,25)
light
title({'\textbf{Rough surface}','\textbf{realization}'}, ...
    'Interpreter','latex','FontSize',25)

% ---- Covariance field inset
ax3 = axes('Position',[xC y0 w h]);
trisurf(tri,P(:,1),P(:,2),P(:,3),field, ...
    'EdgeColor','none','FaceColor','interp');
axis equal tight off
xlim(xl); ylim(yl); zlim(zl)
view(35,25)
colormap(ax3,'hot')
clim([0 1])
title({'\textbf{Normalized}','\textbf{Covariance}'}, ...
    'Interpreter','latex','FontSize',25)

cb = colorbar(ax3);
cb.Ticks = [0 0.5 1];
cb.TickLabels = {'$0$','$0.5$','$1$'};
cb.TickLabelInterpreter = 'latex';
cb.FontSize = 9;


%% FACTOR 10
P_rough = PS(:,:,2,3);
Covariance = C(:,:,3);
alpha = correlation_length(3)^2/2;

C_ref = Covariance(i0,:);           % covariance from node i0 to all nodes
% Geodesic-like distance using mesh edges
TR2 = triangulation(tri,P);
E = edges(TR2);
edge_w = vecnorm(P(E(:,1),:) - P(E(:,2),:),2,2);
Gmesh = graph(E(:,1),E(:,2),edge_w,size(P,1));
dist_geo = distances(Gmesh,i0).';
dist = dist_geo;
%bin
C_num = C_ref(:)/C_ref(i0);
edges_bin = linspace(0,max(dist),nbins+1);
[~,~,binID] = histcounts(dist,edges_bin);
valid = binID > 0;
dist_bin = accumarray(binID(valid),dist(valid),[],@mean);
C_bin    = accumarray(binID(valid),C_num(valid),[],@mean);
% Smooth theoretical curve
dist_smooth = linspace(0,max(dist_bin),10000);
Cth_smooth = exp(-dist_smooth.^2/(8*alpha*0.88));

figure('Color','w','Units','inches','Position',[1 1 10 7])

t = tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
Pall = [P; P_rough];
xl = [min(Pall(:,1)) max(Pall(:,1))];
yl = [min(Pall(:,2)) max(Pall(:,2))];
zl = [min(Pall(:,3)) max(Pall(:,3))];

% ---- Correlation / field plot
             % reference node
field = abs(Covariance(i0,:)./max(abs(Covariance(i0,:))));

% ---- Original mesh
ax1 = nexttile(1);
trisurf(tri,P(:,1),P(:,2),P(:,3), ...
    'EdgeColor','none','LineWidth',0.4,'Facecolor','interp');
title({'\textbf{Nominal}','\textbf{surface}'}, ...
    'Interpreter','latex','FontSize',25)
%title('Nominal Surface','Interpreter','latex','FontSize',20)
axis equal tight off
xlim(xl); ylim(yl); zlim(zl)
view(35,25)
light
colormap(ax1,"white")

%==============================
%==============================
if strcmp(mesh,'almond')
    zmin = min(P(:,1));
    zmax = max(P(:,1));
    height_mesh = zmax - zmin;
x1 = 0.09;
x2 = 0.28;
y  = 0.70;
cap = 0.015;   % cap half-height
% Main horizontal line
annotation('line',[x1 x2],[y y],...
    'Color','k','LineWidth',2.5)
% Left cap
annotation('line',[x1 x1],[y-cap y+cap],...
    'Color','k','LineWidth',2.5)
% Right cap
annotation('line',[x2 x2],[y-cap y+cap],...
    'Color','k','LineWidth',2.5)
% Centered text
annotation('textbox',...
    [(x1+x2)/2-0.04, y-0.05, 0.08, 0.04],...
    'String',sprintf('%.2f m',height_mesh),...
    'Color','k',...
    'EdgeColor','none',...
    'HorizontalAlignment','center',...
    'FontSize',15,...
    'Interpreter','latex')
else 
    zmin = min(P(:,3));
    zmax = max(P(:,3));
    height_mesh = zmax - zmin;
    annotation('line',[0.105 0.105],[0.55 0.8],...
        'Color','k','LineWidth',2.5)
    annotation('line',[0.097 0.113],[0.8 0.8],...
        'Color','k','LineWidth',2.5)
    annotation('line',[0.097 0.113],[0.55 0.55],...
        'Color','k','LineWidth',2.5)
    annotation('textbox',[0.02 0.63 0.08 0.05],...
        'String',sprintf('%.2f m',height_mesh),...
        'Color','k',...
        'EdgeColor','none',...
        'FontSize',15,'Interpreter','latex')
end 
%===============================
%===============================

% ---- Rough mesh
ax2 = nexttile(2);
trisurf(tri,P_rough(:,1),P_rough(:,2),P_rough(:,3), ...
    'EdgeColor','none','LineWidth',0.4,'Facecolor','interp');
title({'\textbf{Rough surface}','\textbf{realization}'}, ...
    'Interpreter','latex','FontSize',25)
%title('Rough Surface Realization','Interpreter','latex','FontSize',20)
axis equal tight off
xlim(xl); ylim(yl); zlim(zl)
view(35,25)
light
colormap(ax2,"white")

% ---- Covariance field
ax3 = nexttile(3);
trisurf(tri,P(:,1),P(:,2),P(:,3),field, ...
    'EdgeColor','none','LineWidth',0.4,'Facecolor','interp');
title({'\textbf{Normalized}','\textbf{Covariance}'}, ...
    'Interpreter','latex','FontSize',25)
axis equal tight off
xlim(xl); ylim(yl); zlim(zl)
view(35,25)
colormap(ax3,"hot")
clim([0 1])
cb = colorbar;
cb.Ticks = [0 0.5 1];
cb.TickLabels = {'$0$','$0.5$','$1$'};

cb.TickLabelInterpreter = 'latex';
cb.FontSize = 12;
cb.Color = 'k';
cb.LineWidth = 1.5;

% ---- Large bottom covariance validation plot
ax4 = nexttile(4,[1 3]);
plot(dist_smooth,Cth_smooth,'k-','LineWidth',2); hold on
plot(dist_bin,abs(C_bin),'ro','MarkerFaceColor','r')
title(sprintf('Correlation length:  $\\ell = %.5f$ m', correlation_length(3)), ...
    'Interpreter','latex', ...
    'FontSize',25)
xlabel('\textbf{Geodesic distance (meters)}','Interpreter','latex','FontSize',20)
ylabel('\textbf{Normalized covariance}','Interpreter','latex','FontSize',20)
legend('Theory','Monte Carlo','Interpreter','latex','FontSize',20)
grid on


%% 10 (2)
figure('Color','w','Units','inches','Position',[1 1 10 4])

% Main covariance plot
ax4 = axes;
plot(dist_smooth,Cth_smooth,'k-','LineWidth',2); hold on
plot(dist_bin,abs(C_bin),'ro','MarkerFaceColor','r')

title(sprintf('Correlation length:  $\\ell = %.5f$ m', correlation_length(1)), ...
    'Interpreter','latex','FontSize',40)

xlabel('\textbf{Geodesic distance (meters)}', ...
    'Interpreter','latex','FontSize',30)

ylabel('\textbf{Normalized covariance}', ...
    'Interpreter','latex','FontSize',30)

legend('Theory','Monte Carlo', ...
    'Interpreter','latex','FontSize',20, ...
    'Location','northeast')

grid on
ylim([0 1])
xlim([0 max(dist_smooth)])

% Save main axis position
mainPos = ax4.Position;

% Inset positions inside the main plot
w = 0.2;
h = 0.3;
y0 = mainPos(2) + 0.1*mainPos(4);

xA = mainPos(1) + 0.2*mainPos(3);
xB = mainPos(1) + 0.5*mainPos(3);
xC = mainPos(1) + 0.8*mainPos(3);

% Common limits
Pall = [P; P_rough];
xl = [min(Pall(:,1)) max(Pall(:,1))];
yl = [min(Pall(:,2)) max(Pall(:,2))];
zl = [min(Pall(:,3)) max(Pall(:,3))];

% ---- Nominal surface inset
ax1 = axes('Position',[xA y0 w h]);
trisurf(tri,P(:,1),P(:,2),P(:,3), ...
    'EdgeColor','none','FaceColor',[0.7 0.7 0.7]);
axis equal tight off
xlim(xl); ylim(yl); zlim(zl)
view(35,25)
light
title({'\textbf{Nominal}','\textbf{surface}'}, ...
    'Interpreter','latex','FontSize',25)

% ---- Rough surface inset
ax2 = axes('Position',[xB y0 w h]);
trisurf(tri,P_rough(:,1),P_rough(:,2),P_rough(:,3), ...
    'EdgeColor','none','FaceColor',[0.7 0.7 0.7]);
axis equal tight off
xlim(xl); ylim(yl); zlim(zl)
view(35,25)
light
title({'\textbf{Rough surface}','\textbf{realization}'}, ...
    'Interpreter','latex','FontSize',25)

% ---- Covariance field inset
ax3 = axes('Position',[xC y0 w h]);
trisurf(tri,P(:,1),P(:,2),P(:,3),field, ...
    'EdgeColor','none','FaceColor','interp');
axis equal tight off
xlim(xl); ylim(yl); zlim(zl)
view(35,25)
colormap(ax3,'hot')
clim([0 1])
title({'\textbf{Normalized}','\textbf{Covariance}'}, ...
    'Interpreter','latex','FontSize',25)

cb = colorbar(ax3);
cb.Ticks = [0 0.5 1];
cb.TickLabels = {'$0$','$0.5$','$1$'};
cb.TickLabelInterpreter = 'latex';
cb.FontSize = 9;

function [P2, tri2] = barycentric_refine(P, tri)

    nV = size(P,1);
    nT = size(tri,1);

    % Barycenter of each triangle
    C = (P(tri(:,1),:) + P(tri(:,2),:) + P(tri(:,3),:)) / 3;

    % Add barycenters as new vertices
    P2 = [P; C];

    % Index of each new barycenter
    cidx = nV + (1:nT)';

    % Each triangle becomes 3 triangles
    tri2 = [
        tri(:,1), tri(:,2), cidx;
        tri(:,2), tri(:,3), cidx;
        tri(:,3), tri(:,1), cidx
    ];

end