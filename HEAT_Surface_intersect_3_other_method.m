% SOLVE HEAT EQUATION OVER A DEFINED SURFACE MESH
clear;
clc;
close all;
%% DEFINE GENERATE THE MESH 

mesh = 'sphere'; %options : 'messi' , 'almond' , 'sphere' 
refine = 0;

switch mesh
    case 'messi'
        messi = 'cabeza lowpoly messi.stl';
        TR = stlread(messi);
        P = TR.Points; P = P/400;
        tri = TR.ConnectivityList;  
    case 'almond'
        almond = 'nasaAlmond.STL';
        TR = stlread(almond);
        P = TR.Points; P = P / 100;            
        tri = TR.ConnectivityList;  
    case 'sphere'
        [ P, tri ] = generateSphereMesh( 4, 'oct' ); P = P';
        TR = triangulation(tri,P);
    otherwise 
        fprintf('Wrong mesh name!!!!!!!!!! \n')
end 

if refine == 1
    [P, tri] = barycentric_refine(P, tri); % BARYCENTRIC REFINE 
    %[P, tri] = barycentric_refine(P, tri); % BARYCENTRIC REFINE 
    TR = triangulation(tri, P);   % new refined triangulation
end 

figure()
trisurf( tri, P(:,1), P(:,2), P(:,3), 'EdgeColor', 'k' );
axis equal 

%% ROUGH SURFACE GENERATION AND MONTE CARLO

E = edges(TR);
edgeLength = vecnorm(P(E(:,1),:) - P(E(:,2),:),2,2);
h_min = min(edgeLength);
h_mean= mean(edgeLength);

% VARIABLES ===============================================================
h = median(edgeLength);          % better than h_mean for irregular meshes
h = h_mean;
lambda = 5*h;
r_target = 0.007;                 % good starting value, between 0.01 and 0.1
t_final = 1; % final time
alpha = (10*lambda/2)^2*(1/t_final);
dt = r_target*h^2/alpha;
nsteps = ceil(t_final/dt);
correlation_length = 2*sqrt(alpha*t_final); % correlation length
epsilon = 1*lambda; % maximum perturbation height after transport map
a = 2; % hyperparameter beta distribution
self_intersection = 1; %1:  for self intersection check 0: for no self intersection check
samples = 2;

% Matrix Building 

[K,M] = FEM_assembly(P,tri); % FEM matrices for solving heat equation
normal_vec = compute_vertex_normals(P,tri);

% CRANK-NICHOLSON ==========

A = M + (dt/2)*alpha*K;
B = M - (dt/2)*alpha*K;
A_solver = decomposition(A,'chol');


% Precompute the face pairs that could intersect for any normal displacement
% RFx in [-epsilon, epsilon]. This broad phase is reused by every sample.
if self_intersection == 1
    candidatePairs0 = buildSweptAABBCandidatePairs(P,tri,normal_vec,epsilon);
    fprintf('Precomputed %d possible self-intersection face pairs \n',size(candidatePairs0,1))
else
    candidatePairs0 = zeros(0,2);
end

%% MONTE CARLO
tic 
U = zeros(length(P),samples);
PS = zeros(length(P),3,samples);
noise = 1;
for s = 1:samples 
    fprintf('Rough surface %d \n',s)
    [PS(:,:,s),U(:,s)] = generate_RF_heatequation_surface2(P,tri,normal_vec,noise,nsteps,M,A,A_solver,B,epsilon,a,self_intersection,candidatePairs0);
end 
toc
P_rough = PS(:,:,1);

%% DISCRETIZATION QUANTITIES

fprintf('epsilon    = %.6g\n',epsilon)
fprintf('correlation lenght  = %.6g\n',correlation_length)
fprintf('alpha   = %.6g\n',alpha)
fprintf('T final   = %.6g\n',t_final)
fprintf('h_min (and time discretization dt)  = %.6g\n',h_min)
fprintf('dt*alpha   = %.6g\n',dt*alpha)

% COVARIANCE CALCULATION AND VERIFICATION
U_sample = U;
C = cov(U_sample');
i0 = 1;                  % reference node
C_ref = C(i0,:);           % covariance from node i0 to all nodes
% Geodesic-like distance using mesh edges
TR2 = triangulation(tri,P);
E = edges(TR2);
edge_w = vecnorm(P(E(:,1),:) - P(E(:,2),:),2,2);
Gmesh = graph(E(:,1),E(:,2),edge_w,size(P,1));
dist_geo = distances(Gmesh,i0).';
dist = dist_geo;
%bin
C_num = C_ref(:)/C_ref(i0);
nbins = 50;
edges_bin = linspace(0,max(dist),nbins+1);
[~,~,binID] = histcounts(dist,edges_bin);
valid = binID > 0;
dist_bin = accumarray(binID(valid),dist(valid),[],@mean);
C_bin    = accumarray(binID(valid),C_num(valid),[],@mean);
% Smooth theoretical curve
dist_smooth = linspace(0,max(dist_bin),10000);
Cth_smooth = exp(-dist_smooth.^2/(8*alpha*t_final));

% figure()
% plot(dist_smooth,Cth_smooth,'k-','LineWidth',2.5); hold on
% plot(dist_bin,C_bin,'ro','MarkerFaceColor','r','MarkerSize',7)
% 
% xlabel('Geodesic distance','Interpreter','latex')
% ylabel('Normalized covariance','Interpreter','latex')
% legend('Theory','Monte Carlo','Interpreter','latex')
% grid on

%%
figure('Color','w','Units','inches','Position',[1 1 10 7])

t = tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
Pall = [P; P_rough];
xl = [min(Pall(:,1)) max(Pall(:,1))];
yl = [min(Pall(:,2)) max(Pall(:,2))];
zl = [min(Pall(:,3)) max(Pall(:,3))];

% ---- Correlation / field plot
i0 = 1;
field = abs(C(i0,:)./max(abs(C(i0,:))));

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
title(sprintf('Correlation length:  $\\ell = %.5f$ m', correlation_length), ...
    'Interpreter','latex', ...
    'FontSize',15)
xlabel('\textbf{Geodesic distance (meters)}','Interpreter','latex','FontSize',15)
ylabel('\textbf{Normalized covariance}','Interpreter','latex','FontSize',15)
legend('Theory','Monte Carlo','Interpreter','latex','FontSize',15)
grid on

%%
P=P_rough;
tria=tri;
save('rougsurfaceTEST.mat','P','tria')

%% TEST A MESH -> DETECT SELF INTERSECTION OF A DATA SET

% % Read mesh from .mat file
% load('101_samples_Rough_mesh_sphere.mat');   % change filename
% 
% % Extract mesh variables
% P   = PS(:,:,2,3);    % or S.P if your file stores original mesh
% 
% % Detect self-intersections
% tol = 1e-12;
% [TF,badNodes,badFaces,badPairs] = detectSelfIntersectionNodesFast(P,tri,tol);
% 
% % Print result
% if TF
%     fprintf('Mesh HAS self-intersections.\n')
%     fprintf('Number of intersecting face pairs: %d\n', size(badPairs,1))
%     fprintf('Number of bad faces: %d\n', numel(badFaces))
%     fprintf('Number of bad nodes: %d\n', numel(badNodes))
% else
%     fprintf('Mesh has NO self-intersections.\n')
% end
%% FUNCTIONS 

function [P_rough,U] = generate_RF_heatequation_surface(P,tri,normal_vec,noise,nsteps,M,A,A_solver,B,epsilon,a,self_intersection,candidatePairs0)
    % DIFFERENT DEFINITIONS OF WHITE NOISE 
    W = randn(size(P,1),1);
    if noise == 1
        mlump = sum(M,2);
        % Define white noise 
        U0 = W ./ sqrt(mlump);
    end 
    % if noise == 2
    %     U0 = M\(W*2*A);
    % end 
    % 
    % if noise == 3
    %     U0 = M\(W*A);
    % end 
    %% CRANK NICHOLSON
    Un = U0;
    for n = 1:nsteps
        U = A_solver\(B*Un);
        Un=U;
    end 
    U = Un;
    %% TRANSPORT MAP
    % TR = triangulation(tri,P);
    % N = vertexNormal(TR);
    % normal_vec = N ./ vecnorm(N,2,2);
    w = mlump / sum(mlump);
    mu = sum(w .* Un);
    sig = sqrt(sum(w .* (Un - mu).^2));
    Un = (Un - mu) / sig;
    Fnx = normcdf(Un, 0, 1); 
    Fbx = betainv(Fnx, a, a); 
    RFx = 2 * epsilon * Fbx - epsilon;
    P_rough = P + (RFx.*normal_vec);
    %% ADD SELF INTERSECTION PREVENTION
    if self_intersection == 1
        iteraciones = 0;
         while iteraciones < 1000
            [~,badNodes,~,~] = detectSelfIntersectionNodesFromPairs(P_rough,tri,candidatePairs0);
            fprintf('Self Intersection detected %d \n',length(badNodes))
            iteraciones = iteraciones + 1;
            if ~isempty(badNodes) 
                XX = zeros(length(P),1);
                XX(badNodes) = 1;
                W20 = randn(length(badNodes),1);
                % DEFINE SHITE NOISE OVER THE SURFACE =========================
                if noise == 1
                    mlump = sum(M,2);
                    U20 = W20 ./ sqrt(mlump(XX==1));
                end 
                % CRANK NICOLSON ==============================================
                A11 = A(XX==0,XX==0);
                A12 = A(XX==0,XX==1);
                A21 = A(XX==1,XX==0);
                A22 = A(XX==1,XX==1);
    
                B11 = B(XX==0,XX==0);
                B12 = B(XX==0,XX==1);
                B21 = B(XX==1,XX==0);
                B22 = B(XX==1,XX==1);
       
                C = (B21 - A21)/(A11 - B11);
                U2n = U20;
                for n = 1:nsteps
                    U2 = (A22- C*A12)\(B22*U2n-C*B12*U2n);
                    U2n = U2;
                end 
                U(XX==1) = U2n;
                % TRANSPORT MAP ===============================================
                w = mlump / sum(mlump);
                mu = sum(w .* U);
                sig = sqrt(sum(w .* (U - mu).^2));
                U = (U - mu) / sig;
                Fnx = normcdf(U, 0, 1); 
                Fbx = betainv(Fnx, a, a); 
                RFx = 2 * epsilon * Fbx - epsilon;
                P_rough = P + (RFx.*normal_vec);
            else 
                fprintf('No self intersections detected \n')
                break;
            end
         end 
    else
        
    end

end 

function [P_rough,U] = generate_RF_heatequation_surface2(P,tri,normal_vec,noise,nsteps,M,A,A_solver,B,epsilon,a,self_intersection,candidatePairs0)
    % DIFFERENT DEFINITIONS OF WHITE NOISE 
    W = randn(size(P,1),1);
    if noise == 1
        mlump = sum(M,2);
        % Define white noise 
        U0 = W ./ sqrt(mlump);
    end 
    % if noise == 2
    %     U0 = M\(W*2*A);
    % end 
    % 
    % if noise == 3
    %     U0 = M\(W*A);
    % end 
    %% CRANK NICHOLSON
    Un = U0;
    for n = 1:nsteps
        U = A_solver\(B*Un);
        Un=U;
    end 
    U = Un;
    %% TRANSPORT MAP
    % TR = triangulation(tri,P);
    % N = vertexNormal(TR);
    % normal_vec = N ./ vecnorm(N,2,2);
    w = mlump / sum(mlump);
    mu = sum(w .* Un);
    sig = sqrt(sum(w .* (Un - mu).^2));
    Un = (Un - mu) / sig;
    Fnx = normcdf(Un, 0, 1); 
    Fbx = betainv(Fnx, a, a); 
    RFx = 2 * epsilon * Fbx - epsilon;
    P_rough = P + (RFx.*normal_vec);
    %% ADD SELF INTERSECTION PREVENTION
    if self_intersection == 1
        iteraciones = 0;
         while iteraciones < 1000
            [~,badNodes,~,~] = detectSelfIntersectionNodesFromPairs(P_rough,tri,candidatePairs0);
            fprintf('Self Intersection detected %d \n',length(badNodes))
            if badNodes == 0
                iteraciones = 1001;
            end
            iteraciones = iteraciones + 1;
            if ~isempty(badNodes) 
                XX = zeros(length(P),1);
                XX(badNodes) = 1;
            %     W20 = randn(length(badNodes),1);
            %     % DEFINE SHITE NOISE OVER THE SURFACE =========================
            %     if noise == 1
            %         mlump = sum(M,2);
            %         U20 = W20 ./ sqrt(mlump(XX==1));
            %     end 
            %     % CRANK NICOLSON ==============================================
            %     A11 = A(XX==0,XX==0);
            %     A12 = A(XX==0,XX==1);
            %     A21 = A(XX==1,XX==0);
            %     A22 = A(XX==1,XX==1);
            % 
            %     B11 = B(XX==0,XX==0);
            %     B12 = B(XX==0,XX==1);
            %     B21 = B(XX==1,XX==0);
            %     B22 = B(XX==1,XX==1);
            % 
            %     C = (B21 - A21)/(A11 - B11);
            %     U2n = U20;
            %     for n = 1:nsteps
            %         U2 = (A22- C*A12)\(B22*U2n-C*B12*U2n);
            %         U2n = U2;
            %     end 
            %     U(XX==1) = U2n;
            %     % TRANSPORT MAP ===============================================
            %     w = mlump / sum(mlump);
            %     mu = sum(w .* U);
            %     sig = sqrt(sum(w .* (U - mu).^2));
            %     U = (U - mu) / sig;
            %     Fnx = normcdf(U, 0, 1); 
            %     Fbx = betainv(Fnx, a, a); 
            %     RFx = 2 * epsilon * Fbx - epsilon;
            %     P_rough = P + (RFx.*normal_vec);
            RFx(XX==1) = RFx(XX==1) -0.1*RFx(XX==1);
            P_rough = P + (RFx.*normal_vec);
            else 
                fprintf('No self intersections detected')
                break;
            end 
         end 
    else
        
    end

end 

function [K,M] = FEM_assembly(P,tri)
    N = size(P,1); %get number of points 
    NT = size(tri,1); %get number of triangles
    I = zeros(9*NT,1);
    J = zeros(9*NT,1);
    Kvals = zeros(9*NT,1); %Complete K matrix values
    Mvals = zeros(9*NT,1); %Complete K matrix values
    
    idx = 1;
    for t = 1:NT %for each triangle in the surface 
        
        nodes = tri(t,:); % connectivity nodes     
        x1 = P(nodes(1),:); % node 1
        x2 = P(nodes(2),:);
        x3 = P(nodes(3),:);

        n = cross(x2-x1, x3-x1);
        A = 0.5*norm(n);
        grad1 = cross(n, x3-x2) / norm(n)^2;
        grad2 = cross(n, x1-x3) / norm(n)^2;
        grad3 = cross(n, x2-x1) / norm(n)^2;
        grads = [grad1; grad2; grad3];
        Kloc = zeros(3,3); %K local (3 points)
        for i = 1:3
            for j = 1:3
                Kloc(i,j) = A * dot(grads(i,:), grads(j,:));
            end
        end
        
        Mloc = (A/12) * [2 1 1; 1 2 1; 1 1 2];  %M local (3 points)
        
        for i = 1:3
            for j = 1:3
                
                I(idx) = nodes(i);
                J(idx) = nodes(j);
                Kvals(idx) = Kloc(i,j);
                Mvals(idx) = Mloc(i,j);
                idx = idx + 1;
            end
        end
    end
    % Make it sparse 
    K = sparse(I,J,Kvals,N,N);
    M = sparse(I,J,Mvals,N,N);
end 

function N = compute_vertex_normals(P,tri)
    Nnodes = size(P,1);
    N = zeros(Nnodes,3);
    for t = 1:size(tri,1)
        nodes = tri(t,:);
        x1 = P(nodes(1),:);
        x2 = P(nodes(2),:);
        x3 = P(nodes(3),:);
        % Non-unit face normal. Its magnitude is 2*area.
        nf = cross(x2-x1, x3-x1);
        % Add area-weighted face normal to each vertex
        N(nodes(1),:) = N(nodes(1),:) + nf;
        N(nodes(2),:) = N(nodes(2),:) + nf;
        N(nodes(3),:) = N(nodes(3),:) + nf;
    end
    % Normalize vertex normals
    N = N ./ vecnorm(N,2,2);
end

%% SELF INTERSECTION 
function candidatePairs = buildSweptAABBCandidatePairs(P,tri,normal_vec,epsilon,tol,leafSize)

if nargin < 5 || isempty(tol)
    tol = 1e-12;
end

if nargin < 6 || isempty(leafSize)
    leafSize = 32;
end

[bmin,bmax] = sweptTriangleAABB(P,tri,normal_vec,epsilon,tol);
tree = buildAABBTree(bmin,bmax,(1:size(tri,1)).',leafSize);
candidatePairs = queryAABBTreeSelfPairs(tree,bmin,bmax,tri,tol);

end

function [bmin,bmax] = sweptTriangleAABB(P,tri,normal_vec,epsilon,tol)

p1m = P(tri(:,1),:) - epsilon*normal_vec(tri(:,1),:);
p1p = P(tri(:,1),:) + epsilon*normal_vec(tri(:,1),:);
p2m = P(tri(:,2),:) - epsilon*normal_vec(tri(:,2),:);
p2p = P(tri(:,2),:) + epsilon*normal_vec(tri(:,2),:);
p3m = P(tri(:,3),:) - epsilon*normal_vec(tri(:,3),:);
p3p = P(tri(:,3),:) + epsilon*normal_vec(tri(:,3),:);

allPts = cat(3,p1m,p1p,p2m,p2p,p3m,p3p);
bmin = min(allPts,[],3) - tol;
bmax = max(allPts,[],3) + tol;

end

function node = buildAABBTree(bmin,bmax,faces,leafSize)

node.bmin = min(bmin(faces,:),[],1);
node.bmax = max(bmax(faces,:),[],1);
node.faces = [];
node.left = [];
node.right = [];

if numel(faces) <= leafSize
    node.faces = faces(:);
    return
end

centers = 0.5*(bmin(faces,:) + bmax(faces,:));
extent = max(centers,[],1) - min(centers,[],1);
[~,axisId] = max(extent);
[~,order] = sort(centers(:,axisId));
faces = faces(order);

mid = floor(numel(faces)/2);
node.left = buildAABBTree(bmin,bmax,faces(1:mid),leafSize);
node.right = buildAABBTree(bmin,bmax,faces(mid+1:end),leafSize);

end

function candidatePairs = queryAABBTreeSelfPairs(tree,bmin,bmax,tri,tol)

candidatePairs = queryNodeSelf(tree,bmin,bmax,tri,tol);

if ~isempty(candidatePairs)
    candidatePairs = unique(sort(candidatePairs,2),'rows');
end

end

function pairs = queryNodeSelf(node,bmin,bmax,tri,tol)

if isLeafNode(node)
    pairs = leafPairs(node.faces,node.faces,bmin,bmax,tri,tol,true);
    return
end

pairsLeft = queryNodeSelf(node.left,bmin,bmax,tri,tol);
pairsMiddle = queryNodePair(node.left,node.right,bmin,bmax,tri,tol);
pairsRight = queryNodeSelf(node.right,bmin,bmax,tri,tol);
pairs = [pairsLeft; pairsMiddle; pairsRight];

end

function pairs = queryNodePair(nodeA,nodeB,bmin,bmax,tri,tol)

if ~boxesOverlap(nodeA.bmin,nodeA.bmax,nodeB.bmin,nodeB.bmax,tol)
    pairs = zeros(0,2);
    return
end

if isLeafNode(nodeA) && isLeafNode(nodeB)
    pairs = leafPairs(nodeA.faces,nodeB.faces,bmin,bmax,tri,tol,false);
    return
end

if isLeafNode(nodeB) || (~isLeafNode(nodeA) && boxVolume(nodeA) >= boxVolume(nodeB))
    pairs = [queryNodePair(nodeA.left,nodeB,bmin,bmax,tri,tol); ...
             queryNodePair(nodeA.right,nodeB,bmin,bmax,tri,tol)];
else
    pairs = [queryNodePair(nodeA,nodeB.left,bmin,bmax,tri,tol); ...
             queryNodePair(nodeA,nodeB.right,bmin,bmax,tri,tol)];
end

end

function pairs = leafPairs(facesA,facesB,bmin,bmax,tri,tol,sameLeaf)

if sameLeaf
    n = numel(facesA);
    if n < 2
        pairs = zeros(0,2);
        return
    end
    [ia,ib] = find(triu(true(n),1));
    pairs = [facesA(ia), facesA(ib)];
else
    [ia,ib] = ndgrid(1:numel(facesA),1:numel(facesB));
    pairs = [facesA(ia(:)), facesB(ib(:))];
end

if isempty(pairs)
    return
end

keep = false(size(pairs,1),1);
for k = 1:size(pairs,1)
    i = pairs(k,1);
    j = pairs(k,2);

    if i == j || any(ismember(tri(i,:),tri(j,:)))
        continue
    end

    keep(k) = boxesOverlap(bmin(i,:),bmax(i,:),bmin(j,:),bmax(j,:),tol);
end

pairs = pairs(keep,:);

end

function tf = boxesOverlap(bminA,bmaxA,bminB,bmaxB,tol)

tf = all(bmaxA >= bminB - tol) && all(bmaxB >= bminA - tol);

end

function tf = isLeafNode(node)

tf = ~isempty(node.faces);

end

function volume = boxVolume(node)

d = max(node.bmax - node.bmin,0);
volume = d(1)*d(2)*d(3);

end

function [TF,badNodes,badFaces,badPairs] = detectSelfIntersectionNodesFromPairs(P,tri,candidatePairs,tol)

if nargin < 4 || isempty(tol)
    tol = 1e-12;
end

if exist('detectSelfIntersectionsFromPairsMex','file') == 3
    [badNodes,badFaces,badPairs] = detectSelfIntersectionsFromPairsMex(P,tri,candidatePairs,tol);
    TF = ~isempty(badPairs);
    return
end

badPairs = zeros(0,2);
badFaces = [];
badNodes = [];

if isempty(candidatePairs)
    TF = false;
    return
end

p1 = P(tri(:,1),:);
p2 = P(tri(:,2),:);
p3 = P(tri(:,3),:);

bmin = min(cat(3,p1,p2,p3),[],3) - tol;
bmax = max(cat(3,p1,p2,p3),[],3) + tol;

a = candidatePairs(:,1);
b = candidatePairs(:,2);
overlap = all(bmax(a,:) >= bmin(b,:) - tol,2) & ...
          all(bmax(b,:) >= bmin(a,:) - tol,2);
pairs = candidatePairs(overlap,:);

if isempty(pairs)
    TF = false;
    return
end

hit = false(size(pairs,1),1);
for k = 1:size(pairs,1)
    i = pairs(k,1);
    j = pairs(k,2);

    T1 = [p1(i,:); p2(i,:); p3(i,:)];
    T2 = [p1(j,:); p2(j,:); p3(j,:)];

    if triTriIntersectSAT(T1,T2,tol)
        hit(k) = true;
    end
end

badPairs = pairs(hit,:);
badFaces = unique(badPairs(:));
badNodes = unique(tri(badFaces,:));
TF = ~isempty(badPairs);

end

function [TF,badNodes,badFaces,badPairs] = detectSelfIntersectionNodesFast(P,tri,tol,cellSize)

if nargin < 3 || isempty(tol)
    tol = 1e-12;
end

F = size(tri,1);
badPairs = zeros(0,2);
badFaces = [];
badNodes = [];

p1 = P(tri(:,1),:);
p2 = P(tri(:,2),:);
p3 = P(tri(:,3),:);

bmin = min(cat(3,p1,p2,p3),[],3) - tol;
bmax = max(cat(3,p1,p2,p3),[],3) + tol;

if nargin < 4 || isempty(cellSize)
    e = [vecnorm(p2-p1,2,2); vecnorm(p3-p2,2,2); vecnorm(p1-p3,2,2)];
    cellSize = 2*median(e(e > 0));
end

origin = min(bmin,[],1);
cmin = floor((bmin - origin)/cellSize) + 1;
cmax = floor((bmax - origin)/cellSize) + 1;
dims = max(cmax,[],1) + 1;

keys = [];
faceIds = [];

for f = 1:F
    xs = cmin(f,1):cmax(f,1);
    ys = cmin(f,2):cmax(f,2);
    zs = cmin(f,3):cmax(f,3);

    [X,Y,Z] = ndgrid(xs,ys,zs);
    k = X(:) + dims(1)*(Y(:)-1) + dims(1)*dims(2)*(Z(:)-1);

    keys = [keys; k];
    faceIds = [faceIds; repmat(f,numel(k),1)];
end

[keys,ord] = sort(keys);
faceIds = faceIds(ord);

candidatePairs = zeros(0,2);
groupStart = [1; find(diff(keys) ~= 0) + 1];
groupEnd = [groupStart(2:end)-1; numel(keys)];

for g = 1:numel(groupStart)
    faces = unique(faceIds(groupStart(g):groupEnd(g)));
    nf = numel(faces);

    if nf < 2
        continue
    end

    [I,J] = find(triu(true(nf),1));
    candidatePairs = [candidatePairs; faces(I) faces(J)];
end

if isempty(candidatePairs)
    TF = false;
    return
end

candidatePairs = unique(sort(candidatePairs,2),'rows');

hit = false(size(candidatePairs,1),1);

for k = 1:size(candidatePairs,1)
    i = candidatePairs(k,1);
    j = candidatePairs(k,2);

    if any(ismember(tri(i,:),tri(j,:)))
        continue
    end

    if any(bmax(i,:) < bmin(j,:) - tol) || any(bmax(j,:) < bmin(i,:) - tol)
        continue
    end

    T1 = [p1(i,:); p2(i,:); p3(i,:)];
    T2 = [p1(j,:); p2(j,:); p3(j,:)];

    if triTriIntersectSAT(T1,T2,tol)
        hit(k) = true;
    end
end

badPairs = candidatePairs(hit,:);
badFaces = unique(badPairs(:));
badNodes = unique(tri(badFaces,:));
TF = ~isempty(badPairs);

end

function flag = triTriIntersectSAT(T1,T2,tol)

E1 = [T1(2,:)-T1(1,:);
      T1(3,:)-T1(2,:);
      T1(1,:)-T1(3,:)];

E2 = [T2(2,:)-T2(1,:);
      T2(3,:)-T2(2,:);
      T2(1,:)-T2(3,:)];

n1 = cross(E1(1,:),E1(2,:));
n2 = cross(E2(1,:),E2(2,:));

axes = [n1; n2];

% Edge cross-product axes
for a = 1:3
    for b = 1:3
        axes = [axes; cross(E1(a,:),E2(b,:))];
    end
end

% Coplanar in-plane separating axes
if norm(cross(n1,n2)) < tol
    for a = 1:3
        axes = [axes; cross(n1,E1(a,:))];
        axes = [axes; cross(n2,E2(a,:))];
    end
end

for k = 1:size(axes,1)

    ax = axes(k,:);
    nax = norm(ax);

    if nax < tol
        continue
    end

    ax = ax/nax;

    pT1 = T1*ax.';
    pT2 = T2*ax.';

    if max(pT1) < min(pT2) - tol || max(pT2) < min(pT1) - tol
        flag = false;
        return
    end
end

flag = true;

end

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