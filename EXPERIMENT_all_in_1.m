% SOLVE HEAT EQUATION OVER A DEFINED SURFACE MESH
clear 
clc
close all
root = 'C:\Users\royes\OneDrive\Escritorio\Lab\rough_surfaces2\CEM-Codes';
addpath(root)
format longe

% root = 'C:\Users\royes\OneDrive\Escritorio\Lab\rough_surfaces2\CEM-Codes';
% addpath(root)
% format longe

%% GENERATE THE ROUGH SURFACES 

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
        P = TR.Points; P = P / 1000;            
        tri = TR.ConnectivityList;  
    case 'sphere'
        [ P, tri ] = generateSphereMesh( 2, 'oct' ); P = P';
        TR = triangulation(tri,P);
    otherwise 
        fprintf('Wrong mesh name!!!!!!!!!! \n')
end 

if refine == 1
    [P, tri] = barycentric_refine(P, tri); % BARYCENTRIC REFINE 
    TR = triangulation(tri, P);   % new refined triangulation
end 
% GET MESH DISCRETIZATION QUANTITIES
E = edges(TR);
edgeLength = vecnorm(P(E(:,1),:) - P(E(:,2),:),2,2);
h_min = min(edgeLength);
h_mean= mean(edgeLength);
normal_vec = compute_vertex_normals(P,tri);

% Constants

lambda = 5*h_mean; % incident wavelength used to define the normalized sweeps
correlation_length_over_lambda = [0.1, 1, 10];
epsilon_over_lambda = [0.01, 0.1, 1];
correlation_length_values = correlation_length_over_lambda*lambda;
epsilon_values = epsilon_over_lambda*lambda;

h = median(edgeLength);
r_target = 0.007;
t_final = 1; % final time
a = 2; % hyperparameter beta distribution
self_intersection = 1; %1: for self intersection check 0: for no self intersection check
samples = 2;
UQps = min(100,samples);
noise = 1;

%% OBTAIN THE SCATTERING UNCERTAINTY QUANTIFICATION FOR EACH correlation_length/lambda and epsilon/lambda parameters

f2p = tri;
k0 = 2*pi/lambda;
eta0 = 377;
E0 = 1;
intord = 1; %quadrature order
theta = linspace(0,pi,181);
nCorr = length(correlation_length_over_lambda);
nEps = length(epsilon_over_lambda);

% Current curve statistics: by default, use a fixed nominal-surface path
% between the minimum-x and maximum-x mesh vertices.
[~,curveStartId] = min(P(:,1));
[~,curveEndId] = max(P(:,1));
curvePoints = P([curveStartId curveEndId],:);
zPlanePenalty = 0;
[sCurve,JCurvePathVertices] = surfaceCurvePath(P,f2p,curvePoints,zPlanePenalty);
nCurve = length(sCurve);

SIGMA = zeros(samples,length(theta),nCorr,nEps);
SIGMA_mean = zeros(length(theta),nCorr,nEps);
SIGMA_std = zeros(length(theta),nCorr,nEps);
JCurveMag = zeros(samples,nCurve,nCorr,nEps);
JCurveMag_mean = zeros(nCurve,nCorr,nEps);
JCurveMag_std = zeros(nCurve,nCorr,nEps);
PS = zeros(length(P),3,UQps,nCorr,nEps);
C = zeros(length(P),length(P),nCorr,nEps);

[K,M] = FEM_assembly(P,tri); % FEM matrices for solving heat equation
normal_vec = compute_vertex_normals(P,tri);

for ie = 1:nEps
    epsilon = epsilon_values(ie); % maximum perturbation height after transport map

    % Precompute the face pairs that could intersect for any normal
    % displacement RFx in [-epsilon, epsilon].
    if self_intersection == 1
        candidatePairs0 = buildSweptAABBCandidatePairs(P,tri,normal_vec,epsilon);
        fprintf('epsilon/lambda = %.1f: precomputed %d possible self-intersection face pairs \n', ...
                epsilon_over_lambda(ie),size(candidatePairs0,1))
    else
        candidatePairs0 = zeros(0,2);
    end

    for ic = 1:nCorr
        correlation_length = correlation_length_values(ic);
        alpha = (correlation_length/2)^2*(1/t_final);
        dt = r_target*h^2/alpha;
        nsteps = round(t_final/dt);
        A = M + (dt/2)*alpha*K;
        B = M - (dt/2)*alpha*K;
        A_solver = decomposition(A,'chol');
        U = zeros(length(P),samples);

        fprintf('Running ell/lambda = %.1f, epsilon/lambda = %.1f \n', ...
                correlation_length_over_lambda(ic),epsilon_over_lambda(ie))

        for s = 1:samples
            fprintf('  Rough surface/scattering sample %d of %d \n',s,samples)
            [p,U(:,s)] = generate_RF_heatequation_surface2(P,tri,normal_vec,noise,nsteps,M,A,A_solver,B,epsilon,a,self_intersection,candidatePairs0);

            if s <= UQps
                PS(:,:,s,ic,ie) = p;
            end

            [SIGMA(s,:,ic,ie),JCurveMag(s,:,ic,ie)] = compute_bistatic_sigma(p,f2p,intord,k0,eta0,E0,theta,JCurvePathVertices);
        end

        C(:,:,ic,ie) = cov(U');
        SIGMA_mean(:,ic,ie) = mean(SIGMA(:,:,ic,ie),1).';
        SIGMA_std(:,ic,ie) = std(SIGMA(:,:,ic,ie),0,1).';
        JCurveMag_mean(:,ic,ie) = mean(JCurveMag(:,:,ic,ie),1).';
        JCurveMag_std(:,ic,ie) = std(JCurveMag(:,:,ic,ie),0,1).';
    end
end

%%
figure('Color','w','Position',[100 100 1400 1000])
tiledlayout(nEps,nCorr,'TileSpacing','compact','Padding','compact')
theta_deg = 180*theta/pi;

for ie = 1:nEps
    for ic = 1:nCorr
        nexttile
        hold on

        SIGMA2 = SIGMA(:,:,ic,ie);
        sigma_mean = mean(SIGMA2,1);
        sigma_std = std(SIGMA2,0,1);
        N = size(SIGMA2,1);
        CI95 = 1.96*sigma_std/sqrt(N);

        sigma_upper = sigma_mean + CI95;
        sigma_lower = sigma_mean - CI95;
        %sigma_lower = max(sigma_mean - CI95, realmin);

        x2 = [theta_deg fliplr(theta_deg)];
        inBetween = [10*log10(sigma_upper), ...
                     fliplr(10*log10(sigma_lower))];

        fill(x2,inBetween,[0.8 0.8 1], ...
            'EdgeColor','none', ...
            'FaceAlpha',0.6)

        plot(theta_deg,10*log10(sigma_mean), ...
            'r--','LineWidth',2.2)

        plot(theta_deg,10*log10(sigma_upper), ...
            'k','LineWidth',1.6)

        plot(theta_deg,10*log10(sigma_lower), ...
            'k','LineWidth',1.6)

        grid on
        box on
        xlim([0 180])
        ylim([-40 40])

        title(sprintf('\\textbf{$\\varepsilon = %.2f\\lambda$ , $\\ell = %.2f\\lambda$}', ...
              epsilon_over_lambda(ie),correlation_length_over_lambda(ic)), ...
              'Interpreter','latex', ...
              'FontSize',18)

        ax = gca;
        ax.FontSize = 14;
        ax.FontWeight = 'bold';
        ax.LineWidth = 1.3;
        ax.TickLabelInterpreter = 'latex';

        if ic == 1 && ie == 2
            ylabel('\textbf{Bistatic RCS [dBsm]}', ...
                'Interpreter','latex', ...
                'FontSize',16)
        end

        if ie == nEps
            xlabel('\textbf{Observation angle $\theta$ [deg]}', ...
                'Interpreter','latex', ...
                'FontSize',16)
        end

        if ie == 1 && ic == 1
            lgd = legend('95\% CI','Monte Carlo mean', ...
                'Interpreter','latex', ...
                'Location','northeast');
            lgd.FontSize = 12;
            lgd.Box = 'on';
        end
    end
end

figure('Color','w','Position',[100 100 1400 1000])
tiledlayout(nEps,nCorr,'TileSpacing','compact','Padding','compact')

JCurveMag_global_min = inf;
JCurveMag_global_max = -inf;
for ie = 1:nEps
    for ic = 1:nCorr
        JCurveMag2 = JCurveMag(:,:,ic,ie);
        JCurveMag_case_mean = mean(JCurveMag2,1);
        JCurveMag_case_std = std(JCurveMag2,0,1);
        N = size(JCurveMag2,1);
        CI95 = 1.96*JCurveMag_case_std/sqrt(N);

        JCurveMag_upper = JCurveMag_case_mean + CI95;
        JCurveMag_lower = max(JCurveMag_case_mean - CI95,0);

        JCurveMag_global_min = min(JCurveMag_global_min,min(JCurveMag_lower));
        JCurveMag_global_max = max(JCurveMag_global_max,max(JCurveMag_upper));
    end
end

JCurveMag_ylim = [JCurveMag_global_min JCurveMag_global_max];
if JCurveMag_ylim(1) == JCurveMag_ylim(2)
    JCurveMag_ylim = JCurveMag_ylim + [-1 1]*max(abs(JCurveMag_ylim(1)),1)*0.05;
else
    yPad = 0.05*diff(JCurveMag_ylim);
    JCurveMag_ylim = JCurveMag_ylim + [-yPad yPad];
end

for ie = 1:nEps
    for ic = 1:nCorr
        nexttile
        hold on

        JCurveMag2 = JCurveMag(:,:,ic,ie)./max(JCurveMag(:,:,ic,ie));
        JCurveMag_case_mean = mean(JCurveMag2,1);
        JCurveMag_case_std = std(JCurveMag2,0,1);
        N = size(JCurveMag2,1);
        CI95 = 1.96*JCurveMag_case_std/sqrt(N);

        JCurveMag_upper = JCurveMag_case_mean + CI95;
        JCurveMag_lower = max(JCurveMag_case_mean - CI95,0);

        x2 = [sCurve.' fliplr(sCurve.')];
        inBetween = [JCurveMag_upper fliplr(JCurveMag_lower)];

        fill(x2,inBetween,[0.8 0.8 1], ...
            'EdgeColor','none', ...
            'FaceAlpha',0.6)

        plot(sCurve,JCurveMag_case_mean, ...
            'r--','LineWidth',2.2)

        plot(sCurve,JCurveMag_upper, ...
            'k','LineWidth',1.6)

        plot(sCurve,JCurveMag_lower, ...
            'k','LineWidth',1.6)

        grid on
        box on
        xlim([sCurve(1) sCurve(end)])
        %ylim(JCurveMag_ylim)

        title(sprintf('\\textbf{$\\varepsilon = %.2f\\lambda$ , $\\ell = %.2f\\lambda$}', ...
              epsilon_over_lambda(ie),correlation_length_over_lambda(ic)), ...
              'Interpreter','latex', ...
              'FontSize',18)

        ax = gca;
        ax.FontSize = 14;
        ax.FontWeight = 'bold';
        ax.LineWidth = 1.3;
        ax.TickLabelInterpreter = 'latex';

        if ic == 1 && ie == 2
            ylabel('\textbf{Current magnitude $|J|$}', ...
                'Interpreter','latex', ...
                'FontSize',16)
        end

        if ie == nEps
            xlabel('\textbf{Distance along surface curve}', ...
                'Interpreter','latex', ...
                'FontSize',16)
        end

        if ie == 1 && ic == 1
            lgd = legend('95\% CI','Monte Carlo mean', ...
                'Interpreter','latex', ...
                'Location','northeast');
            lgd.FontSize = 12;
            lgd.Box = 'on';
        end
    end
end

filename = sprintf('RCS_bell_parameter_sweep_Rough_mesh_%s.mat', mesh);
save(filename, 'SIGMA','SIGMA_mean','SIGMA_std','JCurveMag', ...
               'JCurveMag_mean','JCurveMag_std','sCurve','JCurvePathVertices', ...
               'curvePoints','theta','lambda', ...
               'correlation_length_over_lambda','epsilon_over_lambda', ...
               'correlation_length_values','epsilon_values','PS','C','-v7.3')


%% FUNCTIONS 

function [sigma,JCurveMag] = compute_bistatic_sigma(p,f2p,intord,k0,eta0,E0,theta,pathVertices)
    Eix = @(x,y,z) E0*exp(1i*k0*z);
    Eiy = @(x,y,z) zeros(size(x));
    Eiz = @(x,y,z) zeros(size(x));

    Hix = @(x,y,z) zeros(size(x));
    Hiy = @(x,y,z) exp(1i*k0*z)/eta0;
    Hiz = @(x,y,z) zeros(size(x));

    [xcfie,~,m] = CFIE_densemat(p,f2p,intord,0.5,0.5,k0,eta0, ...
                                Eix,Eiy,Eiz,Hix,Hiy,Hiz);

    J = zeros(m.nf,3);

    for iface = 1:m.nf
        coeff = 1/(2*m.Surf(iface)) * ...
                xcfie(abs(m.f2ed(iface,:))) .* sign(m.f2ed(iface,:)).';

        pcen_face = sum(m.p(m.f2p(iface,:),:),1)/3;

        J(iface,:) = coeff(1)*(pcen_face - m.p(m.f2p(iface,1),:)) + ...
                     coeff(2)*(pcen_face - m.p(m.f2p(iface,2),:)) + ...
                     coeff(3)*(pcen_face - m.p(m.f2p(iface,3),:));
    end

    pcen = (m.p(m.f2p(:,1),:) + m.p(m.f2p(:,2),:) + m.p(m.f2p(:,3),:))/3;
    phi0 = 0;
    sigma = zeros(size(theta));

    for it = 1:length(theta)
        rhat = [sin(theta(it))*cos(phi0), ...
                sin(theta(it))*sin(phi0), ...
                cos(theta(it))];

        phase = exp(1i*k0*(pcen*rhat.'));
        F = sum(J .* (m.Surf(:).*phase), 1);
        Efar_vec = cross(rhat, cross(rhat, F));
        sigma(it) = (k0^2 * eta0^2 / (4*pi*E0^2)) * sum(abs(Efar_vec).^2);
    end

    JCurveMag = currentMagnitudeOnPath(J,m.f2p,size(p,1),pathVertices);
end

function [sCurve,pathVertices] = surfaceCurvePath(p,f2p,curvePoints,zPlanePenalty)
    if size(curvePoints,1) < 2
        error('curvePoints must contain at least a start and an end point.')
    end

    if nargin < 4
        zPlanePenalty = 0;
    end

    curveVertexIds = nearestMeshVertices(p,curvePoints);
    pathVertices = curveVertexIds(1);

    nVertices = size(p,1);
    E = [f2p(:,[1 2]); f2p(:,[2 3]); f2p(:,[3 1])];
    E = unique(sort(E,2),'rows');
    edgeLength = vecnorm(p(E(:,1),:) - p(E(:,2),:),2,2);
    zScale = max(abs(p(:,3)));
    if zScale == 0
        zScale = 1;
    end
    edgeZ = 0.5*(abs(p(E(:,1),3)) + abs(p(E(:,2),3)));
    edgeCost = edgeLength.*(1 + zPlanePenalty*edgeZ/zScale);
    G = graph(E(:,1),E(:,2),edgeCost,nVertices);

    for ip = 1:numel(curveVertexIds)-1
        segmentVertices = shortestpath(G,curveVertexIds(ip),curveVertexIds(ip+1));
        if isempty(segmentVertices)
            error('No connected path found between curve point %d and %d.',ip,ip+1)
        end
        pathVertices = [pathVertices segmentVertices(2:end)]; %#ok<AGROW>
    end

    pathPoints = p(pathVertices,:);
    ds = vecnorm(diff(pathPoints,1,1),2,2);
    sCurve = [0; cumsum(ds)];
end

function JCurveMag = currentMagnitudeOnPath(J,f2p,nVertices,pathVertices)
    JmagFace = sqrt(sum(abs(J).^2,2));
    vertexJMag = accumarray(f2p(:),repmat(JmagFace,3,1),[nVertices 1],@mean,NaN);
    JCurveMag = reshape(vertexJMag(pathVertices),1,[]);
end

function vertexIds = nearestMeshVertices(p,queryPoints)
    vertexIds = zeros(size(queryPoints,1),1);
    for iq = 1:size(queryPoints,1)
        [~,vertexIds(iq)] = min(vecnorm(p - queryPoints(iq,:),2,2));
    end
end

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

% SELF INTERSECTION

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
