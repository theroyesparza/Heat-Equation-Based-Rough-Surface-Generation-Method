clear 
clc
close all
root = 'C:\Users\royes\OneDrive\Escritorio\Lab\rough_surfaces2\CEM-Codes';
addpath(fullfile(root,'peccodes'))
addpath(fullfile(root,'basisfunc'))
addpath(fullfile(root,'meshes'))
addpath(fullfile(root,'singularintegrals'))
savepath
format longe

%% LOAD MESH AND VISUALIZE 

% 1) PERFECT SPHERE FOR VALIDATION
% [ p, f2p ] = generateSphereMesh(4, 'oct' ); p = p';
% p = p + 1;

% 2) NASA ALMOND
% nasa = 'nasaAlmond.STL';
% TR = stlread(nasa);
% p = TR.Points / 1000;     % if STL is in mm 
% f2p = TR.ConnectivityList; 

% messi = 'cabeza lowpoly messi.stl';
% TR = stlread(messi);
% p = TR.Points; p = p/400;
% f2p = TR.ConnectivityList;  

%3) ROUGH MESH
load rougsurfaceTEST.mat
p = P+1;
f2p = tria;

% 2) Rough surface mesh obtained using HEAT_Surface.m 
% load 99_samples_Rough_mesh_sphere.mat
% p = PS(:,:,1,1);
% f2p = tri;

figure()
trisurf(f2p, p(:,1), p(:,2), p(:,3), ...
        'EdgeColor','k');
axis equal;
colorbar;

%================
% E = edges(TR);
% edgeLength = vecnorm(P(E(:,1),:) - P(E(:,2),:),2,2);
% h_min = min(edgeLength);
% h_mean= mean(edgeLength);
% normal_vec = compute_vertex_normals(P,tri);

%% MoM Implementation
lambda = 100; %dimension-less
%lambda = 11.8028/1.19; %from data in GHZ --> in
%lambda = 0.299792458/1.19; % from data in GHZ --> m
% lambda = 0.1*epsilon; % from data in GHZ --> m
k0=2*pi/lambda;
eta0=377;
E0 = 1;
intord=1; %quadrature order
%incident_fields
Eix=@(x,y,z) E0*exp(1i*k0*z); % all the code considers the opposite convention 
Eiy=@(x,y,z) zeros(size(x));
Eiz=@(x,y,z) zeros(size(x));

Hix=@(x,y,z) zeros(size(x));
Hiy=@(x,y,z) exp(1i*k0*z)/eta0;
Hiz=@(x,y,z) zeros(size(x));


[xcfie,TR,m]=...
   CFIE_densemat(p,f2p,intord,0.5,0.5,k0,eta0,Eix,Eiy,Eiz,Hix,Hiy,Hiz);

% [xefie,TR,m]=...
%     EFIE_densemat(p,f2p,intord,k0,eta0,Eix,Eiy,Eiz);
% [xmfie,TR,m]=...
%     MFIE_densemat(p,f2p,intord,k0,eta0,Hix,Hiy,Hiz);

%% CURRENTS CALCULATION AT EACH TRIG. FACE 
% CURRENTS CFIE ===========================================================
J=zeros([m.nf 3]);
for i=1:m.nf
    coeff=1/(2*m.Surf(i))*xcfie(abs(m.f2ed(i,:))).*sign(m.f2ed(i,:)).';
    pcen=sum(m.p(m.f2p(i,:),:),1)/3;
J(i,:)=coeff(1)*(pcen-m.p(m.f2p(i,1),:))+...
       coeff(2)*(pcen-m.p(m.f2p(i,2),:))+...
       coeff(3)*(pcen-m.p(m.f2p(i,3),:));
end

Jmag=sqrt(sum(conj(J).*J,2));

hJsurf = figure();
trisurf(f2p,p(:,1),p(:,2),p(:,3),abs(Jmag(:,1)))
axis equal tight
colorbar
title('$|J|$ over the mesh','Interpreter','latex')

%% CURRENT ALONG A CURVE ON THE SURFACE ==================================
% Define the curve with points on/near the surface. Each point is snapped to
% the nearest mesh vertex, then the curve follows mesh edges from f2p.
meshCenter = 0.5*(min(p,[],1) + max(p,[],1));
curvePoints = [
    min(p(:,1)) meshCenter(2) meshCenter(3)
    max(p(:,1)) meshCenter(2) meshCenter(3)
];
zPlanePenalty = 10;

[sCurve,JCurveMag,JCurve,pathVertices] = plotCurrentOnSurfaceCurve(p,f2p,J,curvePoints,zPlanePenalty);

% CURRENTS EFIE ===========================================================
% Jefie=zeros([m.nf 3]);
% for i=1:m.nf
%     coeff=1/(2*m.Surf(i))*xefie(abs(m.f2ed(i,:))).*sign(m.f2ed(i,:)).';
%     pcen=sum(m.p(m.f2p(i,:),:),1)/3;
% Jefie(i,:)=coeff(1)*(pcen-m.p(m.f2p(i,1),:))+...
%        coeff(2)*(pcen-m.p(m.f2p(i,2),:))+...
%        coeff(3)*(pcen-m.p(m.f2p(i,3),:));
% end
% Jmage=sqrt(sum(conj(Jefie).*Jefie,2));
% subplot(2,2,2),  trisurf(f2p,p(:,1),p(:,2),p(:,3),abs(Jmage(:,1)),'edgealpha',0)


% CURRENTS MFIE ===========================================================
% Jmfie=zeros([m.nf 3]);
% for i=1:m.nf
%     coeff=1/(2*m.Surf(i))*xmfie(abs(m.f2ed(i,:))).*sign(m.f2ed(i,:)).';
%     pcen=sum(m.p(m.f2p(i,:),:),1)/3;
% Jmfie(i,:)=coeff(1)*(pcen-m.p(m.f2p(i,1),:))+...
%        coeff(2)*(pcen-m.p(m.f2p(i,2),:))+...
%        coeff(3)*(pcen-m.p(m.f2p(i,3),:));
% end
% Jmagm=sqrt(sum(conj(Jmfie).*Jmfie,2));
% subplot(2,2,3),  trisurf(f2p,p(:,1),p(:,2),p(:,3),abs(Jmagm(:,1)),'edgealpha',0)


%% VALIDATION FOR SPHERE ===================================================

% MIE series current 
pcen=(m.p(m.f2p(:,1),:)+m.p(m.f2p(:,2),:)+m.p(m.f2p(:,3),:))/3;
 addpath('../')
  [Jth, Jphi,Jx,Jy,Jz] = mie_currents_pec(1, k0, eta0, 10^-5, ...
     pcen(:,1),pcen(:,2),pcen(:,3));

  Jmagmie=sqrt(abs(Jth).^2+abs(Jphi).^2);
 % subplot(2,2,4), trisurf(f2p,p(:,1),p(:,2),p(:,3),Jmagmie,'edgealpha',0)
 % colorbar

 Jmie=[Jx(:) Jy(:) Jz(:)];
 norm(Jmie(:)-conj(J(:)))/norm(Jmie(:)) 
 % norm(Jmie(:)-conj(Jefie(:)))/norm(Jmie(:))
 % norm(Jmie(:)-conj(Jmfie(:)))/norm(Jmie(:))

%% RCS validation MIE CURRENTS
lambda = 1;
a = 1.0*lambda;
phi = 0;

theta_obs = linspace(1e-5, pi-1e-5, 50);
theta_1 = pi - theta_obs; % convention flip
sigma_theory = PEC_sphere_book_RCS(theta_1, phi, a, lambda);

%% BISTATIC SCATTERING CROSS SECTION (ORIGINAL)

% Triangle centroids
pcen = (m.p(m.f2p(:,1),:) + m.p(m.f2p(:,2),:) + m.p(m.f2p(:,3),:))/3;
% Angles
theta = linspace(0,pi,181);      % polar angle
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

%% VALIDATION ==============================================================
% figure()
% plot(theta*180/pi, 10*log10(sigma),'r', 'LineWidth', 1.5)
% hold on 
% plot((pi-theta_1*180/pi)+180, 10*log10(sigma_theory),'ko', 'LineWidth', 1.5)
% grid on
% xlim([0 180])
% xlabel('$\theta$ observation angle [deg]','Interpreter','latex','FontSize',15)
% ylabel('$\sigma$($\theta$,$\phi=0$) [dBsm]','Interpreter','latex','FontSize',15)
% legend('MoM','Analytical','Interpreter','latex','FontSize',15)

% VISUALIZATION OTHER GEOMETRIES 
figure()
plot(theta*180/pi, 10*log10(sigma),'r', 'LineWidth', 1.5)
grid on
xlabel('$\theta$ observation angle [deg]','Interpreter','latex','FontSize',15)
ylabel('$\sigma$($\theta$,$\phi=0$) [dBsm]','Interpreter','latex','FontSize',15)

%% COMPLEX CURRENT 
figure('Color','w')

subplot(2,1,1)
hold on
plot(sCurve,real(JCurve(:,1)),'LineWidth',1.5)
plot(sCurve,real(JCurve(:,2)),'LineWidth',1.5)
plot(sCurve,real(JCurve(:,3)),'LineWidth',1.5)
grid on
xlabel('Distance along surface curve','Interpreter','latex','FontSize',25)
ylabel('$\mathrm{Re}\{\mathbf{J}\}$ [A/m]', ...
    'Interpreter','latex','FontSize',30)
legend('$\mathrm{Re}\{J_x\}$', ...
       '$\mathrm{Re}\{J_y\}$', ...
       '$\mathrm{Re}\{J_z\}$', ...
       'Interpreter','latex','Location','best','FontSize',18)

subplot(2,1,2)
hold on
plot(sCurve,imag(JCurve(:,1)),'LineWidth',1.5)
plot(sCurve,imag(JCurve(:,2)),'LineWidth',1.5)
plot(sCurve,imag(JCurve(:,3)),'LineWidth',1.5)
grid on
xlabel('Distance along surface curve','Interpreter','latex','FontSize',25)
ylabel('$\mathrm{Im}\{\mathbf{J}\}$ [A/m]', ...
    'Interpreter','latex','FontSize',30)
legend('$\mathrm{Im}\{J_x\}$', ...
       '$\mathrm{Im}\{J_y\}$', ...
       '$\mathrm{Im}\{J_z\}$', ...
       'Interpreter','latex','Location','best','FontSize',18)


%% FUNCTIONS 

function sigma = PEC_sphere_book_RCS(theta, phi, a, lambda)

    k = 2*pi/lambda;
    x = k*a;

    theta = theta(:).';
    mu = cos(theta);
    s  = sin(theta);

    Nmax = ceil(x + 4*x^(1/3) + 20);

    Atheta = zeros(size(theta));
    Aphi   = zeros(size(theta));

    for n = 1:Nmax

        alpha_n = (1i)^(-n)*(2*n + 1)/(n*(n + 1));

        % Riccati-Bessel functions
        psi  = riccati_j(n,x);
        xi   = riccati_h2(n,x);

        psip = riccati_j_deriv(n,x);
        xip  = riccati_h2_deriv(n,x);

        bn = -alpha_n * psip/xip;
        cn = -alpha_n * psi/xi;

        % Legendre polynomial and derivative wrt mu = cos(theta)
        [Pn, dPn] = legendreP_and_deriv(n, mu);

        % Book convention: P_n^1(mu) = sin(theta)*P_n'(mu)
        Pn1_over_s = dPn;

        % sin(theta)*d/dmu[P_n^1(mu)]
        sin_dPn1 = mu.*dPn - n*(n + 1)*Pn;

        Atheta = Atheta + (1i)^n * ...
            (bn*sin_dPn1 - cn*Pn1_over_s);

        Aphi = Aphi + (1i)^n * ...
            (bn*Pn1_over_s - cn*sin_dPn1);
    end

    sigma = (lambda^2/pi) * ...
        (cos(phi)^2*abs(Atheta).^2 + sin(phi)^2*abs(Aphi).^2);
end

function y = sphj(n,x)
    y = sqrt(pi/(2*x))*besselj(n+0.5,x);
end

function y = sphh2(n,x)
    y = sqrt(pi/(2*x))*besselh(n+0.5,2,x);
end

function y = riccati_j(n,x)
    y = x*sphj(n,x);
end

function y = riccati_h2(n,x)
    y = x*sphh2(n,x);
end

function y = riccati_j_deriv(n,x)
    y = x*sphj(n-1,x) - n*sphj(n,x);
end

function y = riccati_h2_deriv(n,x)
    y = x*sphh2(n-1,x) - n*sphh2(n,x);
end

function [Pn, dPn] = legendreP_and_deriv(n, mu)

    if n == 0
        Pn = ones(size(mu));
        dPn = zeros(size(mu));
        return;
    end

    Pnm1 = ones(size(mu));    % P_0
    Pn   = mu;                % P_1

    if n == 1
        dPn = ones(size(mu));
        return;
    end

    for m = 2:n
        Pnp1 = ((2*m-1)*mu.*Pn - (m-1)*Pnm1)/m;
        Pnm1 = Pn;
        Pn = Pnp1;
    end

    % Stable derivative formula
    dPn = n*(Pnm1 - mu.*Pn)./(1 - mu.^2);
end

function [sCurve,JCurveMag,JCurve,pathVertices] = plotCurrentOnSurfaceCurve(p,f2p,J,curvePoints,zPlanePenalty)

    if size(curvePoints,1) < 2
        error('curvePoints must contain at least a start and an end point.')
    end

    if nargin < 5
        zPlanePenalty = 0;
    end

    JmagFace = sqrt(sum(abs(J).^2,2));
    nVertices = size(p,1);

    vertexJMag = accumarray(f2p(:),repmat(JmagFace,3,1),[nVertices 1],@mean,NaN);
    vertexJ = zeros(nVertices,3);
    for icomp = 1:3
        vertexJ(:,icomp) = accumarray(f2p(:),repmat(J(:,icomp),3,1),[nVertices 1],@mean,NaN);
    end

    curveVertexIds = nearestMeshVertices(p,curvePoints);
    pathVertices = curveVertexIds(1);

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
    JCurveMag = vertexJMag(pathVertices);
    JCurve = vertexJ(pathVertices,:);

    figure()
    subplot(1,2,1)
    trisurf(f2p,p(:,1),p(:,2),p(:,3),JmagFace, ...
        'EdgeColor','none','FaceAlpha',0.55)
    hold on
    plot3(pathPoints(:,1),pathPoints(:,2),pathPoints(:,3), ...
        'k.-','LineWidth',2,'MarkerSize',14)
    plot3(p(curveVertexIds,1),p(curveVertexIds,2),p(curveVertexIds,3), ...
        'ro','LineWidth',1.5,'MarkerFaceColor','r')
    axis equal tight
    xlabel('$x$','Interpreter','latex')
    ylabel('$y$','Interpreter','latex')
    zlabel('$z$','Interpreter','latex')
    title('Curve on mesh connectivity','Interpreter','latex')
    colorbar
    view(3)

    subplot(1,2,2)
    plot(sCurve,JCurveMag,'k.-','LineWidth',1.5,'MarkerSize',12)
    grid on
    xlabel('Distance along surface curve','Interpreter','latex')
    ylabel('$|J|$','Interpreter','latex')
    title('Current along surface curve','Interpreter','latex')
end

function vertexIds = nearestMeshVertices(p,queryPoints)

    vertexIds = zeros(size(queryPoints,1),1);
    for iq = 1:size(queryPoints,1)
        [~,vertexIds(iq)] = min(vecnorm(p - queryPoints(iq,:),2,2));
    end
end

