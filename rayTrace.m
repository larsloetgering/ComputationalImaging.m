addpath(addpath(genpath('C:\Users\Lars Loetgering\Dropbox\Codes\fracPty')))
% simpleRayTrace

% configuration
set(0, 'DefaultAxesFontSize', 30);
set(0, 'DefaultFigureColor', 'w');
set(0,'defaultAxesFontName', 'serif')

% main program

% free space
z = 10e-2;          % propagation distance
F = [1 z;...        % propagation matrix
     0 1];      
numRays = 50;       % number of rays
maxAngle = 5e-2;    % maximum (paraxial) angle
V = [0*ones(1,numRays); linspace(-maxAngle,maxAngle,numRays)]; % ray bundle at x = 0
V(1,1:2:numRays) = 1.5e-3 * ones(1, numRays/2);     % set half of the positions to x ~= 0

% propagation step: map V to Vp ("Vprime")
Vp = F * V;         % free space propagation (from source to lens)

% plot locations
figure(1), clf      % clf clears previous figure if present
hold on             % hold on allows to draw multiple points
cmap = jet(numRays);% get colormap to draw rays in different colors
for k  = 1:numRays
    plot([0 z], [V(1,k), Vp(1,k)], 'o-','color',cmap(k,:),'MarkerFaceColor',cmap(k,:))
end
hold off

% lens + free-space
f = z/2;            % focal length of lens
L = [1 0; -1/f 1];  % lens matrix
Vpp = F*(L*Vp);     % lens transormation + propagation step 

figure(1)
hold on
for k  = 1:numRays
    plot([z, 2*z], [Vp(1,k), Vpp(1,k)], 'o-','color',cmap(k,:),'MarkerFaceColor',cmap(k,:))
end
% draw ellipse to indicate lens position (not mandatory)
dz = 2e-2;
dx = 1.5e-2;
pos = [z-dz/2 0-dx/2 dz dx]; 
h = rectangle('Position',pos,'Curvature',[1 1],'LineWidth',1);
axis square
grid on
h = gca;
h.LineWidth = 3;
xlabel('z / m'), ylabel('x / m')
hold off

