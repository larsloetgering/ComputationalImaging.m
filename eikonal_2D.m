%% preconfiguration
set(0, 'DefaultAxesFontSize', 28);
set(0, 'DefaultFigureColor', 'w');
set(0,'defaultAxesFontName', 'serif')
set(0, 'DefaultLineLineWidth', 3);
clc

% 2D solver of eikonal equation
wavelength = 500e-9;
Nx = 2000;
Nz = 2000;
chi = zeros(Nx, Nz);
hx = 1e1*wavelength;
oversampling = 600;
hz = hx / oversampling;
Lx = Nx*hx;
Lz = Nz*hz;

% get grid
x = (-Nx/2:Nx/2-1)*hx;
z = (0:Nz-1)*hz;

% boundary condition
% case 1: plane wave
chi(:,1) = ones;

% case 2: grin lens
n = 1 + 1*exp( -( (x.^2/(5e6*hx).^2)' + ...
    ((z-Nz*hz/2)/(25e6*hz)^2).^2)/(2*(10*hx).^2));

figure(1), clf
subplot(1,2,1)
imagesc(z/wavelength, x/wavelength, n, [0.99 2.01])
yticks( [-5000 0 5000] )
axis square
xlabel('z / \lambda')
ylabel('x / \lambda')
colormap gray
colorbar
title('n(x, z)')

%% eikonal solver

for jj = 2:Nz
    max( getGradient(chi(:,jj-1), hx).^2 )
    chi(:,jj) = chi(:,jj-1) + ...
        hz * sqrt( n(:,jj).^2 - abs(getGradient(chi(:,jj-1), hx)).^2 );
end

%% visualization

figure(1)
subplot(1,2,2)
imagesc( z, x, mod( k*(chi - chi(:,1)), 2*pi ) )
axis square off
h = colorbar;
h.Ticks = [0 2*pi-1e-1];
h.TickLabels = {'0', '2\pi'};
h2 = gca;
h2.Colormap = hsv;
title('mod( k \cdot \chi, 2\pi)')
hold on
h3 = rectangle('Position',[Lz/4 -Lx/4 Lz/2 Lx/2],'Curvature',[1,1],'EdgeColor','w','LineWidth',4);
hold off

%% sanity check: from the estimated wavefront, we can retrieve the refractive index
% (by means of the eikonal equation)

n_est = sqrt(( ( chi(2:end,1:end-1)-chi(1:end-1,1:end-1) )/hx).^2 + ...
             ( ( chi(1:end-1,2:end)-chi(1:end-1,1:end-1) )/hz).^2);

figure(4)
subplot(3,1,1)
imagesc(n_est)
axis image off
colorbar
colormap jet
title('forward model')

subplot(3,1,2)
imagesc(n)
axis image off
colorbar
colormap jet
title('ground truth')

subplot(3,1,3)
imagesc( abs(n_est - n(1:end-1,1:end-1)) ./ n(1:end-1,1:end-1) * 100, [0 1] )
axis image off
colorbar
colormap jet
title('rel. error [%]')

%% compare to projection approximation

k = 2*pi/wavelength;
chi_proj = cumsum(n, 2)*hz;
chi_proj = chi_proj(:,end) - chi_proj(:,1);
chi_eik = (chi(:,end)-chi(:,1));

figure(5), clf
plot(x, k*chi_eik,'k')
hold on
plot(x, k*chi_proj,'g--')
grid on
axis square
title('rel. phase shift')

%% compute residuum

residuum = norm(n_est - n(1:end-1,1:end-1), 'fro') / norm(n(1:end-1,1:end-1), 'fro')

function dq = getGradient(q, hx)

dq = 0*q;
dq(2:end) = (q(2:end) - q(1:end-1))/hx;
dq(1) = dq(2);

end
