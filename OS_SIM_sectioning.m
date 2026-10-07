
%   Implements Eqs. (8)-(10) of Neil, Juskaitis & Wilson,
%   Opt. Lett. 22(24), 1905 (1997).
%
%   I_p(u) = |g(2u, nu_tilde)|,  with the Stokseth approximation
%   g(u, nu_tilde) = f(nu_tilde) * 2*J1(x)/x,   x = u*nu_tilde*(1 - nu_tilde/2)
%
%   Usage:  plot_OSSIM_axial_response(1.0)      % single value
%           plot_OSSIM_axial_response([0.4 0.8 1.0])  % compare several
%
%   u is the normalized defocus parameter: u = 8*(pi/lambda)*z*sin^2(alpha/2)

nu_tilde = (0.4:0.2:1);   % default: maximum sectioning strength

% Normalized defocus axis
u = linspace(-20, 20, 2001);

figure(99); hold on; grid on;
% colors = lines(numel(nu_tilde));
colors = gray(5);

for k = 1:numel(nu_tilde)
    nu = nu_tilde(k);
    Ip = osssim_Ip(u, nu);
    plot(u, Ip, 'LineWidth', 1.8, 'Color', colors(k,:), ...
        'DisplayName', sprintf('$\\tilde{f} = %.2f$', nu));
end

xlabel('normalized defocus  $u$', 'Interpreter', 'latex');
ylabel('$I_p(u)$  (axial response)', 'Interpreter', 'latex');
% title('OS-SIM axial response (Neil et al. 1997, Eqs. 8--10)', ...
%     'Interpreter', 'latex');
legend('Interpreter', 'latex', 'Location', 'northeast');
ylim([0 1.05]); xlim([min(u) max(u)]);
axis square tight
set(gca, 'FontSize', 24);

% ---- Core model -------------------------------------------------------
function Ip = osssim_Ip(u, nu_tilde)
% I_p(u) = |g(2u, nu_tilde)|   (Eq. 8), using Stokseth g (Eq. 9).
Ip = abs( g_stokseth(2*u, nu_tilde) );

% Normalize to unity at focus (u = 0) so the FWHM is easy to read off.
Ip = Ip / max(Ip);
% var(Ip(:))
end

function val = g_stokseth(u, nu_tilde)
% Stokseth approximation to g(u, nu_tilde), Eq. (9).
f = 1 - 0.69*nu_tilde + 0.0076*nu_tilde.^2 + 0.043*nu_tilde.^3;  % f(nu~)
x = u * nu_tilde .* (1 - nu_tilde/2);                            % argument

% 2*J1(x)/x, with the removable singularity at x = 0 handled (limit = 1)
jinc = ones(size(x));
nz = (x ~= 0);
jinc(nz) = 2 .* besselj(1, x(nz)) ./ x(nz);

val = f .* jinc;
end
