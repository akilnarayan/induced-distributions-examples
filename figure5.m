% Generates contour plots of cdf's from a 2D gaussian distribution

clear
close all

addpath('induced-distributions');
% Define polynomial family:
alph = 2;
rho = 0;
d = 2;
k = 8;

% Grid for plotting
M = 100;
%xp = linspace(-sqrt(2*k) - 0.5, sqrt(2*k) + 0.5, M).';
xp = linspace(-4.5, 4.5, M).';
[X, Y] = ndgrid(xp, xp);

% Gaussian pdf:
normal_pdf = normpdf(X, 0, 1/sqrt(2)).*normpdf(Y, 0, 1/sqrt(2));

% Equilibrium pdf:
flags0 = (X.^2 + Y.^2) > (2*k);
equil_pdf = zeros(size(X));
equil_pdf = 1/(2*pi*k^2) * (2*k - X.^2 - Y.^2);
equil_pdf(flags0) = 0;

recurrence = @(NN) hermite_recurrence(NN+1, rho);

% Tensor product
[Kx,Ky] = ndgrid((0:k).', (0:k).');
tp_indices = [Kx(:) Ky(:)];
N = size(tp_indices, 1);
V = mpoly_eval([X(:) Y(:)], tp_indices, recurrence);
V = V.^2.*repmat(normal_pdf(:), [1 N]);
%idists = zeros([size(X, 1) size(X, 2) N]);
%fprintf('Idist computation, N = %d, n = ', N)
%for n = 1:N
%  idists(:,:,n) = idist_freud(X, tp_indices(n,1), alph, rho) .* idist_freud(Y, tp_indices(n,2), alph,rho);
%  fprintf('%d, ', n);
%end
%fprintf('\n');

% Tensor product
flags = true([N 1]);
idistpdf_tp = reshape(sum(V(:,flags),2)/sum(flags), size(X));

% Euclidean degree
flags = sqrt(sum(tp_indices.^2, 2)) <= k;
idistpdf_ed = reshape(sum(V(:,flags),2)/sum(flags), size(X));

% Total degree
flags = sum(tp_indices, 2) <= k;
idistpdf_td = reshape(sum(V(:,flags),2)/sum(flags), size(X));

% Hyperbolic cross
flags = prod(tp_indices+1, 2) <= (k+1);
idistpdf_hc = reshape(sum(V(:,flags),2)/sum(flags), size(X));

%levels = linspace(0, 0.3, 21);
levels = sort([linspace(0, max(normal_pdf(:)), 11), ...
               linspace(0, max(equil_pdf(:)), 11), ...
               linspace(0, max(idistpdf_hc(:)), 11)]);
%levels = levels(2:end) - 1/2*(levels(2) - levels(1));

axis_props = {'fontsize', 16, 'fontweight', 'b'};
set_latex = @(hh) set(hh, 'interpreter', 'latex');

figure;
subplot(2,3,1);
[~,h] = contour(X, Y, (normal_pdf), levels);
set(h, 'linewidth', 2);
set(gca, 'ytick', [min(xp), 0, max(xp)]);
set(gca, 'xtick', [min(xp), 0, max(xp)]);
set(gca, axis_props{:});
set_latex(xlabel('$x^{(1)}$'));
set_latex(ylabel('$x^{(2)}$'));
set_latex(title('$w(x,y)$'));
axis equal; axis square;
axis([min(xp) max(xp) min(xp) max(xp)]);
caxis([0 max(levels)]);

subplot(2,3,2);
[~,h] = contour(X, Y, (equil_pdf), levels);
set(h, 'linewidth', 2);
set(gca, 'ytick', [min(xp), 0, max(xp)]);
set(gca, 'xtick', [min(xp), 0, max(xp)]);
set(gca, axis_props{:});
set_latex(xlabel('$x^{(1)}$'));
set_latex(ylabel('$x^{(2)}$'));
set_latex(title('$w_\infty(x,y)$'));
axis equal; axis square;
axis([min(xp) max(xp) min(xp) max(xp)]);
caxis([0 max(levels)]);

subplot(2,3,3);
[~,h] = contour(X, Y, (idistpdf_tp), levels);
set(h, 'linewidth', 2);
set(gca, 'ytick', [min(xp), 0, max(xp)]);
set(gca, 'xtick', [min(xp), 0, max(xp)]);
set(gca, axis_props{:});
set_latex(xlabel('$x^{(1)}$'));
set_latex(ylabel('$x^{(2)}$'));
set_latex(title('$w_{\mathrm{TP}}(x,y)$'));
axis equal; axis square;
axis([min(xp) max(xp) min(xp) max(xp)]);
caxis([0 max(levels)]);

subplot(2,3,4);
[~,h] = contour(X, Y, (idistpdf_ed), levels);
set(h, 'linewidth', 2);
set(gca, 'ytick', [min(xp), 0, max(xp)]);
set(gca, 'xtick', [min(xp), 0, max(xp)]);
set(gca, axis_props{:});
set_latex(xlabel('$x^{(1)}$'));
set_latex(ylabel('$x^{(2)}$'));
set_latex(title('$w_{\mathrm{ED}}(x,y)$'));
axis equal; axis square;
axis([min(xp) max(xp) min(xp) max(xp)]);
caxis([0 max(levels)]);

subplot(2,3,5);
[~,h] = contour(X, Y, (idistpdf_td), levels);
set(h, 'linewidth', 2);
set(gca, 'ytick', [min(xp), 0, max(xp)]);
set(gca, 'xtick', [min(xp), 0, max(xp)]);
set(gca, axis_props{:});
set_latex(xlabel('$x^{(1)}$'));
set_latex(ylabel('$x^{(2)}$'));
set_latex(title('$w_{\mathrm{TD}}(x,y)$'));
axis equal; axis square;
axis([min(xp) max(xp) min(xp) max(xp)]);
caxis([0 max(levels)]);

subplot(2,3,6);
[~,h] = contour(X, Y, (idistpdf_hc), levels);
set(h, 'linewidth', 2);
set(gca, 'ytick', [min(xp), 0, max(xp)]);
set(gca, 'xtick', [min(xp), 0, max(xp)]);
set(gca, axis_props{:});
set_latex(xlabel('$x^{(1)}$'));
set_latex(ylabel('$x^{(2)}$'));
set_latex(title('$w_{\mathrm{HC}}(x,y)$'));
axis equal; axis square;
axis([min(xp) max(xp) min(xp) max(xp)]);
caxis([0 max(levels)]);

hp4 = get(subplot(2,3,6),'Position');
colormap hot
colorbar('Position', [hp4(1)+hp4(3)+0.03  hp4(2)  0.03  hp4(2)+hp4(3)*3.3])
caxis([0 max(levels)]);

set(gcf, 'pos', [1 1 1000 600]);
