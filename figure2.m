% Plots multi-index sets in two dimensions.
 
clear
close all

d = 2;
k = 20;

[X,Y] = ndgrid((0:k).', (0:k).');
alpha_tp = [X(:) Y(:)];

alpha_td = alpha_tp(sum(alpha_tp, 2) <= k, :);
alpha_ed = alpha_tp(sqrt(sum(alpha_tp.^2, 2)) <= k, :);
alpha_hc = alpha_tp(prod(alpha_tp+1, 2) <= k+1, :);

axis_props = {'fontsize', 16, 'fontweight', 'b'};
set_latex = @(hh) set(hh, 'interpreter', 'latex');

figure;
subplot(1,4,1);
plot(alpha_hc(:,1), alpha_hc(:,2), 'sk', 'markerfacecolor', 'k', 'markeredgecolor', 'k');
axis equal; axis square;
set(gca, axis_props{:});
set_latex(xlabel('$\lambda^{(1)}$'));
set_latex(ylabel('$\lambda^{(2)}$'));
set_latex(title(sprintf('$\\Lambda_{\\mathrm{HC}}(%d)$', k)));
set(gca, 'xtick', []);
set(gca, 'ytick', []);
axis([0 k+2 0 k+2]);

subplot(1,4,2);
plot(alpha_td(:,1), alpha_td(:,2), 'sk', 'markerfacecolor', 'k', 'markeredgecolor', 'k');
axis equal; axis square;
set(gca, axis_props{:});
set_latex(xlabel('$\lambda^{(1)}$'));
set_latex(ylabel('$\lambda^{(2)}$'));
set_latex(title(sprintf('$\\Lambda_{\\mathrm{TD}}(%d)$', k)));
set(gca, 'xtick', []);
set(gca, 'ytick', []);
axis([0 k+2 0 k+2]);

subplot(1,4,3);
plot(alpha_ed(:,1), alpha_ed(:,2), 'sk', 'markerfacecolor', 'k', 'markeredgecolor', 'k');
axis equal; axis square;
set(gca, axis_props{:});
set_latex(xlabel('$\lambda^{(1)}$'));
set_latex(ylabel('$\lambda^{(2)}$'));
set_latex(title(sprintf('$\\Lambda_{\\mathrm{ED}}(%d)$', k)));
set(gca, 'xtick', []);
set(gca, 'ytick', []);
axis([0 k+2 0 k+2]);

subplot(1,4,4);
plot(alpha_tp(:,1), alpha_tp(:,2), 'sk', 'markerfacecolor', 'k', 'markeredgecolor', 'k');
axis equal; axis square;
set(gca, axis_props{:});
set_latex(xlabel('$\lambda^{(1)}$'));
set_latex(ylabel('$\lambda^{(2)}$'));
set_latex(title(sprintf('$\\Lambda_{\\mathrm{TP}}(%d)$', k)));
set(gca, 'xtick', []);
set(gca, 'ytick', []);
axis([0 k+2 0 k+2]);
