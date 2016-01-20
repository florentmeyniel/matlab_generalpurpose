% ======================================================
			hgexport (matlab)
% ======================================================

% create a figure
figure;
set(gcf, 'Color', [1 1 1])

% draw some plot
KL = @(p, q) p.*log(p./q) + (1-p).*log((1-p)./(1-q));
q=0.01:0.01:0.99;
plot(KL(0.5, q), 'LineWidth', 2); hold on
plot(KL(0.40, q), 'r', 'LineWidth', 2)
plot([1 length(q)], 1/20*[1 1], 'k')

% adjust font size & type
xlabel('updated parameter: p_0=n_g/n', 'FontName', 'Arial', 'FontSize', 12);
ylabel('Kullback-Leiber divergence', 'FontName', 'Arial', 'FontSize', 12)
xlim([20 80])

% Export image as eps file into the home folder
hgexport(gcf, 'tmp.eps') 

% ======================================================
			export_fig (matlab toolbox)
% ======================================================

export_fig('/home/fmeyniel/test', '-eps', 3);

% NB: inkascape can save in pdf, but not eps. However, pdf can easily be converted to eps: pdftops -eps [file].pdf
