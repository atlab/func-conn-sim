function main
processSim('./data/dataset5/SpikingSim1BWS.mat', 'sim 1', './sim1')
processSim('./data/dataset5/SpikingSim2WS.mat', 'sim 2', './sim2')
processSim('./data/dataset4/SpikingSim3WS.mat', 'sim 3', './sim3')
end


function processSim(path, title_, figurefile)
[Cee, Cei, Cii, See, Sei, Sii] = avgCorrs(path);
subplot 121, plotAverages([Cee Cei Cii]), title corrs
subplot 122, plotAverages([See Sei Sii]), title 'pcorrs', set(gca, 'YTickLabels',{})
suptitle(title_)
set(gcf, 'PaperSize', [1 1]*3, 'PaperPosition', [0 0 1 1]*3)
print('-dpdf', figurefile)
end


function plotAverages(X)
plot(X', 'k-')
hold on
plot(X', 'k^', 'MarkerSize', 4, 'MarkerFaceColor', [.5 .5 .5])
plot([0 4],[0 0])
hold off
box off
xlim([0.5 3.5]) 
set(gca, 'XTick', 1:3, 'XTickLabel', {'E/E' 'E/I' 'I/I'})
yticks = -0.02:0.01:0.06;
ytickLabels = arrayfun(@(x) sprintf('%1.2f', x), yticks, 'uni', false);
set(gca, 'YTick', yticks, 'YTickLabels', ytickLabels) 
ylim([-0.012 0.042])
grid on
end


function [Cee, Cei, Cii, See, Sei, Sii] = avgCorrs(filepath)
[CC, JJ, ne] = load_matrices(filepath);
n = size(CC,1);
fraction = 0.4;  % fraction of the network
ntrials = 11;
Cee = nan(ntrials,1);
Cei = nan(ntrials,1);
Cii = nan(ntrials,1);
See = nan(ntrials,1);
Sei = nan(ntrials,1);
Sii = nan(ntrials,1);
progress_bar = waitbar(0, 'Averaging trials ...');
for trial = 1:ntrials
    waitbar(trial/ntrials)
    ix = rand(n,1) < fraction;
    ne_ = sum(ix(1:ne));
    C = CC(ix, ix);
    J = JJ(ix, ix);
    S = get_partials(C);
    [Cee(trial), Cei(trial), Cii(trial)] = mean_strength(corrcov(C), ne_);
    [See(trial), Sei(trial), Sii(trial)] = mean_strength(-corrcov(S), ne_);
    plot_connectivity(J, C, S);
end
close(progress_bar)
end


function [Ree, Rei, Rii] = mean_strength(R, ne)
n = size(R,1);
[i,j] = meshgrid(1:n,1:n);
Ree = mean(R(i<j & j<=ne));
Rei = mean(R(i<j & j>ne & i<=ne));
Rii = mean(R(i<j & j>ne & i>ne));
end


function K = get_partials(C)
% estimatate partial correlations using lv-glasso
alpha = 0.002;
beta = 0.05;
out = lvglasso(C, alpha, beta);
K = out.S - out.L;
n = size(K,1);
[i,j] = meshgrid(1:n, 1:n);
fprintf('Connectivity: %1.3f     ', 1-mean(~out.S(i<j)));
fprintf('Latent Units: %d / %d\n', sum(logical(out.eigL)), n);
end


function plot_connectivity(truth, cov, icov)
[i, j] = meshgrid(1:size(icov,1), 1:size(icov,1));
subplot 131
imagesc(logical(truth));  axis image
frac = mean(logical(truth(i<j)));
subplot 132
R = corrcov(cov);
thresh = quantile(abs(R(i<j)), 1-frac);
imagesc(abs(R)>thresh & i~=j);
subplot 133
S = -corrcov(icov);
thresh = quantile(abs(S(i<j)), 1-frac);
imagesc(abs(S)>thresh & i~=j);
end


function [C, J, ne] = load_matrices(filename)
s = load(filename);
T = max(max(s.espikes(:,1)), max(s.espikes(:,1)));
ne = max(s.espikes(:,2));   % number of excitatory neurons
ni = max(s.ispikes(:,2));   % number of excitatory interneurons

% compute covariance matrix
binwidth = 100;  % ms
nbins = ceil(T/binwidth);   % number of time bins
ipos = 1;
epos = 1;
C = 0;  % covariances
M = 0;  % means
progress_bar = waitbar(0, 'Computing covariance matrix...');
for bin = 1:nbins
    vector = zeros(ne+ni, 1)+eps;
    while epos <= size(s.espikes, 1) && s.espikes(epos,1)<binwidth*bin
        cellnumber = s.espikes(epos,2);
        vector(cellnumber) = vector(cellnumber)+1;
        epos = epos + 1;
    end
    while ipos <= size(s.ispikes, 1) && s.ispikes(ipos,1)<binwidth*bin
        cellnumber = s.ispikes(ipos,2) + ne;
        vector(cellnumber) = vector(cellnumber)+1;
        ipos = ipos + 1;
    end
    waitbar(bin/nbins)
    C = C + vector*vector'/nbins;
    M = M + vector/nbins;
end
close(progress_bar)
C = C - M*M';

% compute the connectivity matrix
J = [s.Jee s.Jei; s.Jie s.Jii];
end