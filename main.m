function main
processSim('./data/dataset3/SpikingSim1WS.mat', 'simulation 1', './sim1')
processSim('./data/dataset3/SpikingSim2WS.mat', 'simulation 2', './sim2')
processSim('./data/dataset3/SpikingSim3WS.mat', 'simulation 3', './sim3')
end


function processSim(path, title_, figurefile)
[Cee, Cei, Cii, See, Sei, Sii] = avgCorrs(path);
subplot 121, plotAverages([Cee Cei Cii]), title corrs
subplot 122, plotAverages([See Sei Sii]), title 'pcorrs'
suptitle(title_)
set(gcf, 'PaperSize', [4 4], 'PaperPosition', [0 0 4 4])
print('-dpdf', figurefile)
end


function plotAverages(X)
plot(X', 'k.-', 'MarkerSize', 20)
xlim([0.5 3.5]) 
set(gca, 'XTick', 1:3, 'XTickLabel', {'EE' 'EI' 'II'})
yticks = -0.02:0.01:0.03;
ytickLabels = arrayfun(@(x) sprintf('%1.2f', x), yticks, 'uni', false);
set(gca, 'YTick', yticks, 'YTickLabels', ytickLabels) 
ylim([-0.012 0.021])
grid on
end


function [Cee, Cei, Cii, See, Sei, Sii] = avgCorrs(filepath)
[CC, ~, ne] = load_matrices(filepath);
n = size(CC,1);
fraction = 0.5;  % fraction of the network
ntrials = 20;
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
    S = get_partials(C);
    [Cee(trial), Cei(trial), Cii(trial)] = mean_strength(corrcov(C), ne_);
    [See(trial), Sei(trial), Sii(trial)] = mean_strength(-corrcov(S), ne_);
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
K = out.S + out.L;
n = size(K,1);
[i,j] = meshgrid(1:n, 1:n);
fprintf('Connectivity: %1.3f     ', 1-mean(~out.S(i<j)));
fprintf('Latent Units: %d / %d\n', sum(logical(out.eigL)), n);
end


function plot_correlations(R)
subplot 211
n = size(R,1);
imagesc(R(1:10:end,1:10:end), [-1 1]*0.1)
axis image
colormap(doppler)
colorbar
subplot 212
[i,j] = meshgrid(1:n, 1:n);
hist(R(i<j),100);
grid on
hist(R(i<j),100);
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