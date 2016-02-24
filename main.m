function main
[C, ~, ne] = load_matrices('./data/dataset3/SpikingSim1WS.mat');
d = 5;  %downsample factor  -- replace with randsample
ne = floor(ne/d);
C = C(1:d:end,1:d:end);
plot_correlations(corrcov(C));
S = get_partials(C);
[Cee, Cei, Cii] = mean_strength(corrcov(C), ne);
[See, Sei, Sii] = mean_strength(-corrcov(S), ne);
plot([Cee Cei Cii; See Sei Sii]')
legend corrs parcorrs
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
alpha = 0.001;
beta = 0.1;
out = lvglasso(C, alpha, beta);
K = out.S + out.L;
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