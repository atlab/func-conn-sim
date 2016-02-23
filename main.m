% load data
s = load('./data/dataset3/SpikingSim2WS.mat');
dt = 0.2;  %ms
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
progress_bar = waitbar(0,'Computing covariance matrix...');
for bin = 1:nbins
    vector = zeros(ne+ni, 1);
    while s.espikes(epos,1)<bin*dt
        cellnumber = s.espikes(epos,2);
        vector(cellnumber) = vector(cellnumber)+1;
        epos = epos + 1;
        if epos > size(s.espikes, 1)
            break
        end
    end
    while s.ispikes(ipos,1)<bin*dt
        cellnumber = s.ispikes(ipos,2) + ne;
        vector(cellnumber) = vector(cellnumber)+1;
        ipos = ipos + 1;
        if ipos > size(s.ispikes, 1)
            break
        end
    end
    waitbar(bin/nbins)
    C = C + vector*vector'/nbins;
    M = M + vector/nbins;
end
close(progress_bar)

C = C - M*M';
R = corrcov(C);

subplot 211
imagesc(R, [-1 1]*0.05)
axis image
colormap(doppler)
colorbar
subplot 212
[i,j] = meshgrid(1:ne+ni, 1:ne+ni);
hist(R(i<j),100);
ylim([0 200])


% estimatate partial correlations using lv-glasso 