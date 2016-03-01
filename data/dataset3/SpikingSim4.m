% Simulate a recurrent network of integrate-and-fire
% neurons where the exc and inh neurons are driven by an 
% external population of Poisson-spiking neurons.


clear 
close all

% Time to run simulation and time bin size
% All time is in ms
T=5*60*1000;
dt=.2;

% Number of exc and inh neurons
Ne=800;
Ni=200;

% Number of neurons in external
% Poisson-spiking layer
Nx=500;

% Rate of external (Poisson-spiking) population
rx=.005;

% Connection strengths
% PSP size (in mV) is just a bit smaller than
% synaptic weight divided by membrane time
% cnst (in ms)
jee=.5*10;
jei=0.5*10;
jie=1.5*10;
jii=2*10;
jex=1*10;
jix=0.8*10;

% Connection probabilities
pee=.05;
pei=.4;
pie=.05;
pii=.15;
pex=.10;
pix=.10;

% Constant external input bias
Fe=2.2;
%Fi=2.2;%10;
Fi=6;

% % This stuff really only makes sense when N gets big
%wee=jee*pee*Ne;
%wei=jei*pei*Ni;
%wie=jie*pie*Ne;
%wii=jii*pii*Ni;
%disp(sprintf('\nThis list should decrease: %2.2f %2.2f %2.2f',Fe/Fi,wei/wii,wee/wie))
%disp(sprintf('Exc, Inh rates approximately: %2.2fHz %2.2fHz\n',1000*(Fe*wii-Fi*wei)/(wei*wie-wee*wii),1000*(Fe*wie-Fi*wee)/(wei*wie-wee*wii)))


% Synaptic time constants
% Inhibition cannot be slower than 
% excitation, or there will be synchronous
% rate oscillations. Even though this 
% is unrealistic.
tause=5;
tausi=4;
tausx=5;

% Membrane time constants
taume=10;
taumi=10;

% Thresholds (interpreted as distance 
% from reset to thr in mV)
Vthe=10;
Vthi=10;

% Number of time bins
Nt=ceil(T/dt);

% Random initial conditions
% for membrane potentials
Ve=rand(Ne,1)*Vthe;
Vi=rand(Ni,1)*Vthi;

% Zero initial conditions
% for synaptic currents
Se=zeros(Ne,1);
Si=zeros(Ni,1);
Sx=zeros(Nx,1);

% Connection matrices
Jee=sparse(jee*(rand(Ne,Ne)<pee));
Jei=sparse(jei*(rand(Ne,Ni)<pei));
Jie=sparse(jie*(rand(Ni,Ne)<pie));
Jii=sparse(jii*(rand(Ni,Ni)<pii));
Jex=sparse(jex*(rand(Ne,Nx)<pex));
Jix=sparse(jix*(rand(Ni,Nx)<pix));


tic
maxns=T*Ne*.05; % Max number of spikes in all e or i cells
espikes=zeros(maxns,2); % Store spikes here. espikes(:,1) is spike times
ispikes=zeros(maxns,2); % and espikes(:,2) is neuron index
nespikes=0; % Keeps track of number of e
nispikes=0; % and i spikes so far
disp('Percent sim complete   exc-rate   inh-rate ')
for i=1:Nt
        
    % Euler steps for synaptic currents
    Se=Se-dt*Se/tause;
    Si=Si-dt*Si/tausi;
    Sx=Sx-dt*Sx/tausx;
    
    % and for membrane potentials
    Ve=Ve+dt*(-Ve+Fe+Jee*Se-Jei*Si+Jex*Sx)/taume;
    Vi=Vi+dt*(-Vi+Fi+Jie*Se-Jii*Si+Jix*Sx)/taumi;
 
    % Find where spikes occur and propogate them
    temp=find(Ve>Vthe); % Find e spikes    
    if(nespikes+numel(temp)>maxns)
        warning('Max number of e spikes reached, exiting loop.')
        break;
    end
    % Store times and indices 
    espikes(nespikes+1:nespikes+numel(temp),1)=(i-1)*dt;
    espikes(nespikes+1:nespikes+numel(temp),2)=temp;
    Se(temp)=Se(temp)+1/tause; % increment Se
    Ve(temp)=0; % Reset Ve    
    nespikes=nespikes+numel(temp); % Increment nespikes
    
    % Same for i.
    temp=find(Vi>Vthi);
    if(nispikes+numel(temp)>maxns)
        warning('Max number of i spikes reached, exiting loop.')
        break;
    end
    ispikes(nispikes+1:nispikes+numel(temp),1)=(i-1)*dt;
    ispikes(nispikes+1:nispikes+numel(temp),2)=temp;
    Si(temp)=Si(temp)+1/tausi;
    Vi(temp)=0;
    nispikes=nispikes+numel(temp);
    
    
    % Poisson spiking external inputs
    temp=find(rand(Nx,1)<rx*dt);
    Sx(temp)=Sx(temp)+1/tausx;
    
    % Display percent finished every 10%
    if(mod(100*i/Nt,10)==0)
        disp(sprintf('%d                     %.1fHz      %.1fHz',round(100*i/Nt),1000*nespikes/(Ne*i*dt),1000*nispikes/(Ni*i*dt)))
    end
    
end
t0=toc; % time to run sim (in s)

fprintf('Time to run sim: %.1f sec', t0)

% Get rid of padded zeros
espikes=espikes(1:nespikes,:);
ispikes=ispikes(1:nispikes,:);


% Compute rates
erate=1000*nnz(espikes(:,1))/(Ne*T)
irate=1000*nnz(ispikes(:,1))/(Ni*T)

% % Save all data
%save SpikingSim1;

% Save just spikes and weight matrices
save SpikingSim4WS.mat espikes ispikes Jee Jei Jie Jii;

% Set to true or paste the code below
% to plot a raster plot of the first
% 1s of the simulation
if(true)
    
    figure
    subplot(2,1,1)
    Iplot=espikes(:,1)<=1000;
    plot(espikes(Iplot,1),espikes(Iplot,2),'k.')
    ylabel('exc neuron index')
    
    subplot(2,1,2)
    Iplot=ispikes(:,1)<=1000;
    plot(ispikes(Iplot,1),ispikes(Iplot,2),'k.')
    ylabel('inh neuron index')
    xlabel('time (ms)')
    
end
