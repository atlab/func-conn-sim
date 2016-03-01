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
Ne=400%800;
Ni=100%200;

% Number of neurons in external
% Poisson-spiking layer
Nx=400;

% Rate of external (Poisson-spiking) population
rx=.005;

% Connection strengths
% PSP size (in mV) is just a bit smaller than
% synaptic weight divided by membrane time
% cnst (in ms)
jee=.5*10;
jei=2*10;
jie=1.5*10;
jii=2*10;
jex=3*10/2;
jix=1.1*10/2;

% Connection probabilities
pee=.1;
pei=.17;
pie=.08;
pii=.3;
pex=.15;
pix=.7;

% Constant external input bias
Fe=6;
Fi=4;

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
taume=15;
taumi=15;

DTe=2;
DTi=2;
VTe=12;
VTi=12;

% Thresholds (interpreted as distance 
% from reset to thr in mV)
Vthe=20;
Vthi=20;

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
Sxe=zeros(Nx,1);
Sxi=zeros(Nx,1);

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
    Sxe=Sxe-dt*Sxe/tausx;
    Sxi=Sxi-dt*Sxi/tausx;
    
    % and for membrane potentials
    Ve=Ve+dt*(-Ve+DTe*exp((Ve-VTe)/DTe)+Fe+Jee*Se-Jei*Si+Jex*Sxe)/taume;
    Vi=Vi+dt*(-Vi+DTi*exp((Vi-VTi)/DTi)+Fi+Jie*Se-Jii*Si+Jix*Sxi)/taumi;
 
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
    Sxe(temp)=Sxe(temp)+1/tausx;
    
    %temp=find(rand(Nx,1)<rx*dt);
    Sxi(temp)=Sxi(temp)+1/tausx;
    
    % Display percent finished every 10%
    if(mod(100*i/Nt,10)==0)
        disp(sprintf('%d                     %.1fHz      %.1fHz',round(100*i/Nt),1000*nespikes/(Ne*i*dt),1000*nispikes/(Ni*i*dt)))
    end
    
end
t0=toc; % time to run sim (in s)

disp(sprintf('Time to run sim: %.1f sec',t0))

% Get rid of padded zeros
espikes=espikes(1:nespikes,:);
ispikes=ispikes(1:nispikes,:);


% Compute rates
erate=1000*nnz(espikes(:,1))/(Ne*T)
irate=1000*nnz(ispikes(:,1))/(Ni*T)

% % Save all data
%save SpikingSim2;

% Save just spikes and weight matrices
save SpikingSim2WS.mat espikes ispikes Jee Jei Jie Jii;

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


if(true)
    
winsize=100;
Tburn=500;
    

N=Ne+Ni;
ispikes(:,2)=ispikes(:,2)+Ne;
allspikes=[espikes; ispikes];
[~,I]=sort(allspikes(:,1));
allspikes=allspikes(I,:);
rates=StoPN(allspikes(:,1),allspikes(:,2),N,T,winsize);
rates=rates(ceil(Tburn/(2*winsize)):end-ceil(Tburn/(2*winsize)),:);
%rates=rates(:,mean(rates)>.001);
C=corrcoef(rates);
[II,JJ]=meshgrid(1:N,1:N);
C=C(:);II=II(:);JJ=JJ(:);
eicorrs=mean(C(II<=Ne & JJ>Ne & isfinite(C)))
eecorrs=mean(C(II<=Ne & JJ<=Ne & isfinite(C)))
iicorrs=mean(C(II>Ne & JJ>Ne & isfinite(C)))

end
