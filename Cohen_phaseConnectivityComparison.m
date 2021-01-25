%% Effects of time lag and frequency matching on phase-based connectivity
% This code accompanies the paper: 
% Cohen, MX (2014). Effects of time lag and frequency matching on
% phase-based connectivity. Journal of Neuroscience Methods.
%
% mikeXcohen@gmail.com, June 2014

% notes: 
%  1) You will need the following files in the current directory: 
%         lfchans.mat, laplacian_perrinX.m
%  2) eeglab is required for topographical maps, but not for simulations or analyses.
%  3) code for figure 1 is presented at the bottom of this script.

clear

%% easily modifiable parameters

% center frequency in Hz
centfreq = 10;
% mean phase difference between two dipoles (state in ms, converted later to radians)
phasedif_ms = 5; % 5 or 25

% simulationType: 1 (equal stationary frequencies)
%                 2 (unequal stationary frequencies)
%                 3 (nonstationary frequencies)
simulationType = 2;

% Dispersion of frequency nonstationarity, in Hz. In the paper it was 
%  set to 2 or 4. Values between 1 and 5 are reasonable. It has an 
%  effect only when stimulationType is set to 3.
freqDist = 2;

% number of trials
ntrials = 300;

%% initial setup and load in necessary files

% mat file containing leadfield and channel locations
load lfchans

%% other parameters

% sampling rate
srate = 200;

% time for simulation (in seconds)
time  = (-1:1/srate:2)*1000;
phasedif = 2*pi*centfreq*phasedif_ms/1000;

% time points to save from final results, and baseline for power normalization
times2save = -400:50:1600; % in ms
baselinetimerange = [ -500 -200 ];

ntime  = length(time);
nchans = length(chanlocs);

%% pick source dipoles and electrodes

% indices of dipole locations (probably best not to change these)
dipoleOCC =   94;
dipolePFC = 1720;

% use X, Y, or Z oriented dipole (1, 2, or 3, respectively). 
% In the paper, Z was used.
whichOrientation = 3;

% Note: If you want to see the effects using horizontal 
%       diploes, set the orientation above to "2" and use 
%       electrodes F3 and PO3 below.


% in electrode labels; see {chanlocs.labels}
electrodes2use = { 'Fz';'Pz' };
% electrodes2use = { 'F3';'PO3' }; % see comment above for using these electrodes

% convert channel labels to indices
elecs2use = zeros(size(electrodes2use));
for i=1:length(electrodes2use)
    elecs2use(i) = find(strcmpi(electrodes2use{i},{chanlocs.labels}));
end

% To see the scalp dipole projections, use the following code (requires eeglab function topoplot)
%topoplot(squeeze(lf.Gain(:,whichOrientation,dipolePFC)),chanlocs);

%% intialize output matrices

simulatedEEG = zeros(nchans,ntime,ntrials); % simulated electrode data
sourceTimeSeries = zeros(ntime,2,ntrials); % data at source dipoles

% setup time indexing
times2saveidx = dsearchn(time',times2save');
baseidx = dsearchn(time',baselinetimerange');

%% simulate data (noise+sine waves)

for triali=1:ntrials
    
    % data comprise all noise
    data = randn(ntime,size(lf.Gain,3))./20;
    
    % specify frequency
    switch simulationType
        case 1
            % case 1: both dipoles oscillate at centfreq Hz
            trialFreq1 = centfreq;
            trialFreq2 = centfreq;
            [freqmod1,freqmod2,k1,k2] = deal(0); % see case 3
        case 2
            % case 2: one dipole oscillates faster; both oscillators are stationary
            trialFreq1 = centfreq;
            trialFreq2 = centfreq*1.2;
            [freqmod1,freqmod2,k1,k2] = deal(0); % see case 3
        case 3
            % case 3: Both oscillators are nonstationary
            trialFreq1 = centfreq;
            trialFreq2 = centfreq;
            
            % time-varying frequency modulation
            freqmod1 = interp(ceil(freqDist*rand(1,round(ntime/40)))-.5-(freqDist/2),41); freqmod1(ntime+1:end) = [];
            freqmod2 = interp(ceil(freqDist*rand(1,round(ntime/40)))-.5-(freqDist/2),41); freqmod2(ntime+1:end) = [];
            
            % frequency coefficients for generating arbitrary 'chirp' signal
            k1 = (centfreq/srate)*2*pi/trialFreq1;
            k2 = (centfreq/srate)*2*pi/trialFreq2; % was max()
    end
    
    % Gaussian window used for tapering sine waves
    gausWindow = exp( (-((time/1000)-.6).^2)/.1 );
    
    % create signals
    data(:,dipolePFC) = data(:,dipolePFC)' + sin(2*pi.*trialFreq1.*(time/1000) + k1*cumsum(freqmod1) + rand*.1*pi         ) .* gausWindow;
    tempts            = data(:,dipolePFC)' + sin(2*pi.*trialFreq2.*(time/1000) + k2*cumsum(freqmod2) + rand*.1*pi+phasedif) .* gausWindow;
    data(:,dipoleOCC) = data(:,dipoleOCC)' + tempts.*gausWindow;
    
    % simulated EEG data
    simulatedEEG(:,:,triali) = (data*squeeze(lf.Gain(:,whichOrientation,:))')';
    
    % get actual source time series
    sourceTimeSeries(:,:,triali) = data(:,[dipolePFC dipoleOCC]);
end

% also compute laplacian
simulatedLap = laplacian_perrinX(simulatedEEG,[chanlocs.X],[chanlocs.Y],[chanlocs.Z]);
% average reref (leadfield assumes average reference)
simulatedEEG = bsxfun(@minus,simulatedEEG,mean(simulatedEEG,1));

%% end data simulation part of script





%% begin wavelet convolution part of script





%% initial parameters for time-frequency analysis

% seeded phase synchronization (set to empty ("{}") for no analyses)
electrodes4seeded_synch = { chanlocs(elecs2use(1)).labels , chanlocs(elecs2use(2)).labels };

% wavelet parameters
min_freq =  2;
max_freq = 40;
num_frex = 25;
wavelet_cycle_range = [ 3 12 ];

frex = logspace(log10(min_freq),log10(max_freq),num_frex);


% gaussian width and time
s = logspace(log10(wavelet_cycle_range(1)),log10(wavelet_cycle_range(2)),num_frex)./(2*pi.*frex);
t = -2:1/srate:2;

% fft and convolution details
Ltapr  =  length(t);
Ldata  =  prod(ntime*ntrials);
Lconv1 =  Ldata+Ltapr-1;
Lconv  =  pow2(nextpow2( Lconv1 ));

wavelets = zeros(num_frex,length(t));
for fi=1:num_frex
    wavelets(fi,:) = exp(2*1i*pi*frex(fi).*t).*exp(-t.^2./(2*s(fi)^2));
end

% initialize output and hidden-layer TF matrices
tf    = zeros(nchans+2,num_frex,length(times2saveidx),2);
allphasevals = zeros(nchans,num_frex,length(times2save),ntrials,2);
synchOverTrials = zeros(2,length(electrodes4seeded_synch),nchans,num_frex,length(times2saveidx),2);
allAS = zeros(2,num_frex,ntime,ntrials,2);

%% run convolution

% loop around channels
for chani=1:nchans+2
    
    % FFT of data (channel or true dipoles)
    if chani<=nchans
        EEGfft = fft(reshape(simulatedEEG(chani,:,:),1,[]),Lconv);
        Lapfft = fft(reshape(simulatedLap(chani,:,:),1,[]),Lconv);
    else
        [EEGfft,Lapfft] = deal(fft(reshape(sourceTimeSeries(:,mod(chani,nchans),:),1,[]),Lconv));
    end
        
    
    % loop over frequencies and complete convolution
    for fi=1:num_frex
        
        %% Average reference
        % convolve and get analytic signal (as)
        as = ifft(EEGfft.*fft(wavelets(fi,:),Lconv),Lconv);
        as = as(1:Lconv1);
        as = reshape(as(floor((Ltapr-1)/2):end-1-ceil((Ltapr-1)/2)),ntime,ntrials);
        
        % enter into TF matrix
        temppow = mean(abs(as).^2,2);
        tf(chani,fi,:,1) = 10*log10( temppow(times2saveidx)/mean(temppow(baseidx(1):baseidx(2))) );
        
        % save phase values
        if chani<=nchans
            allphasevals(chani,fi,:,:,1) = as(times2saveidx,:);
        end
        
        % all values from all time points
        if chani==elecs2use(1)
            allAS(1,fi,:,:,1) = as;
        elseif chani==elecs2use(2)
            allAS(2,fi,:,:,1) = as;
        end
        
        %% Laplacian
        % convolve and get analytic signal (as)
        as = ifft(Lapfft.*fft(wavelets(fi,:),Lconv),Lconv);
        as = as(1:Lconv1);
        as = reshape(as(floor((Ltapr-1)/2):end-1-ceil((Ltapr-1)/2)),ntime,ntrials);
        
        % enter into TF matrix
        temppow = mean(abs(as).^2,2);
        tf(chani,fi,:,2) = 10*log10( temppow(times2saveidx)/mean(temppow(baseidx(1):baseidx(2))) );
        
        % save phase values
        if chani<=nchans
            allphasevals(chani,fi,:,:,2) = as(times2saveidx,:);
        end
        
        % all values from all time points
        if chani==elecs2use(1)
            allAS(1,fi,:,:,2) = as;
        elseif chani==elecs2use(2)
            allAS(2,fi,:,:,2) = as;
        end
        
    end % end frequency loop
end % end channel loop

%% compute phase connectivity over trials

for chanx=1:length(electrodes4seeded_synch)
    
    % ISPC (average reference and laplacian)
    synchOverTrials(1,chanx,:,:,:,1) = mean(exp(1i* bsxfun(@minus,angle(allphasevals(strcmpi(electrodes4seeded_synch{chanx},{chanlocs.labels}),:,:,:,1)),angle(allphasevals(:,:,:,:,1))) ),4);
    synchOverTrials(1,chanx,:,:,:,2) = mean(exp(1i* bsxfun(@minus,angle(allphasevals(strcmpi(electrodes4seeded_synch{chanx},{chanlocs.labels}),:,:,:,2)),angle(allphasevals(:,:,:,:,2))) ),4);
    
    
    % wPLI (average reference and laplacian)
    cdd = bsxfun(@times,allphasevals(strcmpi(electrodes4seeded_synch{chanx},{chanlocs.labels}),:,:,:,1),conj(allphasevals(:,:,:,:,1)));
    cdi = imag(cdd);
    synchOverTrials(2,chanx,:,:,:,1) = mean( abs( mean( abs(cdi).*sign(cdi) ,4) )./mean(abs(cdi),4) ,4);
    
    cdd = bsxfun(@times,allphasevals(strcmpi(electrodes4seeded_synch{chanx},{chanlocs.labels}),:,:,:,2),conj(allphasevals(:,:,:,:,2)));
    cdi = imag(cdd);
    synchOverTrials(2,chanx,:,:,:,2) = mean( abs( mean( abs(cdi).*sign(cdi) ,4) )./mean(abs(cdi),4) ,4);
end

% NaN's cause electrode data shifts in the eeglab topoplot function
synchOverTrials(isnan(synchOverTrials))=0;

%% plot data

times2plot = dsearchn(times2save',[200 800]');
freq2plot  = dsearchn(frex',centfreq);

clim = .8;

figure(1), clf, colormap hot
subplot(221)
topoplot(squeeze(mean(mean(abs(synchOverTrials(1,1,:,freq2plot-1:freq2plot+1,times2plot(1):times2plot(2),1)),5),4)),chanlocs,'maplimits',[0 clim],'plotrad',.63,'numcontour',0,'style','map','electrodes','off','emarker2',{elecs2use '.' 'g' 8 1});
title('seeded connectivity, average reference')

subplot(222)
topoplot(squeeze(mean(mean(abs(synchOverTrials(1,1,:,freq2plot-1:freq2plot+1,times2plot(1):times2plot(2),2)),5),4)),chanlocs,'maplimits',[0 clim],'plotrad',.63,'numcontour',0,'style','map','electrodes','off','emarker2',{elecs2use '.' 'g' 8 1});
title('seeded connectivity, Laplacian')

subplot(223)
contourf(times2save,frex,abs(squeeze(synchOverTrials(1,1,elecs2use(2),:,:,1))),40,'linecolor','none')
set(gca,'clim',[0 clim])
hold on
toplot = squeeze(mean(abs(synchOverTrials(1,1,elecs2use(2),freq2plot-1:freq2plot+1,:,1)),4));
plot(times2save,toplot*30+10,'w','linew',2)
plot(get(gca,'xlim'),[10 10],'w--')

subplot(224)
contourf(times2save,frex,abs(squeeze(synchOverTrials(1,1,elecs2use(2),:,:,2))),40,'linecolor','none')
set(gca,'clim',[0 clim])
hold on
toplot = squeeze(mean(abs(synchOverTrials(1,1,elecs2use(2),freq2plot-1:freq2plot+1,:,2)),4));
plot(times2save,toplot*30+10,'w','linew',2)
plot(get(gca,'xlim'),[10 10],'w--')
xlabel('Time (ms)'), ylabel('Frequency (Hz)')


figure(2), clf, colormap hot
subplot(221)
topoplot(squeeze(mean(mean(abs(synchOverTrials(2,1,:,freq2plot-1:freq2plot+1,times2plot(1):times2plot(2),1)),5),4)),chanlocs,'maplimits',[0 clim],'plotrad',.63,'numcontour',0,'style','map','electrodes','off','emarker2',{elecs2use '.' 'g' 8 1});
title('seeded connectivity, average reference')

subplot(222)
topoplot(squeeze(mean(mean(abs(synchOverTrials(2,1,:,freq2plot-1:freq2plot+1,times2plot(1):times2plot(2),2)),5),4)),chanlocs,'maplimits',[0 clim],'plotrad',.63,'numcontour',0,'style','map','electrodes','off','emarker2',{elecs2use '.' 'g' 8 1});
title('seeded connectivity, Laplacian')

subplot(223)
contourf(times2save,frex,abs(squeeze(synchOverTrials(2,1,elecs2use(2),:,:,1))),40,'linecolor','none')
hold on
toplot = squeeze(mean(abs(synchOverTrials(2,1,elecs2use(2),freq2plot-1:freq2plot+1,:,1)),4));
plot(times2save,toplot*30+10,'w','linew',2)
plot(get(gca,'xlim'),[10 10],'w--')
set(gca,'clim',[0 clim])

subplot(224)
contourf(times2save,frex,abs(squeeze(synchOverTrials(2,1,elecs2use(2),:,:,2))),40,'linecolor','none')
hold on
toplot = squeeze(mean(abs(synchOverTrials(2,1,elecs2use(2),freq2plot-1:freq2plot+1,:,1)),4));
plot(times2save,toplot*30+10,'w','linew',2)
plot(get(gca,'xlim'),[10 10],'w--')
set(gca,'clim',[0 clim])
xlabel('Time (ms)'), ylabel('Frequency (Hz)')

%% end data simulation/analyses part of script





%% code for figure 1





%% panel A

% number of time steps
n = 100;

t = [17 86]; % time points where polar distributions are drawn

figure(3), clf

phasedat = rand(n,1)*pi*(pi/10);
connres  = zeros(n,4);

for i=1:n
    phasedat = mod(phasedat + pi/20,2*pi);
    connres(i,1) = abs(mean(exp(1i*phasedat)));
    connres(i,2) = angle(mean(exp(1i*phasedat)));
    
    cdd = exp(1i*phasedat);
    cdi = imag(cdd);
    connres(i,3) = squeeze(mean( abs( mean( abs(cdi).*sign(cdi)) )./mean(abs(cdi))))';
    connres(i,4) = abs(mean(sign(imag(cdd))));
    
    
    if i==t(1)
        subplot(223)
        polar([zeros(n,1) phasedat]',[zeros(n,1) ones(n,1)]','k')
    elseif i==t(2)
        subplot(224)
        polar([zeros(n,1) phasedat]',[zeros(n,1) ones(n,1)]','k')
    end
end

subplot(211)
plot(connres(:,[1 3 4]))
xlabel('Time (a.u.)'), ylabel('ISPC or dwPLI')
set(gca,'ylim',[0 1.05],'ytick',0:.2:1)
legend({'ISPC';'wPLI';'PLI'})

hold on
plot([t(1) t(1)],get(gca,'ylim'),'k')
plot([t(2) t(2)],get(gca,'ylim'),'k')

%% panel B

t = [25 90]; % time points where polar distributions are drawn

figure(4), clf
connres2  = zeros(n,4);

for i=1:n
    
    phasedat = pi/2 + rand(n,1)*pi*(i/n) - (pi*(i/n))/2;

    connres2(i,1) = abs(mean(exp(1i*phasedat)));
    connres2(i,2) = angle(mean(exp(1i*phasedat)));
    
    cdd = exp(1i*phasedat);
    cdi = imag(cdd);
    connres2(i,3) = squeeze(mean( abs( mean( abs(cdi).*sign(cdi)) )./mean(abs(cdi))))';
    connres2(i,4) = abs(mean(sign(imag(cdd))));
    
    
    if i==t(1)
        subplot(223)
        polar([zeros(n,1) phasedat]',[zeros(n,1) ones(n,1)]','k')
    elseif i==t(2)
        subplot(224)
        polar([zeros(n,1) phasedat]',[zeros(n,1) ones(n,1)]','k')
    end
end

subplot(211)
plot(connres2(:,[1 3 4]))
xlabel('Time (a.u.)'), ylabel('ISPC or dwPLI')
set(gca,'ylim',[0 1.05],'ytick',0:.2:1)
legend({'ISPC';'wPLI';'PLI'})

hold on
plot([t(1) t(1)],get(gca,'ylim'),'k')
plot([t(2) t(2)],get(gca,'ylim'),'k')

%% end of script
 