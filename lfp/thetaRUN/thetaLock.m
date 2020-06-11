% thetaLock(session,mapfile,ch)
%
% This code calculates theta phase locking during RUN sessions for individual
% units. The main processing is mostly from AxEEGlock.m.
%
% OUTPUT
%   phlock.
%       nspkT:          number of spikes
%       rateT:          firing rate [Hz]
%       resultVectLenT: resultant vector length
%       prefPhaseT:     preferred phase [deg]
%       pRayleighT:     p-value of Rayleigh test
%       zRayleighT:     z-value of Rayleigh test
%       normSpkCountT:  normalized spike count
%       p:              params
%
% Takuma Kitanishi, OCU, 2019/022/19

clear;

% mapfile = 'Buzsaki256.mat';
% refch   = [6 15];
% sessList = {
%     'tk0056-171014-02'
%     'tk0056-171014-04'
%     'tk0056-171014-06'
%     'tk0056-171014-08'};
% mainFunc(sessList,mapfile,refch)

% mapfile = 'Buzsaki256.mat';
% refch   = [6 15];
% sessList = {
%     'tk0056-171015-02'
%     'tk0056-171015-04'
%     'tk0056-171015-06'
%     'tk0056-171015-08'};
% mainFunc(sessList,mapfile,refch)

% mapfile = 'Buzsaki256.mat';
% refch   = [2 11];
% sessList = {
%     'tk0062-171129-02'
%     'tk0062-171129-04'
%     'tk0062-171129-06'
%     'tk0062-171129-08'};
% mainFunc(sessList,mapfile,refch)

% mapfile = 'Buzsaki256.mat';
% refch   = [2 11];
% sessList = {
%     'tk0062-171130-02'
%     'tk0062-171130-04'
%     'tk0062-171130-06'
%     'tk0062-171130-08'};
% mainFunc(sessList,mapfile,refch)
%  
% mapfile = 'Buzsaki256.mat';
% refch   = [5 14];
% sessList = {
%     'tk0064-180206-02'
%     'tk0064-180206-04'
%     'tk0064-180206-06'
%     'tk0064-180206-08'};
% mainFunc(sessList,mapfile,refch)

% mapfile = 'Buzsaki256.mat';
% refch   = [5 14];
% sessList = {
%     'tk0064-180207-02'
%     'tk0064-180207-04'
%     'tk0064-180207-06'
%     'tk0064-180207-08'};
% mainFunc(sessList,mapfile,refch)

% mapfile = 'Buzsaki256.mat';
% refch = [5 15];    % not the highest power refchannel
% sessList = {
%     'tk0067-180426-02'
%     'tk0067-180426-04'
%     'tk0067-180426-07'
%     'tk0067-180426-09'};
% mainFunc(sessList,mapfile,refch)
% 
% mapfile = 'Buzsaki256.mat';
% refch = [5 15];    % not the highest power refchannel
% sessList = {
%     'tk0067-180427-02'
%     'tk0067-180427-04'
%     'tk0067-180427-07'
%     'tk0067-180427-09'};
% mainFunc(sessList,mapfile,refch)

mapfile = 'A8x32Edge.mat';
refch   = [5 17];     % <-- irregular theta shape 
% refch   = [8 10];
sessList = {
    'tk0068-180530-02'
    'tk0068-180530-04'
    'tk0068-180530-07'
    'tk0068-180530-09'};
mainFunc(sessList,mapfile,refch)

mapfile = 'A8x32Edge.mat';
refch   = [5 17];     % <-- irregular theta shape 
% refch   = [8 10];
sessList = {
    'tk0068-180531-02'
    'tk0068-180531-04'
    'tk0068-180531-07'
    'tk0068-180531-09'};
mainFunc(sessList,mapfile,refch)

% mapfile = 'A8x32_5mm_35_300_160.mat';
% refch = [5 18];
% sessList = {
%     'tk0069-180726-02'
%     'tk0069-180726-04'
%     'tk0069-180726-07'
%     'tk0069-180726-09'};
% mainFunc(sessList,mapfile,refch)

% mapfile = 'A8x32_5mm_35_300_160.mat';
% refch = [5 18];
% sessList = {
%     'tk0069-180727-02'
%     'tk0069-180727-04'
%     'tk0069-180727-07'
%     'tk0069-180727-09'};
% mainFunc(sessList,mapfile,refch)

% mapfile = 'A8x32_5mm_35_300_160.mat';
% refch = [3 14];    % not the highest power refchannel
% sessList = {
%     'tk0070-180829-02'
%     'tk0070-180829-04'
%     'tk0070-180829-07'
%     'tk0070-180829-09'};
% mainFunc(sessList,mapfile,refch)

% mapfile = 'A8x32_5mm_35_300_160.mat';
% refch = [3 14];    % not the highest power refchannel
% sessList = {
%     'tk0070-180830-02'
%     'tk0070-180830-04'
%     'tk0070-180830-07'
%     'tk0070-180830-09'};
% mainFunc(sessList,mapfile,refch)
% 
% 
% mapfile = 'A8x32Edge.mat';
% refch = [5 14];
% sessList = {
%     'tk0072-181025-02'
%     'tk0072-181025-04'
%     'tk0072-181025-07'
%     'tk0072-181025-09'};
% mainFunc(sessList,mapfile,refch)
% 
% mapfile = 'A8x32Edge.mat';
% refch = [5 14];
% sessList = {
%     'tk0072-181026-02'
%     'tk0072-181026-04'
%     'tk0072-181026-07'
%     'tk0072-181026-09'};
% mainFunc(sessList,mapfile,refch)
% 
% mapfile = 'A8x32Edge.mat';
% refch = [1 19];    % not the highest power refchannel
% sessList = {
%     'tk0074-181219-02'
%     'tk0074-181219-04'
%     'tk0074-181219-07'
%     'tk0074-181219-09'};
% mainFunc(sessList,mapfile,refch)

% mapfile = 'A8x32Edge.mat';
% refch = [1 19];    % not the highest power refchannel
% sessList = {
%     'tk0074-181220-02'
%     'tk0074-181220-04'
%     'tk0074-181220-07'
%     'tk0074-181220-09'};
% mainFunc(sessList,mapfile,refch)
% 
% 
% mapfile = 'A8x32Edge.mat';
% refch = [4 14];
% sessList = {
%     'tk0075-190326-02'
%     'tk0075-190326-04'
%     'tk0075-190326-07'
%     'tk0075-190326-09'};
% mainFunc(sessList,mapfile,refch)
% 
% 
% mapfile = 'A8x32Edge.mat';
% refch = [4 14]; 
% sessList = {
%     'tk0075-190327-02'
%     'tk0075-190327-04'
%     'tk0075-190327-07'
%     'tk0075-190327-09'};
% mainFunc(sessList,mapfile,refch)

% mapfile = 'A8x32Edge.mat';
% refch = [6 8];
% sessList = {
%     'tk0076-190424-02'
%     'tk0076-190424-04'
%     'tk0076-190424-07'
%     'tk0076-190424-09'};
% mainFunc(sessList,mapfile,refch)

% mapfile = 'A8x32Edge.mat';
% refch = [6 8];
% sessList = {
%     'tk0076-190425-02'
%     'tk0076-190425-04'
%     'tk0076-190425-07'
%     'tk0076-190425-09'};
% mainFunc(sessList,mapfile,refch)


function mainFunc(sessList,mapfile,refch)


% Parameters_______________________________________________________________

p.savePath = 'D:\data\rec\lfp\thetaRUN';

% LFP sampling rate [Hz]
p.fs = 1250;

% theta range [Hz]
p.f1 = 4;   % stop band
p.f2 = 5;   % pass band
p.f3 = 10;  % pass band
p.f4 = 11;  % stop band

% angle bin size [deg]
p.angbin=30;
p.angaxis = p.angbin/2:p.angbin:360;

% plot max value (for visualization only)
p.maxPolar = 0.5;
p.maxHist = 0.3;
p.maxIm   = 0.25;

% make variables (removing p.)
pList = fieldnames(p);
for n=1:length(pList)
    eval(sprintf('%s=p.%s;',pList{n},pList{n}))
end
%__________________________________________________________________________


%% stack spike theta phases during RUN across session

spkPhaseAll = [];
cluAll = [];
DurationAll = 0;

for ii=1:length(sessList)
    session = sessList{ii};
    idx = strfind(session,'-');
    rat      = session(1:idx(1)-1);
    day      = session(idx(1)+1:idx(2)-1);
    sess     = session(idx(2)+1:end);
    
    lfpPath  = fullfile('D:\data\rec\',rat,[day '-' sess]);
    spkPath  = fullfile('D:\data\rec\',rat,[day '_sorted']);
    
    
    %% load lfp
    [lfp,t,n] = loadLfp(lfpPath,[],mapfile);
    
    % theta band
    lfpt = fftbandpass(lfp((refch(1)-1)*32+refch(2),:),p.fs,p.f1,p.f2,p.f3,p.f4);
    
    % theta phase
    [phaseT,~] = lfpphase(lfpt);  % [degree,~]
    
    
    %% get RUN period (TheStateEditor, state 1)
    load(fullfile(lfpPath,[session '.SleepState.states.mat']),'SleepState');
    
    % Select satate 1 of the TheStateEditor
    stateName = [SleepState.idx.statenames{1} 'state'];
    
    if isfield(SleepState.ints,stateName)
        % RUN period time stamp
        ints = SleepState.ints.(stateName);
        

        % load spikes
        [ts,clu,cluOri,Fs,roi,sh,ch]=loadSpkHist(spkPath,[],session);
        % session duration
        hdr = loadMeta(fullfile(lfpPath,[session '.meta']));
        Duration = hdr.filelength_sec;  % [sec]
        
        % get spikes in RUN periods
        disp([num2str(size(ints,1)) ' RUN periods in ' session])
        includes = false(size(ts));
        for jj=1:size(ints,1)
            includes(ints(jj,1)<ts & ts<ints(jj,2)) = true;
            DurationAll = DurationAll + (ints(jj,2)-ints(jj,1));
        end
        tsRUN = ts(includes);
        cluRUN = clu(includes);
        
        % spike theta phase
        spkPhaseRUN = spikePhase(tsRUN,phaseT,p.fs);
        
        nanIdx = isnan(spkPhaseRUN);
        spkPhaseRUN(nanIdx) = [];
        cluRUN(nanIdx) = [];
    
        % stack
        if isempty(spkPhaseAll)
            spkPhaseAll = spkPhaseRUN;
            cluAll      = cluRUN;
        else
            spkPhaseAll = [spkPhaseAll; spkPhaseRUN];
            cluAll      = [cluAll; cluRUN];
        end
    else
        error(['ERROR: No RUN period in: ' session])
    end
end


if length(spkPhaseAll)~=length(cluAll)
    error('ERROR: data length does not match!')
end


%% calculate phase locking statistics

for ii=1:max(clu)
    fprintf('  %s: processing unit %u\n',datestr(now),ii)
    %     tsTmp = ts(clu==ii);
    %
    %     spikePhaseT = spikePhase(tsTmp,phaseT,p.fs);
    %     spikePhaseT(isnan(spikePhaseT)) = [];
    
    spkPhaseTmp = spkPhaseAll(cluAll==ii);
    
    nspkT(ii) = length(spkPhaseTmp);
    rateT(ii) = nspkT(ii)./DurationAll;
    
    % no spikes
    if nspkT(ii)<1
        disp(['# of spikes to be analyzed is 0 for cluster: ',num2str(ii)]);
        rateT(ii)= 0;
        rT(ii)   = nan;
        angT(ii) = nan;
        pT(ii)   = nan;
        zT(ii)   = nan;
        nT(ii,:) = nan(1,length(angaxis));
        continue;
    end
    
    % Degree -> Radian spike phase
    spkRadPhaseT = circ_ang2rad(spkPhaseTmp);
    
    % Resultant vector length
    rT(ii) = circ_r(spkRadPhaseT);
    % angular deviation
    % stdT(ii) = circ_std(spkRadPhaseT);
    % Pairwise phase consistency
    ppcT(ii) = circ_ppc(spkRadPhaseT);
    % Resultant vector phase [deg]
    angT(ii) = mod(circ_rad2ang(circ_mean(spkRadPhaseT))+360,360);
    % Rayleigh test (circ_rtest)
    [pT(ii),zT(ii)] = circ_rtest(spkRadPhaseT);
    
    % calculate PPC (pairwise phase consistency)  -----------------?
    
    % normalized spike count histogram
    nT(ii,:) = hist(spkPhaseTmp,angaxis);
    nT(ii,:) = nT(ii,:)./nspkT(ii);
    
    % plot
    figure(1); clf
    subplot(121)
    [tout,rout] = rose(spkRadPhaseT,360/angbin);
    polar_lim(tout,rout./nspkT(ii),maxPolar); hold on
    h = compass_lim(rT(ii)*cos(circ_ang2rad(angT(ii))),rT(ii)*sin(circ_ang2rad(angT(ii))),maxPolar);
    set(h,'Color',[0 0 0],'LineWidth',2)
    title('Spike theta phase')
    
    subplot(122)
    if pT(ii)<0.05
        bar([angaxis angaxis+360],[nT(ii,:) nT(ii,:)],1,'FaceColor',[1 .5 .5])
        text(30,0.27,'p<0.05')
    else
        bar([angaxis angaxis+360],[nT(ii,:) nT(ii,:)],1,'FaceColor',[.5 .5 .5])
    end
    yticks([0:0.1:maxHist])
    axis([0 720 0 maxHist]);
    xlabel('Phase (deg)');
    ylabel('Normalized spike counts')
    title([session, ', clu ', num2str(ii)])
    figSetting
    
    set(gcf,'Position',[100 900 700 320]);
    drawnow
end


%% save
nclu = max(clu);

phlock.rat           = repmat({rat},1,nclu);
phlock.day           = repmat({day},1,nclu);
phlock.clu           = 1:nclu;
% phlock.goodidx       = goodidx;
phlock.roi           = roi;
phlock.sh            = sh;
phlock.ch            = ch;
phlock.nspkT         = nspkT;
phlock.rateT         = rateT;   % mean rate [Hz]
phlock.resultVectLenT= rT;
phlock.ppcT          = ppcT;
phlock.prefPhaseT    = angT;    % [deg]
phlock.pRayleighT    = pT;
phlock.zRayleighT    = zT;
phlock.normSpkCountT = nT;
phlock.p = p;

save(fullfile(p.savePath,[rat '-' day '.mat']),'phlock')


%% plot for this rat

% load unit classification
classifyFile = fullfile('D:\data\rec\classifyUnits',[rat '-' day '.mat']);
load(classifyFile,'Tclassify')


% select units
% (SUBpyr or CA1pyr) & locked & rate>0.1
unitIdx = (Tclassify.SUBp'==1 | Tclassify.CA1p'==1) & pT<0.05 & rateT>0.1;

roiTmp = roi(unitIdx);
nTtmp  = nT(unitIdx,:);
rTtmp  = rT(unitIdx);
angTtmp= angT(unitIdx);

nT1  = nTtmp(roiTmp==3,:);  % CA1 locked good unit
nT2  = nTtmp(roiTmp==1,:);  % SUB locked good unit
rT1  = rTtmp(roiTmp==3);
rT2  = rTtmp(roiTmp==1);
angT1= angTtmp(roiTmp==3);
angT2= angTtmp(roiTmp==1);
% locked1 = sum(roiTmp==3)/sum(roi==3 & goodidx & rateT>0.1)*100;
% locked2 = sum(roiTmp==1)/sum(roi==1 & goodidx & rateT>0.1)*100;

% sort nT
% [~,sortIdx1]=sort(sum(nT1(:,0<=angaxis & angaxis<180),2),'descend');
% [~,sortIdx2]=sort(sum(nT2(:,0<=angaxis & angaxis<180),2),'descend');

angT1tmp = angT1; angT1tmp(angT1tmp<180) = angT1tmp(angT1tmp<180)+360;
angT2tmp = angT2; angT2tmp(angT2tmp<180) = angT2tmp(angT2tmp<180)+360;
[~,sortIdx1]=sort(angT1tmp,'descend');
[~,sortIdx2]=sort(angT2tmp,'descend');

% for visulization
x = 0:1:360;
y = -cos(x/180*pi())/2+.5;

% plot
figure(2); clf
subplot(2,3,[1 2])
imagesc([angaxis angaxis+360],[],[nT1(sortIdx1,:) nT1(sortIdx1,:)])
caxis([0 maxIm])
xticks([0:180:720])
colorbar; colormap jet
ylabel('CA1 unit (p<0.05 only, sorted)')
title(['RUN theta locking: ' rat '-' day])
% title([session  ' (' num2str(round(locked1)) '% units locked)'])

subplot(2,3,[4 5])
imagesc([angaxis angaxis+360],[],[nT2(sortIdx2,:) nT2(sortIdx2,:)])
caxis([0 maxIm])
xticks([0:180:720])
colorbar; colormap jet
xlabel('Theta phase (deg, 0=trough)')
ylabel('SUB unit (p<0.05 only, sorted)')
% title([' (' num2str(round(locked2)) '% units locked)'])

subplot(2,3,3)
raxis = 0.05:0.1:1;
r1 = hist(rT1,raxis);
r2 = hist(rT2,raxis);
plot(raxis,[r1/sum(r1); r2/sum(r2)],'.-')
box off
xlim([0 1])
xlabel('Resultant vector length')
ylabel('Normalized count')
legend({'CA1','SUB'})
title('p<0.05 only')

subplot(2,3,6)
ang1 = hist(angT1,angaxis);
ang2 = hist(angT2,angaxis);
plot([angaxis angaxis+360],[ang1/sum(ang1) ang1/sum(ang1); ang2/sum(ang2) ang2/sum(ang2)],'.-'); hold on
plot([x x+360],[y y]/20,'k')
box off
xlim([0 720])
xticks([0:180:720])
xlabel('Preferred phase (deg)')
ylabel('Normalized count')
title('p<0.05 only')


set(gcf,'Position',[100 400 700 500])


print(figure(2),fullfile(p.savePath,[rat '-' day]),'-dpng')
% print(figure(2),fullfile(p.savepath,p.session),'-depsc')


end



%% Functions -----------------------------------------------------

% Seting figure parameters
function figSetting

set(gca,'Box','Off')
set(gca, 'TickDir', 'Out');
% set(gca, 'TickLength', [0.02 0.02]);
% set(gca,'FontSize',8);
set(gca,'XTick',[0 180 360 540 720])

end

% Extract EEG phase from the band pass filtered EEG signal by Hilbert
% transform. Troughs are set to 0/360 deg. Output is in deg
function [Phase, Amp] = lfpphase(Eeg)

% Hilbert transform
Z = hilbert(Eeg);

% Wave amplitude
Amp = abs(Z);

% EEG phase in rad
Phase = angle(Z);

% Rad to Degree (-180 to +180)
Phase = Phase / pi *180;

% Degree (0 to +360)
Phase = Phase + 180;

end

function spikehase = spikePhase(ts, eegPhase, Fs)
% spikePhase calculates phase of each spike according the eeg phase data.
% ts,        spike timing [sec]
% eegPhase,  eeg phase data generated by funtion thetaPhase. Length of
%            eegPhase is equal to that of EEG recording data, and its
%            values are between 0-360 deg.
% Fs,        sampling rate of eeg (Hz)
% spikePhase,spike phase of each spike. The length is equal to that of ts.

% Check spike number
if isempty(ts)
    spikehase = NaN;
    return;
end

% spike timing --> eeg index conversion
idx = round(ts./(1/Fs));

% find out-of-range index
idxStart = 1;
while 1
    if idx(idxStart) < 1
        idxStart = idxStart+1;
    else
        break;
    end
end
idxStop = length(idx);
while 1
    if idx(idxStop) > length(eegPhase)
        idxStop = idxStop-1;
    else
        break;
    end
end
idxRange = idxStart:idxStop;

% spike phase
spikehase = nan(length(ts),1);
spikehase(idxRange) = eegPhase(idx(idxRange));

end


function h = polar_lim(theta, rho, max_lim, varargin)
%
% 例題:
%    t = 0:.01:2*pi;
%    polar_lim(t,sin(2*t).*cos(2*t),1.5,'--r')

% max_lim のサイズで表示
x_fake=[0 max_lim 0 -max_lim];
y_fake=[max_lim 0 -max_lim 0];
h_fake=polar(x_fake,y_fake);
set(h_fake,'Visible','off')

% 通常のプロットを重ね書き
hold on
h = polar(theta,rho,varargin{:});
hold off

end


function h = compass_lim(x, y, max_lim)
%
% 例題:
%   x = eig(randn(20,20));
%   compass_lim(x,10);
%   %% or
%   compass_lim(real(x),imag(x),10);

if nargin == 2;
    max_lim = y;
    y = imag(x);
    x = real(x);
end

% max_lim のサイズで表示
x_fake=[0 max_lim 0 -max_lim];
y_fake=[max_lim 0 -max_lim 0];
h_fake=compass(x_fake,y_fake);
set(h_fake,'Visible','off')

% 通常のプロットを重ね書き
hold on
h = compass(x,y);
hold off

end


%% Circular functions --------------------------------------------

function r = circ_r(alpha, w, d)
% r = circ_r(alpha, w, d)
%   Computes mean resultant vector length for circular data.
%
%   Input:
%     alpha	sample of angles in radians
%     [w		number of incidences in case of binned angle data]
%     [d    spacing of bin centers for binned data, if supplied
%           correction factor is used to correct for bias in
%           estimation of r, in radians (!)]
%
%   Output:
%     r		mean resultant length
%
% PHB 7/6/2008
%
% References:
%   Statistical analysis of circular data, N.I. Fisher
%   Topics in circular statistics, S.R. Jammalamadaka et al.
%   Biostatistical Analysis, J. H. Zar
%
% Circular Statistics Toolbox for Matlab

% By Philipp Berens, 2009
% berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens/circStat.html

% check vector size
if size(alpha,2) > size(alpha,1)
    alpha = alpha';
end

if nargin<2
    % if no specific weighting has been specified
    % assume no binning has taken place
    w = ones(size(alpha));
else
    if size(w,2) > size(w,1)
        w = w';
    end
end

if nargin<3
    % per default do not apply correct for binned data
    d = 0;
end

% compute weighted sum of cos and sin of angles
r = w'*exp(1i*alpha);

% obtain length
r = abs(r)/sum(w);

% for data with known spacing, apply correction factor to correct for bias
% in the estimation of r (see Zar, p. 601, equ. 26.16)
if d ~= 0
    c = d/2/sin(d/2);
    r = c*r;
end

end


function [mu ul ll] = circ_mean(alpha, w)
%
% mu = circ_mean(alpha, w)
%   Computes the mean direction for circular data.
%
%   Input:
%     alpha	sample of angles in radians
%     [w		weightings in case of binned angle data]
%
%   Output:
%     mu		mean direction
%     ul    upper 95% confidence limit
%     ll    lower 95% confidence limit
%
% PHB 7/6/2008
%
% References:
%   Statistical analysis of circular data, N. I. Fisher
%   Topics in circular statistics, S. R. Jammalamadaka et al.
%   Biostatistical Analysis, J. H. Zar
%
% Circular Statistics Toolbox for Matlab

% By Philipp Berens, 2009
% berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens/circStat.html

% check vector size
if size(alpha,2) > size(alpha,1)
    alpha = alpha';
end

if nargin<2
    % if no specific weighting has been specified
    % assume no binning has taken place
    w = ones(size(alpha));
else
    if size(w,2) > size(w,1)
        w = w';
    end
end

% compute weighted sum of cos and sin of angles
r = w'*exp(1i*alpha);

% obtain mean by
mu = angle(r);

% confidence limits if desired
if nargout > 1
    t = circ_confmean(alpha,0.05,w);
    ul = mu + t;
    ll = mu - t;
end

end


function alpha = circ_ang2rad(alpha)

% alpha = circ_ang2rad(alpha)
%   converts values in degree to radians
%
% Circular Statistics Toolbox for Matlab

% By Philipp Berens, 2009
% berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens/circStat.html

alpha = alpha * pi /180;

end


function alpha = circ_rad2ang(alpha)

% alpha = circ-rad2ang(alpha)
%   converts values in radians to values in degree
%
% Circular Statistics Toolbox for Matlab

% By Philipp Berens, 2009
% berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens/circStat.html

alpha = alpha / pi *180;

end


function [pval z] = circ_rtest(alpha, w, d)
%
% [pval, z] = circ_rtest(alpha,w)
%   Computes Rayleigh test for non-uniformity of circular data.
%   H0: the population is uniformly distributed around the circle
%   HA: the populatoin is not distributed uniformly around the circle
%   Assumption: the distribution has maximally one mode and the data is
%   sampled from a von Mises distribution!
%
%   Input:
%     alpha	sample of angles in radians
%     [w		number of incidences in case of binned angle data]
%     [d    spacing of bin centers for binned data, if supplied
%           correction factor is used to correct for bias in
%           estimation of r, in radians (!)]
%
%   Output:
%     pval  p-value of Rayleigh's test
%     z     value of the z-statistic
%
% PHB 7/6/2008
%
% References:
%   Statistical analysis of circular data, N. I. Fisher
%   Topics in circular statistics, S. R. Jammalamadaka et al.
%   Biostatistical Analysis, J. H. Zar
%
% Circular Statistics Toolbox for Matlab

% By Philipp Berens, 2009
% berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens/circStat.html

if size(alpha,2) > size(alpha,1)
    alpha = alpha';
end

if nargin < 2
    r =  circ_r(alpha);
    n = length(alpha);
else
    if length(alpha)~=length(w)
        error('Input dimensions do not match.')
    end
    if nargin < 3
        d = 0;
    end
    r =  circ_r(alpha,w(:),d);
    n = sum(w);
end

% compute Rayleigh's R (equ. 27.1)
R = n*r;

% compute Rayleigh's z (equ. 27.2)
z = R^2 / n;

% compute p value using approxation in Zar, p. 617
pval = exp(sqrt(1+4*n+4*(n^2-R^2))-(1+2*n));

end


function [pval m] = circ_otest(alpha, sz, w)
%
% [pval, m] = circ_otest(alpha,sz,w)
%   Computes Omnibus or Hodges-Ajne test for non-uniformity of circular data.
%   H0: the population is uniformly distributed around the circle
%   HA: the population is not distributed uniformly around the circle
%
%   Alternative to the Rayleigh and Rao's test. Works well for unimodal,
%   bimodal or multimodal data. If requirements of the Rayleigh test are
%   met, the latter is more powerful.
%
%   Input:
%     alpha	sample of angles in radians
%     [sz   step size for evaluating distribution, default 1 degree
%     [w		number of incidences in case of binned angle data]

%   Output:
%     pval  p-value
%     m     minimum number of samples falling in one half of the circle
%
% PHB 3/16/2009
%
% References:
%   Biostatistical Analysis, J. H. Zar
%   A bivariate sign test, J. L. Hodges et al., 1955
%   A simple test for uniformity of a circular distribution, B. Ajne, 1968
%
% Circular Statistics Toolbox for Matlab

% By Philipp Berens, 2009
% berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens/circStat.html

if size(alpha,2) > size(alpha,1)
    alpha = alpha';
end

if nargin < 2 || isempty(sz)
    sz = circ_ang2rad(1);
end

if nargin < 3
    w = ones(size(alpha));
else
    if length(alpha)~=length(w)
        error('Input length does not match.')
    end
    w =w(:);
end

alpha = mod(alpha,2*pi);
n = sum(w);
dg = 0:sz:pi;

m = zeros(size(dg));
for i=1:length(dg)
    m(i) = sum((alpha > dg(i) & alpha < pi + dg(i)).*w);
end
m = min(m);

if n > 50
    % approximation by Ajne (1968)
    A = pi*sqrt(n) / 2 / (n-2*m);
    pval = sqrt(2*pi) / A * exp(-pi^2/8/A^2);
else
    % exact formula by Hodges (1955)
    pval = 2^(1-n) * (n-2*m) * nchoosek(n,m);
end

end


function [p U UC] = circ_raotest(alpha)

% [p U UC] = circ_raotest(alpha)
%   Calculates Rao's spacing test by comparing distances between points on
%   a circle to those expected from a uniform distribution.
%
%   H0: Data is distributed uniformly around the circle.
%   H1: Data is not uniformly distributed around the circle.
%
%   Alternative to the Rayleigh test and the Omnibus test. Less powerful
%   than the Rayleigh test when the distribution is unimodal on a global
%   scale but uniform locally.
%
%   Due to the complexity of the distributioin of the test statistic, we
%   resort to the tables published by
%       Russell, Gerald S. and Levitin, Daniel J.(1995)
%       'An expanded table of probability values for rao's spacing test'
%       Communications in Statistics - Simulation and Computation
%   Therefore the reported p-value is the smallest alpha level at which the
%   test would still be significant. If the test is not significant at the
%   alpha=0.1 level, we return the critical value for alpha = 0.05 and p =
%   0.5.
%
%   Input:
%     alpha     sample of angles
%
%   Output:
%     p         smallest p-value at which test would be significant
%     U         computed value of the test-statistic u
%     UC        critical value of the test statistic at sig-level
%
%
%   References:
%     Batschelet, 1981, Sec 4.6
%
% Circular Statistics Toolbox for Matlab

% By Philipp Berens, 2009
% berens@tuebingen.mpg.de

alpha = alpha(:);

% for the purpose of the test, convert to angles
alpha = circ_rad2ang(alpha);
n = length(alpha);
alpha = sort(alpha);

% compute test statistic
U = 0;
lambda = 360/n;
for j = 1:n-1
    ti = alpha(j+1) - alpha(j);
    U = U + abs(ti - lambda);
end

tn = (360 - alpha(n) + alpha(1));
U = U + abs(tn-lambda);

U = (1/2)*U;

% get critical value from table
[p UC] = getVal(n,U);

end


function [p UC] = getVal(N, U)

% Table II from Russel and Levitin, 1995

alpha = [0.001, .01, .05, .10];
table = [ 4   247.32, 221.14, 186.45, 168.02;
    5   245.19, 211.93, 183.44, 168.66;
    6   236.81, 206.79, 180.65, 166.30;
    7   229.46, 202.55, 177.83, 165.05;
    8   224.41, 198.46, 175.68, 163.56;
    9   219.52, 195.27, 173.68, 162.36;
    10  215.44, 192.37, 171.98, 161.23;
    11  211.87, 189.88, 170.45, 160.24;
    12  208.69, 187.66, 169.09, 159.33;
    13  205.87, 185.68, 167.87, 158.50;
    14  203.33, 183.90, 166.76, 157.75;
    15  201.04, 182.28, 165.75, 157.06;
    16  198.96, 180.81, 164.83, 156.43;
    17  197.05, 179.46, 163.98, 155.84;
    18  195.29, 178.22, 163.20, 155.29;
    19  193.67, 177.08, 162.47, 154.78;
    20  192.17, 176.01, 161.79, 154.31;
    21  190.78, 175.02, 161.16, 153.86;
    22  189.47, 174.10, 160.56, 153.44;
    23  188.25, 173.23, 160.01, 153.05;
    24  187.11, 172.41, 159.48, 152.68;
    25  186.03, 171.64, 158.99, 152.32;
    26  185.01, 170.92, 158.52, 151.99;
    27  184.05, 170.23, 158.07, 151.67;
    28  183.14, 169.58, 157.65, 151.37;
    29  182.28, 168.96, 157.25, 151.08;
    30  181.45, 168.38, 156.87, 150.80;
    35  177.88, 165.81, 155.19, 149.59;
    40  174.99, 163.73, 153.82, 148.60;
    45  172.58, 162.00, 152.68, 147.76;
    50  170.54, 160.53, 151.70, 147.05;
    75  163.60, 155.49, 148.34, 144.56;
    100 159.45, 152.46, 146.29, 143.03;
    150 154.51, 148.84, 143.83, 141.18;
    200 151.56, 146.67, 142.35, 140.06;
    300 148.06, 144.09, 140.57, 138.71;
    400 145.96, 142.54, 139.50, 137.89;
    500 144.54, 141.48, 138.77, 137.33;
    600 143.48, 140.70, 138.23, 136.91;
    700 142.66, 140.09, 137.80, 136.59;
    800 142.00, 139.60, 137.46, 136.33;
    900 141.45, 139.19, 137.18, 136.11;
    1000  140.99, 138.84, 136.94, 135.92  ];

ridx = find(table(:,1)>=N,1);
cidx = find(table(ridx,2:end)<U,1);

if ~isempty(cidx)
    UC = table(ridx,cidx+1);
    p = alpha(cidx);
else
    UC = table(ridx,end-1);
    p = .5;
end

end
