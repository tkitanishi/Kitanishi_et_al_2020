% Detect SPW-Rs and save time stamp.
% 
% INPUT
%   session:    rat-day-sess name (ex. 'tk0062-171129-09')
%   mapfile:    channel map file (ex. 'Buzsaki256.mat')
%   ch:         channel used for ripple detection ([Shank# ChannelFromTopOfShank#])
%   pixel2cm:   pixel to cm conversion factor [cm/pixels]
% 
% OUTPUT
%   .rpl.evt:   ripple timing event file (for Neuroscope)
%   .mat:       ripple time stamp [sec]
% 
% Takum Kitanishi, Osaka City Univ, 2018/11/01
% 
% mostly taken from rippleDetection.m by Hiro Miyawaki


function rippleDetect(session,mapfile,ch,pixel2cm)

% test input
% session = 'tk0062-171129-09'
% mapfile = 'Buzsaki256.mat'
% ch   = [3 8]; % <- shank 3 channel 8 (from top of shank 3)
% pixel2cm = 0.5217;

% session = 'tk0067-180427-07'
% mapfile = 'Buzsaki256.mat';
% ch = [4 21];
% pixel2cm = 0.5217;


% Parameters_______________________________________________________________
idx = strfind(session,'-');

p.session  = session;
p.rat      = session(1:idx(1)-1);
p.day      = session(idx(1)+1:idx(2)-1);
p.sess     = session(idx(2)+1:end);
p.mapfile  = mapfile;
p.ch       = ch;
p.pixel2cm = pixel2cm;
p.lfppath  = fullfile('D:\data\rec\',p.rat,[p.day '-' p.sess]);
p.spkpath  = fullfile('D:\data\rec\',p.rat,[p.day '_sorted']);

% LFP sampling rate [Hz]
p.fs = 1250;

% ripple range [Hz]
p.f1 = 130; % stop band
p.f2 = 140; % pass band
p.f3 = 230; % pass band
p.f4 = 240; % stop band

p.windowLength = 11;     % ripple smoothing window length [samples]
p.lowTh    = 3;          % lower threshold for detecting ripples
p.highTh   = 7;          % higher threshold for detecting ripples
p.minInterRplInterval= 30;  % [ms]
p.minRplDuration     = 15;  % [ms]
p.maxRplDuration     = 300; % [ms]
p.halfWinLength = 0.200; % half window size [sec] (for visualization only)

p.selectState   = true;  % whether selecting brain state defined by TheStateEditor
p.state         = 3;     % TheStateEditor state (1=RUN,2=WAKE,3=SWS,4=REM)

p.selectImmobile= false; % whether selecting immobile periods   
p.vTh = 5;               % running speed more than this speed was excluded [cm/s]
p.winLength = 11;        % position smoothing window length [samples]

% save
p.evtFile  = fullfile('D:\data\rec\',p.rat,[p.day '-' p.sess],[session '.rpl.evt']);
p.savepath = fullfile('D:\data\rec\lfp\rippleSWS\imRefCA1');
p.matFile  = fullfile('D:\data\rec\lfp\rippleSWS\imRefCA1',[session '.mat']);

%__________________________________________________________________________

% make variables (removing p.)
pList = fieldnames(p);
for n=1:length(pList)
    eval(sprintf('%s=p.%s;',pList{n},pList{n}))
end

%% load lfp
[lfp,t,n] = loadLfp(lfppath,[],mapfile);

% ripple band
signal = fftbandpass(lfp((ch(1)-1)*32+ch(2),:),fs,f1,f2,f3,f4);


%% (optional) selecting brain state defined by TheStateEditor

if selectState
    fprintf('%s loading results from TheStateEditor\n',datestr(now));
    
    % get state start/stop [sec]
    load(fullfile(lfppath,[session '.SleepState.states.mat']),'SleepState');
    stateName = SleepState.idx.statenames{p.state};
    if isfield(SleepState.ints,[stateName 'state'])
        ints = SleepState.ints.([stateName 'state']);

        % analyze only in the selected state
        includes = false(size(signal));
        for ii=1:size(ints,1)
            includes(ints(ii,1)<t & t<ints(ii,2)) = true;
        end
    else
        disp('Indicated state does not exist in this session!')
        includes = false(size(signal));
        
        rpl.timestamps      = nan;
        rpl.peaks           = nan;
        rpl.peakNormedPower = nan;
        rpl.stdev           = nan;
        save(matFile,'rpl','p','-v7.3');
        
        return;
    end
else 
    fprintf('%s entire session included\n',datestr(now));
    
    % analyze the entire session
    includes = true(size(signal));
end


%% (optional) selecting immobile periods

if selectImmobile
    fprintf('%s speed threshold applied\n',datestr(now));
    
    % position
    [post,posx1,posy1,posx2,posy2,posx3,posy3] = loadWhl(lfppath,pixel2cm);
    fspos = 1/mean(diff(post));
    
    win = ones(p.winLength,1)./p.winLength;
    posx = conv(posx1,win,'same');
    posy = conv(posy1,win,'same');
    
    % speed (cm/s)
    v = sqrt(gradient(posx).^2+gradient(posy).^2)*fspos;
    v = conv(v,win,'same');
    
    % start/end time points of running periods
    vStart = find(diff(v>p.vTh)>0);
    vStop  = find(diff(v>p.vTh)<0);
    
    % Include 1st run if it is incomplete
    if vStart(1)>vStop(1)
        vStart = [1; vStart];
    end
    % Include last run period if it is incomplete
    if vStart(end)>vStop(end)
        vStop = [vStop; length(v)];
    end
    
    % exclude run periods
    includes2 = true(size(t));
    for ii=1:length(vStart)
        includes2(post(vStart(ii))<t & t<post(vStop(ii))) = false;
    end
else
    fprintf('%s speed threshold not applied\n',datestr(now));
    includes2 = true(size(signal));
end

% unify
includes = [includes & includes2];


%% Normalized ripple power
fprintf('%s nomalizing lfp\n',datestr(now));

squaredSignal = signal'.^2;
squaredSignal(~includes) = 0;

window = ones(windowLength,1)/windowLength;
smoothed = Filter0(window,squaredSignal);

sd  = std(smoothed(includes));
avg = mean(smoothed(includes));
smoothed(~includes) = 0;
normalized = (smoothed-avg)/sd;


%% Detect candidate ripples
fprintf('%s first detection lfp\n',datestr(now));

% Detect ripple periods by thresholding normalized squared signal
thresholded = normalized > lowTh;
start = find(diff(thresholded)>0);
stop = find(diff(thresholded)<0);

% Exclude last ripple if it is incomplete
if length(stop) == length(start)-1
	start = start(1:end-1);
end

% Exclude first ripple if it is incomplete
if length(stop)-1 == length(start)
    stop = stop(2:end);
end

% Correct special case when both first and last ripples are incomplete
if start(1) > stop(1)
	stop(1) = [];
	start(end) = [];
end
firstPass = [start,stop];
if isempty(firstPass)
	disp('  Detection by thresholding failed');
	return
else
	disp(['  After detection by thresholding: ' num2str(length(firstPass)) ' events.']);
end

%% Merge ripples if inter-ripple period is too short

% fprintf('%s marging detected ripples separated with <%0.1f ms intervals\n',datestr(now),minInterRplInterval);

minInterRplSamples = minInterRplInterval/1000*fs;
secondPass = [];
ripple = firstPass(1,:);
for i = 2:size(firstPass,1)
	if firstPass(i,1) - ripple(2) < minInterRplSamples
		% Merge
		ripple = [ripple(1) firstPass(i,2)];
	else
		secondPass = [secondPass ; ripple];
		ripple = firstPass(i,:);
	end
end
secondPass = [secondPass ; ripple];
if isempty(secondPass)
	disp('  Rpl merge failed');
	return
else
	disp(['  After ripple merge: ' num2str(length(secondPass)) ' events.']);
end


%% Discard ripples with a peak power < highTh
thirdPass = [];
peakNormalizedPower = [];
for i = 1:size(secondPass,1)
	[maxValue,~] = max(normalized(secondPass(i,1):secondPass(i,2)));
	if maxValue > highTh
		thirdPass = [thirdPass ; secondPass(i,:)];
		peakNormalizedPower = [peakNormalizedPower ; maxValue];
	end
end
if isempty(thirdPass)
	disp('Peak thresholding failed.');
	return
else
	disp(['  After peak thresholding: ' num2str(length(thirdPass)) ' events.']);
end


%% Detect negative peak position for each ripple
peakPosition = zeros(size(thirdPass,1),1);
for i=1:size(thirdPass,1)
	[~,minIndex] = min(signal(thirdPass(i,1):thirdPass(i,2)));
	peakPosition(i) = minIndex + thirdPass(i,1) - 1;
end


%% Discard ripples that are too long or too short
duration=diff(thirdPass,1,2)/fs*1000;
idx = duration<minRplDuration | maxRplDuration<duration;
thirdPass(idx,:) = [];
peakPosition(idx,:) = [];
peakNormalizedPower(idx,:) = [];
disp(['  After duration test: ' num2str(size(thirdPass,1)) ' events.']);


%% Save detection result

rpl.timestamps      = thirdPass/fs;
rpl.peaks           = peakPosition/fs;
rpl.peakNormedPower = peakNormalizedPower;
rpl.stdev           = sd;

save(matFile,'rpl','p','-v7.3');


%% Making an evt file

if ~isempty(evtFile)
    evtList = sortrows(...
    [rpl.timestamps(:,1),1*ones(size(rpl.timestamps,1),1);
     rpl.timestamps(:,2),2*ones(size(rpl.timestamps,1),1);
     rpl.peaks,3*ones(size(rpl.timestamps,1),1)]);
    disp([datestr(now) ' making evt file for ripple'])
    
    fid = fopen(evtFile,'w'); 
    evtType={'onset','offset','peak'};
    for n=1:size(evtList,1)
        fprintf(fid,'%f %s\n',...
        evtList(n,1)*1e3,evtType{evtList(n,2)});
    end
    fclose(fid);
end


%% Show ripple triggered average

% ripple-triggered average of lfp 
nWin = round(p.halfWinLength*p.fs); % half window size [samples]
taxis = (-nWin:nWin)/p.fs;          % [sec]

lfpTrig = zeros(size(lfp,1),nWin*2+1);
for ii=1:length(rpl.peaks)
    lfpTrig = lfpTrig + lfp(:,round(rpl.peaks(ii)*fs)-nWin:round(rpl.peaks(ii)*fs)+nWin);
end
lfpTrig = lfpTrig/length(rpl.peaks);

lfpTrig = lfpTrig-mean(lfpTrig,2); % baseline correction

% plot
figure(1); clf
scale = 80;
Pad = repmat( (0:-1:-size(lfp,1)+1)'*scale, 1,nWin*2+1);
for sh=1:8
    plot(taxis,lfpTrig(32*(sh-1)+1:32*sh,:)+Pad(32*(sh-1)+1:32*sh,:),'Color',hsv2rgb([(sh-1)/8 1 .7])); hold on
end
% plot selected channel
plot(taxis,lfpTrig(32*(ch(1)-1)+ch(2),:)+Pad(32*(ch(1)-1)+ch(2),:),'Color',[0 0 0],'LineWidth',1)

plot([0 0],[min(min(Pad))-500 max(max(Pad))+500],'k:')
title(['Ripple-triggered average of LFP:  ',p.session])
xlabel('Time from ripple negative peak (sec)')
xlim([-p.halfWinLength p.halfWinLength])
ylim([[min(min(Pad))-500 max(max(Pad))+500]])

set(gcf,'Position',[100 100 600 1200])
drawnow

print(figure(1),fullfile(p.savepath,p.session),'-dpng')
print(figure(1),fullfile(p.savepath,p.session),'-dsvg')

% figure(2); clf
% % plot selected channel
% plot(taxis,lfpTrig(32*(ch(1)-1)+ch(2),:),'Color',[0 0 0],'LineWidth',1)
% title(['Ripple-triggered average of LFP:  ',p.session])
% xlabel('Time from ripple negative peak (sec)')
% xlim([-p.halfWinLength p.halfWinLength])
% box off
% set(gcf,'Position',[100 100 600 300])
% print(figure(2),fullfile(p.savepath,p.session),'-depsc')


%% Show ripple position

if selectImmobile
    % find position index
    posInd = [];
    for ii=1:length(rpl.peaks)
        [~,posInd(ii)] = min(abs(post-rpl.peaks(ii)));
    end
    
    % plot
    figure(2); clf
    plot(posx(1<post & post<max(post)-1),posy(1<post & post<max(post)-1),'Color',[.5 .5 .5]); hold on
    plot(posx(posInd),posy(posInd),'r.')
    axis equal ij; box off
    xlabel('X position (cm)');
    ylabel('Y position (cm)');
    title(['Ripple position: ' p.session ' (' num2str(length(rpl.peaks)) ' events)'])
    set(gcf,'Position',[100 900 560 420])
    
    print(figure(2),fullfile(p.savepath,[p.session,'_pos']),'-dpng')
end


%% Functions
function y = Filter0(b,x)

if size(x,1) == 1
	x = x(:);
end

if mod(length(b),2)~=1
	error('filter order should be odd');
end

shift = (length(b)-1)/2;

[y0 z] = filter(b,1,x);

y = [y0(shift+1:end,:) ; z(1:shift,:)];

function [U,stdA] = unity(A,sd,restrict)

if ~isempty(restrict),
	meanA = mean(A(restrict));
	stdA = std(A(restrict));
else
	meanA = mean(A);
	stdA = std(A);
end
if ~isempty(sd),
	stdA = sd;
end

U = (A - meanA)/stdA;
