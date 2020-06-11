% ampxSpace(basepath,trange)
% 
% The code is modified from AxInfo2.m/AxSpace.m
% 
% INPUT
%   path: path to the folder including npy/whl files.
%   trange: time range for data analysis.
% 
% OUTPUT
%   Grid score (g):
%   Mean vector length (v): if positions are tracked with 2 LEDs
%   Border score (b):
%   Theta modulation (t):   theta modulated if t>5 (Boccara NN 2010)
% 
% *** Input arguments ***
%
% smoothingChoice   Set this to one (1) if you want the program to calculate
%                   the optimal smoothing factor for each cell and use this
%                   for the ratemap. Set it to zero (0) to use the
%                   smoothing value set in the beginning of the code in the
%                   program.
%
% Takuma Kitanishi, Osaka Univ 2015 / Osaka City Univ 2018

function spaceT()

fclose all;

% rec only
io('tk0056-171014-06',[0 1200]+20);
io('tk0056-171015-06',[0 1200]+15);

io('tk0062-171129-06',[0 1200]+6);
io('tk0062-171130-06',[0 1200]+5);

% with optoid
io('tk0064-180206-06',[0 1200]+5);
io('tk0064-180207-06',[0 1200]+2);

io('tk0067-180426-07',[0 1200]+5);
io('tk0067-180427-07',[0 1200]+7);

io('tk0068-180530-07',[0 1200]+5);
io('tk0068-180531-07',[0 1200]+5);

io('tk0069-180726-07',[0 1200]+8);
io('tk0069-180727-07',[0 1200]+5);

io('tk0070-180829-07',[0 1200]+5);
io('tk0070-180830-07',[0 1200]+2);

io('tk0072-181025-07',[0 1200]+2);
io('tk0072-181026-07',[0 1200]+5);

io('tk0074-181219-07',[0 1200]+2);
io('tk0074-181220-07',[0 1200]+5);

io('tk0075-190326-07',[0 1200]+5);
io('tk0075-190327-07',[0 1200]+5);

io('tk0076-190424-07',[0 1200]+5);
io('tk0076-190425-07',[0 1200]+5);

spaceT_stack()


function io(session,trange)

% % test input
% session = 'tk0067-180427-07';
% trange = [10 1210];


% setting
idx = strfind(session,'-');
rat = session(1:idx(1)-1);
day = session(idx(1)+1:idx(2)-1);
sess= session(idx(2)+1:end);

% Parameters_______________________________________________________________

p.basepath = fullfile('D:\data\rec\',rat,[day '-' sess]);
p.spkpath = fullfile('D:\data\rec\',rat,[day '_sorted']);
p.savepath = 'D:\data\rec\06T\space';

% pixel to cm conversion (cm/pixels)
p.pixel2cm = 0.5217;

% Bin width for the ratemap (Fyhn Hippocampus 2008 -> 1.5)
p.binWidth = 2.5;   % [cm]

% Smoothing factor for rate maps (Fyhn Hippocampus 2008 -> 3.0)
p.smoothing = 5;    % [cm]

% box size
p.boxsize = 120;     % [cm]

% Minimum number of bins in a placefield. Fields with fewer bins than this
% treshold are not considered as placefields (Fyhn Hippocampus 2008)
p.minNumBins = 10;

% Bins with rate at 0.2 * peak rate and higher will be considered as part
% of a place field (Fyhn Hippocampus 2008)
p.fieldTreshold = 0.2;

% Lowest field rate in Hz.
p.lowestFieldRate = 1;

% video sampling rate (Hz)
p.fsv = 39.0625;

% save ratemap images (=1) or not (=0)
p.storeImage = 1;

% show histological roi name (=1) or not (=0)
p.roi = 1;
%__________________________________________________________________________


% path to the data folder
fprintf('%s%s\n','Session name: ',session);

% chnage recursive limit
set(0,'RecursionLimit',1156)

% Open the output file and write header
fid = fopen(fullfile(p.savepath,[session '.xls']),'w');

fprintf(fid,'%s%s\n','spaceS.m ', datestr(now)');
fprintf(fid,'%s\t','Session');
fprintf(fid,'%s\t','Unit');
fprintf(fid,'%s\t','ROI');
fprintf(fid,'%s\t','ProxDist');
fprintf(fid,'%s\t','SupDeep');
fprintf(fid,'%s\t','Shank');
fprintf(fid,'%s\t','Channel');
fprintf(fid,'%s\t','NumSpikes');
fprintf(fid,'%s\t','MeanRateSession');
fprintf(fid,'%s\t','PeakRate_Spatial');
fprintf(fid,'%s\t','Field');
fprintf(fid,'%s\t','FieldX');
fprintf(fid,'%s\t','FieldY');
fprintf(fid,'%s\t','FieldSize_cm2');
fprintf(fid,'%s\t','AvgRateField');
fprintf(fid,'%s\t','PeakRateField');
fprintf(fid,'%s\t','BurstPrct');
fprintf(fid,'%s\t','Sparseness');
fprintf(fid,'%s\t','Information_BitsPerSpk');
fprintf(fid,'%s\t','Selectivity');
fprintf(fid,'%s\t','Zcoherence');
fprintf(fid,'%s\t','GridScore');
fprintf(fid,'%s\t','DirPeakRate');
fprintf(fid,'%s\t','MeanVectLength');
fprintf(fid,'%s\t','MeanPhase');
fprintf(fid,'%s\t','DirInfo_BitsPerSpk');
fprintf(fid,'%s\t','BorderScore');
fprintf(fid,'%s\t','SpeedScore');
fprintf(fid,'%s\t','ThetaModulation');
fprintf(fid,'%s\n','ModulationFreq');

% load position
[t,x1,y1,x2,y2] = loadWhl(p.basepath,p.pixel2cm);
% flip y data
y1 = -y1;
y2 = -y2;
% trim pos data
if trange(1)<t(1) || t(end)<trange(2)
    error 'Error: trange setting is out of range!'
end
idx = trange(1)<=t & t<trange(2);
x1(~idx)=[];
y1(~idx)=[];
x2(~idx)=[];
y2(~idx)=[];
t(~idx)=[];

% Use the better-tracking LED
% (For head direction analysis, t,x1,y1,x2,y2 are used)
if sum(x1==-1)>sum(x2==-1)
    posx = x2;
    posy = y2;
else
    posx = x1;
    posy = y1;
end
post = t;

% % Treshold for how far a rat can move (100cm/s), in one sample (sampFreq =
% % 39.0625 Hz)
% treshold = 100/p.fsv;
% % Remove bad tracking from the position samples
% [posx,posy,post] = remBadTrack(posx,posy,post,treshold);
% [x1,y1,x2,y2,t] = remBadTrack2Diode(x1,y1,x2,y2,t,treshold);

% Remove NaN from the position samples if present
indNaN = isnan(posx);
% Copy of the postion timestamps before illegal times are removed
cPost = post;
posx(indNaN) = [];
posy(indNaN) = [];
post(indNaN) = [];

% Remove NaN from the position samples if present
indNaN = isnan(x1) | isnan(x2);
% Copy of the postion timestamps before illegal times are removed
ct = t;
x1(indNaN) = [];
y1(indNaN) = [];
x2(indNaN) = [];
y2(indNaN) = [];
t(indNaN) = [];

% Smooth the position samples using a moving mean filter
[posx,posy] = posMeanFilter(posx,posy);
[x1,y1] = posMeanFilter(x1,y1);
[x2,y2] = posMeanFilter(x2,y2);

% Centre the box in the coordinate system
centre = centreBox(posx,posy);
posx = posx - centre(1);
posy = posy - centre(2);
x1 = x1 - centre(1);
y1 = y1 - centre(2);
x2 = x2 - centre(1);
y2 = y2 - centre(2);

% Observed diameter/sidelength of the box
obsLength = max(max(posx)-min(posx),max(posy)-min(posy));
bins = ceil(obsLength/p.binWidth);
sLength = p.binWidth*bins;
mapAxis = (-sLength/2+p.binWidth/2):p.binWidth:(sLength/2-p.binWidth/2);

% Calulate what areas of the box that have been visited
visited = visitedBins(posx,posy,mapAxis);

% load spike timing and cluster index
if p.roi
    [ts,clu,~,~,roi,sh,ch] = loadSpkHist(p.spkpath,[],session);
    proxdist = loadProxDist(rat,1);
    supdeep = loadSupDeep(rat,1);
else
    [ts,clu]=loadSpk(p.spkpath,session);
end
idx = trange(1)<=ts & ts<trange(2);
ts(~idx)=[];
clu(~idx)=[];

% call the main func for each cluster
u = struct;
for unit=1:max(clu)
    if p.roi
        roiName = getRoiName(roi(unit));
        if strcmp(roiName,'SUB')
            ProxDist = proxdist(sh(unit));
            SupDeep = supdeep(ch(unit),sh(unit));
        else
            ProxDist = '-';
            SupDeep = '-';
        end
        u = mainFunc(posx,posy,post,cPost,x1,y1,x2,y2,t,ct,ts(clu==unit),roiName,ProxDist,SupDeep,sh(unit),ch(unit),session,unit,mapAxis,fid,visited,p,u);
    else
        u = mainFunc(posx,posy,post,cPost,x1,y1,x2,y2,t,ct,ts(clu==unit),'-','-','-',nan,nan,session,unit,mapAxis,fid,visited,p,u);
    end
end

% save
save(fullfile(p.savepath,session),'session','u','mapAxis','visited','p');

% Close the files
fclose(fid);
disp('Done!')

%__________________________________________________________________________
%
%                           Main function
%__________________________________________________________________________

function u = mainFunc(posx,posy,post,cPost,x1,y1,x2,y2,t,ct,ts,roiName,proxdist,supdeep,sh,ch,session,unit,mapAxis,fid2,visited,p,u)

fprintf('%s%u\n','Analyzing unit ',unit);
u(unit).roiName = roiName;
u(unit).proxdist = proxdist;
u(unit).supdeep = supdeep;
u(unit).sh = sh;
u(unit).ch = ch;

% number of bins
bins = length(mapAxis);

% Construct image file names
imname = fullfile(p.savepath,sprintf('%s%s%u',session,'_u',unit));
psname = fullfile(p.savepath,session);

if isempty(ts) % No spikes for this cell
    % Make an empty map
    map = zeros(bins,bins);
    map(visited==0) = NaN;
    u(unit).map = map;
    
    nFields = 0;
    meanRate = 0;
    peakRate = 0;
    bursting = NaN;
    sparseness = NaN;
    information = NaN;
    selectivity = NaN;
    zCoherence = NaN;
    g = NaN;
    hdPeakRate = 0;
    v = NaN;
    ph = NaN;
    dirInfo = NaN;
    b = NaN;
    s = NaN;
    tm = NaN;
    maxF = NaN;
    
else
    
    % Plot empty map and store it to file
    %     figure(2)
    % %     drawfield(map,mapAxis,'jet',max(max(map)));
    %     drawfield(map,mapAxis,'jet',max(max(map)),p.binWidth,p.smoothing);
    %     axis off;
    %     f = getframe(gcf);
    %     [pic, cmap] = frame2im(f);
    %     if p.storeImage
    %         imwrite(pic,bmpImage,'bmp');
    %     end
    %     print('-depsc2',epsImage);
    %     print('-depsc', '-tiff',epsImage);
    
    
    %% ----- General spatial params -----
    
    % Average rate for the whole session
    meanRate = length(ts)/(max(cPost)-min(cPost));
    % spike positions, removing spikes from areas where the tracking is bad
    [spkx,spky,ts,spkInd] = spikePos(ts,posx,posy,post,cPost);
    
    % Calculate the spike firing map
    figure(1); clf; subplot(241)
    plot(posx,posy,'k',spkx,spky,'.r');
    set(gcf,'color',[1 1 1]);
    set(gca,'TickDir','out')
    axis equal
    box off
    title(sprintf('%s%s%u%s%s%s',session,'u',unit,' (',roiName,')'))
    ylabel('Distance (cm)')
    
    % Calculate the rate map and position probability density function
    [map,pospdf] = ratemap(spkx,spky,posx,posy,post,p.smoothing,mapAxis,mapAxis);
    
    % Set unvisited parts of the box to NaN
    map(visited==0) = NaN;
    
    %     % field properties
    [nFields,fieldProp] = placefield(map,p,mapAxis);
    
    % z-coherence
    zCoherence = fieldcohere(map);
    
    % Peak rate of rate map
    peakRate = nanmax(nanmax(map));
    
    pospdf(visited==0) = 0;
    % Calculate the information (shannon), sparseness and selectivity
    [information,sparseness,selectivity] = mapstat(map,pospdf);
    
    % Calculate the percentage of bursting
    [bursts,singleSpikes] = burstfinder(ts,0.010);  % 10 ms criterion
    bursting = length(bursts)/(length(bursts)+length(singleSpikes)) *100;
    
    % draw ratemap
    figure(1); subplot(245)
    drawfield(map,mapAxis,'jet',max(map(:)));
    title(sprintf('%s%3.1f%s%3.1f%s%u%s','Mean=',meanRate,'Hz Peak=',peakRate,'Hz (0-',ceil(max(map(:))),'Hz)'))
    box off
    ylabel('Distance (cm)')
    set(gca,'TickDir','out')
    
    % register
    u(unit).map = map;
    
    
    %% ----- head direction (Langston Science 2010) -----
    % based on AxHeadCellClass.m
    
    % head direction (0 to 360 deg)
    direct = mod( atan2(y2-y1,x2-x1)*180/pi + 360, 360);
    direct(direct==360) = 0;
    
    % Map spike time stamps to postion time stamps
    spkIndHd = getSpkInd(ts,t,ct);
    
    % Find the direction of the rat at the spike times
    spkDir = direct(spkIndHd);
    clear spkIndHd;
    
    % spike phase distribution (0.5 deg bins)
    hdAxis = 0.25:0.5:360;
    posDirBin = hist(direct,hdAxis);
    spkDirBin = hist(spkDir,hdAxis);
    % Convert posDirBin from number of samples to time
    posDirBin = posDirBin ./ p.fsv;
    
    % circular smoothing (14.5 degrees, 14 bins on each side)
    posDirBin = hdsmooth(posDirBin);
    spkDirBin = hdsmooth(spkDirBin);
    
    % Calculate the raw rate of each bin
    rateBin = spkDirBin ./ posDirBin;
    
    % directional peak rate [Hz]
    hdPeakRate = max(rateBin);
    
    % mean vector length (v)
    v = circ_r( hdAxis(:)*pi/180,rateBin );
    fprintf('%s%3.3f\n','  Mean vector length: ',v);
    
    % mean firing phase
    ph = mod(circ_mean(hdAxis(:)*pi/180,rateBin)*180/pi +360, 360);
    
    % position direction probability density function (sum == 1)
    posDirPdf = posDirBin / sum(posDirBin);
    
    % directional information score (Langston Science 2010)
    dirInfo = 0;
    for ii = 1:length(rateBin)
        if rateBin(ii) > 0
            dirInfo = dirInfo + posDirPdf(ii) * rateBin(ii)/meanRate * log2(rateBin(ii)/meanRate);
        end
    end
    
    % directional stability
    % Langston Science 2010
    
    % plot
    figure(1); subplot(242)
    h = polar_lim(hdAxis(:)*pi/180, rateBin, ceil(hdPeakRate)); hold on
    plot([0 cos(ph/180*pi)*ceil(hdPeakRate)],[0 sin(ph/180*pi)*ceil(hdPeakRate)],'r')
    set(h,'Color',[0 0 0],'LineWidth',2)
    title(sprintf( '%s%3.1f%s%4.1f%s%1.3f','Peak=',hdPeakRate,'Hz Ph=',ph,'deg VecLen=',v ))
    
    % register
    u(unit).hdRateBin = rateBin;
    u(unit).hdAxis = hdAxis;
    
    
    %% ----- grid socre (Langston Science 2010) -----
    
    % autocorrelation for the rate map
    Rxx = correlation(map,map);
    % set the axis for the correlation map
    corrAxis = p.binWidth * linspace(-((size(Rxx,1)-1)/2),((size(Rxx,1)-1)/2),size(Rxx,1));
    
    % determine the radius of center peak [cm] (Stensola & Moser Nature 2012)
    centreRadius = centreradius(Rxx,corrAxis,p.binWidth,p.boxsize);
    
    % calculate grid score by changin the radius (Langstron Science 2010)
    g = gridscore(Rxx,corrAxis,p.binWidth,centreRadius+10,p.boxsize-10);
    fprintf('%s%3.3f\n','  Grid score: ',g);
    
    % Draw the Rxx
    subplot(246)
    drawfield(Rxx+1,corrAxis,'jet',max(max(Rxx+1)));
    set(gca,'TickDir','out')
    box off
    title(sprintf('%s%2.2f','Grid score = ',g));
    ylabel('Distance (cm)')
    
    % register
    u(unit).Rxx = Rxx;
    u(unit).corrAxis = corrAxis;
    
    
    %% ----- border score (Solstad, Science, 2008) -----
    
    % putative border field
    % nBorder:              number of putative border fields
    % xy(ii).x, xy(ii).y:   bin coordinates of ii-th fields
    [nBorder,xy] = borderfield(map,mapAxis);
    
    if nBorder
        % cm, maximum coverage of any single single field over any of the four
        % walls of the environment
        cm = cmCalc(map,nBorder,xy);
        % dm, mean firing distance
        dm = dmCalc(map,nBorder,xy);
        
        % border score
        b = (cm-dm)/(cm+dm);
    else
        b = NaN;
    end
    fprintf('%s%3.3f\n','  Border score: ',b);
    
    
    %% ----- Speed tuning (Kropff Nature 2015) -----
    
    % 250-ms gaussian window smoothing
    w = gausswin(round(0.250*p.fsv));
    w = w./sum(w);
    
    % instantaneous speed
    spInst = sqrt(diff(posx).^2 + diff(posy).^2) ./ diff(post);
    spInst = [spInst(1); spInst];
    spInst = conv(spInst,w,'same');
    
    % instantaneous firing rate
    rateInst = histcounts(spkInd,0.5:length(post)+0.5) * p.fsv;
    rateInst = conv(rateInst,w,'same');
    
    % speed tuning curve
    spEdges = 0:5:50;
    spAxis = (spEdges(1:end-1)+spEdges(2:end))/2;
    for ii=1:length(spAxis)
        idx = spEdges(ii)<=spInst & spInst<spEdges(ii+1);
        spRate(ii) = mean(rateInst(idx));
    end
    
    % speed socre (calculated from points 2-50cm/s)
    idx = 2<spInst & spInst<50;
    s = corrcoef(spInst(idx),rateInst(idx));
    s = s(2);
    fprintf('%s%3.3f\n','  Speed score: ',s);
    
    % plot
    subplot(244)
    plot(spAxis,spRate,'-')
    xlabel('Speed (cm/s)'); ylabel('Firing rate (Hz)')
    title(sprintf('%s%.2f','SpeedScore=',s))
    xlim([min(spEdges) max(spEdges)]);
    ylim([0 ceil(max(spRate))+1]);
    box off
    
    % register
    u(unit).spAxis = spAxis;
    u(unit).spRate = spRate;
    
    
    %% ----- Theta modulation (Langston Science 2010) -----
    
    % sort spike timing into 2ms bins
    spkHist = hist( ts, 0.001:0.002:(max(cPost)-min(cPost)) );
    % autocorrelogram of spike timing (0-500 ms)
    [ac,lags] = xcorr(spkHist,500/2);
    ac(lags==0) = max(ac(lags~=0));
    ac = ac(251:end);
    
    % plot auto correlogram
    subplot(247)
    bar(0:0.002:0.5, ac, 1)
    set(gcf,'color',[1 1 1]);
    xlim([0 0.5]); ylim([0 Inf])
    box off
    xlabel('Time (sec)')
    ylabel('Spike counts')
    title('Spike autocorrelogram')
    
    % normalize autocorrelogram
    ac = ac - mean(ac);
    
    % aplly Hamming window
    w = hamming(length(ac));
    ac = ac(:) .* w(:);
    
    % FFT of autocorrelogram
    Fs = 500;                           % Sampling frequency [Hz]
    L = length(ac);                     % Length of signal
    nfft = 2^16;                        % length of FFT
    f = Fs/2*linspace(0,1,nfft/2+1);    % freqency vector [Hz]
    FFT = fft(ac,nfft)/L;
    fftPow = abs(FFT).^2;
    fftPow = fftPow(1:nfft/2+1);        % 0-250 Hz
    % trim to 0-125 Hz
    idx = 0<=f & f<=125;
    fftPow(~idx) = [];
    f(~idx)      = [];
    
    % peak frequency in the theta (4-11 Hz) frequency range
    idx = 4<f & f<11;
    [maxPow,maxIdx] = max(fftPow(idx));
    fTmp = f(idx);
    maxF = fTmp(maxIdx);
    
    % theta modulation
    % Giocomo Cell 2011 (theta modulated = t>3)
    % Boccara NN 2010   (theta modulated = t>5)
    idx = maxF-1<f & f<maxF+1;
    tm = mean(fftPow(idx)) / mean(fftPow);
    fprintf('%s%3.3f\n','  Theta modulation: ',tm);
    
    % Plot single-sided amplitude spectrum
    subplot(248)
    plot(f,fftPow); hold on
    xlim([0 20])
    box off
    plot(maxF,maxPow,'ro')
    title(sprintf('%s%4.2f','Theta modulation = ',tm))
    xlabel('Frequency (Hz)')
    ylabel('FFT power')
    
    
    %% ----- Theta skipping (Brandon & Hasselmo NN 2013) -----
    
    %     % Collect all inter-spike intervals <=400ms
    %     % and construct 10-ms-bin autocorrelogram
    %     acAxis2 = (5:10:400)/1000;
    %     ac2 = zeros(size(acAxis2));
    %     for ii=1:length(ts)-1
    %         idx = find(ts-ts(ii)<=0.400,1,'last');
    %         ac2 = ac2 + hist(ts(ii+1:idx)-ts(ii),acAxis2);
    %     end
    %     ac2(1) = min([ac2(1) max(ac2(2:end))]);
    % %     bar(acAxis2,ac2)
    
    %     curve fitting
    %     http://jp.mathworks.com/help/optim/examples/nonlinear-data-fitting.html
    
    
    %% save image
    
    drawnow
    if p.storeImage
        set(gcf,'Position',[80 918 1100 420])
        print(gcf,imname,'-dpng');
        print(gcf,imname,'-dsvg');
    end

end

%% Write results to file
fprintf(fid2,'%s\t', session);
fprintf(fid2,'%u\t', unit);
fprintf(fid2,'%s\t', roiName);
fprintf(fid2,'%u\t', proxdist);
fprintf(fid2,'%u\t', supdeep);
fprintf(fid2,'%u\t', sh);
fprintf(fid2,'%u\t', ch);
fprintf(fid2,'%u\t', length(ts));   u(unit).nSpike = length(ts);
fprintf(fid2,'%3.4f\t',meanRate);   u(unit).meanRate = meanRate;
fprintf(fid2,'%3.4f\t',peakRate);   u(unit).peakRate = peakRate;
if nFields==0
    fprintf(fid2,'%u\t', 0);        u(unit).field = 0;
    fprintf(fid2,'%3.4f\t',NaN);    u(unit).field = nan;
    fprintf(fid2,'%3.4f\t',NaN);    u(unit).field = nan;
    fprintf(fid2,'%5.4f\t',0);      u(unit).fieldSize = 0;
    fprintf(fid2,'%3.4f\t',NaN);    u(unit).fieldAvgRate = nan;
    fprintf(fid2,'%3.4f\t',NaN);    u(unit).fieldPeakRate= nan;
else
    ii=1;
    fprintf(fid2,'%u\t', ii);                       u(unit).field = ii;
    fprintf(fid2,'%3.4f\t',fieldProp(ii).x);        u(unit).fieldX = fieldProp(ii).x;
    fprintf(fid2,'%3.4f\t',fieldProp(ii).y);        u(unit).fieldY = fieldProp(ii).y;
    fprintf(fid2,'%5.4f\t',fieldProp(ii).size);     u(unit).fieldSize = fieldProp(ii).size;
    fprintf(fid2,'%3.4f\t',fieldProp(ii).avgRate);  u(unit).fieldAvgRate = fieldProp(ii).avgRate;
    fprintf(fid2,'%3.4f\t',fieldProp(ii).peakRate); u(unit).fieldPeakRate= fieldProp(ii).peakRate;
end
fprintf(fid2,'%3.4f\t',bursting);   u(unit).burstPrct = bursting;
fprintf(fid2,'%3.4f\t',sparseness); u(unit).sparseness = sparseness;
fprintf(fid2,'%3.4f\t',information);u(unit).information = information;
fprintf(fid2,'%3.4f\t',selectivity);u(unit).selectivity = selectivity;
fprintf(fid2,'%3.4f\t',zCoherence); u(unit).zCoherence = zCoherence;
fprintf(fid2,'%3.4f\t',g);          u(unit).gridScore = g;
fprintf(fid2,'%3.4f\t',hdPeakRate); u(unit).hdPeakRate = hdPeakRate;
fprintf(fid2,'%3.4f\t',v);          u(unit).hdVectLen = v;
fprintf(fid2,'%3.4f\t',ph);         u(unit).hdMeanPhase = ph;
fprintf(fid2,'%3.4f\t',dirInfo);    u(unit).hdDirInfo = dirInfo;
fprintf(fid2,'%3.4f\t',b);          u(unit).borderScore = b;
fprintf(fid2,'%3.4f\t',s);          u(unit).speedScore = s;
fprintf(fid2,'%3.4f\t',tm);         u(unit).thetaMod = tm;
fprintf(fid2,'%3.4f\n',maxF);       u(unit).thetaModFreq = maxF;


%__________________________________________________________________________
%
%                  Functions
%__________________________________________________________________________

function roiName = getRoiName(roi)

switch roi
    case -1; roiName = 'ND';
    case 0; roiName = 'Other';
    case 1; roiName = 'SUB';
    case 2; roiName = 'SUBmol';
    case 3; roiName = 'CA1';
    case 4; roiName = 'DG';
    case 5; roiName = 'CA3';
    otherwise; roiName = 'ND';
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


function r = circ_r(alpha, w, d, dim)
% r = circ_r(alpha, w, d)
%   Computes mean resultant vector length for circular data.
%
%   Input:
%     alpha	sample of angles in radians
%     [w		number of incidences in case of binned angle data]
%     [d    spacing of bin centers for binned data, if supplied 
%           correction factor is used to correct for bias in 
%           estimation of r, in radians (!)]
%     [dim  compute along this dimension, default is 1]
%
%     If dim argument is specified, all other optional arguments can be
%     left empty: circ_r(alpha, [], [], dim)
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

if nargin < 4
  dim = 1;
end

if nargin < 2 || isempty(w) 
  % if no specific weighting has been specified
  % assume no binning has taken place
	w = ones(size(alpha));
else
  if size(w,2) ~= size(alpha,2) || size(w,1) ~= size(alpha,1) 
    error('Input dimensions do not match');
  end 
end

if nargin < 3 || isempty(d)
  % per default do not apply correct for binned data
  d = 0;
end

% compute weighted sum of cos and sin of angles
r = sum(w.*exp(1i*alpha),dim);

% obtain length 
r = abs(r)./sum(w,dim);

% for data with known spacing, apply correction factor to correct for bias
% in the estimation of r (see Zar, p. 601, equ. 26.16)
if d ~= 0
  c = d/2/sin(d/2);
  r = c*r;
end


function [mu ul ll] = circ_mean(alpha, w, dim)
%
% mu = circ_mean(alpha, w)
%   Computes the mean direction for circular data.
%
%   Input:
%     alpha	sample of angles in radians
%     [w		weightings in case of binned angle data]
%     [dim  compute along this dimension, default is 1]
%
%     If dim argument is specified, all other optional arguments can be
%     left empty: circ_mean(alpha, [], dim)
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

if nargin < 3
  dim = 1;
end

if nargin < 2 || isempty(w)
  % if no specific weighting has been specified
  % assume no binning has taken place
	w = ones(size(alpha));
else
  if size(w,2) ~= size(alpha,2) || size(w,1) ~= size(alpha,1) 
    error('Input dimensions do not match');
  end 
end

% compute weighted sum of cos and sin of angles
r = sum(w.*exp(1i*alpha),dim);

% obtain mean by
mu = angle(r);

% confidence limits if desired
if nargout > 1
  t = circ_confmean(alpha,0.05,w,[],dim);
  ul = mu + t;
  ll = mu - t;
end


% circular smoothing (mean window filter, 14 bins on each side)
% Langston et al., Science 2010
function r = hdsmooth(r)

r = r(:);

% number of bins
n = length(r);

% duplicate the data
r2 = [r;r;r];

% window function
w = ones(29,1)/29;

% smoothing
r3 = conv(r2,w,'same');

% take the center of smoothed data
r = r3(n+1:2*n);



% Finds the position timestamp indexes for spike timestamps
function spkInd = getSpkInd(ts,post,cPost)
% Number of spikes
N = length(ts);
spkInd = zeros(N,1);

count = 0;
for ii = 1:N
    tdiff = (post-ts(ii)).^2;
    tdiff2 = (cPost-ts(ii)).^2;
    [m,ind] = min(tdiff);
    [m2,ind2] = min(tdiff2);
    % Check if spike is in legal time sone
    if m == m2
        count = count + 1;
        spkInd(count) = ind(1);
    end
end
spkInd = spkInd(1:count);



% Smooth samples with a mean filter over 15 samples
function [posx,posy] = posMeanFilter(posx,posy)

n = 1;
w = ones(2*n+1,1)/(2*n+1);
posx = conv(posx,w,'valid');
posy = conv(posy,w,'valid');
posx = [ones(n,1)*posx(1); posx; ones(n,1)*posx(end)];
posy = [ones(n,1)*posy(1); posy; ones(n,1)*posy(end)];


% Removes position "jumps", i.e position samples that imply that the rat is
% moving quicker than physical possible.
function [x1,y1,x2,y2,t] = remBadTrack2Diode(x1,y1,x2,y2,t,treshold)

% Indexes to position samples that are to be removed
remInd = [];

diffX = diff(x1);
diffY = diff(y1);
diffR = sqrt(diffX.^2 + diffY.^2);
ind = find(diffR > treshold);

if ind(end) == length(x1)
    offset = 2;
else
    offset = 1;
end

for ii = 1:length(ind)-offset
    if ind(ii+1) == ind(ii)+1
        % A single sample position jump, tracker jumps out one sample and
        % then jumps back to path on the next sample. Remove bad sample.
        remInd = [remInd; ind(ii)+1];
        ii = ii+1;
        continue
    else
        % Not a single jump. 2 possibilities:
        % 1. Tracker jumps out, and stay out at the same place for several
        % samples and then jumps back.
        % 2. Tracker just has a small jump before path continues as normal,
        % unknown reason for this. In latter case the samples are left
        % untouched.
        idx = find(x1(ind(ii)+1:ind(ii+1)+1)==x1(ind(ii)+1));
        if length(idx) == length(x1(ind(ii)+1:ind(ii+1)+1))
            remInd = [remInd; (ind(ii)+1:ind(ii+1)+1)'];
        end
    end
end
% Remove the samples
x1(remInd) = [];
y1(remInd) = [];
x2(remInd) = [];
y2(remInd) = [];
t(remInd) = [];


% calculate the parameter dm for border score (Solstad, Science, 2008)
function dm = dmCalc(map,nBorder,xy)

% map size
[nY,nX] = size(map);

% list all bin coordinates
binsX =[];
binsY =[];
for ii=1:nBorder
    binsX = [binsX; xy(ii).x];
    binsY = [binsY; xy(ii).y];
end
linearIdx = sub2ind(size(map),binsY,binsX);

% firing rate at each bin
rate = map(linearIdx);
% normalize
rate = rate./sum(rate);

% distance to nearest wall [bins]
d = min( [binsY-1 nY-binsY binsX-1 nX-binsX],[],2 );
% weighted average distance
d = sum(d .* rate);

% normalize by the half length of open field
dm = d/(min([nX nY])/2);

% calculate the parameter cm for border score (Solstad, Science, 2008)
function cm = cmCalc(map,nBorder,xy)

% map size
[nY,nX] = size(map);

% coverage of walls (expcted to be a rectangular open field)
coverage =[];
for ii=1:nBorder
    % define bins in wall 1
    Y = 1:nY;
    X = ones(1,nY);
    linearIdx = sub2ind(size(map),Y,X);
    % select the 2nd-inside bin if the outermost bin is nonvisited
    X(isnan(map(linearIdx))) = 2;
    % count the number of visited/occupied bins
    coverage = [coverage; countbins(map,xy(ii),X,Y)];
    
    % define bins in wall 2
    Y = 1:nY;
    X = ones(1,nY)*nX;
    linearIdx = sub2ind(size(map),Y,X);
    % select the 2nd-inside bin if the outermost bin is nonvisited
    X(isnan(map(linearIdx))) = nX-1;
    % count the number of occupied bins
    coverage = [coverage; countbins(map,xy(ii),X,Y)];
    
    % define bins in wall 3
    Y = ones(1,nX);
    X = 1:nX;
    linearIdx = sub2ind(size(map),Y,X);
    % select the 2nd-inside bin if the outermost bin is nonvisited
    Y(isnan(map(linearIdx))) = 2;
    % count the number of visited/occupied bins
    coverage = [coverage; countbins(map,xy(ii),X,Y)];
    
    % define bins in wall 4
    Y = ones(1,nX)*nY;
    X = 1:nX;
    linearIdx = sub2ind(size(map),Y,X);
    % select the 2nd-inside bin if the outermost bin is nonvisited
    Y(isnan(map(linearIdx))) = nY-1;
    % count the number of visited/occupied bins
    coverage = [coverage; countbins(map,xy(ii),X,Y)];
end
% cm (maximum coverage of any single field over)
cm = max(coverage);

% count bins for border calculation
function coverage = countbins(map,xy,X,Y)

% number of bins to be tested
n = length(X);

% count bins
nBins = 0;
nOccupied = 0;
for jj=1:n
    % visited bins

    if ~isnan(map(Y(jj),X(jj))); nBins = nBins + 1; end
    % occupied bins
    if sum(Y(jj)==xy.y & X(jj)==xy.x); nOccupied = nOccupied+1; end
end
% count the number of occupied bins
coverage = nOccupied/nBins;


%  calculate grid score by changing the radius (Langston Science 2010)
function g = gridscore(Rxx,corrAxis,binWidth,inRadius,outRadius)
    
% rotation step [deg]
rotStep = 30;   

% Number of iterations
N = floor(180/rotStep);

% start rotation
g = [];
for r = inRadius:binWidth:outRadius
    
    % rotate in steps of 30 deg
    corrVal = nan(N,1);
    for ii = 1:N
        % rotate Rxx
        RxxR = imrotate(Rxx,rotStep*(ii-1),'bilinear','crop');
        
        % Set the part of the map outside the radius to NaN
        RxxTrim  = adjustMap(Rxx, r+binWidth,r,corrAxis);
        RxxRTrim = adjustMap(RxxR,r+binWidth,r,corrAxis);
        
%         figure(99);
%         subplot(221); imagesc(RxxTrim);  axis image;
%         subplot(222); imagesc(RxxRTrim); axis image;
%         drawnow
        
        % Calculate the crosscorrelation, only the central value.
        Rxy = pointCorr(RxxTrim,RxxRTrim,0,0,size(Rxx,1));
        corrVal(ii) = Rxy;
    end
    % grid score with the current radius
    g = [g; min([corrVal(3),corrVal(5)]) - max([corrVal(2),corrVal(4),corrVal(6)])];
end
% subplot(223);
% plot(g)
% ylim([-2 2])
% xlabel('radius idx')
% ylabel('grid score')
% drawnow

% grid score
g = max(g);


% Determine the radius of center peak of autocorrelation map based on 
% Stensola & Moser Nature 2012.
function centreRadius = centreradius(Rxx,corrAxis,binWidth,boxsize)

% the range of radius to be examined
rRange = binWidth:binWidth:boxsize-10;
n = length(rRange);
meanRxx = nan(n,1);
minRxx  = nan(n,1);

% for each radius
for ii = 1:n
    % Set the part of the map outside the radius to NaN
    RxxTmp  = adjustMap(Rxx, rRange(ii)+binWidth,rRange(ii),corrAxis);
    % 2 methods
    meanRxx(ii) = nanmean(RxxTmp(:));
    minRxx(ii)  = min(RxxTmp(:));
end

% method1: local minimum of average correlation values
[~,idx1] = findpeaks(-meanRxx);
idx1 = min(idx1);
% method2: the first index with correlation <0.2
idx2 = find(minRxx<0.2,1);

% inner radius (radius of central peak of Rxx)
idx = min([idx1 idx2]);
centreRadius = idx * binWidth;


% Rotated the position angles according to the angle tAngle [radians]
function [newX,newY] = rotatePath(x,y,tAngle)

newX = x * cos(tAngle) - y * sin(tAngle); % Tilted x-coord
newY = x * sin(tAngle) + y * cos(tAngle); % Tilted y-coord


% Locates the postion to the 6 closest fields to the centre field, and sets
% the radius of the autocorrelation map equal to the distance from oriogo
% to the bin, of the bins constructing the 6 fields, farthers from the
% centre.
function [radius,cBinsX,cBinsY] = setCorrelationArea(Rxx,binWidth,aTreshold)

% Number of bins in the map
[M,N] = size(Rxx);

% Set the map axis
mapAxis = binWidth * linspace(-((size(Rxx,1)-1)/2),((size(Rxx,1)-1)/2),size(Rxx,1));


% Allocate memory for the array
visited = zeros(M,N);
% Set bins with rate below treshold to visited.
 snitt=nanmean(nanmean(Rxx));
 standardavvik=mean(nanstd(Rxx));
 terskel = snitt + (aTreshold*standardavvik);
% ind = find(Rxx < aTreshold);
ind = find(Rxx < terskel);
sprintf('%s%3.2f%s%3.2f%s%i','Snitt-korr:_',snitt,'   Terskel:_',terskel,'   Antall bins under terskel: ',length(ind));

visited(ind) = 1;
visited(isnan(Rxx)) = 1;

nFields = 0;
% Minimum number of bins in a field
pBins = 30;
disp(strcat('Minimum number of bins in a field, pBin = ',num2str(pBins)));

% Array that will hold the bins for each correlation field
fieldBins = cell(100,2);
fieldPos = [];

% Go as long as there are unvisited parts of the map left
while ~prod(prod(visited))
     % Array that will contain the bin positions to the current placefield
    binsX = [];
    binsY = [];
    % Find the unvisited bins in the map
    [I,J] = find(visited==0);
    % The first unvisited bin are set as the starting point
    legalI = I(1);
    legalJ = J(1);
    % Go as long as there still are bins left in the current placefield
    
    while 1
        % Add current bin to the bin position arrays
        binsX = [binsX; legalI(1)];
        binsY = [binsY; legalJ(1)];
        % Set the current bin to visited
        visited(legalI(1),legalJ(1)) = 1;
        
        % Find which of the 4 neigbour bins that are part of the placefield
        [legalI,legalJ] = getLegals(visited,legalI,legalJ);
        
        % Remove current bin from the array containing bins to be added
        legalI(1) = [];
        legalJ(1) = [];
        % Check if we are finished with this placefield
        if isempty(legalI)
            break;
        end
    end
    if length(binsX)>=pBins % Minimum size of a placefield
        nFields = nFields + 1;
        % Store the bin indexes for the field
        fieldBins(nFields,1) = {binsX};
        fieldBins(nFields,2) = {binsY};
        % Find centre of mass (com)
        comX = 0;
        comY = 0;
        % Total rate
        R = 0;
        for ii = 1:length(binsX)
            R = R + Rxx(binsX(ii),binsY(ii));
            comX = comX + Rxx(binsX(ii),binsY(ii))*mapAxis(binsX(ii));
            comY = comY + Rxx(binsX(ii),binsY(ii))*mapAxis(binsY(ii));
        end
        fieldPos = [fieldPos; [comX/R, comY/R]];
    end
end

% ** Should alway be at least one field (centre field) **
% Locate the field closest to origo.
if isempty(fieldPos)
    disp('No fields was found');
    cBinsX = [];
    cBinsY = [];
    radius = 1000;
    return
end
y = fieldPos(:,1);
x = fieldPos(:,2);
dist = sqrt(x.^2 + y.^2);
[minDist, minInd] = min(dist);
% Position for centre field
cField = [x(minInd), y(minInd)];
cBinsX = fieldBins{minInd,1};
cBinsY = fieldBins{minInd,2};

if length(dist) > 1
    dist = zeros(nFields,1);
    for ii  = 1:nFields
        if ii ~= minInd
            % Distance between this field and the centre field
            dist(ii) = sqrt((x(ii)-cField(1))^2 + (y(ii)-cField(2))^2);
        else
            % Set the distance to the centre field itself to infinity
            dist(ii) = inf;
        end
    end

    % Find the 6 closest fields to the centre field if enough fields
    fieldDist = zeros(6,2);
    if length(dist)>6
        % Keep the 6 closest fields
        for ii=1:6
            [m,ind] = min(dist);
            fieldDist(ii,1) = m;
            fieldDist(ii,2) = ind;
            dist(ind) = inf;
        end
    else
        % Keep all the fields
        for ii=1:length(dist)
            [m,ind] = min(dist);
            if m~=inf
                fieldDist(ii,1) = m;
                fieldDist(ii,2) = ind;
                dist(ind) = inf;
            end
        end
    end
nFields;
    % Find what bin of the bins that construct the 6 (or fewer) fields are
    % fartherst off the centre field
    if length(dist) < 6
        N = length(dist);
    else
        N = 6;
    end
    sDist = 0;
    sInd = 0;

    for ii = 1:N
        d = fieldDist(ii,1);
        if d>sDist && ~isnan(d)
            sDist = d;
            sInd = ii;
        end
    end

    if sDist ~= 0
        % Set the radius
        r = sqrt(mapAxis(fieldBins{fieldDist(sInd,2),1}).^2 + mapAxis(fieldBins{fieldDist(sInd,2),2}).^2);
        radius = max(r);
    else
        disp('Only one field was found');
        radius = 1000;
    end
elseif length(dist) == 1
    disp('Only one field was found');
    radius = 1000;
else
    disp('No fields was found, treshold must be adjusted')
    cBinsX = [];
    cBinsY = [];
    radius = 1000;
end
if N < 6
    disp('Less than 6 peaks was found and no radius could be calculated');
    radius = 1000;
end


function [legalI,legalJ] = getLegals(visited,legalI,legalJ)
% Current bin
cI = legalI(1);
cJ = legalJ(1);
% Neigbour bins
leftI   = cI-1;
rightI  = cI+1;
upI     = cI;
downI   = cI;
leftJ   = cJ;
rightJ  = cJ;
upJ     = cJ-1;
downJ   = cJ+1;

% Check left
if leftI >= 1 % Inside map
    if visited(leftI,leftJ)==0 % Unvisited bin
        if ~(~isempty(find(legalI==leftI, 1)) && ~isempty(find(legalJ==leftJ, 1))) % Not part of the array yet
            % Left bin is part of placefield and must be added
            legalI = [legalI;leftI];
            legalJ = [legalJ;leftJ];
        end
    end
end
% Check rigth
if rightI <= size(visited,2) % Inside map
    if visited(rightI,rightJ)==0 % Unvisited bin
        if ~(~isempty(find(legalI==rightI, 1)) && ~isempty(find(legalJ==rightJ, 1))) % Not part of the array yet
            % Right bin is part of placefield and must be added
            legalI = [legalI;rightI];
            legalJ = [legalJ;rightJ];
        end
    end
end
% Check up
if upJ >= 1 % Inside map
    if visited(upI,upJ)==0 % Unvisited bin
        if ~(~isempty(find(legalI==upI, 1)) && ~isempty(find(legalJ==upJ, 1))) % Not part of the array yet
            % Up bin is part of placefield and must be added
            legalI = [legalI;upI];
            legalJ = [legalJ;upJ];
        end
    end
end
% Check down
if downJ <= size(visited,1) % Inside map
    if visited(downI,downJ)==0 % Unvisited bin
        if ~(~isempty(find(legalI==downI, 1)) && ~isempty(find(legalJ==downJ, 1))) % Not part of the array yet
            % Right bin is part of placefield and must be added
            legalI = [legalI;downI];
            legalJ = [legalJ;downJ];
        end
    end
end


% Finds the radius of the centre field
function centreRadius = removeCentreField(cBinsX,cBinsY,mapAxis)

% Set the radius
if ~isempty(cBinsX)
    r = sqrt(mapAxis(cBinsX).^2 + mapAxis(cBinsY).^2);
    centreRadius = max(r);
else
    centreRadius = -1;
end


% Sets the bins of the map outside the radius to NaN
function Rxx = adjustMap(Rxx,radius,centreRadius,mapAxis)

% Size, in number of bins, of the map
[N,M] = size(Rxx);
centreX = (N+1)/2;
centreY = (M+1)/2;

binWidth = mapAxis(2) - mapAxis(1);

% Calculate the distance from origo for each bin in the map
oDist = zeros(N,M);
for ii = 1:N
    for jj = 1:M
        oDist(ii,jj) = sqrt(mapAxis(ii)^2 + mapAxis(jj)^2);
    end
end
Rxx(oDist>radius) = NaN;
Rxx(oDist<=centreRadius) = NaN;


% Calculates the correlation between map1 and map 2. By using the same map
% for map1 and map2, the autocorrelation is calculated for that map. It is
% assumed that the size of the 2 maps are equal.
% Author: Raymond Skjerpeng
function Rxy = correlation(map1,map2)

map1(isnan(map2)) = NaN;
map2(isnan(map1)) = NaN;

bins = size(map1,1);
% N = bins + round(0.8*bins);
N = bins + round(1*bins);
if ~mod(N,2)
    N = N - 1;
end
% Centre bin
cb = (N+1)/2;
Rxy = zeros(N);
for ii = 1:N
    rowOff = ii-cb;
    for jj = 1:N
        colOff = jj-cb;
        Rxy(ii,jj) = pointCorr(map1,map2,rowOff,colOff,bins);
    end
end

% Calculates the correlation for a point in the autocorrelogram. It is
% using the Pearsons correlation method (Sargolini Science 2006).
function Rxy = pointCorr(map1,map2,rowOff,colOff,N)

% Number of rows in the correlation for this lag
numRows = N - abs(rowOff);
% Number of columns in the correlation for this lag
numCol = N - abs(colOff);

% Set the start and the stop indexes for the maps
if rowOff > 0
    rSt1 = 1+abs(rowOff)-1;
    rSt2 = 0;
else
    rSt1 = 0;
    rSt2 = abs(rowOff);
end
if colOff > 0
    cSt1 = abs(colOff);
    cSt2 = 0;
else
    cSt1 = 0;
    cSt2 = abs(colOff);
end

sumXY = 0;
sumX = 0;
sumY = 0;
sumX2 = 0;
sumY2 = 0;
NB = 0;
for ii = 1:numRows
    for jj = 1:numCol
        if ~isnan(map1(rSt1+ii,cSt1+jj)) && ~isnan(map2(rSt2+ii,cSt2+jj))
            NB = NB + 1;
            sumX = sumX + map1(rSt1+ii,cSt1+jj);
            sumY = sumY + map2(rSt2+ii,cSt2+jj);
            sumXY = sumXY + map1(rSt1+ii,cSt1+jj) * map2(rSt2+ii,cSt2+jj);
            sumX2 = sumX2 + map1(rSt1+ii,cSt1+jj)^2;
            sumY2 = sumY2 + map2(rSt2+ii,cSt2+jj)^2;
        end
    end
end

if NB >= 20
    sumx2 = sumX2 - sumX^2/NB;
    sumy2 = sumY2 - sumY^2/NB;
    sumxy = sumXY - sumX*sumY/NB;
    if (sumx2<=0 && sumy2>=0) || (sumx2>=0 && sumy2<=0)
        Rxy = NaN;
    else
        Rxy = sumxy/sqrt(sumx2*sumy2);
    end
else
    Rxy = NaN;
end


%__________________________________________________________________________
%
%                       Statistic function
%__________________________________________________________________________

% Shannon information, sparseness, and selectivity  
function [information,sparsity,selectivity] = mapstat(map,posPDF)

% Sparseness
n = size(map,1);
meanrate = nansum(nansum( map .* posPDF ));
meansquarerate = nansum(nansum( (map.^2) .* posPDF ));
if meansquarerate == 0
    sparsity = NaN;
else
    sparsity = meanrate^2 / meansquarerate;
end

% Selectivity
maxrate = max(max(map));
if meanrate == 0
   selectivity = NaN;
else
   selectivity = maxrate/meanrate;
end

% Shannon information
% Note: current formula is "information density (bits/spike)".
% Note: The following formula is "information rate (bits/sec)":
%   akksum = akksum + posPDF(ii1,ii2) * map(ii1,ii2) * log2( map(ii1,ii2) / meanrate );,
%
[i1, i2] = find( (map>0) & (posPDF>0) );  % the limit of x*log(x) as x->0 is 0 
if ~isempty(i1)
    akksum = 0;
    for i = 1:length(i1);
        ii1 = i1(i);
        ii2 = i2(i);
        % bits/spike
        akksum = akksum + posPDF(ii1,ii2) * (map(ii1,ii2)/meanrate) * log2( map(ii1,ii2) / meanrate ); 
    end
    information = akksum;
else
    information = NaN;
end


function [bursts,singlespikes] = burstfinder(ts,maxisi)
bursts = [];
singlespikes = [];
isi = diff(ts);
n = length(ts);
if numel(isi)>0  % <--modified by Takuma 100423
    if isi(1) <= maxisi
       bursts = 1;
    else
       singlespikes = 1;
    end
    for t = 2:n-1;
       if (isi(t-1)>maxisi) && (isi(t)<=maxisi)
          bursts = [bursts; t];
       elseif (isi(t-1)>maxisi) && (isi(t)>maxisi)
          singlespikes = [singlespikes; t];      
       end
    end   
    if (isi(n-1)>maxisi) 
        singlespikes = [singlespikes; n];      
    end
else
    bursts=0;
    singlespikes=0;
end


function z = fieldcohere(map)
[n,m] = size(map);
tmp = zeros(n*m,2);
k=0;
for y = 1:n
    for x = 1:m
        k = k + 1;
        xstart = max([1,x-1]);
        ystart = max([1,y-1]);
        xend = min([m x+1]);
        yend = min([n y+1]);
        nn = sum(sum(isfinite(map(ystart:yend,xstart:xend)))) - isfinite(map(y,x));
        if (nn > 0)
            tmp(k,1) = map(y,x);
            tmp(k,2) = nansum([ nansum(nansum(map(ystart:yend,xstart:xend))) , -map(y,x) ]) / nn;
        else
            tmp(k,:) = [NaN,NaN];    
        end
    end
end
index = find( isfinite(tmp(:,1)) & isfinite(tmp(:,2)) );
if length(index) > 3
    cc = corrcoef(tmp(index,:));
    z = atanh(cc(2,1));
else
    z = NaN;
end

%__________________________________________________________________________
%
%           Function for setting the axis for drawing
%__________________________________________________________________________

function mapAxis = setMapAxis(posx,posy,mapAxis,binWidth)

% Check for asymmetri in the path. If so correct acount for it in
% mapAxis
minX = min(posx);
maxX = max(posx);
minY = min(posy);
maxY = max(posy);
if minX < mapAxis(1)
    nXtra = ceil(abs(minX-mapAxis(1))/binWidth);
else
    nXtra = 0;
end
if maxX > mapAxis(end)
    pXtra = ceil(abs(maxX-mapAxis(end))/binWidth);
else
    pXtra = 0;
end
if nXtra
    for nn =1:nXtra
        tmp = mapAxis(1) - binWidth;
        mapAxis = [tmp; mapAxis'];
        mapAxis = mapAxis';
    end
end
if pXtra
    tmp = mapAxis(end) + binWidth;
    mapAxis = [mapAxis'; tmp];
    mapAxis = mapAxis';
end

if minY < mapAxis(1)
    nXtra = ceil(abs(minX-mapAxis(1))/binWidth);
else
    nXtra = 0;
end
if maxY > mapAxis(end)
    pXtra = ceil(abs(maxX-mapAxis(end))/binWidth);
else
    pXtra = 0;
end
if nXtra
    for nn =1:nXtra
        tmp = mapAxis(1) - binWidth;
        mapAxis = [tmp; mapAxis'];
        mapAxis = mapAxis';
    end
end
if pXtra
    tmp = mapAxis(end) + binWidth;
    mapAxis = [mapAxis'; tmp];
    mapAxis = mapAxis';
end

% Put on 3 extra cm on each side.
mapAxis = [mapAxis(1)-1.5;mapAxis'];
mapAxis = [mapAxis; mapAxis(end)+1.5];
mapAxis = mapAxis';
mapAxis = [mapAxis(1)-1.5;mapAxis'];
mapAxis = [mapAxis; mapAxis(end)+1.5];
mapAxis = mapAxis';


%__________________________________________________________________________
%
%                           Field functions
%__________________________________________________________________________


% Calculates the rate map and the position probability density function (pospdf).
function [map,pospdf] = ratemap(spkx,spky,posx,posy,post,h,yAxis,xAxis)
invh = 1/h;
map = zeros(length(xAxis),length(yAxis));
pospdf = zeros(length(xAxis),length(yAxis));

yy = 0;
for y = yAxis
    yy = yy + 1;
    xx = 0;
    for x = xAxis
        xx = xx + 1;
        [map(yy,xx),pospdf(yy,xx)] = rate_estimator(spkx,spky,x,y,invh,posx,posy,post);
    end
end

% Normalize the pospdf (Should actually normalize the integral of the
% pospdf by dividing by total area too, but the functions making use of the 
% pospdf assume that this is not done (see mapstat()).
pospdf = pospdf ./ sum(sum(pospdf)); % Position Probability Density Function


% Calculate the rate for one position value
function [r,edge_corrector] = rate_estimator(spkx,spky,x,y,invh,posx,posy,post)
% edge-corrected kernel density estimator

conv_sum = sum(gaussian_kernel(((spkx-x)*invh),((spky-y)*invh)));
edge_corrector =  trapz(post,gaussian_kernel(((posx-x)*invh),((posy-y)*invh)));
%edge_corrector(edge_corrector<0.15) = NaN;
r = (conv_sum / (edge_corrector + 0.0001)) + 0.0001; % regularised firing rate for "wellbehavedness"
                                                       % i.e. no division by zero or log of zero
% Gaussian kernel for the rate calculation
function r = gaussian_kernel(x,y)
% k(u) = ((2*pi)^(-length(u)/2)) * exp(u'*u)
r = 0.15915494309190 * exp(-0.5*(x.*x + y.*y));


% Identifies border fields with >0.3*peakrate, >200cm2 (Solstad, Science,
% 2008). Modified from furnction placefield.
%
% map           Rate map
% mapAxis       The map axis
function [nBorder,xy]= borderfield(map,mapAxis)

binWidth = mapAxis(2) - mapAxis(1);

% Counter for the number of fields
nBorder = 0;
% Field properties will be stored in this struct array
xy = [];

% Allocate memory to the arrays
[M,N] = size(map);
% Array that contain the bins of the map this algorithm has visited
visited = zeros(M,N);
nanInd = isnan(map);
visited(nanInd) = 1;
visited2 = visited;

% maximum rate of the entire map
peak = max(map(:));

% Go as long as there are unvisited parts of the map left
while ~prod(prod(visited))
    % Array that will contain the bin positions to the current placefield
    binsX = [];
    binsY = [];
    
    % Find the current maximum
    [peak2,r] = max(map);
    [peak2,pCol] = max(peak2);
    pCol = pCol(1);
    pRow = r(pCol);
    
    % Check if peak rate is high enough
    if peak2 <= 0.3 * peak
        break;
    end
    
    % Find firing field with peak rate larger than 0.3*peak
    visited2( map <= 0.3*peak ) = 1;
    % Find the bins that construct the peak field
    [binsX,binsY,visited2] = recursiveBins(map,visited2,[],[],pRow,pCol,N,M);

    if length(binsX) >= 200/binWidth^2  % Min size of border field (200 cm2)
        nBorder = nBorder + 1;
        % Put the field properties in the struct array
        xy = [xy; struct('x',binsY,'y',binsX)];
    end
    visited(binsX,binsY) = 1;
    map(visited2==1) = 0;
end


% placefield identifies the placefields in the firing map. It returns the
% number of placefields and the location of the peak within each
% placefield.
%
% map           Rate map
% pTreshold     Field treshold
% pBins         Minimum number of bins in a field
% mapAxis       The map axis
function [nFields,fieldProp] = placefield(map,p,mapAxis)

binWidth = mapAxis(2) - mapAxis(1);


% Counter for the number of fields
nFields = 0;
% Field properties will be stored in this struct array
fieldProp = [];

% Allocate memory to the arrays
[M,N] = size(map);
% Array that contain the bins of the map this algorithm has visited
visited = zeros(M,N);
nanInd = isnan(map);
visited(nanInd) = 1;
visited2 = visited;

% Go as long as there are unvisited parts of the map left
while ~prod(prod(visited))
    % Array that will contain the bin positions to the current placefield
    binsX = [];
    binsY = [];
    
    % Find the current maximum
    [peak,r] = max(map);
    [peak,pCol] = max(peak);
    pCol = pCol(1);
    pRow = r(pCol);
    
    % Check if peak rate is high enough
    if peak < p.lowestFieldRate
        break;
    end
    
    visited2(map<p.fieldTreshold*peak) = 1;
    % Find the bins that construct the peak field
    [binsX,binsY,visited2] = recursiveBins(map,visited2,[],[],pRow,pCol,N,M);

    if length(binsX)>=p.minNumBins % Minimum size of a placefield
        nFields = nFields + 1;
        % Find centre of mass (com)
        comX = 0;
        comY = 0;
        % Total rate
        R = 0;
        for ii = 1:length(binsX)
            R = R + map(binsX(ii),binsY(ii));
            comX = comX + map(binsX(ii),binsY(ii))*mapAxis(binsX(ii));
            comY = comY + map(binsX(ii),binsY(ii))*mapAxis(binsY(ii));
        end
        % Average rate in field
        avgRate = nanmean(nanmean(map(binsX,binsY)));
        % Peak rate in field
        peakRate = nanmax(nanmax(map(binsX,binsY)));
        % Size of field
        fieldSize = length(binsX) * binWidth^2;
        % Put the field properties in the struct array
        fieldProp = [fieldProp; struct('x',comY/R,'y',comX/R,'avgRate',avgRate,'peakRate',peakRate,'size',fieldSize)];
    end
    visited(binsX,binsY) = 1;
    map(visited2==1) = 0;
end


function [binsX,binsY,visited] = recursiveBins(map,visited,binsX,binsY,ii,jj,N,M)
% If outside boundaries of map -> return.
if ii<1 || ii>N || jj<1 || jj>M
    return;
end
% If all bins are visited -> return.
if prod(prod(visited))
    return;
end
if visited(ii,jj) % This bin has been visited before
    return;
else
    binsX = [binsX;ii];
    binsY = [binsY;jj];
    visited(ii,jj) = 1;
    % Call this function again in each of the 4 neighbour bins
    [binsX,binsY,visited] = recursiveBins(map,visited,binsX,binsY,ii,jj-1,N,M);
    [binsX,binsY,visited] = recursiveBins(map,visited,binsX,binsY,ii-1,jj,N,M);
    [binsX,binsY,visited] = recursiveBins(map,visited,binsX,binsY,ii,jj+1,N,M);
    [binsX,binsY,visited] = recursiveBins(map,visited,binsX,binsY,ii+1,jj,N,M);
end


% function [legalI,legalJ] = getLegals(visited,legalI,legalJ)
% % Current bin
% cI = legalI(1);
% cJ = legalJ(1);
% % Neigbour bins
% leftI   = cI-1;
% rightI  = cI+1;
% upI     = cI;
% downI   = cI;
% leftJ   = cJ;
% rightJ  = cJ;
% upJ     = cJ-1;
% downJ   = cJ+1;
% 
% % Check left
% if leftI >= 1 % Inside map
%     if visited(leftI,leftJ)==0 % Unvisited bin
%         if ~(length(find(legalI==leftI)) & length(find(legalJ==leftJ))) % Not part of the array yet
%             % Left bin is part of placefield and must be added
%             legalI = [legalI;leftI];
%             legalJ = [legalJ;leftJ];
%         end
%     end
% end
% % Check rigth
% if rightI <= size(visited,2) % Inside map
%     if visited(rightI,rightJ)==0 % Unvisited bin
%         if ~(length(find(legalI==rightI)) & length(find(legalJ==rightJ))) % Not part of the array yet
%             % Right bin is part of placefield and must be added
%             legalI = [legalI;rightI];
%             legalJ = [legalJ;rightJ];
%         end
%     end
% end
% % Check up
% if upJ >= 1 % Inside map
%     if visited(upI,upJ)==0 % Unvisited bin
%         if ~(length(find(legalI==upI)) & length(find(legalJ==upJ))) % Not part of the array yet
%             % Up bin is part of placefield and must be added
%             legalI = [legalI;upI];
%             legalJ = [legalJ;upJ];
%         end
%     end
% end
% % Check down
% if downJ <= size(visited,1) % Inside map
%     if visited(downI,downJ)==0 % Unvisited bin
%         if ~(length(find(legalI==downI)) & length(find(legalJ==downJ))) % Not part of the array yet
%             % Right bin is part of placefield and must be added
%             legalI = [legalI;downI];
%             legalJ = [legalJ;downJ];
%         end
%     end
% end


% Finds the position to the spikes
function [spkx,spky,newTs,spkInd] = spikePos(ts,posx,posy,post,cPost)
N = length(ts);
spkx = zeros(N,1);
spky = zeros(N,1);
newTs = zeros(N,1);
spkInd = zeros(N,1);
count = 0;
for ii = 1:N
    tdiff = (post-ts(ii)).^2;
    tdiff2 = (cPost-ts(ii)).^2;
    [m,ind] = min(tdiff);
    [m2,ind2] = min(tdiff2);
    % Check if spike is in legal time sone
    if m == m2
        count = count + 1;
        spkx(count) = posx(ind(1));
        spky(count) = posy(ind(1));
        newTs(count) = ts(ii);
        spkInd(count) = ind(1);
    end
end
spkx = spkx(1:count);
spky = spky(1:count);
newTs = newTs(1:count);
spkInd = spkInd(1:count);

% Calculates what area of the map that has been visited by the rat
function visited = visitedBins(posx,posy,mapAxis)

% binWidth = mapAxis(2)-mapAxis(1);

% Takuma modification, 3 cm by (Hippocampus,2008,18:1230)
binWidth = (mapAxis(2)-mapAxis(1)); 

% Number of bins in each direction of the map
N = length(mapAxis);
visited = zeros(N);

for ii = 1:N
    for jj = 1:N
        px = mapAxis(ii);
        py = mapAxis(jj);
        distance = sqrt( (px-posx).^2 + (py-posy).^2 );
        
        if min(distance) <= binWidth
            visited(jj,ii) = 1;
        end
    end
end

%__________________________________________________________________________
%
%                   Function for modifying the path
%__________________________________________________________________________

% Removes position "jumps", i.e position samples that imply that the rat is
% moving quicker than physical possible.
function [x,y,t] = remBadTrack(x,y,t,treshold)

% Indexes to position samples that are to be removed
remInd = [];

diffX = diff(x);
diffY = diff(y);
diffR = sqrt(diffX.^2 + diffY.^2);
ind = find(diffR > treshold);

if ind(end) == length(x)
    offset = 2;
else
    offset = 1;
end

for ii = 1:length(ind)-offset
    if ind(ii+1) == ind(ii)+1
        % A single sample position jump, tracker jumps out one sample and
        % then jumps back to path on the next sample. Remove bad sample.
        remInd = [remInd; ind(ii)+1];
        ii = ii+1;
        continue
    else
        % Not a single jump. 2 possibilities:
        % 1. Tracker jumps out, and stay out at the same place for several
        % samples and then jumps back.
        % 2. Tracker just has a small jump before path continues as normal,
        % unknown reason for this. In latter case the samples are left
        % untouched.
        idx = find(x(ind(ii)+1:ind(ii+1)+1)==x(ind(ii)+1));
        if length(idx) == length(x(ind(ii)+1:ind(ii+1)+1));
            remInd = [remInd; (ind(ii)+1:ind(ii+1)+1)'];
        end
    end
end
% Remove the samples
x(remInd) = [];
y(remInd) = [];
t(remInd) = [];


% If 1-5 positon samples are missing, this function will interpolate
% the missing samples. If more than 5 samples are missing in a row they are
% left as NaN.
function [x1,y1] = interporPos(x1,y1)

% Find the indexes to the missing samples, for the red tracking diode
ind = isnan(x1);
ind = find(ind==1);
N = length(ind);

if N == 0
    % No samples missing and we return
    return
end

% Set start and stop points for the loop
if ind(N) >= length(x1)-5
    endLoop = N-5;
else
    endLoop = N;
end
if ind(1) == 1
    startLoop = 2;
else
    startLoop = 1;
end

for ii = startLoop:endLoop
    if isempty(find(ind==ind(ii)+1, 1)) && isempty(find(ind==ind(ii)-1, 1))
        % Only one missing sample in a row
        x1(ind(ii)) = (x1(ind(ii)-1)+x1(ind(ii)+1))/2;
        y1(ind(ii)) = (y1(ind(ii)-1)+y1(ind(ii)+1))/2;
    else
        if isempty(find(ind==ind(ii)+2, 1)) && isempty(find(ind==ind(ii)-1, 1))
            % 2 missing samples in a row
            xDist = abs(x1(ind(ii)-1)-x1(ind(ii)+2));
            yDist = abs(y1(ind(ii)-1)-y1(ind(ii)+2));
            x1(ind(ii)) = x1(ind(ii)-1) + 1/3*xDist;
            y1(ind(ii)) = y1(ind(ii)-1) + 1/3*yDist;
            x1(ind(ii)+1) = x1(ind(ii)-1) + 2/3*xDist;
            y1(ind(ii)+1) = y1(ind(ii)-1) + 2/3*yDist;
            ii = ii+1;
        else
            if isempty(find(ind==ind(ii)+3, 1)) && isempty(find(ind==ind(ii)-1, 1))
                % 3 missing samples in a row
                xDist = abs(x1(ind(ii)-1)-x1(ind(ii)+3));
                yDist = abs(y1(ind(ii)-1)-y1(ind(ii)+3));
                x1(ind(ii)) = x1(ind(ii)-1) + 1/4*xDist;
                y1(ind(ii)) = y1(ind(ii)-1) + 1/4*yDist;
                x1(ind(ii)+1) = x1(ind(ii)-1) + 1/2*xDist;
                y1(ind(ii)+1) = y1(ind(ii)-1) + 1/2*yDist;
                x1(ind(ii)+2) = x1(ind(ii)-1) + 3/4*xDist;
                y1(ind(ii)+2) = y1(ind(ii)-1) + 3/4*yDist;
                ii = ii+2;
            else
                if isempty(find(ind==ind(ii)+4, 1)) && isempty(find(ind==ind(ii)-1, 1))
                    % 4 missing samples in a row
                    xDist = abs(x1(ind(ii)-1)-x1(ind(ii)+4));
                    yDist = abs(y1(ind(ii)-1)-y1(ind(ii)+4));
                    x1(ind(ii)) = x1(ind(ii)-1) + 1/5*xDist;
                    y1(ind(ii)) = y1(ind(ii)-1) + 1/5*yDist;
                    x1(ind(ii)+1) = x1(ind(ii)-1) + 2/5*xDist;
                    y1(ind(ii)+1) = y1(ind(ii)-1) + 2/5*yDist;
                    x1(ind(ii)+2) = x1(ind(ii)-1) + 3/5*xDist;
                    y1(ind(ii)+2) = y1(ind(ii)-1) + 3/5*yDist;
                    x1(ind(ii)+3) = x1(ind(ii)-1) + 4/5*xDist;
                    y1(ind(ii)+3) = y1(ind(ii)-1) + 4/5*yDist;
                    ii = ii+3;
                else
                    if isempty(find(ind==ind(ii)+5, 1)) && isempty(find(ind==ind(ii)-1, 1))
                        % 5 missing samples in a row
                        xDist = abs(x1(ind(ii)-1)-x1(ind(ii)+5));
                        yDist = abs(y1(ind(ii)-1)-y1(ind(ii)+5));
                        x1(ind(ii)) = x1(ind(ii)-1) + 1/6*xDist;
                        y1(ind(ii)) = y1(ind(ii)-1) + 1/6*yDist;
                        x1(ind(ii)+1) = x1(ind(ii)-1) + 2/6*xDist;
                        y1(ind(ii)+1) = y1(ind(ii)-1) + 2/6*yDist;
                        x1(ind(ii)+2) = x1(ind(ii)-1) + 3/6*xDist;
                        y1(ind(ii)+2) = y1(ind(ii)-1) + 3/6*yDist;
                        x1(ind(ii)+3) = x1(ind(ii)-1) + 4/6*xDist;
                        y1(ind(ii)+3) = y1(ind(ii)-1) + 4/6*yDist;
                        x1(ind(ii)+4) = x1(ind(ii)-1) + 5/6*xDist;
                        y1(ind(ii)+4) = y1(ind(ii)-1) + 5/6*yDist;
                        ii = ii+4;
                    end
                end
            end
        end
    end
end


% Find the centre of the box
function centre = centreBox(posx,posy)
% Find border values for path and box
maxX = max(posx);
minX = min(posx);
maxY = max(posy);
minY = min(posy);

% Set the corners of the reference box
NE = [maxX, maxY];
NW = [minX, maxY];
SW = [minX, minY];
SE = [maxX, minY];

% Get the centre coordinates of the box
centre = findCentre(NE,NW,SW,SE);

% Calculates the centre of the box from the corner coordinates
function centre = findCentre(NE,NW,SW,SE)

% The centre will be at the point of interception by the corner diagonals
a = (NE(2)-SW(2))/(NE(1)-SW(1)); % Slope for the NE-SW diagonal
b = (SE(2)-NW(2))/(SE(1)-NW(1)); % Slope for the SE-NW diagonal
c = SW(2);
d = NW(2);
x = (d-c+a*SW(1)-b*NW(1))/(a-b); % X-coord of centre
y = a*(x-SW(1))+c; % Y-coord of centre
centre = [x,y];



%__________________________________________________________________________
%
%           Function for fixing the position timestamps
%__________________________________________________________________________

function [didFix,fixedPost] = fixTimestamps(post)

% First time stamp in file
first = post(1);
% Number of timestamps
N = length(post);
uniqePost = unique(post);

if length(uniqePost)~=N
    disp('Position timestamps are corrected');
    didFix = 1;
    numZeros = 0;
    % Find the number of zeros at the end of the file
    while 1
        if post(end-numZeros)==0
            numZeros = numZeros + 1;
        else
            break;
        end
    end
    
    last = first + (N-1-numZeros) *0.02;
    fixedPost = first:0.02:last;
    fixedPost = fixedPost';
else
    didFix = 0;
    fixedPost = [];
end


%__________________________________________________________________________
%
%                           Import rutines
%__________________________________________________________________________

% Read 2 Diode Axona position data. This is the case when two diodes have
% been use for recording, but must not be confused with the option in the
% Axona recording system to do a two spot recording, because this is
% something else.
function [x1,y1,x2,y2,t] = get2DiodePos(posfile,colour1,colour2,arena)

% Call Sturla Moldens import rutine to get the position data
[tracker,trackerparam] = importvideotracker(posfile);
  
t = zeros(trackerparam.num_pos_samples,1);


% Find the size of the coordinate array. It should be size 4, but can be
% size 2 if the 2 spot recording option has been chosen in the Axona system
% while recording.
N = size(tracker(1).xcoord,2);
if N == 2 % 2 spot chosen
    temp = zeros(trackerparam.num_pos_samples,4);
else 
    temp = zeros(trackerparam.num_pos_samples,8);
end
for ii = 1:trackerparam.num_pos_samples
    t(ii) = tracker(ii).timestamp;
    temp(ii,:) = [tracker(ii).xcoord tracker(ii).ycoord];
end

[didFix,fixedPost] = fixTimestamps(t);
if didFix
    t = fixedPost;
    disp('Continue to read data');
end

if N == 4
    switch colour1 % Red, Green, Blue, Black on White
        case {'red LED'}
            x1 = temp(:,1) + trackerparam.window_min_x;
            y1 = temp(:,5) + trackerparam.window_min_y;
        case {'green LED'}
            x1 = temp(:,2) + trackerparam.window_min_x;
            y1 = temp(:,6) + trackerparam.window_min_y;
        case {'blue LED'}
            x1 = temp(:,3) + trackerparam.window_min_x;
            y1 = temp(:,7) + trackerparam.window_min_y;
        case {'black on white'}
            x1 = temp(:,4) + trackerparam.window_min_x;
            y1 = temp(:,8) + trackerparam.window_min_y;
        otherwise
            error(sprintf('%s%s','unknown colour',colour1));
    end
    switch colour2
        case {'red LED'}
            x2 = temp(:,1) + trackerparam.window_min_x;
            y2 = temp(:,5) + trackerparam.window_min_y;
        case {'green LED'}
            x2 = temp(:,2) + trackerparam.window_min_x;
            y2 = temp(:,6) + trackerparam.window_min_y;
        case {'blue LED'}
            x2 = temp(:,3) + trackerparam.window_min_x;
            y2 = temp(:,7) + trackerparam.window_min_y;
        case {'black on white'}
            x2 = temp(:,4) + trackerparam.window_min_x;
            y2 = temp(:,8) + trackerparam.window_min_y;
        otherwise
            error(sprintf('%s%s','unknown colour',colour2));
    end
end
if N == 2
    x1 = temp(:,1) + trackerparam.window_min_x;
    y1 = temp(:,3) + trackerparam.window_min_y;
    x2 = temp(:,2) + trackerparam.window_min_x;
    y2 = temp(:,4) + trackerparam.window_min_y;
end

numPos = length(x1);
numPost = length(t);
if numPos ~= numPost
    x1 = x1(1:numPost);
    y1 = y1(1:numPost);
    x2 = x2(1:numPost);
    y2 = y2(1:numPost);
end

% index = find ( (posx==0) & (posy==511) ); % this is internal to CBM's equipment.
% posx(index) = NaN;                        % 
% posy(index) = NaN;                        % 


if (nargin > 3)
    % Adjust the coordinates according to the room in use
    [x1, y1] = arena_config(x1,y1,arena);
    [x2, y2] = arena_config(x2,y2,arena);
end
% Force the time stamps to start on zero
t = t - t(1);

% Remove NaN in position if these appear at the end of the file
[x1,y1,x2,y2,t] = removeNaN(x1,y1,x2,y2,t);

function [posx, posy] = arena_config(posx,posy,arena)
switch arena
    case 'Kyoto'            % Kyoto Univ room H101
        centre = [364, 264];
        conversion = 656;   % pixels/m
    case 'room0'
        centre = [433, 256];
        conversion = 387;
    case 'room1'
        centre = [421, 274];
        conversion = 403;
    case 'room2'
        centre = [404, 288];
        conversion = 400;
    case 'room3'
        centre = [356, 289]; % 2. Nov 2004
        conversion = 395;
    case 'room4'
        centre = [418, 186]; % 2. Nov 2004
        conversion = 313;
    case 'room5'
        centre = [406, 268]; % 2. Nov 2004
        conversion = 317;   % 52 cm abow floor, Color camera
    case 'room5bw'
        centre = [454,192]; % 2. Nov 2004
        conversion = 314;   % 52 cm abow floor, Black/White camera
    case 'room5floor'
        centre = [422,266]; % 2. Nov 2004
        conversion = 249;   % On floor, Color camera
    case 'room5bwfloor'
        centre = [471,185]; % 2. Nov 2004
        conversion = 249;   % On floor, Black/White camera
    case 'room6'
        centre = [402,292]; % 1. Nov 2004
        conversion = 398;
    case 'room7'
        centre = [360,230]; % 11. Oct 2004
        conversion = 351;
    case 'room8'
        centre = [390,259]; % Color camera
        conversion = 475;   % January 2006, 52 cm abow floor, 
    case 'room8bw'
        centre = [383,273]; % 11. Oct 2004
        conversion = 391;   % 66 cm abow floor, Black/White camera
    case 'room9'
        centre = [393,290];
        conversion = 246;
    case 'room10'
        centre = [385,331];
        conversion = 247;
    case 'room11'   % New Lab (Jamie and Ros)
        centre = [350, 325]; % 26. August 2008
        conversion = 420;    
    case 'BoxA'
        centre = [230, 245];
        conversion = 350;
    case 'BoxB' % Masato
        centre = [230, 245];
        conversion = 350;
    case 'Linear'
        centre = [339, 120];
        conversion = 295;
    otherwise
        disp('Uknown room info.')
        return
end
posx = 100 * (posx - centre(1))/conversion;
posy = 100 * (centre(2) - posy)/conversion;


% Removes NaN from the position samples if these occure at the end of the
% position file
function [x1,y1,x2,y2,t] = removeNaN(x1,y1,x2,y2,t)

while 1
    if isnan(x1(end)) | isnan(y1(end)) | isnan(x1(end)) | isnan(y2(end))
        x1(end) = [];
        y1(end) = [];
        x2(end) = [];
        y2(end) = [];
        t(end) = [];
    else
        break;
    end
end

              
%__________________________________________________________________________
%
%                   Additional graphics functions
%__________________________________________________________________________

function drawfield(map,ax,cmap,maxrate)
maxrate = ceil(maxrate);
if maxrate < 1
    maxrate = 1;
end    
n = size(map,1);
plotmap = ones(n,n,3);
for jj = 1:n
   for ii = 1:n
      if isnan(map(jj,ii))
         plotmap(jj,ii,1) = 1;
         plotmap(jj,ii,2) = 1;
		 plotmap(jj,ii,3) = 1;
      else
         rgb = pixelcolour(map(jj,ii),maxrate,cmap);
         plotmap(jj,ii,1) = rgb(1);
         plotmap(jj,ii,2) = rgb(2);
		 plotmap(jj,ii,3) = rgb(3);
      end   
   end
end   
image(ax,ax,plotmap);
set(gca,'YDir','Normal');
s = sprintf('%s%u%s','Peak rate ',maxrate,' Hz');
title(s);

axis image;
set(gcf,'color',[1 1 1]);
drawnow


function adjustaxis(minx, maxx, miny, maxy)
axis equal;
axis([minx maxx miny maxy]);


function rgb = pixelcolour(map,maxrate,cmap)
cmap1 = ...
    [    0         0    0.5625; ...
         0         0    0.6875; ...
         0         0    0.8125; ...
         0         0    0.9375; ...
         0    0.0625    1.0000; ...
         0    0.1875    1.0000; ...
         0    0.3125    1.0000; ...
         0    0.4375    1.0000; ...
         0    0.5625    1.0000; ...
         0    0.6875    1.0000; ...
         0    0.8125    1.0000; ...
         0    0.9375    1.0000; ...
    0.0625    1.0000    1.0000; ...
    0.1875    1.0000    0.8750; ...
    0.3125    1.0000    0.7500; ...
    0.4375    1.0000    0.6250; ...
    0.5625    1.0000    0.5000; ...
    0.6875    1.0000    0.3750; ...
    0.8125    1.0000    0.2500; ...
    0.9375    1.0000    0.1250; ...
    1.0000    1.0000         0; ...
    1.0000    0.8750         0; ...
    1.0000    0.7500         0; ...
    1.0000    0.6250         0; ...
    1.0000    0.5000         0; ...
    1.0000    0.3750         0; ...
    1.0000    0.2500         0; ...
    1.0000    0.1250         0; ...
    1.0000         0         0; ...
    0.8750         0         0; ...
    0.7500         0         0; ...
    0.6250         0         0 ];

cmap2 = ...
   [0.0417         0         0; ...
    0.1250         0         0; ...
    0.2083         0         0; ...
    0.2917         0         0; ...
    0.3750         0         0; ...
    0.4583         0         0; ...
    0.5417         0         0; ...
    0.6250         0         0; ...
    0.7083         0         0; ...
    0.7917         0         0; ...
    0.8750         0         0; ...
    0.9583         0         0; ...
    1.0000    0.0417         0; ...
    1.0000    0.1250         0; ...
    1.0000    0.2083         0; ...
    1.0000    0.2917         0; ...
    1.0000    0.3750         0; ...
    1.0000    0.4583         0; ...
    1.0000    0.5417         0; ...
    1.0000    0.6250         0; ...
    1.0000    0.7083         0; ...
    1.0000    0.7917         0; ...
    1.0000    0.8750         0; ...
    1.0000    0.9583         0; ...
    1.0000    1.0000    0.0625; ...
    1.0000    1.0000    0.1875; ...
    1.0000    1.0000    0.3125; ...
    1.0000    1.0000    0.4375; ...
    1.0000    1.0000    0.5625; ...
    1.0000    1.0000    0.6875; ...
    1.0000    1.0000    0.8125; ...
    1.0000    1.0000    0.9375];
if strcmp(cmap,'jet')
   steps = (31*(map/maxrate))+1;
   steps = round(steps);
   if steps>32; steps = 32; end
   if steps<1; steps = 1; end
   rgb = cmap1(steps,:);
else
   steps = (31*(map/maxrate))+1;
   steps = round(steps);
   if steps>32; steps = 32; end
   if steps<1; steps = 1; end
   rgb = cmap2(steps,:);
end    