% Speed cell analysis (Kropff et al., 2015; Iwase et al., in prep)
% 
% The code is modified from spaceS.m
% 
% INPUT
%   path: path to the folder including npy/whl files.
%   trange: time range for data analysis.
% 
% OUTPUT
%
% Takuma Kitanishi, Osaka Univ 2015 / Osaka City Univ 2018
%   last modified 2019/12/06


% function spaceS()
% 
trange = [10 1210];

% io('tk0056-171014-02',trange);
% io('tk0056-171015-02',trange);
% 
% io('tk0062-171129-02',trange);
% io('tk0062-171130-02',trange);
% 
% io('tk0064-180206-02',trange);
% io('tk0064-180207-02',trange);
% 
% io('tk0067-180426-02',trange);
% io('tk0067-180427-02',trange);
% 
% io('tk0068-180530-02',trange);
% io('tk0068-180531-02',trange);
% 
% io('tk0069-180726-02',trange);
% io('tk0069-180727-02',trange);
% 
% io('tk0070-180829-02',trange);
% io('tk0070-180830-02',trange);
% 
% io('tk0072-181025-02',trange);
% io('tk0072-181026-02',trange);
% 
% io('tk0074-181219-02',trange);
% io('tk0074-181220-02',trange);
% 
% io('tk0075-190326-02',trange);
% io('tk0075-190327-02',trange);
% 
% io('tk0076-190424-02',trange);
% io('tk0076-190425-02',trange);

speedCells_stack('D:\data\rec\02S\speed\sd256ms\tk*.mat','D:\data\rec\02S\speed\sd256ms\')

run('speedCells_sum.m')


function io(session,trange)

idx = strfind(session,'-');
rat = session(1:idx(1)-1);
day = session(idx(1)+1:idx(2)-1);
sess= session(idx(2)+1:end);

% Parameters_______________________________________________________________
fprintf('%s%s\n','Session name: ',session);

p.basepath = fullfile('D:\data\rec\',rat,[day '-' sess]);
p.spkpath = fullfile('D:\data\rec\',rat,[day '_sorted']);

% video pixel to cm conversion (cm/pixels)
p.pixel2cm = 0.5217;

% video sampling rate (Hz, /s)
p.fsv = 39.0625;

% speed (cm/s), threshold and bin size
p.vTh1 = 2;
p.vTh2 = 50;
p.sbin = 4;

% tangential acceleration (cm/s^2), threshold and bin size
p.atTh1 = -40;
p.atTh2 = 40;
p.atbin = 10;

% normal acceleration (cm/s^2), threshold and bin size
p.anTh1 = -30;
p.anbin = 5;
p.anTh2 = 30;

% firing rate bin size (Hz) (for calculating speed information)
p.rbin = 1;

% SD of gaussian for position/firing rate smoothing (sec)
p.sd = 0.256;

p.savepath = ['D:\data\rec\02S\speed\sd' num2str(p.sd*1000) 'ms'];

% temporal shift max (sec)
p.tShiftMax = 1.536; 

% number of shuffles (for temporal shift)
p.nShuffle = 100;

% number of shuffles (for calculating information etc)
p.nShuffle2 = 1000;

% minimum time shift for shuffuling (sec)
p.minShift = 30;

%__________________________________________________________________________

%% position

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
t = t-trange(1);

% Use the better-tracking LED
if sum(x1==-1)>sum(x2==-1)
    posx = x2;
    posy = y2;
else
    posx = x1;
    posy = y1;
end
post = t;

% Centre the box in the coordinate system
centre = centreBox(posx,posy);
posx = posx - centre(1);
posy = posy - centre(2);

% Smooth pos first (then calculate velocity)
w = normpdf(-p.sd*3:1/p.fsv:p.sd*3,0,p.sd);
w = w./sum(w);
posxSm = conv(posx,w,'same');
posySm = conv(posy,w,'same');

% instanteneous speed, smoothed (cm/s)
vxSm = diff(posxSm)./diff(post); vxSm = [vxSm(1); vxSm];
vySm = diff(posySm)./diff(post); vySm = [vySm(1); vySm];
vSm = sqrt( vxSm.^2 + vySm.^2 );

% instanteneous acceleration (cm/s^2)
axSm = diff(vxSm)./diff(post); axSm = [axSm(1); axSm];
aySm = diff(vySm)./diff(post); aySm = [aySm(1); aySm];

% Úü‰Á‘¬“x (tangential acceleration)
at = dot([axSm aySm],[vxSm vySm],2) ./ vSm;
% at = diff(v)./diff(post);
% at = [a(1); a];

% –@ü‰Á‘¬“x (Normal acceleration, centripetal acceleration)
anx = axSm - dot([axSm aySm],[vxSm vySm],2)./vSm.^2 .* vxSm;
any = aySm - dot([axSm aySm],[vxSm vySm],2)./vSm.^2 .* vySm;
an = sqrt(anx.^2 + any.^2);
% whether the Normal aceleration is leftward (negative) or rightward (positive)
idx = vxSm.*aySm - vySm.*axSm>0;
an(idx) = -an(idx);


%% load spikes

[ts,clu,cluOri,Fs,roi,sh,ch] = loadSpkHist(p.spkpath,[],session);
idx = trange(1)<=ts & ts<trange(2);
ts(~idx)=[];
clu(~idx)=[];
ts = ts-trange(1);


%% call main func for each cluster
u = struct;
for ii=1:max(clu)
    % roi name 
    roiName = getRoiName(roi(ii));
    
    u = mainFunc(trange,post,vSm,at,an,ts(clu==ii),session,ii,roiName,p,u);
end

% save
save(fullfile(p.savepath,session),'session','post','vSm','u','p');

end


function u = mainFunc(trange,post,v,at,an,ts,session,unit,roiName,p,u)

fprintf('%s%u\n','  Analyzing unit ',unit);
imname = fullfile(p.savepath,sprintf('%s%s%u%s',session,'_u',unit));


%% Speed tuning analysis (Kropff Nature 2015; Iwase et al, submitted)

% instantaneous firing rate (Hz)
edges = [post(1)-1/p.fsv/2; (post(1:end-1)+post(2:end))./2; post(end)+1/p.fsv/2];
r = histcounts(ts,edges) * p.fsv;

% Smooth instantaneous rate 
w = normpdf(-p.sd*3:1/p.fsv:p.sd*3,0,p.sd);
w = w./sum(w);

r = conv(r,w,'same');
r = r(:);

% remove too low/too high speed periods
remIdx = v<p.vTh1 | p.vTh2<v;
vTrim = v(~remIdx);
rTrim = r(~remIdx);

% speed score (real)
spScore = corrcoef(vTrim,rTrim);
spScore = spScore(2);

% speed score (first/second halves)
ind = round(length(v)/2);

idx1 = remIdx; idx1(ind+1:end) = true;
spScore1 = corrcoef(v(~idx1),r(~idx1));
spScore1 = spScore1(2);

idx2 = remIdx; idx2(1:ind) = true;
spScore2 = corrcoef(v(~idx2),r(~idx2));
spScore2 = spScore2(2);

% speed score (temporal shift)
spScoreShift = nan(1,p.nShuffle);
shiftAxis = -p.tShiftMax:1/p.fsv:p.tShiftMax;
nShift = length(shiftAxis);

for ii=1:nShift
    % shifting
    tShift = shiftAxis(ii);
    indShift = round(tShift * p.fsv);
    rShift = circshift(r,indShift);
    
    % remove too low/too high speed periods
    rTrim2 = rShift(~remIdx);
    
    % speed score (real)
    spScoreTmp = corrcoef(vTrim,rTrim2);
    spScoreShift(ii) = spScoreTmp(2);
end
% normalized correlation
normCorr = spScoreShift./spScoreShift(shiftAxis==0);
% preferred temporal shift (sec)
[~,idxShift] = max(normCorr);
prefShift = shiftAxis(idxShift); 

% speed tuning curve
vEdges = p.vTh1:p.sbin:p.vTh2;
vAxis = (vEdges(1:end-1)+vEdges(2:end))/2;
for ii=1:length(vAxis)
    idx = vEdges(ii)<=vTrim & vTrim<vEdges(ii+1);
    spRate(ii) = mean(rTrim(idx));
end

% speed information (Iwase et al., submitted)
pdfV = histcounts(vTrim,vEdges);
pdfV = pdfV./sum(pdfV);

Ispk = sum( pdfV .* spRate/mean(rTrim) .* log2(spRate/mean(rTrim)) );
Isec = sum( pdfV .* spRate .* log2(spRate/mean(rTrim)) );

% speed score (shuffled)
spScoreShuff = nan(1,p.nShuffle2);
IspkShuff = nan(1,p.nShuffle2);
IsecShuff = nan(1,p.nShuffle2);
for ii=1:p.nShuffle2
    % shuffling
    tShift = rand()*(max(post)-p.minShift*2)+p.minShift;
    indShift = round(tShift * p.fsv);
    rShuffle = circshift(r,indShift);
    
    % remove too low/too high speed periods
    rTrim2 = rShuffle(~remIdx);
    
    % speed score (shuffled)
    spScoreTmp = corrcoef(vTrim,rTrim2);
    spScoreShuff(ii) = spScoreTmp(2);
        
    % speed tuning curve
    spRateTmp = [];
    for jj=1:length(vAxis)
        idx = vEdges(jj)<=vTrim & vTrim<vEdges(jj+1);
        spRateTmp(jj) = mean(rTrim2(idx));
    end
    % speed information (Shuffled) (Iwase et al., submitted)
    IspkShuff(ii) = sum( pdfV .* spRateTmp/mean(rTrim2) .* log2(spRateTmp/mean(rTrim2)) );
    IsecShuff(ii) = sum( pdfV .* spRateTmp .* log2(spRateTmp/mean(rTrim2)) );
end
% surrogate-normalized speed information
normIspk = (Ispk - mean(IspkShuff))./std(IspkShuff);
normIsec = (Isec - mean(IsecShuff))./std(IsecShuff);

% speed slope
tmp = polyfit(vTrim,rTrim,1);
spSlope = tmp(1);
spIntrcept = tmp(2);

% speed tuning curve
vEdges = p.vTh1:p.sbin:p.vTh2;
vAxis = (vEdges(1:end-1)+vEdges(2:end))/2;
for ii=1:length(vAxis)
    idx = vEdges(ii)<=vTrim & vTrim<vEdges(ii+1);
    spRate(ii) = mean(rTrim(idx));
end

% speed information (Iwase et al., submitted)
pdfV = histcounts(vTrim,vEdges);
pdfV = pdfV./sum(pdfV);

Ispk = sum( pdfV .* spRate/mean(rTrim) .* log2(spRate/mean(rTrim)) );
Isec = sum( pdfV .* spRate .* log2(spRate/mean(rTrim)) );

% z score function
zscore = @(x) (x-mean(x))/std(x);

% plot
figure(1); clf
set(gcf,'Position',[80 200 800 640])

subplot(2,4,[1 2])
idx = 60<post & post<90;
plot(post(idx),zscore(v(idx)),'Color',[0 .4470 .7410]); hold on
plot(post(idx),zscore(r(idx)),'Color',[.0 .0 .0])
box off
xlabel('Time (sec)')
ylabel('Speed (blue), firing rate (black)(z-score)')
title(sprintf('%s unit %u (%s)',session,unit,roiName))

subplot(2,4,3)
plot(vAxis,spRate)
xlim([0 p.vTh2])
ylim([0 max([ceil(max(spRate)) 1])])

box off
xlabel('Speed (cm/s)')
ylabel('Firing rate (Hz)')
title(sprintf('Speed score = %.4f',spScore))

subplot(2,4,4)
plot([0 0],[0 2],'k:',shiftAxis,normCorr,'',prefShift,normCorr(idxShift),'+')
xlim([-p.tShiftMax p.tShiftMax])
ylim([0 2])
box off
xlabel('Time (sec)')
ylabel('Normalized correlation')
title(sprintf('Preferred shift (ms) = %.0f',prefShift*1000))


%% Acceleration analysis

subplot(2,4,[5 6])
idx = 60<post & post<90;
plot(post(idx),zscore(at(idx)),'Color',[0.9290 0.6940 0.1250]); hold on
plot(post(idx),zscore(an(idx)),'Color',[0.4940 0.1840 0.5560])
box off
xlabel('Time (sec)')
ylabel('Tang acc (yellow), norm acc (purple)(z)')
title(sprintf('%s unit %u (%s)',session,unit,roiName))


% tangential accel score (real)
remIdx = isnan(at) | at<p.atTh1 | p.atTh2<at;
atTrim = at(~remIdx);
rTrim  = r(~remIdx);

atScore = corrcoef(atTrim,rTrim);
atScore = atScore(2);

% tuning curve
atEdges = p.atTh1:p.atbin:p.atTh2;
atAxis = (atEdges(1:end-1)+atEdges(2:end))/2;
for ii=1:length(atAxis)
    idx = atEdges(ii)<=at & at<atEdges(ii+1);
    atRate(ii) = mean(r(idx));
end

subplot(2,4,7)
plot(atAxis,atRate,'Color',[0.9290 0.6940 0.1250])
xlim([p.atTh1 p.atTh2])
ylim([0 max([1 ceil(max(atRate))])])
title(sprintf('at score = %.4f',atScore))
box off
xlabel('Tangential accel (cm/s^2)')
ylabel('Firing rate (Hz)')


% normal accel score (real)
remIdx = isnan(an) | an<p.anTh1 | p.anTh2<an;
anTrim = an(~remIdx);
rTrim  = r(~remIdx);

anScore = corrcoef(anTrim,rTrim);
anScore = anScore(2);

% accel tuning curve
anEdges = p.anTh1:p.anbin:p.anTh2;
anAxis = (anEdges(1:end-1)+anEdges(2:end))/2;
for ii=1:length(anAxis)
    idx = anEdges(ii)<=an & an<anEdges(ii+1);
    anRate(ii) = mean(r(idx));
end

subplot(2,4,8)
plot(anAxis,anRate,'Color',[0.4940 0.1840 0.5560])
xlim([p.anTh1 p.anTh2])
ylim([0 max([1 ceil(max(anRate))])])
title(sprintf('an score = %.4f',anScore))
box off
xlabel('Normal accel (cm/s^2)')
ylabel('Firing rate (Hz)')


%% save image
print(gcf,imname,'-dpng')
print(gcf,imname,'-dsvg')

% register
% u(unit).session = session;
u(unit).roiName = roiName;
u(unit).rate    = r;    % firing rate along time
% u(unit).meanRate= length(ts)/(max(post)-min(post));
u(unit).meanRate= length(ts)/(max(trange)-min(trange)); % 200310 debug
u(unit).spAxis  = vAxis;
u(unit).spRate  = spRate;
u(unit).spScore = spScore;
u(unit).spScore1 = spScore1;    % 1st half 
u(unit).spScore2 = spScore2;    % 2nd half 
u(unit).spScoreShuffle = spScoreShuff;
u(unit).spSlope = spSlope;
u(unit).spIntrcept = spIntrcept;
u(unit).Ispk = Ispk;
u(unit).Isec = Isec;
u(unit).normIspk    = normIspk;
u(unit).normIsec    = normIsec;
u(unit).IspkShuffle = IspkShuff;
u(unit).IsecShuffle = IsecShuff;
u(unit).Isec = Isec;
u(unit).shiftAxis = shiftAxis;
u(unit).shiftPref = prefShift;
u(unit).shiftNormCorr = normCorr;
u(unit).atAxis  = atAxis;
u(unit).atRate  = atRate;
u(unit).atScore = atScore;
u(unit).anAxis  = anAxis;
u(unit).anRate  = anRate;
u(unit).anScore = anScore;


end


function speedCells_stack(matpath,savepath)

% savepath = 'D:\data\rec\02S\speed';

matList = dir(matpath);

uall = [];
for ii=1:length(matList)
    fprintf('Stacking: %s\n',matList(ii).name)
    
    load(fullfile(matList(ii).folder,matList(ii).name),'u','session')
    
    % add session and unit num
    for jj=1:length(u)
        u(jj).session = session;
        u(jj).unit = jj;
    end
    
    % reorder field names and stack fields
    nField = length(fieldnames(u));
    perim = [nField-1 nField 1:nField-2];
    
    uall =[uall, orderfields(u,perim)];
end

% save
save(fullfile(savepath,'speedCells.mat'),'uall')
    

end


%__________________________________________________________________________
%
%                  Functions
%__________________________________________________________________________


function H = getEntropy(pdf)

% linearlize
pdf = pdf(:);

% normalize (just in case)
pdf = pdf./sum(pdf);

% remove 0
pdf(pdf==0) = [];

% entropy 
H = sum(-pdf.*log2(pdf));

end

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

end

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

end