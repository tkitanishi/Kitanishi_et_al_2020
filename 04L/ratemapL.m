% Rate map/theta phase locking analysis of linear track
% 
% OUPUT
%   u.
%     roiName
%     nspk1         spike number map (trial-by-trial) (1=L2R, 2=R2L trials)
%     rate1         rate map (trial-by-trial)
%     r1, r2        rate map, mean
%     e1,r2         rate map, sem
%     meanRate1     mean rate
%     infoPerSec1   spatial information
%     infoPerSpk1   spatial information
%     sparsity1     sparsity
%     selectivity1  spatial selectivity
%     peakRate1     peak rate
%     phcount1      spike theta phase histogram (counts)
%     cPrefPh1      preferred spike theta phase (deg)
%     cR1           resultant vector length of spike theta phase
%     cPPC1         pairwise phase consistency of spike theta phase
%     cSTD1         circular standard deviation of spike theta phase
%     cPval         Rayleigh test p-value for spike theta phase
%     nspkPhmap1    pos-theta phase spike map
%     ratePhmap1    pos-theta phase rate map
% 
%   tMap1,2:    time occupancy map (sec) (1=left trial, 2=right trial)
%   vMap1,2:    velocity map (cm/s) (1=left trial, 2=right trial)
% 
% Takuma Kitanishi, OCU, 2018/07/18
%   modifed on 2019/07/10, bug fix of tMap1

clear; 
fclose all;

mainFunc('tk0056-171014-04','Buzsaki256.mat',[6 15],[0 1200]+20)
mainFunc('tk0056-171015-04','Buzsaki256.mat',[6 15],[0 1200]+15)

mainFunc('tk0062-171129-04','Buzsaki256.mat',[2 11],[0 1200]+10)
mainFunc('tk0062-171130-04','Buzsaki256.mat',[2 11],[0 1200]+30)

mainFunc('tk0064-180206-04','Buzsaki256.mat',[5 14],[0 1200]+15)
mainFunc('tk0064-180207-04','Buzsaki256.mat',[5 14],[0 1200]+15)

mainFunc('tk0067-180426-04','Buzsaki256.mat',[5 15],[0 1200]+10)
mainFunc('tk0067-180427-04','Buzsaki256.mat',[5 15],[0 1200]+20)

mainFunc('tk0068-180530-04','A8x32Edge.mat', [5 17],[0 1200]+20)
mainFunc('tk0068-180531-04','A8x32Edge.mat', [5 17],[0 1200]+15)

mainFunc('tk0069-180726-04','A8x32_5mm_35_300_160.mat',[5 18],[0 1200]+10)
mainFunc('tk0069-180727-04','A8x32_5mm_35_300_160.mat',[5 18],[0 1200]+40)

mainFunc('tk0070-180829-04','A8x32_5mm_35_300_160.mat',[3 14],[0 1200]+12)
mainFunc('tk0070-180830-04','A8x32_5mm_35_300_160.mat',[3 14],[0 1200]+10)

mainFunc('tk0072-181025-04','A8x32Edge.mat',[5 14],[0 1200]+10)
mainFunc('tk0072-181026-04','A8x32Edge.mat',[5 14],[0 1200]+12)

mainFunc('tk0074-181219-04','A8x32Edge.mat',[1 19],[0 1200]+10)
mainFunc('tk0074-181220-04','A8x32Edge.mat',[1 19],[0 1200]+20)

mainFunc('tk0075-190326-04','A8x32Edge.mat',[4 14],[0 1200]+10)
mainFunc('tk0075-190327-04','A8x32Edge.mat',[4 14],[0 1200]+15)

mainFunc('tk0076-190424-04','A8x32Edge.mat',[6 8],[0 1200]+10)
mainFunc('tk0076-190425-04','A8x32Edge.mat',[6 8],[0 1200]+15)


function mainFunc(session,mapfile,refch,trange)

%% config 
p.session = session;
idx = strfind(session,'-');
p.rat = session(1:idx(1)-1);
p.day = session(idx(1)+1:idx(2)-1);
p.sess= session(idx(2)+1:end);
p.mapfile = mapfile;
p.refch = refch;
p.trange = trange;

% file path
p.savePath = 'D:\data\rec\04L\ratemap';
p.basePath = fullfile('D:\data\rec',p.rat,[p.day '-' p.sess]);
p.spkPath  = fullfile('D:\data\rec',p.rat,[p.day '_sorted']);

% position scaling
p.pixel2cm = 0.4130;
p.centerX = 422*p.pixel2cm;
p.centerY = 230*p.pixel2cm;

% linear track half length (cm) used for analysis
% (outside this length is omitted from the analysis)
p.halfLen = 100; 

% position bin (cm)
p.posBin = 2;
p.nbin = p.halfLen*2/p.posBin;
p.posbins = -p.halfLen:p.posBin:p.halfLen;
p.posaxis = -p.halfLen+p.posBin/2:p.posBin:p.halfLen-p.posBin/2;

% postion sampling rate after upsampling (Hz)
p.fsBehav = 100;

% phase bin (deg)
p.phbins = 0:15:360;
p.phaxis = (p.phbins(1:end-1)+p.phbins(2:end))/2;

% LFP sampling rate [Hz]
p.fs = 1250;

% theta range [Hz]
p.f1 = 4;   % stop band
p.f2 = 5;   % pass band
p.f3 = 10;  % pass band
p.f4 = 11;  % stop band

% smoothing window for time occupancy map & spike number map (bins)
p.w1 = gaussian(15,2)/sum(gaussian(15,2));

% smoothing window for position-phase space
p.w2 = customgauss([9 9], 2, 2, 0, 0, 1, [0 0]);
p.w2 = p.w2./sum(p.w2(:));

% smoothing window for mutual window calculation
p.sd = 0.128;   % [s]
p.w3 = normpdf(-p.sd*3:1/p.fsBehav:p.sd*3,0,p.sd);
p.w3 = p.w3./sum(p.w3);

% firing rate bin number for mutual information (4)
p.nqtile = 4;

% minimum rate (for visualization only)
p.minRate = 5;

% SD of gaussian for position/firing rate smoothing (sec)
p.sd = 0.128;

% number of shuffling to calculate chance level of spatial information
p.nshuff = 1000;


%% load position

% load pos (cm)
[postTmp,posxTmp,posyTmp] = loadWhl(p.basePath,p.pixel2cm);

% interporate and upsampling position
post = postTmp(1):1/p.fsBehav:postTmp(end);
posx = interp1(postTmp,posxTmp,post);
posy = interp1(postTmp,posyTmp,post);

% post edge
dt = 1/p.fsBehav;
tedge = [post(1)-dt/2 post+dt/2];

% velocity (cm/s)
% w1 = gausswin(11)/sum(gausswin(11));
% v = conv( sqrt(gradient(posx).^2+gradient(posy).^2)*p.fsBehav,w1,'same');

% Smooth pos first (then calculate velocity)
w = normpdf(-p.sd*3:1/p.fsBehav:p.sd*3,0,p.sd); w = w./sum(w);
posxSm = conv(posx,w,'same');
posySm = conv(posy,w,'same');
% instanteneous speed, smoothed (cm/s)
vxSm = diff(posxSm)./diff(post); vxSm = [vxSm(1) vxSm];
vySm = diff(posySm)./diff(post); vySm = [vySm(1) vySm];
v = sqrt( vxSm.^2 + vySm.^2 );

% trim and center data
idx = trange(1)<post & post<trange(2);
post = post(idx)-trange(1);
posx = posx(idx)-p.centerX;
posy = -(posy(idx)-p.centerY);
v = v(idx);

% linearlized pos
pos = posx;


%% trial definition

[state1,trial1,n(1)] = defineTrials(pos,p.halfLen,1);
[state2,trial2,n(2)] = defineTrials(-pos,p.halfLen,2);
state = max([state1;state2]);
trial = max([trial1;trial2]);

% plot
figure(1); clf
subplot(3,2,[1 5])
plot(pos,post,'k'); hold on
plot(pos(state==1),post(state==1),'b.')
plot(pos(state==2),post(state==2),'r.')
plot([-p.halfLen -p.halfLen],[0 range(trange)],'k:')
plot([p.halfLen p.halfLen],[0 range(trange)],'k:')
title(session)
xlabel('Position (cm)'); ylabel('Time (sec)')
set(gca,'YDir','reverse')
box off


%% Speed map (cm/s)

vMap1 = vmap(v,post,state,1,trial,n(1),pos,p.halfLen,p.posBin);
vMap2 = vmap(v,post,state,2,trial,n(2),pos,p.halfLen,p.posBin);

% plot
clims = [ 0 min([200 ceil(max([vMap1(:);vMap2(:)]))+10 ])];

subplot(3,2,2)
imagesc(p.posaxis,1:n(1),vMap1,clims)
colorbar; colormap jet;
ylabel('L-to-R trials')
title('Speed map (cm/s)')

subplot(3,2,4)
imagesc(p.posaxis,1:n(2),vMap2,clims)
colorbar; colormap jet;
ylabel('R-to-L trials')

subplot(3,2,6)
v1 = mean(vMap1,'omitnan'); 
v2 = mean(vMap2,'omitnan'); 
ve1 = std(vMap1,'omitnan')/sqrt(n(1));
ve2 = std(vMap2,'omitnan')/sqrt(n(2));
plot(p.posaxis,v1,'b',p.posaxis,v2,'r'); hold on
plot(p.posaxis,v1-ve1,'b:',p.posaxis,v1+ve1,'b:',p.posaxis,v2-ve2,'r:',p.posaxis,v2+ve2,'r:');
% ylim([0 150])
box off
% colorbar
xlabel('Position (cm)')
ylabel('Speed (cm/s)')
legend({'LtoR','RtoL'},'Location','eastoutside')

% save
set(gcf,'Position',[80 700 800 600])
print(figure(1),fullfile(p.savePath,session),'-dpng')
print(figure(1),fullfile(p.savePath,session),'-dsvg')


%% time occupancy map (trial-by-trial)

tMap1 = tmap(post,state,1,trial,n(1),pos,p.halfLen,p.posBin,p.fsBehav);
tMap2 = tmap(post,state,2,trial,n(2),pos,p.halfLen,p.posBin,p.fsBehav);

% smoothed
tMap1sm = conv2(tMap1,p.w1,'same');
% remove 0
nbin = size(tMap1sm,2);
[rmap,c]=find(tMap1sm==0);
for ii=1:length(rmap)
    if c==1
        tMap1sm(rmap,c:c+1) = repmat( mean(tMap1sm(rmap,c:c+1)) ,1,2 );
    elseif c==nbin
        tMap1sm(rmap,c-1:c) = repmat( mean(tMap1sm(rmap,c-1:c)) ,1,2 );
    else
        tMap1sm(rmap,c-1:c+1) = repmat( mean(tMap1sm(rmap,c-1:c+1)) ,1,3 );
    end
end

% smoothed
tMap2sm = conv2(tMap2,p.w1,'same');
% remove 0
nbin = size(tMap2sm,2);
[rmap,c]=find(tMap2sm==0);
for ii=1:length(rmap)
    if c==1
        tMap2sm(rmap,c:c+1) = repmat( mean(tMap2sm(rmap,c:c+1)) ,1,2 );
    elseif c==nbin
        tMap2sm(rmap,c-1:c) = repmat( mean(tMap2sm(rmap,c-1:c)) ,1,2 );
    else
        tMap2sm(rmap,c-1:c+1) = repmat( mean(tMap2sm(rmap,c-1:c+1)) ,1,3 );
    end
end


%% load spikes

% load spike timing (sec)
[ts,clu,cluOri,Fs,roi,sh,ch]=loadSpkHist(p.spkPath,[],session);
idx = trange(1)<ts & ts<trange(2);
ts = ts(idx)-trange(1);
clu = clu(idx);
cluOri = cluOri(idx);

% spike position (cm) [1 x length(ts)] 
idx = round((ts-post(1))*p.fsBehav)+1;
idx(idx==0) = 1;
idx(idx>length(pos)) = length(pos);
spkpos = pos(idx);


%% load 1250-Hz lfp

% load lfp
[lfp,t,~] = loadLfp(p.basePath,[],mapfile);

idx = trange(1)<=t & t<=trange(2);
lfp = lfp(:,idx);
t   = t(idx)-trange(1);

% theta band
lfpt = fftbandpass(lfp((refch(1)-1)*32+refch(2),:),p.fs,p.f1,p.f2,p.f3,p.f4);

% theta phase
[phaseT,~] = lfpphase(lfpt);  % [degree,~]

% find spike theta phase (deg) [1 x length(ts)] 
idx = round((ts-t(1))*p.fs+1);
idx(idx==0) = 1;
spkph = phaseT(idx);

% find position theta phase (deg) [1 x length(post)]
idx = round((post-t(1))*p.fs+1);
idx(idx==0) = 1;
% idx(idx>length(phaseT)) = [];
posph = phaseT(idx);


%% position sample map in position-phase space

% left trial
for j=1:n(1)
    % trial start/end timing (sec)
    t1 = min(post(state==1 & trial==j));
    t2 = max(post(state==1 & trial==j));
    % cumulative spike pos
    if j==1
        pospos1 = pos(t1<post & post<t2);
        posph1 = posph(t1<post & post<t2);
    else
        pospos1 = [pospos1 pos(t1<post & post<t2)];
        posph1  = [posph1 posph(t1<post & post<t2)];
    end
end
% right trial
for j=1:n(2)
    % trial start/end timing (sec)
    t1 = min(post(state==2 & trial==j));
    t2 = max(post(state==2 & trial==j));
    % cumulative spike pos
    if j==1
        pospos2 = pos(t1<post & post<t2);
        posph2 = posph(t1<post & post<t2);
    else
        pospos2 = [pospos2 pos(t1<post & post<t2)];
        posph2  = [posph2 posph(t1<post & post<t2)];
    end
end

% occupied time map in position-phase space
tPhmap1 = zeros(length(p.phaxis),length(p.posaxis));
for ii=1:size(tPhmap1,1)
    for jj=1:size(tPhmap1,2)
        tPhmap1(ii,jj) = sum(p.phbins(ii)<posph1 & posph1<=p.phbins(ii+1) & p.posbins(jj)<pospos1 & pospos1<=p.posbins(jj+1))./p.fsBehav;
    end
end
tPhmap2 = zeros(length(p.phaxis),length(p.posaxis));
for ii=1:size(tPhmap2,1)
    for jj=1:size(tPhmap2,2)
        tPhmap2(ii,jj) = sum(p.phbins(ii)<posph2 & posph2<=p.phbins(ii+1) & p.posbins(jj)<pospos2 & pospos2<=p.posbins(jj+1))./p.fsBehav;
    end
end


%% for each unit
disp('  Calculating rate maps...')

for i=1:max(clu)
    switch roi(i)
        case 1
            u(i).roiName= 'SUB';
        case 3 
            u(i).roiName= 'CA1';
        case 4
            u(i).roiName= 'DG';
        otherwise
            u(i).roiName= 'na';
    end
end

% for each unit
posTrim1 = pos(state==1);     % for mutual information
posTrim2 = pos(state==2);     % for mutual information

for i=1:max(clu)
    fprintf('  Unit %u\n',i)
    spkposTmp = spkpos(clu==i);
    spkphTmp = spkph(clu==i);
    tsTmp = ts(clu==i);

    % ----- mutual inforamtion -----
    
    % time-binned rate
    fr = histcounts(tsTmp,tedge)./diff(tedge);
    % smoothing
    fr = conv(fr,p.w3,'same');
    % remove track edge data points
    frTrim1  = fr(state==1);
    frTrim2  = fr(state==2);
    
    % firing rate bins
    fredge1 = prctile(frTrim1(frTrim1~=0), 0:100/p.nqtile:100); fredge1(1) = 0;
    fredge2 = prctile(frTrim2(frTrim2~=0), 0:100/p.nqtile:100); fredge2(1) = 0;
    
    % mutual information
    if any(isnan(fredge1))
        u(i).MI1 = nan;
    else
        u(i).MI1 = mutualinfo(posTrim1,frTrim1,p.posbins,fredge1);
    end
    if any(isnan(fredge2))
        u(i).MI2 = nan;
    else
        u(i).MI2 = mutualinfo(posTrim2,frTrim2,p.posbins,fredge2);
    end

    % ----- Other params -----
    
    % spike number map (trial-by-trial)
    u(i).nspk1 = spkmap(tsTmp,spkposTmp,post,state,1,trial,n(1),p.halfLen,p.posBin);
    u(i).nspk2 = spkmap(tsTmp,spkposTmp,post,state,2,trial,n(2),p.halfLen,p.posBin);
    
    % ratemap (trial-by-trial)
    u(i).rate1 = conv2(u(i).nspk1,p.w1,'same')./tMap1sm;
    u(i).rate2 = conv2(u(i).nspk2,p.w1,'same')./tMap2sm;
    
    % mean & sem of ratemap along lineartrack
    u(i).r1 = mean(u(i).rate1); 
    u(i).r2 = mean(u(i).rate2); 
    u(i).e1 = sem(u(i).rate1);
    u(i).e2 = sem(u(i).rate2);
    
    % map statistics
    pospdf1 = mean(tMap1) ./ sum(mean(tMap1));  % probability density map
    pospdf2 = mean(tMap2) ./ sum(mean(tMap2));  % probability density map

    [u(i).meanRate1,u(i).infoPerSec1,u(i).infoPerSpk1,u(i).sparsity1,u(i).selectivity1] = mapstat(u(i).r1,pospdf1);
    [u(i).meanRate2,u(i).infoPerSec2,u(i).infoPerSpk2,u(i).sparsity2,u(i).selectivity2] = mapstat(u(i).r2,pospdf2);
    if isinf(u(i).meanRate1); warning('meanRate1 is Inf!'); end
    if isinf(u(i).meanRate2); warning('meanRate2 is Inf!'); end
    
    % shuffle spatial information
    for jj=1:p.nshuff
        [u(i).infoPerSec1shuff(jj),u(i).infoPerSpk1shuff(jj)] = mapstat2(mean(getShuffledMap(u(i).rate1)),pospdf1);
        [u(i).infoPerSec2shuff(jj),u(i).infoPerSpk2shuff(jj)] = mapstat2(mean(getShuffledMap(u(i).rate2)),pospdf2);
    end
    
    % peakrate 
    u(i).peakRate1 = max(u(i).r1);
    u(i).peakRate2 = max(u(i).r2);
    
    
    % ----- Plot-----
    
    % plot ratemap
    clims = [ 0 max([p.minRate ceil(prctile([5;u(i).rate1(:);u(i).rate2(:)],99))+1]) ];
    
    figure(2); clf
    subplot(4,2,1)
    imagesc(p.posaxis,[],u(i).rate1,clims); 
    colorbar
    ylabel('LtoR trial')
    title(sprintf('Rate map: %s unit %u (%s)',session,i,u(i).roiName))
    
    subplot(4,2,3)
    imagesc(p.posaxis,[],u(i).rate2,clims); 
    colorbar
    ylabel('RtoL trial')
    
    subplot(4,2,5)
    plotshaded(p.posaxis,[u(i).r1-u(i).e1; u(i).r1; u(i).r1+u(i).e1],'b'); hold on
    plotshaded(p.posaxis,[u(i).r2-u(i).e2; u(i).r2; u(i).r2+u(i).e2],'r'); hold on
%     plot(p.posaxis,u(i).r1,'b',p.posaxis,u(i).r2,'r'); hold on
%     plot(p.posaxis,u(i).r1-u(i).e1,'b:',p.posaxis,u(i).r1+u(i).e1,'b:',p.posaxis,u(i).r2-u(i).e2,'r:',p.posaxis,u(i).r2+u(i).e2,'r:');
    title(sprintf('Isec(bits/s)=%.2f, %.2f, MI(bits)=%.2f, %.2f (LR, RL)',u(i).infoPerSec1,u(i).infoPerSec2,u(i).MI1,u(i).MI2))
    yMax = ceil(max([4 u(i).r1+u(i).e1 u(i).r2+u(i).e2]))+1;
    ylim([0 yMax])
    ylabel('Firing rate (Hz)')
    colorbar
    box off
   
    
    %% theta phase precession, pos-spike map
    
    % left trial
    for j=1:n(1)
        % trial start/end timing (sec)
        t1 = min(post(state==1 & trial==j));
        t2 = max(post(state==1 & trial==j));
        % cumulative spike pos
        if j==1
            spkpos1 = spkposTmp(t1<tsTmp & tsTmp<t2);
            spkph1 = spkphTmp(t1<tsTmp & tsTmp<t2);
        else
            spkpos1 = [spkpos1 spkposTmp(t1<tsTmp & tsTmp<t2)];
            spkph1 = [spkph1 spkphTmp(t1<tsTmp & tsTmp<t2)];
        end
    end
    % right trial
    for j=1:n(2)
        % trial start/end timing (sec)
        t1 = min(post(state==2 & trial==j));
        t2 = max(post(state==2 & trial==j));
        % cumulative spike pos
        if j==1
            spkpos2 = spkposTmp(t1<tsTmp & tsTmp<t2);
            spkph2 = spkphTmp(t1<tsTmp & tsTmp<t2);
        else
            spkpos2 = [spkpos2 spkposTmp(t1<tsTmp & tsTmp<t2)];
            spkph2 = [spkph2 spkphTmp(t1<tsTmp & tsTmp<t2)];
        end
    end
    % phase histogram
    u(i).phcount1 = hist(spkph1,p.phaxis);
    u(i).phcount2 = hist(spkph2,p.phaxis);
    maxcount = max([u(i).phcount1(:); u(i).phcount2(:)]);
    
    % circular statistics
    spkph1rad = spkph1(:)/180*pi();
    spkph2rad = spkph2(:)/180*pi();
    u(i).cPrefPh1 = circ_mean(spkph1rad) /pi()*180;  % deg
    u(i).cR1      = circ_r(spkph1rad);               % resultant vector length
    u(i).cPPC1    = circ_ppc(spkph1rad);
    u(i).cSTD1    = sqrt(-2*log(u(i).cR1));
    u(i).cPval1   = circ_rtest(spkph1rad);
    
    u(i).cPrefPh2 = circ_mean(spkph2rad) /pi()*180;  % deg
    u(i).cR2      = circ_r(spkph2rad);               % resultant vector length
    u(i).cPPC2    = circ_ppc(spkph2rad);
    u(i).cSTD2    = sqrt(-2*log(u(i).cR2));
    u(i).cPval2   = circ_rtest(spkph2rad);
    
    % plot
    subplot(4,2,2)
    plot([spkpos1; spkpos1],[spkph1; spkph1+360],'b.'); hold on
    colorbar
    xlim([-p.halfLen p.halfLen])
    ylim([0 720]); 
    yticks(0:180:720)
    ylabel('Theta phase')
    title('Position-phase relationship')
    
    subplot(4,2,4)
    plot([spkpos2; spkpos2],[spkph2; spkph2+360],'r.'); hold on
    colorbar
    xlim([-p.halfLen p.halfLen])
    ylim([0 720]); 
    yticks(0:180:720)
    ylabel('Theta phase')
    set(gca,'XDir','reverse')
    
    subplot(4,2,7)
    plot([p.phaxis p.phaxis+360],[u(i).phcount1 u(i).phcount1],'b'); hold on
    plot([p.phaxis p.phaxis+360],[u(i).phcount2 u(i).phcount2],'r')
    YLim = get(gca,'YLim');
    refTheta = (-cos(p.phaxis/180*pi())+1)*YLim(2)/20;
    plot([p.phaxis p.phaxis+360],[refTheta refTheta],'k:')
    xlim([0 720])
    xticks(0:180:720)
    box off
    xlabel('Theta phase (deg)')
    ylabel('Spike count')
    legend({'L2R','R2L'},'location','eastoutside')
    
    
    %% theta phase precession, colormap
    
    % map for spike number and occupied time in position-phase space
    nspkPhmap1 = zeros(length(p.phaxis),length(p.posaxis));
    for ii=1:size(nspkPhmap1,1)
        for jj=1:size(nspkPhmap1,2)
            nspkPhmap1(ii,jj) = sum(p.phbins(ii)<spkph1 & spkph1<=p.phbins(ii+1) & p.posbins(jj)<spkpos1 & spkpos1<=p.posbins(jj+1));
        end
    end
    
    nspkPhmap2 = zeros(length(p.phaxis),length(p.posaxis));
    for ii=1:size(nspkPhmap2,1)
        for jj=1:size(nspkPhmap2,2)
            nspkPhmap2(ii,jj) = sum(p.phbins(ii)<spkph2 & spkph2<=p.phbins(ii+1) & p.posbins(jj)<spkpos2 & spkpos2<=p.posbins(jj+1));
        end
    end
    
    % ratemap
    ratePhmap1 = conv2(repmat(nspkPhmap1,3,1),p.w2,'same')./conv2(repmat(tPhmap1,3,1),p.w2,'same');
    ratePhmap1 = ratePhmap1(length(p.phaxis)+1:2*length(p.phaxis),:);
    
    ratePhmap2 = conv2(repmat(nspkPhmap2,3,1),p.w2,'same')./conv2(repmat(tPhmap2,3,1),p.w2,'same');
    ratePhmap2 = ratePhmap2(length(p.phaxis)+1:2*length(p.phaxis),:);
    
    % plot
    subplot(4,2,6)
    imagesc(p.posaxis,[p.phaxis p.phaxis+360],[ratePhmap1; ratePhmap1])
    colorbar
    colormap jet
    yticks(0:180:720)
    axis xy
    clim([0 max([p.minRate ceil(max(ratePhmap1(:)))])])
    ylabel('Theta phase')
    
    subplot(4,2,8)
    imagesc(p.posaxis,[p.phaxis p.phaxis+360],[ratePhmap2; ratePhmap2])
    colorbar
    colormap jet
    yticks(0:180:720)
    axis xy
    clim([0 max([p.minRate ceil(max(ratePhmap2(:)))])])
    set(gca,'XDir','reverse')
    xlabel('Position (flipped for R2L) (cm)')
    ylabel('Theta phase')
    
    % save
    set(gcf,'Position',[100 300 900 800])
    print(figure(2),fullfile(p.savePath,[session '_' num2str(i)]),'-dpng')
    print(figure(2),fullfile(p.savePath,[session '_' num2str(i)]),'-dsvg')

    
    % register
    u(i).nspkPhmap1 = nspkPhmap1;
    u(i).nspkPhmap2 = nspkPhmap2;
    u(i).ratePhmap1 = ratePhmap1;
    u(i).ratePhmap2 = ratePhmap2;
end

save(fullfile(p.savePath,session),'u','vMap1','vMap2','tMap1','tMap2','tPhmap1','tPhmap2','session','p')

end


%% ---------------------------------------------------------------
% functions
% -------------------------------------------------------------------------

% get position-shuffled rate map
function rmap = getShuffledMap(rmap)

shiftind = randi(size(rmap,2),[size(rmap,1),1]);

for jj=1:size(rmap,1)
    rmap(jj,:) = circshift(rmap(jj,:),shiftind(jj));
end

end


% CUSTOMGAUSS    Generate a custom 2D gaussian
%
%    gauss = customgauss(gsize, sigmax, sigmay, theta, offset, factor, center)
%
%          gsize     Size of the output 'gauss', should be a 1x2 vector
%          sigmax    Std. dev. in the X direction
%          sigmay    Std. dev. in the Y direction
%          theta     Rotation in degrees
%          offset    Minimum value in output
%          factor    Related to maximum value of output, should be 
%                    different from zero
%          center    The center position of the gaussian, should be a
%                    1x2 vector                     
function ret = customgauss(gsize, sigmax, sigmay, theta, offset, factor, center)

ret     = zeros(gsize);
rbegin  = -round(gsize(1) / 2);
cbegin  = -round(gsize(2) / 2);
for r=1:gsize(1)
    for c=1:gsize(2)
        ret(r,c) = rotgauss(rbegin+r,cbegin+c, theta, sigmax, sigmay, offset, factor, center);
    end
end
end
function val = rotgauss(x, y, theta, sigmax, sigmay, offset, factor, center)
xc      = center(1);
yc      = center(2);
theta   = (theta/180)*pi;
xm      = (x-xc)*cos(theta) - (y-yc)*sin(theta);
ym      = (x-xc)*sin(theta) + (y-yc)*cos(theta);
u       = (xm/sigmax)^2 + (ym/sigmay)^2;
val     = offset + factor*exp(-u/2);
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

function [state,trial,n] = defineTrials(pos,halfLen,stateNum)

state = zeros(size(pos));
trial = zeros(size(pos));
n = 0;
idx = 2;
trialFlag = false;

while idx<length(pos)
    % trial start
    if pos(idx-1)<-halfLen && pos(idx)>=-halfLen
        trialFlag = true;
        idx1 = idx;
        idx2 = idx+1;
        while idx2<length(pos)
            % going back to left edge
            if pos(idx2)<-halfLen
                trialFlag = false;
                idx = idx2;
                break;
            end
            if pos(idx2)>halfLen
                trialFlag = true;
                idx = idx2;
                idx2 = idx2-1;
                break;
            end
            idx2 = idx2+1;
        end
        if idx2>=length(pos)
            break;
        end
        if trialFlag
            n = n+1;
            state(idx1:idx2) = stateNum;
            trial(idx1:idx2) = n;
        end
    else
        idx = idx+1;
    end
end

end


% Shannon information, sparseness, and selectivity
function [meanrate,infoPerSec,infoPerSpk,sparsity,selectivity] = mapstat(map,pospdf)

% peakrate = nanmax(map(:));
meanrate = nansum(nansum( map .* pospdf ));
meansquarerate = nansum(nansum( (map.^2) .* pospdf ));

% sparsity
if meansquarerate == 0
    sparsity = NaN;
else
    sparsity = meanrate^2 / meansquarerate;
end

% selectivity
maxrate = max(max(map));
if meanrate == 0
    selectivity = NaN;
else
    selectivity = maxrate/meanrate;
end

% Shannon information
%   infoPerSec, information (bits/sec)
%   infoPerSpk, information (bits/spike)
% Advances in Neural Information Processing Systems 1993. An information-theoretic approach to deciphering the Hippocampal code
[i1, i2] = find( (map>0) & (pospdf>0) );  % the limit of x*log(x) as x->0 is 0
if ~isempty(i1)
    akksum = 0;
    for i = 1:length(i1)
        ii1 = i1(i);
        ii2 = i2(i);
        akksum = akksum + pospdf(ii1,ii2) * (map(ii1,ii2)) * log2( map(ii1,ii2) / meanrate );
    end
    infoPerSec= akksum;
else
    infoPerSec = NaN;
end

if meanrate == 0
    infoPerSpk = NaN;
else
    infoPerSpk = infoPerSec./meanrate;
end

end


% Shannon information only
function [infoPerSec,infoPerSpk] = mapstat2(map,pospdf)

% peakrate = nanmax(map(:));
meanrate = nansum(nansum( map .* pospdf ));

% Shannon information
%   infoPerSec, information (bits/sec)
%   infoPerSpk, information (bits/spike)
% Advances in Neural Information Processing Systems 1993. An information-theoretic approach to deciphering the Hippocampal code
[i1, i2] = find( (map>0) & (pospdf>0) );  % the limit of x*log(x) as x->0 is 0
if ~isempty(i1)
    akksum = 0;
    for i = 1:length(i1)
        ii1 = i1(i);
        ii2 = i2(i);
        akksum = akksum + pospdf(ii1,ii2) * (map(ii1,ii2)) * log2( map(ii1,ii2) / meanrate );
    end
    infoPerSec= akksum;
else
    infoPerSec = NaN;
end

if meanrate == 0
    infoPerSpk = NaN;
else
    infoPerSpk = infoPerSec./meanrate;
end

end


% spike number map (sec)
function nSpk = spkmap(tsTmp,spkposTmp,post,state,stateIdx,trial,n,halfLen,posBin)

nbin = round(2*halfLen/posBin);
posAxis = (1:nbin)*posBin-mean((1:nbin)*posBin);

nSpk = nan(n,nbin);    % spike number map [trial bin]
for jj=1:n
    % trial start/end timing (sec)
    t1 = min(post(state==stateIdx & trial==jj));
    t2 = max(post(state==stateIdx & trial==jj));
    % spike pos
    spkposTmp2 = spkposTmp(t1<tsTmp & tsTmp<t2);
    for kk=1:nbin
        nSpk(jj,kk) = sum(posAxis(kk)-posBin/2<spkposTmp2 & spkposTmp2<=posAxis(kk)+posBin/2);
    end
end

end


% time occupancy map (sec)
function tMap = tmap(post,state,stateIdx,trial,n,pos,halfLen,posBin,fsBehav)

nbin = round(2*halfLen/posBin);
posAxis = (1:nbin)*posBin-mean((1:nbin)*posBin);

tMap = nan(n,nbin); 

% for each trial
for jj=1:n
    % trial start/end timing (sec)
    t1 = min(post(state==stateIdx & trial==jj));
    t2 = max(post(state==stateIdx & trial==jj));
    % animal pos
    posTmp = pos(t1<=post & post<t2);
    for kk=1:nbin
        tMap(jj,kk) = sum(posAxis(kk)-posBin/2<posTmp & posTmp<=posAxis(kk)+posBin/2) /fsBehav;
    end
end

end


% velocity map

function vMap = vmap(v,post,state,stateIdx,trial,n,pos,halfLen,posBin)

nbin = round(2*halfLen/posBin);
posAxis = (1:nbin)*posBin-mean((1:nbin)*posBin);

vMap = nan(n,nbin);

% for each trial
for ii=1:n
    % trial start/end timing (sec)
    t1 = min(post(state==stateIdx & trial==ii));
    t2 = max(post(state==stateIdx & trial==ii));
    for kk=1:nbin
        vMap(ii,kk) = mean(v( (t1<=post & post<t2) & (posAxis(kk)-posBin/2<pos & pos<=posAxis(kk)+posBin/2) ));
    end
end

end

% MI = mutualinfo(x,y,xedge,yedge)
% 
% This code calculates mutual information between two varibales x and y
% (such as position and rate).
% 
% INPUT: 
%   x,      vector variable (such as animal position)
%   y,      vector variable (such as firing rate)
%   xedge,  edges for x binning
%   yedge,  edges for y binning
% 
% OUTPUT:
%   MI,     mutual informtion between x and y (bits)
% 
% REFERENCE:
%   Souza et al., Neuroscience, 375: 62-73 (2018), On information metrics
%   for spatial coding.
% 
% Takuma Kitanishi, OCU, 2020/02/01

function MI = mutualinfo(x,y,xedge,yedge)

% input check
if ~isvector(x) || ~isvector(y) || ~isvector(xedge) || ~isvector(yedge)
    error('All inputs shuould be a vector!');
end

% probability distribution funcitons
nx = histcounts(x,xedge);   
ny = histcounts(y,yedge);
nxy= histcounts2(x,y,xedge,yedge);

px = nx./sum(nx);
py = ny./sum(ny);
pxy= nxy./sum(sum(nxy));

% mutual information
MI = 0;
for i=1:length(px)
    for j=1:length(py)
        mi = pxy(i,j) * log2(pxy(i,j)/(px(i)*py(j)));
        if ~isnan(mi) && ~isinf(mi)
            MI = MI + mi;
        end
    end
end

end