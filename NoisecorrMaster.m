%% Noise correlation script
% call Noisecorr function
% load the Noisecorr files the fucntion made. 
Scriptsdir = 'C:\Users\gillissen\Documents\GitHub\Mouse' %Direction of scripts
storepath = 'C:\Users\gillissen\Desktop\InternshipCédric\FGmainanylsis' %Main storage --> Will contain all processed data & model of the brain
DataDirectory = '\\VC2NIN\wbimaging'%Raw images.
miceopt = {'Marsellus'} %options for mice
Stim2Check = 'FGTask'%Name of the stimulus as written in the LOG-file
fijiLoc = 'C:\Users\gillissen\Desktop\InternshipCédric\Fiji.app\scripts' %FijiLocation
tempstorage = 'C:\Users\gillissen\Desktop\InternshipCédric\tempstorage'%Temporary storage of 
baselinemethod = 1 %1: trial by trial baseline (no detrending necessary per se), 2: filtered baseline (F0), no trial specific baseline (basically F), 3: average baseline (averaged over all trials per condition and then averaged over conditions, after detrending)
takeequalsample = 0;
smoothfact = 2;
RedoAll = 0;
trialtype = 'FG';
timelim = [50 200];
%ACHTUNG: you need to install: https://sites.google.com/site/qingzongtseng/template-matching-ij-plugin#downloads

addpath(genpath(Scriptsdir))
global UserQuestions %Defines whether gui's show up to ask what to do with existing datasets. If 0, existing datasets are used and not overwritten
UserQuestions =0;
  
%% Make info file
% info = BuildInfo(miceopt,DataDirectory,storepath,Stim2Check); %Get all logs of WM-imaging sessions
% % save(fullfile(storepath,'sessionstruct'),'info')
load('C:\Users\gillissen\Desktop\InternshipCédric\FGmainanylsis\sessionstructmarsellus')
logs = info.logs(2:3);
paths = info.paths(2:3);
% logs = info.logs;
% paths = info.paths;
info.logs = info.logs(2:3);
info.paths = info.paths(2:3);

% NoiseCorrelations(info,miceopt,storepath,DataDirectory,Stim2Check,1,trialtype,timelim,2,0,0);
%% EvokedActivity ROI selection
% Use timewindow of 80-120 to get clean figure representation in the brain
% 
% ROIselectionforFIG(info,miceopt,storepath,DataDirectory,Stim2Check,baselinemethod,[-300 1500],{'FG','GREY'},[-300 500],takeequalsample); %Last one: plotlim
% DifferencesCheck(info,miceopt,storepath,DataDirectory,Stim2Check,baselinemethod,[-300 1500],{'FG','GREY'},[-300 500],takeequalsample); %Last one: plotlim

       
%% USER INPUT

trialtypes = {'FG'};
StorePath = storepath;
originallim = [50 200];

nback = 20;
createvideosandfigurespermouse =1;
latencyana = 0;
contrastROIs = 1; %Contrast between conditions to determine ROI?
EvokedActivityROI = 1;
removeErrors = 0 %If you don't want to include Errors, make it 1
removeHits = 0 %if you don't want to include hits, make it 1
removeMiss = 0 %If you don't want to include misses, make it 1
fgtw = [120 250]; %FOR FG     

%% STOP USER INPUT
%Timelimit: Don't need data from time after this.
%Make colormaps
%Make colormaps
posmap = fliplr([linspace(1,1,128);linspace(0,1,128);zeros(1,128)]);
% blackmap = fliplr([linspace(0.2,0.40,12);linspace(0.2,0.40,12);linspace(0.2,0.40,12)]);
negmap = fliplr([zeros(1,128);linspace(1,1,128);fliplr(linspace(0,1,128))]);
PSCOREMAP = fliplr(cat(2,posmap,negmap))';
blackval = round(0.95*size(PSCOREMAP,1)/2);
blackrange = (size(PSCOREMAP,1)/2)-blackval:(size(PSCOREMAP,1)/2)+(blackval-1);
blackmap = [fliplr(linspace(0,0.6,blackval)),linspace(0,0.6,blackval)]; %make 0.6 or sth instead of 1 to have more 'abrupt' black to color
PSCOREMAP(blackrange,:) = PSCOREMAP(blackrange,:).*repmat(blackmap,[3,1])';
for i = 1:3
    PSCOREMAP(:,i) = smooth(PSCOREMAP(:,i),5);
end
x = 1:256;
y = 1:256;
X = meshgrid(x,y);
figure; imagesc(X)
colormap(PSCOREMAP)
ActSupColorMap = fliplr(cat(2,posmap,negmap))';
%Mix in black in the middle
blackval = 60;
blackrange = (size(ActSupColorMap,1)/2)-blackval:(size(ActSupColorMap,1)/2)+(blackval-1);
blackmap = [fliplr(linspace(0,1,blackval)),linspace(0,1,blackval)]; %make 0.6 or sth instead of 1 to have more 'abrupt' black to color
%Smooth
ActSupColorMap(blackrange,:) = ActSupColorMap(blackrange,:).*repmat(blackmap,[3,1])';
for i = 1:3
    ActSupColorMap(:,i) = smooth(ActSupColorMap(:,i),5);
end
x = 1:256;
y = 1:256;
X = meshgrid(x,y);
figure; imagesc(X)
colormap(ActSupColorMap)
% Make line map
%Green for hit, red for erros, black for misses
greenmap = [zeros(1,5);linspace(0.5,1,5);zeros(1,5)];
redmap = [linspace(0.5,1,5);zeros(1,5);zeros(1,5)];
blackmap = [linspace(0,0.5,5);linspace(0,0.5,5);linspace(0,0.5,5)];
yellowmap = [linspace(0.5,1,5);linspace(0.5,1,5);zeros(1,5)];
LineMap = cat(3,fliplr(redmap),fliplr(greenmap),fliplr(blackmap),fliplr(yellowmap));
micechosen = zeros(1,length(miceopt));

%% Gather all data
nrMouse = 1:length(miceopt);
mousecount = 0;
for midx = 1:nrMouse %For this mouse
    if sum(~cellfun(@isempty, {logs{midx,:,:}})) < 1 %If not recorded that day, skip
        continue
    end
    
    SUMSQALL = cell(1,length(trialtypes));
    zeesc = SUMSQALL;
    nrt = zeesc;
    nrtperpix = zeesc;
    trialavg = zeesc;
    
    
%     for c = 1:length(trialtypes)
%         zeesc{c}(1) = zeros(800,800);
%         zeesc{c}(2) = zeros(800,800);
%     end
%     
    
    newroiscount = 1;
    mouse = miceopt{midx};
    mousecount = mousecount+1
    referenceimage = uint8(imread(fullfile(StorePath,mouse,'\RefFile.bmp')));
    sessioncount = 0;
    
    %Load Alan Brain model
    BrainModel{midx} = load(fullfile(StorePath,mouse,'brainareamodel.mat'))
    FigPos(midx) = load(fullfile(StorePath,mouse,'figpos.mat'));
    
      for didx = 1:size(logs,2) %Loop over days
        if sum(~cellfun(@isempty, {logs{midx,didx,:}})) < 1 %If not recorded that day, skip
            continue
        end
        for sidx = 1:size(logs,3) %If no xth session, continue
            if sum(~cellfun(@isempty,{logs{midx,didx,sidx}})<1)
                continue
            end
            micechosen(midx) = 1;
            sessioncount = sessioncount+1
            clear tosave;
            clear LOG
            clear this
            tmppath = paths{midx,didx,sidx};
            date = strsplit(tmppath,mouse);
            date = date{3}(1:end-1) %Find date
            expnr = strsplit(tmppath,mouse);
            expnr = str2num(expnr{end});%find session nr
            
            disp(['Loading data ' mouse ', day ' date ', session ' num2str(expnr)])
            
            %% LOG
            load([tmppath '\' mouse num2str(expnr) '.mat']);
            
            if exist('tosave','var')
                try
                    LOG=tosave.LOG;
                catch
                    LOG = tosave.Log;
                end
            end
            if exist('Log','var')
                LOG = Log;
                clear Log;
            end
            
            if strcmp(Stim2Check,'DelayedOriTuningSound')
                %Make this.log.Orientation longer with nans
                LOG.Orientation(end:length(LOG.Reaction)) = 500; %Only goes till 360
                
                while ~isfield(LOG,'correctReaction') %Check whether reactions were registered okay
                    CheckReactions([this.folder this.expname '\' this.mouse this.expnum '.mat'])
                    tmp = load([folder expname '\' mouse expnum '.mat']);
                    if isfield(tmp,'tosave')
                        tmp = tmp.tosave;
                    end
                end
                LOG.Reaction = LOG.correctReaction; %Change the reactions into checked reactions
            end
            
            OriOpttmp = unique(LOG.Orientation);
            if isfield(LOG,'Side')
                SideOpttmp = unique(LOG.Side);
            else
                SideOpttmp = 1;
            end
            if ~iscell(SideOpttmp)
                SideOpttmp = {num2str(SideOpttmp)};
            end
            if isfield(LOG,'Reactions') || isfield(LOG,'Reaction')
                ReactionOpttmp = {'Miss','Hit','Error','Too Early','TooFast'};
                LOG.Condition = zeros(length(LOG.Reaction), 1);
            end
            
            if isfield(LOG,'figon') %MAKE FIGON SIDE == 'none' (or 3)
                FigureOpt = unique(LOG.figon);
                SideOpttmp = [SideOpttmp 'none'];
                LOG.Side(LOG.figon==0) = {'none'};
            end
            
            count = 0;
            for oidx = 1:length(OriOpttmp)
                for soidx = 1:length(SideOpttmp)
                    if isfield(LOG,'Reactions') | isfield(LOG,'Reaction') %active
                        for ridx = 1:length(ReactionOpttmp)
                            count = count + 1;
                            LOG.Condition(strcmp(LOG.Reaction,ReactionOpttmp{ridx})& LOG.Orientation == OriOpttmp(oidx) & ...
                                strcmp(LOG.Side,SideOpttmp{soidx})) = count;
                            ConditionNamestmp{count} = [ReactionOpttmp{ridx} ' Ori' num2str(OriOpttmp(oidx)) ' Side ' SideOpttmp{soidx}];
                        end
                    else %Passive
                        count = count + 1;
                        LOG.Condition(LOG.Orientation == OriOpttmp(oidx) & ...
                            strcmp(LOG.Side,SideOpttmp{soidx})) = count;
                        ConditionNamestmp{count} = ['Ori' num2str(OriOpttmp(oidx)) ' Side ' SideOpttmp{soidx}];
                    end
                    
                end
            end
            
            LOG.Conditions = unique(LOG.Condition);
            cvec = LOG.Conditions;
            if size(cvec,1) > size(cvec,2)
                cvec = cvec'
            end
            cvec(cvec==0) = [];
            
            ConditionNamestmp(cvec)
            %% Load Noisecorrss data and other datafiles
            iddone = zeros(1,length(trialtypes));
            for id = 1:length(trialtypes)
                if ~exist(fullfile(StorePath,mouse,[mouse date],[mouse num2str(expnr)],['Baseline' num2str(baselinemethod) ,'_' trialtypes{id}, '_eqsample' num2str(takeequalsample)],'NoiseCorr.mat'))
                    disp(fullfile(StorePath,mouse,[mouse date],[mouse num2str(expnr)],['Baseline' num2str(baselinemethod) ,'_' trialtypes{id}, '_eqsample' num2str(takeequalsample)],'NoiseCorr.mat'))
                    disp('Skipping this session, cause no NoiseCorr data detected. First run WMspecificAnalysis with this baselinemethod')
                    continue
                end
                
                disp('Loading data...')
                load(fullfile(StorePath,mouse,[mouse date],[mouse num2str(expnr)],['Baseline' num2str(baselinemethod) ,'_' trialtypes{id}, '_eqsample' num2str(takeequalsample)],'NoiseCorr.mat'))
                if ~exist('NoiseCorr','var')
                    continue
                    disp('Skipping this session, cause no NoiseCorr data detected. There were no trials of this type in this session')                    
                else
                    tmp = ~cellfun(@isempty,NoiseCorr.nrt);
                    
                    if ~any(tmp(:))
                        disp('Skipping this session, cause no NoiseCorr data detected. There were no trials of this type in this session')
                        
                        continue
                    end
                end 
                
                rawdatfiles = dir(fullfile(StorePath,mouse,[mouse date],[mouse num2str(expnr)],[mouse num2str(expnr) '_RawData*']));
                iddone(id) = 1;
                    if ~isempty(zeesc{id})
                        try
                            for rr = 1:size(zeesc{id},1)
                                for cc = 1:size(zeesc{id},2)
                                    zeesc{id}{rr,cc} = cat(3,zeesc{id}{rr,cc},NoiseCorr.zeesc{rr,cc});
                                    trialavg{id}{rr,cc} = cat(3,trialavg{id}{rr,cc},NoiseCorr.Trialavg{rr,cc});
                                    nrtperpix{id} = cellfun(@(X,Y) X+Y, NoiseCorr.nrtPerPix,nrtperpix{id},'UniformOutput',0);
                                end
                            end
                            
                    catch ME
                        disp(ME)
                        keyboard
                    end
                    else
                       zeesc{id} = NoiseCorr.zeesc; 
                       nrtperpix{id} = NoiseCorr.nrtPerPix;
                       trialavg = NoiseCorr.Trialavg;
                    end
             
            
                if ~isempty(nrt{id})
                    nrt{id} = cellfun(@(X,Y) X+Y,nrt{id},NoiseCorr.nrt,'UniformOutput',0);
                else
                    nrt{id} = NoiseCorr.nrt(:,~ismember(NoiseCorr.ReactionOpt,'TooEarly'));
                end
                meanRT{sessioncount,id} = NoiseCorr.meanRT;
                stdRT{sessioncount,id} = NoiseCorr.stdRT;
                SideOpt{sessioncount,id} = NoiseCorr.SideOpt;
                ReactionOpt{sessioncount,id} = NoiseCorr.ReactionOpt;
                clear NoiseCorr
                load(fullfile(StorePath,mouse,[mouse date],[mouse num2str(expnr)],rawdatfiles(strcmp({rawdatfiles(:).name},[mouse num2str(expnr) '_RawData_C' num2str(length(cvec)) '.mat'])).name));
                clear conddata
                timevectmp{mousecount,sessioncount,id} = timeline(timeline>=originallim(1)&timeline<=originallim(2));
                
            end
        end
    end
    
     %% Define conditions
    maxnrReactOpt = max(max(cellfun(@length,ReactionOpt)));
    AllReactionOpt = unique([ReactionOpt{:}],'stable');
    maxnrReactOpt =maxnrReactOpt- sum(strcmp(AllReactionOpt,'TooEarly'));
    AllReactionOpt(strcmp(AllReactionOpt,'TooEarly')) = [];
    if removeErrors
        AllReactionOpt(strcmp(AllReactionOpt,'Error')) = [];
    end
    if removeHits
        AllReactionOpt(strcmp(AllReactionOpt,'Hit')) = [];
    end
    if removeMiss
        AllReactionOpt(strcmp(AllReactionOpt,'Miss')) = [];
    end

    maxnrSideOpt = max(max(cellfun(@length,SideOpt)));
    AllSideOpt = unique([SideOpt{:}],'stable');
    fgopt = find(iddone);
    [xpix ypix ~] = size(zeesc{fgopt(1)}{1,1});
    
    %Check timevecs: if not the same; warning!
    timevec = cat(1,timevectmp{:});
    timevec = unique(timevec,'rows');
    if size(timevec,1) > 1
        error('Multiple versions of timevec...')
    end
    
    if ~exist(fullfile(StorePath,'Figures',['Baseline' num2str(baselinemethod) '_eqsample' num2str(takeequalsample)]))
        mkdir(fullfile(StorePath,'Figures',['Baseline' num2str(baselinemethod) '_eqsample' num2str(takeequalsample)]))
    end
    
end

    
    %% Z scored NC per seed pixel
%     SideOpt = SideOpt(1,1);
%     ReactionOpt = ReactionOpt(1,1);
    zseedpixelcorr = cell(size(1:length(trialtypes)));
    for trialidx = 1:length(trialtypes)
        zseedpixelcorr{trialidx} = cell(length(maxnrSideOpt),length(maxnrReactOpt));
    end
    model = BrainModel{midx};
    seed = FigPos(midx).figmap;
    halfidx = round(model.Model.Lambda(2));
    seedtmp = cell(2,1);
    seedtmp{1} = seed;
    seedtmp{2} = seed;
    seedtmp{2}(1:halfidx,:) = 0; % if stimulus is left, right v1 as seed
    seedtmp{1}(halfidx:end,:) = 0; 
    clipping = 2;
     for id = 1:length(trialtypes)
           if ~isempty(zeesc{id})
               for sideidx = 1:length(SideOpt)
                for ridx = 1:length(ReactionOpt)
                    tmp = zeesc{sideidx,ridx};
                    tmp(tmp<-clipping) = -clipping;
                    tmp(tmp>clipping) =clipping;
                    tmp2 = tmp;
                    seedtempmat = repmat(seedtmp{sideidx},[1,1,size(tmp,3)]);
                    tmp2(~seedtempmat) = nan;
                    seedtempmat = nanmean(reshape(tmp2,[size(tmp2,1)*size(tmp2,2),size(tmp2,3)]),1);
                    %                    seedtempmat = repmat(seedtmp{sideidx},[1,1,size(tmp,3)]); %repmat the seed in 3rd dim
                    %                    tmp2(~seedtempmat) = nan; % nan all nonrelevant pixels
                    %                    select = ~isnan(tmp2); %find index numbers in the 3d mat that belong to seed
                    %                    seedidx = find(seedtmp{sideidx}); % find idx in 2d seed so you know how to reshape
                    %                    tmp2 = tmp2(select);
                    %                    tmp2 = reshape(tmp2,[size(seedidx,1),size(tmp,3)]);
                    
                    %                  seedtempmat = nanmean(tmp2,1); % Don't average noise corrs
                    %                    seedtempmat = squeeze(nanmean(seedtempmat,2));
                    corrvec = corr(seedtempmat',reshape(tmp,[size(tmp,1)*size(tmp,2),size(tmp,3)])');
                    zseedpixelcorr{sideidx,ridx} = reshape(corrvec,[size(tmp,1),size(tmp,2)]);
                    figure; imagesc(zseedpixelcorr{sideidx,ridx},[0.5 1])
                    colormap(ActSupColorMap)
                end
               end
            else zseedpixelcorr{id} = [];
          end
     end
     
%% LEFT vs RIGHT correlationmaps

%% Amplitude pRF
TW = [150 500];
twidx = find(timeline>=TW(1)&timeline<=TW(2));
for midx = 1:length(miceopt)
    %load pRF
    %Check with mapping if available
    load(fullfile(storepath,miceopt{midx},'brainareamodel.mat'))
    
    pRFfile = fullfile(DataDirectory,'pRF Results',miceopt{midx},'4Mapping','pRFmaps');
    pRFbrain = fullfile(DataDirectory,'pRF Results',miceopt{midx},'4Mapping','RefImg');
    
    %Load up teh reference image from the pRF maps
    load(pRFfile)
    load(pRFbrain)
    
    if ~exist('FGAmplPerRegio','var')
        FGAmplPerRegio = nan(length(Model.Regions),length(miceopt));
    end
    
    mask = abs(out.SIGNMAPt);
    mask = cat(2,zeros([size(mask,1),size(Model.Regions{1},2)-size(mask,2)]),mask);
    for regidx = 1:length(Model.Regions)
        
        mask2 = false(size(mask,1),size(mask,2));
        mask2(Model.Regions{regidx}' == 1 & mask == 1) = true;
        
        if sum(mask2(:)) <= 20
            continue
        end
        figure; subplot(2,3,1); imagesc(Model.Regions{regidx}');
        axis square
        subplot(2,3,2)
        imagesc(mask);         axis square
        subplot(2,3,3);
        imagesc(mask2);
        axis square
        
        title(Model.Rnames{regidx})
        mask2 = repmat(mask2,[1,1,37]);
        %         tmp = abs(AllMicedFF{midx}{1}{1}-AllMicedFF{midx}{1}{2}); %FG-Mod
        %         tmp(~mask2) = nan;
        tmp1=AllMicedFF{midx}{1}{1};
       
        %Half the cortex is where lambda y coord is is: 
        halfidx = round(Model.Lambda(2));
        tmp1(~mask2) = nan;
        tmp2 = AllMicedFF{midx}{1}{2};
        tmp2(~mask2)= nan;
        
        tmp = nan(size(tmp1));
        %Right side of the cortex (L-R)
        tmp(1:halfidx,:) = (tmp1(1:halfidx,:) - tmp2(1:halfidx,:));
        %Left Side of the cortex (R-L)
        tmp(halfidx+1:end,:) = (tmp2(halfidx+1:end,:) - tmp1(halfidx+1:end,:));
        
        %Right side traces
        tmpdFFR = squeeze(nanmean(nanmean(tmp(1:halfidx,:,:),1),2));
        tmpdFF1R = squeeze(nanmean(nanmean(tmp1(1:halfidx,:,:),1),2));
        tmpdFF2R = squeeze(nanmean(nanmean(tmp2(1:halfidx,:,:),1),2));
        subplot(2,3,4)
        plot(timeline(timeline>=-300&timeline<=1500),tmpdFF1R,'g'); hold on; plot(timeline(timeline>=-300&timeline<=1500),tmpdFF2R,'r');
        plot(timeline(timeline>=-300&timeline<=1500),tmpdFFR,'b')
        title('Traces Right Side cortex (L-R)')
        xlim([-200 800])
        ylim([-3*10^-3 17*10^-3])
        
        %Left Side traces
        tmpdFFL = squeeze(nanmean(nanmean(tmp(halfidx:end,:,:),1),2));
        tmpdFF1L = squeeze(nanmean(nanmean(tmp2(halfidx:end,:,:),1),2));
        tmpdFF2L = squeeze(nanmean(nanmean(tmp1(halfidx:end,:,:),1),2));
        subplot(2,3,5)
        plot(timeline(timeline>=-300&timeline<=1500),tmpdFF1L,'g'); hold on; plot(timeline(timeline>=-300&timeline<=1500),tmpdFF2L,'r');
        plot(timeline(timeline>=-300&timeline<=1500),tmpdFFL,'b')
        title('Traces Left Side cortex (R-L)')
        xlim([-200 800])
        ylim([-3*10^-3 17*10^-3])
        
        subtmp = subplot(2,3,6);
        posleg = get(subtmp,'position');
        delete(subtmp)
    
        legend({'Figure','Ground','Modulation'},'Position',posleg)
        FGAmplPerRegio(regidx,midx) = nanmean(nanmean([tmpdFFR(twidx),tmpdFFL(twidx)],2),1);
        
        saveas(gcf,fullfile(storefigures,[miceopt{midx} '_' Model.Rnames{regidx} '_.bmp']))
        saveas(gcf,fullfile(storefigures,[miceopt{midx} '_' Model.Rnames{regidx} '_.fig']))
        saveas(gcf,fullfile(storefigures,[miceopt{midx} '_' Model.Rnames{regidx} '_.eps']))      

    end
    
    %Plotting extra Rois:
    figure; imagesc(brain)
    colormap gray
    hold on
    scatter(Model.AllX,Model.AllY,'y.')
    title(['Press k for enough ROIs or any other key to draw more'])
    roicount = 1;
    stop = 0;
    while ~stop
        BW{roicount} = roipoly;
        roidrawn(roicount) = bwboundaries(BW{roicount});

        scatter(roidrawn{roicount}(:,2),roidrawn{roicount}(:,1),'r.')
        text(nanmean(roidrawn{roicount}(:,2)),nanmean(roidrawn{roicount}(:,1)),num2str(roicount))
        disp('Press keyboard!') 
        waitforbuttonpress
         key = get(gcf,'CurrentCharacter');
         if strcmp(key,'k')
             stop = 1;
         end     
         roicount = roicount+1;
    end
    
      for regidx = 1:length(roidrawn)
        
     
        figure; subplot(2,2,1); imagesc(BW{regidx});
        axis square
      
        
        title(['roi ' num2str(regidx)])
      
        tmp1=AllMicedFF{midx}{1}{1};
        mask2 = repmat(BW{regidx},[1,1,size(tmp1,3)]);
        %Half the cortex is where lambda y coord is is: 
        halfidx = round(Model.Lambda(2));
        tmp1(~mask2) = nan;
        tmp2 = AllMicedFF{midx}{1}{2};
        tmp2(~mask2)= nan;
        
        tmp = nan(size(tmp1));
        %Right side of the cortex (L-R)
        tmp(1:halfidx,:) = (tmp1(1:halfidx,:) - tmp2(1:halfidx,:));
        %Left Side of the cortex (R-L)
        tmp(halfidx+1:end,:) = (tmp2(halfidx+1:end,:) - tmp1(halfidx+1:end,:));
        
        %Right side traces
        tmpdFFR = squeeze(nanmean(nanmean(tmp(1:halfidx,:,:),1),2));
        tmpdFF1R = squeeze(nanmean(nanmean(tmp1(1:halfidx,:,:),1),2));
        tmpdFF2R = squeeze(nanmean(nanmean(tmp2(1:halfidx,:,:),1),2));
        if any(~isnan((unique(tmpdFFR))))
        
        subplot(2,2,3)
        plot(timeline(timeline>=-300&timeline<=1500),tmpdFF1R,'g'); hold on; plot(timeline(timeline>=-300&timeline<=1500),tmpdFF2R,'r');
        plot(timeline(timeline>=-300&timeline<=1500),tmpdFFR,'b')
        title('Traces Right Side cortex (L-R)')
        xlim([-200 800])
        ylim([-3*10^-3 17*10^-3])
        end
        
        %Left Side traces
        tmpdFFL = squeeze(nanmean(nanmean(tmp(halfidx:end,:,:),1),2));
        tmpdFF1L = squeeze(nanmean(nanmean(tmp2(halfidx:end,:,:),1),2));
        tmpdFF2L = squeeze(nanmean(nanmean(tmp1(halfidx:end,:,:),1),2));
        if any(~isnan((unique(tmpdFFL))))
            subplot(2,2,4)
            plot(timeline(timeline>=-300&timeline<=1500),tmpdFF1L,'g'); hold on; plot(timeline(timeline>=-300&timeline<=1500),tmpdFF2L,'r');
            plot(timeline(timeline>=-300&timeline<=1500),tmpdFFL,'b')
            title('Traces Left Side cortex (R-L)')
            xlim([-200 800])
            ylim([-3*10^-3 17*10^-3])
        end
        
        subtmp = subplot(2,2,2);
        posleg = get(subtmp,'position');
        delete(subtmp)
    
        legend({'Figure','Ground','Modulation'},'Position',posleg)
        FGAmplPerRegio(regidx,midx) = nanmean(nanmean([tmpdFFR(twidx),tmpdFFL(twidx)],2),1);
        
        saveas(gcf,fullfile(storefigures,[miceopt{midx} '_ExtraROI' num2str(regidx)  '_.bmp']))
        saveas(gcf,fullfile(storefigures,[miceopt{midx} '_ExtraROI' num2str(regidx) '_.fig']))
        saveas(gcf,fullfile(storefigures,[miceopt{midx} '_ExtraROI' num2str(regidx) '_.eps']))      

    end
    
    
end
 

%% Plot Results


      % plot results
      for ridx = 1:length(ReactionOpt)
          for sideidx = 1:length(SideOpt)
              if ~isempty(zseedpixelcorr{sideidx,ridx})
                  figure 
                  h = imagesc(zseedpixelcorr{sideidx,ridx});
                  colormap(ActSupColorMap)
%                   set(h,'AlphaData',brainmask~=0)
                  str = sprintf('%s %s n = %1.0f',ReactionOpt{ridx},SideOpt{sideidx});
                  title(str);
                  colorbar;
                  
              end
          end
      end
      
      
      
      
      
      
      
      