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
                                  
                    try
                        zeesc{id} = cellfun(@(X,Y) X+Y, NoiseCorr.zeesc,zeesc(id),'UniformOutput',0);
                        nrtperpix{id} = cellfun(@(X,Y) X+Y, NoiseCorr.nrtPerPix,nrtperpix{id},'UniformOutput',0);
                    catch ME
                        disp(ME)
                        keyboard
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
    
     seed = FigPos(midx);
      for ridx = 1:length(ReactionOpt)
        for sideidx = 1:length(SideOpt)
           if ~isempty(zeesc{sideidx,ridx})
            tmp = zeesc{sideidx,ridx};
            tmp(tmp>2|tmp<-2) = 2;
            tmp2 = tmp;
            seedtempmat = repmat(seedtemp,[1,1,size(tmp,3)]);
            tmp2(~seedtempmat) = nan;
            seedtempmat = nanmean(tmp2,1);
            seedtempmat = squeeze(nanmean(seedtempmat,2));
            corrvec = corr(seedtempmat,reshape(tmp,[size(tmp,1)*size(tmp,2),size(tmp,3)])');
            zseedpixelcorr{sideidx,ridx} = reshape(corrvec,[size(tmp,1),size(tmp,2)]);
           else; zseedpixelcorr{sideidx,ridx} = [];
           end
        end
      end
      

 

%% Plot Results


      % plot results
      for ridx = 1:length(ReactionOpt)
          for sideidx = 1:length(SideOpt)
              if ~isempty(zseedpixelcorr{sideidx,ridx})
                  figure 
                  h = imagesc(zseedpixelcorr{sideidx,ridx},imrange);
                  colormap(ActSupColorMap)
                  set(h,'AlphaData',brainmask~=0)
                  str = sprintf('%s %s n = %1.0f',ReactionOpt{ridx},SideOpt{sideidx},nrt{sideidx,ridx});
                  title(str);
                  colorbar;
                  
              end
          end
      end
      
      
      
      
      
      
      
      