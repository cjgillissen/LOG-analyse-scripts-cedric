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
trialtype = 1500;
timelim = [0 500];
%ACHTUNG: you need to install: https://sites.google.com/site/qingzongtseng/template-matching-ij-plugin#downloads

addpath(genpath(Scriptsdir))
global UserQuestions %Defines whether gui's show up to ask what to do with existing datasets. If 0, existing datasets are used and not overwritten
UserQuestions =0;
  
%% Make info file
% info = BuildInfo(miceopt,DataDirectory,storepath,Stim2Check); %Get all logs of WM-imaging sessions
% save(fullfile(storepath,'sessionstruct'),'info')
load('C:\Users\gillissen\Desktop\InternshipCédric\FGmainanylsis\sessionstruct')
info.logs = info.logs(4,7,1); % session Frey20161121
info.paths = info.paths(4,7,1);

%% EvokedActivity ROI selection
% Use timewindow of 80-120 to get clean figure representation in the brain
% Make masks for each side of the brain, mask ~=V1 

load('C:\Users\gillissen\Desktop\InternshipCédric\FGmainanylsis\sessionstruct')
info.paths = info.paths(2,:,:);
info.logs = info.logs(2,:,:);

ROIselectionforFIG(info,miceopt,storepath,DataDirectory,Stim2Check,baselinemethod,[-300 1500],{'FG','GREY'},[-300 500],takeequalsample); %Last one: plotlim
DifferencesCheck(info,miceopt,storepath,DataDirectory,Stim2Check,baselinemethod,[-300 1500],{'FG','GREY'},[-300 500],takeequalsample); %Last one: plotlim


%% Noise Correlations


    
    
NoiseCorrelations(info,miceopt,storepath,DataDirectory,Stim2Check,baselinemethod,'FG',timelim,smoothfact,takeequalsample,RedoAll)                   
 %% Load the NoiseCorr structs... 
 % Generalize code
 % If only interested in visual processing than concatenate the 1500 and 0
 % trialtypes
 

% load('C:\Users\gillissen\Desktop\InternshipCédric\MainAnaStorage\Frey\Frey20161121\Frey1\Frey1_RawData_C1','timeline');
% timelim = [-300 2500];
% timevec = timeline(timeline>=timelim(1)&timeline<=timelim(2));
%in real version here the task epochs are also defined. 0-500 visual. 500-1450 delay. 1500- onwards response period. 
 
load('C:\Users\gillissen\Desktop\InternshipCédric\MainAnaStorage\Frey\brainareamodel');
brainmask = zeros(800,800);

   %% Colormap
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
      

%% Brainmask
            % create cell with the proper area logicals
            areas = Model.Regions;
            brainmask = zeros(800,800);
            for i = 1:length(Model.Regions)
                Borders = Model.Boundaries{i};
                for j = 1:length(Borders)
                    tmp = poly2mask(Borders{j}(:,1),Borders{j}(:,2),800,800);
                    tmp = imfill(tmp,'holes');
                    brainmask(tmp) = 0+i;
                end
            end
            brainmask = imfill(brainmask,'holes');
            
            %% Load variables
            
                dFFav = NoiseCorr.dFFav;
                nrt = NoiseCorr.nrt;
                meanRT = NoiseCorr.meanRT;
                stdRT = NoiseCorr.stdRT;
                SideOpt = NoiseCorr.SideOpt;
                ReactionOpt = NoiseCorr.ReactionOpt;
                Visualavg = NoiseCorr.Trialavg;
                zeesc = NoiseCorr.zeesc;
                
                ConditionNames =  NoiseCorr.ConditionNames;
                
                conditionparts = cellfun(@(X) strsplit(X,' '),ConditionNames,'UniformOutput',0);
                reaction = cellfun(@(X) X{1},conditionparts,'UniformOutput',0); %Reaction
                orientation = cellfun(@(X) X{2},conditionparts,'UniformOutput',0); %orientations
                side = cellfun(@(X) X{4},conditionparts,'UniformOutput',0); %SIdes
                
                clear NoiseCorr
  
      ConAreaAv = cell(length(SideOpt),length(ReactionOpt));
      corrFCM1 = cell(length(SideOpt),length(ReactionOpt));
      cohereFCM1 = cell(length(SideOpt),length(ReactionOpt));
      Seedpixelcorr = cell(length(SideOpt),length(ReactionOpt));
      zpixelcorr = cell(length(SideOpt),length(ReactionOpt));

%% Noise Correlation per pixel
 
load('C:\Users\gillissen\Desktop\InternshipCédric\MainAnaStorage\Frey\Frey20161121\Frey1\Baseline4_1500_eqsample0\NoiseCorr');


     fig = imagesc(brainmask);
     seed = roipoly;
     seedtemp = seed;
     
      for ridx = 1:length(ReactionOpt)
        for sideidx = 1:length(SideOpt)
            tmp = Visualavg{sideidx,ridx};
            tmp2 = tmp;
            seedtemp = repmat(seedtemp,[1,1,size(tmp,3)]);
            tmp2(~seedtemp) = nan;
            seedtemp = nanmean(tmp2,1);
            seedtemp = squeeze(nanmean(seedtemp,2));
            
            FC = nan(size(tmp)); % initiliaze vector with correlation measures     
                for pixelx = 1:size(tmp,1)
                    for pixely = 1:size(tmp,2) 
                     
                    Seedpixelcorr{sideidx,ridx}(pixely,pixelx) = corr(seedtemp,squeeze(tmp(pixelx,pixely,:)));
                    end
                end
        end
      end
      
   
      
      %% Z scored NC per seed pixel
     imrange =  [-0.8,0.95];
     fig = imagesc(brainmask);
     seed = roipoly;
     seedtemp = seed;
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
            zpixelcorr{sideidx,ridx} = reshape(corrvec,[size(tmp,1),size(tmp,2)]);
           else; zpixelcorr{sideidx,ridx} = [];
           end
        end
      end
      
      % plot results
      for ridx = 1:length(ReactionOpt)
          for sideidx = 1:length(SideOpt)
              if ~isempty(zpixelcorr{sideidx,ridx})
                  figure 
                  h = imagesc(zpixelcorr{sideidx,ridx},imrange);
                  colormap(ActSupColorMap)
                  set(h,'AlphaData',brainmask~=0)
                  str = sprintf('%s %s n = %1.0f',ReactionOpt{ridx},SideOpt{sideidx},nrt{sideidx,ridx});
                  title(str);
                  colorbar;
                  
              end
          end
      end
      
      
      
      
      %% Pairwise Noise correlations between pixels of different areas
      % For every pair of pixels correlate noise correlations per condition
      % Per area otherwise too many pairs??
     
      CorrCond = cell(length(SideOpt),length(ReactionOpt));
      
      
      selectseed = listdlg('PromptString','Select one seed area:',...
                'SelectionMode','single',...
                'ListString',Model.Rnames);
      selecttarget = listdlg('PromptString','Select target areas:',...
                'SelectionMode','multiple',...
                'ListString',Model.Rnames);
        
         
        corrcell = cell(1:length(selecttarget)); % intiliaze noise correlation cell where all the correlation matrices go in
        targetim = cell(1:length(selecttarget));
        seedim = brainmask==selectseed;
        for i = 1:length(targetim)
            targetim{i} = brainmask==selecttarget(i);
        end
         for ridx = 1:length(ReactionOpt)
          for sideidx = 1:length(SideOpt)
            if ~isempty(zeesc{sideidx,ridx})
              seedtmp = zeesc{sideidx,ridx};
              seedtmp(~repmat(seedim,[1,1,size(seedtmp,3)])) = nan;
              seedtmp = reshape(seedtmp,[size(seedtmp,1)*size(seedtmp,2),size(seedtmp,3)]);
              seedtmp=seedtmp(~isnan(seedtmp(:,2)),:);
                  for targetidx = 1:length(targetim)
                      targettmp = zeesc{sideidx,ridx};
                      targettmp(~repmat(targetim{targetidx},[1,1,size(targettmp,3)])) = nan;
                      targettmp = reshape(targettmp,[size(targettmp,1)*size(targettmp,2),size(targettmp,3)]);
                      targettmp=targettmp(~isnan(targettmp(:,2)),:);
                      corrcell{targetidx} = corr(targettmp',seedtmp');
                  end
                  CorrCond{sideidx,ridx} = corrcell;
                 else 
                  CorrCond{sideidx,ridx} = [];
            end
          end
         end
                          
            %% USER INPUT
nback = 20;
createvideosandfigurespermouse =1;
latencyana = 0;
contrastROIs = 1; %Contrast between conditions to determine ROI?
EvokedActivityROI = 1;
removeE rrors = 0 %If you don't want to include Errors, make it 1
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
mousecount = 0;
for midx = 1:nrMouse %For this mouse
    if sum(~cellfun(@isempty, {logs{midx,:,:}})) < 1 %If not recorded that day, skip
        continue
    end
    
    SUMSQALL = cell(1,length(trialtypes));
    SumdFF = SUMSQALL;
    nrt = SumdFF;
    nrtperpix = SumdFF;
    
    newroiscount = 1;
    mouse = miceopt{midx};
    mousecount = mousecount+1
    referenceimage = uint8(imread(fullfile(StorePath,mouse,'\RefFile.bmp')));
    sessioncount = 0;
    
    %Load Alan Brain model
    BrainModel{midx} = load(fullfile(StorePath,mouse,'brainareamodel.mat'))
    
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
                
                if isfield(NoiseCorr,'SUMSqr')
                    %replace empty matrices with zeros
                    [r,c] = find(cell2mat(cellfun(@isempty,NoiseCorr.SUMSqr,'UniformOutput',0)));
                    [rtempl,ctempl] = find(~cell2mat(cellfun(@isempty,NoiseCorr.SUMSqr,'UniformOutput',0)),1);
                    if ~isempty(r)
                        for idx = 1:length(r)
                            try
                                NoiseCorr.SUMSqr{r(idx),c(idx)} = zeros(size(NoiseCorr.SUMSqr{rtempl,ctempl}));
                            catch ME
                                disp(ME)
                                keyboard
                            end
                            
                        end
                    end
                    %Make 0 of nans, cause summation with 0's is going
                    %wrong.
                    for rr = 1:size(NoiseCorr.SUMSqr,1)
                        for cc = 1:size(NoiseCorr.SUMSqr,2)
                            NoiseCorr.SUMSqr{rr,cc}(isnan(NoiseCorr.SUMSqr{rr,cc})) = 0;
                        end
                    end
                    if ~isempty(SUMSQALL{id})
                        if size(NoiseCorr.SUMSqr,2) < size(SUMSQALL{id},2)
                            c = find(~ismember(unique([ReactionOpt{:}]),NoiseCorr.ReactionOpt));
                            for idx = 1:length(c)
                                for rowidx = 1:size(SUMSQALL{id},1)
                                    NoiseCorr.SUMSqr{rowidx,c(idx)} = zeros(size(SUMSQALL{id}{rowidx,c(idx)}));
                                end
                            end
                        end
                        try
                            SUMSQALL{id} = cellfun(@(X,Y) X+Y,SUMSQALL{id},NoiseCorr.SUMSqr(:,~ismember(NoiseCorr.ReactionOpt,'TooEarly')),'UniformOutput',0);
                        catch ME
                            disp(ME)
                            keyboard
                        end
                        
                    else
                        SUMSQALL{id} = NoiseCorr.SUMSqr(:,~ismember(NoiseCorr.ReactionOpt,'TooEarly'));
                    end
                else
                    warning(['No SUM of squared dFF variable, redo FGSpecificAnalysis for this session!'])
                    disp('Skipping session')
                    continue
                end
                % SUMdFF, SUMSQdFF & NRT need to be saved out. throw away
                % NoiseCorr after that
                if ~isfield(NoiseCorr,'SumdFF')
                    error(['Redo ConditionAveraginScript for session ' fullfile(StorePath,mouse,[mouse date],[mouse num2str(expnr)])])
                end
                
                %replace empty matrices with zeros
                [r,c] = find(cell2mat(cellfun(@isempty,NoiseCorr.SumdFF,'UniformOutput',0)));
                [rtempl,ctempl] = find(~cell2mat(cellfun(@isempty,NoiseCorr.SumdFF,'UniformOutput',0)),1);
                if ~isempty(r)
                    for idx = 1:length(r)
                        NoiseCorr.SumdFF{r(idx),c(idx)} = zeros(size(NoiseCorr.SumdFF{rtempl,ctempl}));
                        NoiseCorr.nrtPerPix{r(idx),c(idx)} = zeros(size(NoiseCorr.nrtPerPix{rtempl,ctempl}));
                    end
                end
                
                %replace empty with 0
                [r,c] = find(cellfun(@isempty,NoiseCorr.nrt));
                if ~isempty(r)
                    for idx = 1:length(r)
                        NoiseCorr.nrt{r(idx),c(idx)}=0;
                    end
                end
                %Make 0 of nans, cause summation with nans is going
                %wrong.
                for rr = 1:size(NoiseCorr.dFFav,1)
                    for cc = 1:size(NoiseCorr.dFFav,2)
                        NoiseCorr.SumdFF{rr,cc}(isnan(NoiseCorr.SumdFF{rr,cc})) = 0;
                    end
                end
                
                if ~isempty(SumdFF{id})
                    if size(NoiseCorr.SumdFF,2) < size(SumdFF{id},2)
                        c = find(~ismember(unique([ReactionOpt{:}]),NoiseCorr.ReactionOpt));
                        for idx = 1:length(c)
                            for rowidx = 1:size(SumdFF{id},1)
                                NoiseCorr.SumdFF{rowidx,c(idx)} = zeros(size(SumdFF{id}{rowidx,c(idx)}));
                                NoiseCorr.nrt{rowidx,c(idx)} = 0;
                                NoiseCorr.nrtPerPix{rowidx,c(idx)} = zeros(size(SumdFF{id}{rowidx,c(idx)}));
                            end
                        end
                    end
                    
                    
                    try
                        SumdFF{id} = cellfun(@(X,Y) X+Y, NoiseCorr.SumdFF(:,~ismember(NoiseCorr.ReactionOpt,'TooEarly')),SumdFF{id},'UniformOutput',0);
                        nrtperpix{id} = cellfun(@(X,Y) X+Y, NoiseCorr.nrtPerPix(:,~ismember(NoiseCorr.ReactionOpt,'TooEarly')),nrtperpix{id},'UniformOutput',0);
                    catch ME
                        disp(ME)
                        keyboard
                    end
                else
                    SumdFF{id} = NoiseCorr.SumdFF(:,~ismember(NoiseCorr.ReactionOpt,'TooEarly'));
                    nrtperpix{id} = NoiseCorr.nrtPerPix(:,~ismember(NoiseCorr.ReactionOpt,'TooEarly'));
                end
                
                
                
                if ~isempty(nrt{id})
                    nrt{id} = cellfun(@(X,Y) X+Y,nrt{id},NoiseCorr.nrt(:,~ismember(NoiseCorr.ReactionOpt,'TooEarly')),'UniformOutput',0);
                else
                    nrt{id} = NoiseCorr.nrt(:,~ismember(NoiseCorr.ReactionOpt,'TooEarly'));
                end
                meanRT{sessioncount,id} = NoiseCorr.meanRT;
                stdRT{sessioncount,id} = NoiseCorr.stdRT;
                SideOpt{sessioncount,id} = NoiseCorr.SideOpt;
                ReactionOpt{sessioncount,id} = NoiseCorr.ReactionOpt;
                zeesc{sessioncount,id} = NoiseCorr.zeesc;
                clear NoiseCorr
                load(fullfile(StorePath,mouse,[mouse date],[mouse num2str(expnr)],rawdatfiles(strcmp({rawdatfiles(:).name},[mouse num2str(expnr) '_RawData_C' num2str(length(cvec)) '.mat'])).name));
                clear conddata
                timevectmp{mousecount,sessioncount,id} = timeline(timeline>=originallim(1)&timeline<=originallim(2));
                
            end
        end
    end
    
    %% Differences
      
      
      
      
      
      
      
      
      
      
      