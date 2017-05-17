%% Noise correlation script
% call Noisecorr function
% load the Noisecorr files the fucntion made. 
Scriptsdir = 'C:\Users\gillissen\Documents\GitHub\Mouse' %Direction of scripts
storepath = 'C:\Users\gillissen\Desktop\InternshipCédric\MainAnaStorage' %Main storage --> Will contain all processed data & model of the brain
DataDirectory = 'C:\Users\gillissen\Desktop\InternshipCédric\Datadirectproxy'%Raw images.
miceopt = {'Frey'} %options for mice
Stim2Check = 'DelayedOriTuningSound'%Name of the stimulus as written in the LOG-file
fijiLoc = 'C:\Users\gillissen\Desktop\InternshipCédric\Fiji.app\scripts' %FijiLocation
tempstorage = 'C:\Users\gillissen\Desktop\InternshipCédric\tempstorage'%Temporary storage of 
baselinemethod = 4 %1: trial by trial baseline (no detrending necessary per se), 2: filtered baseline (F0), no trial specific baseline (basically F), 3: average baseline (averaged over all trials per condition and then averaged over conditions, after detrending)
takeequalsample = 0;
smoothfact = 2;
RedoAll = 0;
%ACHTUNG: you need to install: https://sites.google.com/site/qingzongtseng/template-matching-ij-plugin#downloads

addpath(genpath(Scriptsdir))
global UserQuestions %Defines whether gui's show up to ask what to do with existing datasets. If 0, existing datasets are used and not overwritten
UserQuestions =0;

%% Make info file
% info = BuildInfo(miceopt,DataDirectory,storepath,Stim2Check); %Get all logs of WM-imaging sessions
% save(fullfile(storepath,'sessionstruct'),'info')
load('C:\Users\gillissen\Desktop\InternshipCédric\MainAnaStorage\sessionstruct')
info.logs = info.logs(4,7,1); % session Frey20161121
info.paths = info.paths(4,7,1);

%% Noise Correlations

trialtype = 1500;
timelim = [0 600];
NoiseCorrelations(info,miceopt,storepath,DataDirectory,Stim2Check,baselinemethod,trialtype,timelim,smoothfact,takeequalsample,RedoAll)

 %% Load the NoiseCorr structs... 
 % Generalize code
 % If only interested in visual processing than concatenate the 1500 and 0
 % trialtypes
 

% load('C:\Users\gillissen\Desktop\InternshipCédric\MainAnaStorage\Frey\Frey20161121\Frey1\Frey1_RawData_C1','timeline');
% timelim = [-300 2500];
% timevec = timeline(timeline>=timelim(1)&timeline<=timelim(2));
%in real version here the task epochs are also defined. 0-500 visual. 500-1450 delay. 1500- onwards response period. 
 

brainmask = zeros(800,800);

%% Brainmask
            % create cell with the proper area logicals
            areas = Model.Regions;
            brainmask = zeros(800,800);
            for i = 1:length(Model.Regions)
                Borders = Model.Boundaries{i};
                for j = 1:length(Borders)
                    tmp = poly2mask(Borders{j}(:,1),Borders{j}(:,2),800,800);
                    tmp = imfill(tmp,'holes');
                    brainmask(tmp) = 1+i;
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
                Visualavg = NoiseCorr.TrialavgVisual;
                
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

%% Noise Correlation per pixel
     
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