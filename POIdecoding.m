%% First steps in decoding the brain
%All preprocessing steps have been completed until condition averaging. 


%% Initiliaze paths

Scriptsdir = 'C:\Users\gillissen\Documents\Github\Mouse\ImagingAnalysis' %Direction of scripts
storepath = '\\vcnin\mouse_working_memory\Cédric\MainStorage' %Main storage --> Will contain all processed data & model of the brain
DataDirectory = '\\vcnin\mouse_working_memory\Imaging\'%Raw images.
miceopt = {'Chief'} %options for mice
Stim2Check = 'DelayedOriTuningSound'%Name of the stimulus as written in the LOG-file
fijiLoc = 'C:\Users\gillissen\Desktop\InternshipCédric\Fiji.app' %FijiLocation
tempstorage = 'C:\Users\gillissen\Desktop\InternshipCédric\TEMPstorage'%Temporary storage of 
baselinemethod = 3 %1: trial by trial baseline (no detrending necessary), 2: filtered baseline (F0), no trial specific baseline (basically F), 3: average baseline (averaged over all trials per condition and then averaged over conditions, after detrending)
smoothfact = 2;
RedoAll = 0;
%ACHTUNG: you need to install: https://sites.google.com/site/qingzongtseng/template-matching-ij-plugin#downloads
addpath('C:\Users\gillissen\Desktop\InternshipCédric');
addpath(genpath('C:\Users\gillissen\Desktop\InternshipCédric\Scripts\Mouse'));
mouse = 'Chief';
date = 20170103;
expnr =1;
%% Load LOG file and create condition names



            load('\\vcnin\mouse_working_memory\Cédric\MainStorage\Chief\Chief20170103\Chief_20170103_B1');

            
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
                    CheckReactions('\\vcnin\mouse_working_memory\Cédric\MainStorage\Chief\Chief20170103\Chief_20170103_B1')
                    tmp = load('\\vcnin\mouse_working_memory\Cédric\MainStorage\Chief\Chief20170103\Chief_20170103_B1');
                    if isfield(tmp,'tosave')
                        tmp = tmp.tosave;
                    end
                end
                LOG.Reaction = LOG.correctReaction; %Change the reactions into checked reactions
            end
            
            OriOpt = unique(LOG.Orientation);
            if isfield(LOG,'Side')
                SideOpt = unique(LOG.Side);
            else
                SideOpt = 1;
            end
            if ~iscell(SideOpt)
                SideOpt = {num2str(SideOpt)};
            end
            if isfield(LOG,'Reactions') || isfield(LOG,'Reaction')
                ReactionOpt = {'Miss','Hit','Error','Too Early','TooFast'};
                LOG.Condition = zeros(length(LOG.Reaction), 1);
            end
            
            count = 0;
            for oidx = 1:length(OriOpt)
                for soidx = 1:length(SideOpt)
                    if isfield(LOG,'Reactions') | isfield(LOG,'Reaction') %active
                        for ridx = 1:length(ReactionOpt)
                            count = count + 1;
                            LOG.Condition(strcmp(LOG.Reaction,ReactionOpt{ridx})& LOG.Orientation == OriOpt(oidx) & ...
                                strcmp(LOG.Side,SideOpt{soidx})) = count;
                            ConditionNames{count} = [ReactionOpt{ridx} ' Ori' num2str(OriOpt(oidx)) ' Side ' SideOpt{soidx}];
                        end
                    else %Passive
                        count = count + 1;
                        LOG.Condition(LOG.Orientation == OriOpt(oidx) & ...
                            strcmp(LOG.Side,SideOpt{soidx})) = count;
                        ConditionNames{count} = ['Ori' num2str(OriOpt(oidx)) ' Side ' SideOpt{soidx}];
                    end
                    
                end
            end
            
            LOG.Conditions = unique(LOG.Condition);
            cvec = LOG.Conditions;
            if size(cvec,1) > size(cvec,2)
                cvec = cvec'
            end
            cvec(cvec==0) = [];
            
            ConditionNames = ConditionNames(cvec);


%% Load raw condition data

cd('\\vcnin\mouse_working_memory\Cédric\MainStorage\Chief\Chief20170103\Chief1')
load(fullfile(storepath,mouse,'brainareamodel.mat'))
load(fullfile(storepath,mouse,[mouse date],[mouse num2str(expnr)],'BASELINEMAT.mat'))

% select some random pixels  inside V1 of left hemisphere
% might want to include t-score based POI selection
x = find(brainmask);
x = x(randperm(length(x),10));
[x,y] = ind2sub([800 800],x(1:10));
POI = [x y];

for j = 1:length(ConditionNames)

load(fullfile(['Chief1_RawData_C' num2str(j)]))


%% Apply brainmask

V1 = single(conddata);
clear conddata
V1(V1==0)=nan; 
brainmask = zeros(800,800);
brainmask(Model.Regions{11}) = 1;

for i = 1:size(V1,4) %loop over trials 
    
    tmp = V1(:,:,:,i);
    tmp(~repmat(brainmask,[1,1,size(tmp,3)]))=nan;
    V1(:,:,:,i) = imgaussfilt(tmp,smoothfact); %Gaussian filter
                
end

clear tmp

%% Slow trent correction

trialidx = single(ctrials{j});

for i = 1:100:size(V1,1)
     QQ = single(V1(i:i+99,:,:,:));
     QQ = QQ ./ permute(repmat(BASELINEMAT(i:i+99,:,trialidx),[1,1,1,size(QQ,3)]),[1,2,4,3]);
     V1(i:i+99,:,:,:) = QQ;
end

clear QQ

V1 = V1(:,:,:,find(removeidx(j,:))); %Delete movement trials.



 %% Decodign attempts of just a couple of pixels
    % located in V1 of left hemisphere.
    % Select some random pixels of interest. 
    tmp = nan(length(POI),length(POI),size(V1,3),size(V1,4));
        for i = 1:length(POI) % only use POI datapoints
            tmp(i,i,:,:) = QQ(POI(i,1),POI(i,2),:,:);             
        end
        
    
                
end

