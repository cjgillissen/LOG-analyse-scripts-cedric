%% First steps in decoding the brain
%All preprocessing steps have been completed until condition averaging. 


%% Initiliaze paths

Scriptsdir = 'C:\Users\gillissen\Documents\Github\Mouse\ImagingAnalysis' %Direction of scripts
storepath = '\\vcnin\mouse_working_memory\C�dric\MainStorage' %Main storage --> Will contain all processed data & model of the brain
DataDirectory = '\\vcnin\mouse_working_memory\Imaging\'%Raw images.
miceopt = {'Chief'} %options for mice
Stim2Check = 'DelayedOriTuningSound'%Name of the stimulus as written in the LOG-file
fijiLoc = 'C:\Users\gillissen\Desktop\InternshipC�dric\Fiji.app' %FijiLocation
tempstorage = 'C:\Users\gillissen\Desktop\InternshipC�dric\TEMPstorage'%Temporary storage of 
baselinemethod = 3 %1: trial by trial baseline (no detrending necessary), 2: filtered baseline (F0), no trial specific baseline (basically F), 3: average baseline (averaged over all trials per condition and then averaged over conditions, after detrending)
smoothfact = 2;
RedoAll = 0;
%ACHTUNG: you need to install: https://sites.google.com/site/qingzongtseng/template-matching-ij-plugin#downloads
addpath('C:\Users\gillissen\Desktop\InternshipC�dric');
addpath(genpath('C:\Users\gillissen\Desktop\InternshipC�dric\Scripts\Mouse'));
mouse = 'Chief';
date = 20170103;
expnr =1;
%% Load LOG file and create condition names



            load('\\vcnin\mouse_working_memory\C�dric\MainStorage\Chief\Chief20170103\Chief_20170103_B1');

            
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
                    CheckReactions('\\vcnin\mouse_working_memory\C�dric\MainStorage\Chief\Chief20170103\Chief_20170103_B1')
                    tmp = load('\\vcnin\mouse_working_memory\C�dric\MainStorage\Chief\Chief20170103\Chief_20170103_B1');
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

cd('\\vcnin\mouse_working_memory\C�dric\MainStorage\Chief\Chief20170103\Chief1')
load(fullfile(storepath,mouse,'brainareamodel.mat'))
load(fullfile(storepath,mouse,[mouse date],[mouse num2str(expnr)],'BASELINEMAT.mat'))

% for j = 1:length(ConditionNames)

load(fullfile(['Chief1_RawData_C' num2str(j)]))


%% ROI and downsample data

V1 = single(conddata);
clear conddata
V1(V1==0)=nan; 
% brainmask = zeros(800,800);
% brainmask(Model.Regions{11}) = 1;
[x,y] = find(Model.Regions{11,1});
maxy = max(y);
V1 = V1(round(Model.Bregma(2)):maxy,100:400,3:end,:);
BASELINEMAT = BASELINEMAT(round(Model.Bregma(2)):maxy,100:400,:);

clear tmp

%% Slow trent correction

% trialidx = single(ctrials{j});
trialidx = ctrials{2};
tmp = V1(:,:,:,find(~removeidx(2,:))); %Delete movement trials.



% tmp  = tmp ./permute(repmat(BASELINEMAT(:,:,trialidx(find(removeidx(2,:)))),[1,1,1,size(tmp,3)]),[1,2,4,3]);
% 
% bleach correction to harsh for optical flow... values around 1 

base = single(squeeze(nanmean(tmp(:,:,timeline>=-350 & timeline<0,:),3))); %Baseline
base(base==0) = nan; %Remove 0 and make nan; cannot divide by 0
tmp = single(tmp); %F
tmp = (tmp - permute(repmat(base,[1,1,1,size(tmp,3)]),[1,2,4,3]))./permute(repmat(base,[1,1,1,size(tmp,3)]),[1,2,4,3]); %dF/F



mask = ones(size(tmp));
mask(isnan(tmp)) = 0;
x = mask(:,:,5);
save('logcialmask','x');



V1 = double(squeeze(nanmean(tmp,4)));
save('chiefC2processed.mat','tmp');


