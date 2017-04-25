function FunctionalConnectivity(info,miceopt,StorePath,DataDirectory,Stim2Check,baselinemethod,originallim,trialtypes,plotlim)
global UserQuestions
paths = info.paths;
logs = info.logs;
nrMouse = size(paths,1);
%Timelimit: Don't need data from time after this.
%Make colormaps
%Make colormaps
posmap = fliplr([linspace(1,1,128);linspace(0,1,128);zeros(1,128)]);
% blackmap = fliplr([linspace(0.2,0.40,12);linspace(0.2,0.40,12);linspace(0.2,0.40,12)]);
negmap = fliplr([zeros(1,128);linspace(1,1,128);fliplr(linspace(0,1,128))]);
ActSupColorMap = fliplr(cat(2,posmap,negmap))';

%Mix in black in the middle
blackval = 40;
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
    mouse = miceopt{midx};
    mousecount = mousecount+1
    referenceimage = uint8(imread(fullfile(StorePath,mouse,'\RefFile.bmp')));
    sessioncount = 0;
    
    %Load Alan Brain model
    BrainModel{mousecount} = load(fullfile(StorePath,mouse,'brainareamodel.mat'))
       
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
            %% Load LEFTvsRight data and other datafiles
            iddone = zeros(1,length(trialtypes));
            for id = 1:length(trialtypes)
                if ~exist(fullfile(StorePath,mouse,[mouse date],[mouse num2str(expnr)],['Baseline' num2str(baselinemethod) ,'_' trialtypes{id}],'LEFTVSRIGHT.mat'))
                    disp(fullfile(StorePath,mouse,[mouse date],[mouse num2str(expnr)],['Baseline' num2str(baselinemethod) ,'_' trialtypes{id}],'LEFTVSRIGHT.mat'))
                    disp('Skipping this session, cause no LEFTvsRIGHT data detected. First run WMspecificAnalysis with this baselinemethod')
                    continue
                end
                
                iddone(id) = 1;
                disp('Loading data...')
                load(fullfile(StorePath,mouse,[mouse date],[mouse num2str(expnr)],['Baseline' num2str(baselinemethod) ,'_' trialtypes{id}],'LEFTVSRIGHT.mat'))
                rawdatfiles = dir(fullfile(StorePath,mouse,[mouse date],[mouse num2str(expnr)],[mouse num2str(expnr) '_RawData*']));
                if isfield(LEFTVSRIGHT,'SUMSqr')
                    SUMSQALL{mousecount,sessioncount,id} = LEFTVSRIGHT.SUMSqr;
                else
                    warning(['No SUM of squared dFF variable, redo FGSpecificAnalysis for this session!'])
                    disp('Skipping session')
                    continue
                end
                dFFav{mousecount,sessioncount,id} = LEFTVSRIGHT.dFFav;
                nrt{mousecount,sessioncount,id} = LEFTVSRIGHT.nrt;
                meanRT{mousecount,sessioncount,id} = LEFTVSRIGHT.meanRT;
                stdRT{mousecount,sessioncount,id} = LEFTVSRIGHT.stdRT;
                SideOpt{mousecount,sessioncount,id} = LEFTVSRIGHT.SideOpt;
                ReactionOpt{mousecount,sessioncount,id} = LEFTVSRIGHT.ReactionOpt;
               
                    
                load(fullfile(StorePath,mouse,[mouse date],[mouse num2str(expnr)],rawdatfiles(strcmp({rawdatfiles(:).name},[mouse num2str(expnr) '_RawData_C' num2str(length(cvec)) '.mat'])).name));
                clear conddata
                timevec{mousecount,sessioncount,id} = timeline(timeline>=originallim(1)&timeline<=originallim(2));
                try
                    
                    if isfield(LEFTVSRIGHT,'ConditionNames')
                        ConditionNames{mousecount,sessioncount} =  LEFTVSRIGHT.ConditionNames;
                        conditionparts{mousecount,sessioncount} = cellfun(@(X) strsplit(X,' '),ConditionNames{mousecount,sessioncount},'UniformOutput',0);
                        reaction{mousecount,sessioncount} = cellfun(@(X) X{1},conditionparts{mousecount,sessioncount},'UniformOutput',0); %Reaction
                        orientation{mousecount,sessioncount} = cellfun(@(X) X{2},conditionparts{mousecount,sessioncount},'UniformOutput',0); %orientations
                        side{mousecount,sessioncount} = cellfun(@(X) X{4},conditionparts{mousecount,sessioncount},'UniformOutput',0); %SIdes
                    else
                        ConditionNames{mousecount,sessioncount} = ConditionNamestmp(cvec);
                        conditionparts{mousecount,sessioncount} = cellfun(@(X) strsplit(X,' '),ConditionNames{mousecount,sessioncount},'UniformOutput',0);
                        reaction{mousecount,sessioncount} = cellfun(@(X) X{1},conditionparts{mousecount,sessioncount},'UniformOutput',0); %Reaction
                        orientation{mousecount,sessioncount} = cellfun(@(X) X{2},conditionparts{mousecount,sessioncount},'UniformOutput',0); %orientations
                        side{mousecount,sessioncount} = cellfun(@(X) X{4},conditionparts{mousecount,sessioncount},'UniformOutput',0); %SIdes
                    end
                catch ME
                    disp(ME)
                    keyboard
                end
            end
        end
    end        
end

clear LEFTVSRIGHT
            
            
        
            
            

            