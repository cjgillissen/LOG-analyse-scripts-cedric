function NoiseCorrelations(info,miceopt,StorePath,DataDirectory,Stim2Check,baselinemethod,trialtype,timelim,smoothfact,takeequalsample,Redo)
global UserQuestions
%Make colormaps
redColorMap = [zeros(1,132),linspace(0,1,124)];
blueColorMap =  [linspace(1,0,124),zeros(1,132)];
greenColorMap = [zeros(1,256)];
ActSupColorMap = [redColorMap;greenColorMap;blueColorMap];

% Make line map
%Green for hit, red for erros, black for misses
greenmap = [zeros(1,5);linspace(0.5,1,5);zeros(1,5)];
redmap = [linspace(0.5,1,5);zeros(1,5);zeros(1,5)];
blackmap = [linspace(0,0.5,5);linspace(0,0.5,5);linspace(0,0.5,5)];
LineMap = cat(3,fliplr(redmap),fliplr(greenmap),fliplr(blackmap));

if isa(trialtype,'char')
    TRIALTYPE = trialtype;
else
    TRIALTYPE = num2str(trialtype);
end

paths = info.paths;
logs = info.logs;
nrMouse = size(paths,1);
for midx = 1:nrMouse %For this mouse
    if sum(~cellfun(@isempty, {logs{midx,:,:}})) < 1 %If not recorded that day, skip
        continue
    end
    mouse = miceopt{midx};
    referenceimage = uint8(imread(fullfile(StorePath,mouse,'\RefFile.bmp')));
    for didx = 1:size(logs,2) %Loop over days
        if sum(~cellfun(@isempty, {logs{midx,didx,:}})) < 1 %If not recorded that day, skip
            continue
        end
        for sidx = 1:size(logs,3) %If no xth session, continue
            if sum(~cellfun(@isempty,{logs{midx,didx,sidx}})<1)
                continue
            end
            clear ctrials
            clear tosave;
            clear LOG
            clear this
            tmppath = paths{midx,didx,sidx};
            date = strsplit(tmppath,mouse);
            date = date{3}(1:end-1) %Find date
            expnr = strsplit(tmppath,mouse);
            expnr = str2num(expnr{end});%find session nr
            
            disp(['Loading data ' mouse ', day ' date ', session ' num2str(expnr)])
            
            %Load 'drift correction'
            try
                load(fullfile(StorePath,mouse,[mouse date],[mouse num2str(expnr)],'BASELINEMAT.mat'))
            catch ME
                disp(ME)
                keyboard
            end
            
            %% Log file
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
            
            if isfield(LOG,'figon') %MAKE FIGON SIDE == 'none' (or 3)
                FigureOpt = unique(LOG.figon);
                SideOpt = [SideOpt 'none'];
                LOG.Side(LOG.figon==0) = {'none'};
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
            
            ConditionNames(cvec)
            rawdatfiles = dir(fullfile(StorePath,mouse,[mouse date],[mouse num2str(expnr)],[mouse num2str(expnr) '_RawData*']));
            
            %% Extract Alan Brain Model
            load(fullfile(StorePath,mouse,'brainareamodel.mat'))
            
            %% Make direction
            if ~exist(fullfile(StorePath,mouse,[mouse date],[mouse num2str(expnr)],['Baseline' num2str(baselinemethod) '_' TRIALTYPE, '_eqsample' num2str(takeequalsample)]))
                mkdir(fullfile(StorePath,mouse,[mouse date],[mouse num2str(expnr)],['Baseline' num2str(baselinemethod) '_' TRIALTYPE, '_eqsample' num2str(takeequalsample)]))
            end
            
            %% Load different conditions
            if  ~Redo && exist(fullfile(StorePath,mouse,[mouse date],[mouse num2str(expnr)],['Baseline' num2str(baselinemethod) '_' TRIALTYPE, '_eqsample' num2str(takeequalsample)],'LEFTVSRIGHT.mat'))
                if UserQuestions
                    button = questdlg(['LEFTVSRIGHT.mat detected, want to redo?' fullfile(StorePath, mouse, [mouse date], [mouse num2str(expnr)])],'Data Redo', 'Redo','Keep current','Keep current');
                else
                    button = 'Keep current';
                end
            else
                button = 'Redo';
            end
            if ~exist(fullfile(StorePath,mouse,[mouse date],[mouse num2str(expnr)],['Baseline' num2str(baselinemethod) '_' TRIALTYPE, '_eqsample' num2str(takeequalsample)],'LEFTVSRIGHT.mat')) || strcmp(button,'Redo')
                clear RawData
                countermin = 0;
                for cidx = 1:length(cvec)
                    disp(['Loading data condition ' num2str(cidx) ' of ' num2str(length(cvec))])
                    clear conddata
                    load(fullfile(StorePath,mouse,[mouse date],[mouse num2str(expnr)],rawdatfiles(strcmp({rawdatfiles(:).name},[mouse num2str(expnr) '_RawData_C' num2str(cidx) '.mat'])).name));
                    try
                        RawData{cvec(cidx)-countermin} = conddata(:,:,timeline<=timelim(2),:); %+
                    catch
                        keyboard
                    end
                    %Remove conditions that did not happen
                    removevec = cellfun(@isempty,RawData);
                    if any(removevec==1)
                        ConditionNames(removevec) = []
                        RawData(removevec) = [];
                        countermin = countermin+sum(removevec);
                    end
                end
                ConditionNames(sum(removevec==0)+1:end) = [];
                timeline(timeline>timelim(2))=[];
                
                idx = find(~cellfun(@isempty,strfind(ConditionNames,'Too Early')));
                ConditionNames(idx) = cellfun(@(X) strrep(X,X(strfind(X,'Too Early'):9),'TooEarly'),ConditionNames(idx),'UniformOutput',0);
                
                %Average over orientations
                conditionparts = cellfun(@(X) strsplit(X,' '),ConditionNames,'UniformOutput',0);
                
                %Find all reactions
                reaction = cellfun(@(X) X{1},conditionparts,'UniformOutput',0); %Reaction
                ReactionOpt = unique(reaction);
                orientation = cellfun(@(X) X{2},conditionparts,'UniformOutput',0); %orientations
                
                side = cellfun(@(X) X{4},conditionparts,'UniformOutput',0); %SIdes
                SideOpt = unique(side,'stable');
                if isempty(timeline)
                    timeline = LOG.Start*1000:LOG.Exposure:round(LOG.Duration*1000)-abs(LOG.Start*1000);
                end
                
                timevec = timeline(timeline>=timelim(1)&timeline<=timelim(2));
                
                %% Get rid of motion trials
                if UserQuestions
                    mbutton = questdlg('Which motion recording method was used?','Motiontype','FG_Ulf','FG_Chris','WM_Chris','FG_Chris');
                elseif ~UserQuestions && strcmp(Stim2Check,'FGTask')
                    mbutton = 'FG_Chris';
                else
                    mbutton = 'WM_Chris';
                end
                
                if strcmp(mbutton,'FG_Ulf')
                    movepath = fullfile(DataDirectory,'Motionlogs (motion ONLY_Ulf''s method)',mouse);
                    newdateform = date;
                    tmpdir = dir(fullfile(movepath,'*.mat'));
                elseif strcmp(mbutton,'FG_Chris')
                    movepath = fullfile(DataDirectory,'Eye + Motion Data (CvdT method)',mouse);
                    newdateform = [date(1:4) '_' date(5:6) '_' date(7:8)];
                    tmpdir = dir(fullfile(movepath,'*.dat'));
                elseif strcmp(mbutton,'WM_Chris')
                    movepath = fullfile(DataDirectory,'Eye_and_Motion',mouse);
                    movepath2 = strsplit(movepath,'Imaging\');
                    movepath = [movepath2{:}];
                    newdateform = [date(1:4) '_' date(5:6) '_' date(7:8)];
                    tmpdir = dir(fullfile(movepath,'*.dat'));
                end
                
                ids = strfind({tmpdir(:).name},newdateform);
                
                %If multiple, search the one with the matching session nr.
                if isfield(LOG,'subsess')
                    if ischar(LOG.subsess)
                        tmpses = num2cell(LOG.subsess);
                        LOG.subsess = cellfun(@str2num,tmpses(~cellfun(@(X) ismember(X,''),tmpses)),'UniformOutput',1);
                    end
                    sess = unique(LOG.subsess);
                    if sum(~cell2mat(cellfun(@isempty,ids,'UniformOutput',0)))>1
                        ids2 = strfind({tmpdir(:).name},['B' num2str(sess(1))]);
                        ids(cellfun(@isempty,ids2)) = {[]};
                    end
                else
                    sess = [];
                    if sum(~cell2mat(cellfun(@isempty,ids,'UniformOutput',0)))>1
                        ids2 = strfind({tmpdir(:).name},['B' num2str(expnr)]);
                        ids(cellfun(@isempty,ids2)) = {[]};
                    end
                end
                try
                    [MoveMat,PupilMat]= extractmotion(fullfile(movepath,tmpdir(~cellfun(@isempty, ids)).name),mbutton,LOG,timeline,ctrials);
                catch ME
                    disp(ME)
                    keyboard
                end
                DifMotion = abs(cat(1,zeros(1,size(MoveMat,2),size(MoveMat,3)),diff(MoveMat,[],1)));
                clear zpupil
                for tmpi = 1:size(PupilMat,4)
                    tmp = PupilMat(:,:,:,tmpi);
                    zpupil(:,:,:,tmpi) = (tmp - repmat(nanmean(tmp(:)),[size(tmp,1),size(tmp,2),size(tmp,3)]))./repmat(nanstd(tmp(:)),[size(tmp,1),size(tmp,2),size(tmp,3)]);
                end
            else
                %                 names = dir(fullfile(StorePath,mouse,[mouse date],[mouse num2str(expnr)],[mouse num2str(expnr) '_RawData_*.mat']));
                %                 tmp2 = cellfun(@(X) strsplit(X,'.mat'),{names(:).name},'UniformOutput',0);
                %                 tmp2 = cellfun(@(X) strsplit(X{1},'C'),tmp2,'UniformOutput',0);
                %                 [tmp2 idx]= max(cell2mat(cellfun(@(X) str2num(X{2}),tmp2,'UniformOutput',0)));
                %
                %                 load(fullfile(StorePath,mouse,[mouse date],[mouse num2str(expnr)],names(idx).name))
            end
            
            %% Motion split
            if ~Redo && exist(fullfile(StorePath,mouse,[mouse date],[mouse num2str(expnr)],'ThrowAwayIdx.mat'))
                if UserQuestions
                    button2 = questdlg(['Motion Thresholds detected. Want to do it again?' fullfile(StorePath, mouse, [mouse date], [mouse num2str(expnr)])],'Motion data', 'Redo','Keep current','Keep current');
                else
                    button2 = 'Keep current';
                end
            else
                button2 = 'Redo';
            end
            
            if strcmp(button2,'Redo')
                if isfield(LOG,'permissiontime')
                    timeconst = 500; %in ms
                    threshold = 30;
                elseif strcmp(Stim2Check,'DelayedOriTuningSound')
                    timeconst = min(unique(LOG.Stimdur))+trialtype%in ms
                    threshold = 200;
                else
                    timeconst = 300; %in ms
                    threshold = 50;
                end
                done = 0;
                %Choose timeconst
                disp('Choose timepoint from where movement doesn''t matter, and threshold for throwing movement out (in terms of value)')
                disp('Press ''k'' for okay')
                HH = figure;
                tp = 1;
                while ~done
                    
                    if tp == 0 && (sld.Value ~= timeconst || sld2.Value ~= threshold)
                        delete(sh)
                    end
                    if tp == 0
                        timeconst = sld.Value
                        threshold = sld2.Value
                    end
                    
                    %Plots of movement over time
                    sh(1) = subplot(2,2,1);
                    plot(timeline,reshape(DifMotion,[size(DifMotion,1),size(DifMotion,2)*size(DifMotion,3)]))
                    set(gca,'xlim',[min(timeline) max(timeline)]);
                    hbar = line([timeconst,timeconst],get(gca,'ylim'),'color',[0 0 0],'LineWidth',2);
                    hold on
                    xtmplim = get(gca,'xlim');
                    title('Movement all trials')
                    
                    if tp ==1
                        sld = uicontrol('Style','slider','Min',min(xtmplim),'Max',max(xtmplim),'Value',timeconst,'Position', [200 500 120 20]);
                    end
                    
                    %Histogram of maximum movement before timepoint
                    tmp = nanmax(reshape(DifMotion(timeline<timeconst,:,:),[sum(timeline<timeconst),size(DifMotion,2)*size(DifMotion,3)]),[],1);
                    tmp(tmp==0) = []; %No trials
                    sh(2) = subplot(2,2,2);
                    histogram(tmp)
                    title(['Movement before ' num2str(timeconst), 'threshold ' num2str(threshold)])
                    hbar2 = line([threshold,threshold],get(gca,'ylim'),'color',[0 0 0],'LineWidth',2);
                    xtmplim = get(gca,'xlim');
                    if tp == 1
                        sld2 = uicontrol('Style','slider','Min',min(xtmplim),'Max',max(xtmplim),'Value',threshold,'Position', [800 500 120 20]);
                    end
                    
                    
                    %Histogram of maximum movement after timepoint
                    tmp2 = nanmax(reshape(DifMotion(timeline>timeconst,:,:),[sum(timeline>timeconst),size(DifMotion,2)*size(DifMotion,3)]),[],1);
                    tmp2(tmp2==0) = []; %No trials
                    sh(3) = subplot(2,2,4);
                    histogram(tmp2)
                    title(['Movement after ' num2str(timeconst), 'threshold ' num2str(threshold)])
                    xlabel('Difference in motion')
                    
                    %Calculate which trials should be removed
                    removeidx =squeeze((nanmax(DifMotion(timeline<timeconst,:,:),[],1)>threshold));
                    
                    sh(4) = subplot(2,2,3);
                    tmp3 = reshape(DifMotion,[size(DifMotion,1),size(DifMotion,2)*size(DifMotion,3)]);
                    plot(timeline,tmp3(:,~reshape(removeidx,[1,size(removeidx,1)*size(removeidx,2)])))
                    hbar = line([timeconst,timeconst],get(gca,'ylim'),'color',[0 0 0],'LineWidth',2);
                    set(gca,'xlim',[min(timeline) max(timeline)]);
                    title(['Movement of selected trials, n = ' num2str(sum(cellfun(@length,ctrials))-sum(removeidx(:)))])
                    xlabel('Time (ms)')
                    
                    if get(HH,'CurrentCharacter') == 'k'
                        done = 1;
                    end
                    
                    if tp == 1
                        tp = 0;
                    end
                    
                    if UserQuestions
                        waitforbuttonpress %Wait for input
                    else
                        done = 1;
                    end
                end
                
                save(fullfile(StorePath,mouse,[mouse date],[mouse num2str(expnr)],'ThrowAwayIdx'),'removeidx','threshold','timeconst')
            else
                load(fullfile(StorePath,mouse,[mouse date],[mouse num2str(expnr)],'ThrowAwayIdx.mat'))
            end
            
            %% Brainmask
            brainmask = zeros(800,800);
            for i = 1:length(Model.Regions)
                Borders = Model.Boundaries{i};
                for j = 1:length(Borders)
                    tmp = poly2mask(Borders{j}(:,1),Borders{j}(:,2),800,800);
                    tmp = imfill(tmp,'holes');
                    brainmask(tmp) = 1;
                end
            end
            brainmask = imfill(brainmask,'holes');