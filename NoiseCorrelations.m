function [corrmat] = NoiseCorrelations(info,miceopt,StorePath,DataDirectory,Stim2Check,baselinemethod,trialtype,timelim,smoothfact,takeequalsample,Redo,seed)
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
            % create cell with the proper area logicals
            brainmask = zeros(800,800);
            areamask = zeros(800,800);
            for i = 1:length(Model.Regions)
                Borders = Model.Boundaries{i};
                for j = 1:length(Borders)
                    tmp = poly2mask(Borders{j}(:,1),Borders{j}(:,2),800,800);
                    tmp = imfill(tmp,'holes');
                    areamask(tmp) = 1+i;
                    brainmask(tmp) = 1;
                end
            end
            brainmask = imfill(brainmask,'holes');
        
            
             %% Extract raw calcium traces per trialtype per condition
             % Don't apply drift correction, just trial specific baseline
             
            if strcmp(Stim2Check,'DelayedOriTuningSound')
                fullfgtr = find(LOG.currentdelay==trialtype & LOG.Gavepassive(LOG.RealTask==1) == 0 & LOG.Ezbox == 0 & LOG.TotalStimDur == 500);
            else
                if strcmp(trialtype,'FG')
                    fullfgtr = find(LOG.BGContrast==1 & LOG.Gavepassive==0&LOG.Ezbox==0);
                else
                    fullfgtr = find(LOG.BGContrast==0 & LOG.Gavepassive==0&LOG.Ezbox==0);
                end
            end
            
            
            
            if ~exist(fullfile(StorePath,mouse,[mouse date],[mouse num2str(expnr)],['Baseline' num2str(baselinemethod) '_' TRIALTYPE, '_eqsample' num2str(takeequalsample)],'LEFTVSRIGHT.mat'))|| strcmp(button,'Redo')
                
                nrtr = zeros(length(ReactionOpt),length(SideOpt));
                %Calculate minimum nr of trials for all conditoins
                for ridx = 1:length(ReactionOpt)
                    for stidx = 1:length(SideOpt)
                        cidx = find(ismember(reaction,ReactionOpt{ridx})&ismember(side,SideOpt{stidx})&~ismember(orientation,'Ori500'));
                        for ccidx = 1:length(cidx)
                            tmp = ctrials{cidx(ccidx)};
                            try
                            tmp = tmp(~removeidx(1:length(tmp),cidx(ccidx))' & ismember(tmp,fullfgtr));
                            catch ME
                                disp(ME)
                                keyboard
                            end
                            nrtr(ridx,stidx) = nrtr(ridx,stidx) + length(tmp);
                        end
                    end
                end
                nrtrials2take = min(nrtr(:));
                if isempty(fullfgtr)
                    disp('This session does not have trials with requested properties...')
                    save(fullfile(StorePath,mouse,[mouse date],[mouse num2str(expnr)],['Baseline' num2str(baselinemethod) '_' TRIALTYPE, '_eqsample' num2str(takeequalsample)],'LEFTVSRIGHT.mat'),'fullfgtr')
                    continue
                else
                    if isa(trialtype,'char')
                        disp(['Only taking ' trialtype ' trials into account'])
                    else
                        disp(['Only taking ' num2str(trialtype) 'ms delay trials into account'])
                    end
                end
                
                %initialize
                dFFav = cell(length(SideOpt),length(ReactionOpt));
                nrt = dFFav;
                meanRT = dFFav;
                stdRT = dFFav;
                SUMSqr = dFFav;
                nrtPerPix = dFFav;
                SumdFF = dFFav;
                
                
                
                disp('Calculating dFF for different conditions')
                %                 makeplots = 0;
                %                 if makeplots
                %                 %Look at raw F
                %                 FRAW = figure;
                %
                %                 %Look at individual TC vs Pixel TC
                %                 FIV = figure;
                %
                %                 %Check properties of raw signal
                %                 load(fullfile(StorePath,mouse,'Averages','ROIs.mat'))
                %                 RV1mask = poly2mask(rois{1}.xi,rois{1}.yi,800,800);
                %                 %Shrink to not have border effects
                %                 RV1mask = bwmorph(RV1mask,'shrink',1);
                %                 legendname = {};
                %                 end
                
                %Calculate baseline
                if baselinemethod==2
                    tracesthiscond = uint8([]); %find traces
                    rm2idx = logical([]);
                    trialidx = [];
                    for ridx = 1:length(ReactionOpt)
                        for stidx = 1:length(SideOpt)
                            ccidx = find(ismember(reaction,ReactionOpt{ridx}) & ismember(side,SideOpt{stidx}) & ~ismember(orientation,'Ori500')); %Find condition index
                            if isempty(ccidx)
                                continue
                            end
                            trialidxtmp = [];
                            for indx = 1:length(ccidx)
                                trialidxtmp = [trialidxtmp ctrials{ccidx}]; %trialindex
                                rmtmp = removeidx(1:length(ctrials{ccidx(indx)}),ccidx(indx))';
                                rm2tmp = ~ismember(ctrials{ccidx(indx)},fullfgtr);
                                rm3tmp = (rmtmp==1 | rm2tmp==1);
                                try
                                    if takeequalsample
                                        tmpvec = find(rm3tmp == 0);
                                        tmprand = randperm(length(tmpvec));
                                        rm3tmp(tmpvec(tmprand(1:length(tmpvec)-nrtrials2take))) = 1;
                                        savetmprand{ridx,stidx} = tmprand;
                                    end
                                    tracesthiscond = cat(4,tracesthiscond,RawData{ccidx(indx)}(:,:,timeline>=-300 & timeline<0,~rm3tmp));
                                    trialidx = [trialidx trialidxtmp(~rm3tmp)];
                                catch ME
                                    disp(ME)
                                    if strcmp(ME.identifier,'MATLAB:nomem')
                                        keyboard
                                    end
                                end
                            end
                        end
                    end
                    tmp =  single(tracesthiscond);
                    tmp(tmp==0)=nan;
                    base = tmp./permute(repmat(BASELINEMAT(:,:,trialidx),[1,1,1,size(tmp,3)]),[1,2,4,3]); %apply drift correction
                    base = nanmean(nanmean(base,3),4);
                end
                for ridx = 1:length(ReactionOpt)
                    for stidx = 1:length(SideOpt)
                        ccidx = find(ismember(reaction,ReactionOpt{ridx}) & ismember(side,SideOpt{stidx}) & ~ismember(orientation,'Ori500')); %Find condition index
                        
                        if isempty(ccidx)
                            continue
                        end
                        try
                            trialidx = [ctrials{ccidx}];%trialindex
                        catch ME
                            disp(ME)
                        end
                        rm2idx = logical([]);
                        tracesthiscond = uint8([]); %find traces
                        for indx = 1:length(ccidx)
                            rmtmp = removeidx(1:length(ctrials{ccidx(indx)}),ccidx(indx))';
                            rm2tmp = ~ismember(ctrials{ccidx(indx)},fullfgtr);
                            rm3tmp = (rmtmp==1 | rm2tmp==1);
                            if takeequalsample
                                tmpvec = find(rm3tmp == 0);
                                if baselinemethod==2
                                    tmprand=savetmprand{ridx,stidx};
                                else
                                    tmprand = randperm(length(tmpvec));
                                end
                                rm3tmp(tmpvec(tmprand(1:length(tmpvec)-nrtrials2take))) = 1;
                            end
                            try
                                tracesthiscond = cat(4,tracesthiscond,RawData{ccidx(indx)}(:,:,:,~rm3tmp));
                                rm2idx = [rm2idx, rm3tmp];
                            catch ME
                                disp(ME)
%                                 if strcmp(ME.identifier,'MATLAB:nomem')
                                    keyboard
%                                 end
                            end
                            
                        end
                        trialidx = trialidx(~rm2idx);
                        
                        if isempty(trialidx)
                            continue
                        end
                        
                        %
                        %                         if makeplots
                        %                             tmp2 = single(tracesthiscond);
                        %                             tmp2(tmp2==0) = nan;
                        %                             tmp2(~repmat(RV1mask,[1,1,size(tmp2,3),size(tmp2,4)]))=nan;
                        %
                        %                             [sortedid idx] = sort(trialidx);
                        %                             figure(FRAW)
                        %                             %Take all baselines
                        %                             subplot(2,2,1)
                        %                             tmp = squeeze(nanmean(reshape(tmp2,[800*800,size(tmp2,3),size(tmp2,4)]),1));
                        %                             base = squeeze(nanmean(tmp(timeline>=-300 & timeline<=0,:),1));
                        %                             plot(sortedid,nanmean(base(:,idx),1),'-*','Color',LineMap(:,stidx*2,ridx),'MarkerSize',5,'LineWidth',2)
                        %                             hold on
                        %                             title('F0')
                        %
                        %                             subplot(2,2,2)
                        %                             tmp = squeeze(nanmean(tmp(timeline>=80 & timeline<=200,:),1)) - base; %F-F0
                        %                             plot(sortedid,nanmean(tmp(:,idx),1),'-*','Color',LineMap(:,stidx*2,ridx),'MarkerSize',5,'LineWidth',2)
                        %                             title('F - F0')
                        %                             hold on
                        %
                        %                             subplot(2,2,3)
                        %                             tmp = tmp./base;
                        %                             plot(sortedid,nanmean(tmp(:,idx),1),'-*','Color',LineMap(:,stidx*2,ridx),'MarkerSize',5,'LineWidth',2)
                        %                             title('(F-F0)/F0')
                        %                             hold on
                        %
                        %                             figure(FIV)
                        %                             tmp = reshape(tmp2,[800*800,size(tmp2,3),size(tmp2,4)]);
                        %                             [r] = find(~isnan(tmp(:,10,1)));
                        %
                        %                             %Pick 5 neighbouring pixels
                        %                             randchoose = randsample(r,1);
                        %
                        %                             subplot(2,2,1)
                        %                             tmp = tmp(randchoose-2:randchoose+2,timeline>=-300 & timeline<=0,1);
                        %                             plot(repmat(timeline(timeline>=-300 & timeline<=0),[size(tmp,1),1])',tmp')
                        %                             xlabel('Time')
                        %                             ylabel('F')
                        %                             tmp = reshape(tmp2,[800*800,size(tmp2,3),size(tmp2,4)]);
                        %                             randchoose = randsample(r(50:end-50),1);
                        %                             tmp = reshape(tmp(randchoose-50:randchoose+50,timeline>=-300 & timeline<=0,:),[101,length(timeline(timeline>=-300 & timeline<=0))*size(tmp,3)]);
                        %
                        %                             cormat = cov(tmp');
                        %                             subplot(2,2,2); imagesc(cormat,[-15 60])
                        %                             title('Covariance over time between pixels')
                        %                             xlabel('Distance between pixels')
                        %                             colorbar
                        %
                        %                             tmp = reshape(tmp2,[800*800,size(tmp2,3),size(tmp2,4)]);
                        %
                        %                             subplot(2,2,3)
                        %                             tmp = tmp(randchoose-2:randchoose+2,timeline>=80 & timeline<=250,1);
                        %                             plot(repmat(timeline(timeline>=80 & timeline<=250),[size(tmp,1),1])',tmp')
                        %                             xlabel('Time')
                        %                             ylabel('F')
                        %                             tmp = reshape(tmp2,[800*800,size(tmp2,3),size(tmp2,4)]);
                        %                             randchoose = randsample(r(50:end-50),1);
                        %                             tmp = reshape(tmp(randchoose-50:randchoose+50,timeline>=80 & timeline<=250,:),[101,length(timeline(timeline>=80 & timeline<=250))*size(tmp,3)]);
                        %
                        %                             cormat = cov(tmp');
                        %                             subplot(2,2,4); imagesc(cormat,[-15 60])
                        %                             title('Covariance over time between pixels')
                        %                             xlabel('Distance between pixels')
                        %                             colorbar
                        %
                        %                         end
                        try
                            tracesthiscond = single(tracesthiscond); %make single
                            tracesthiscond(tracesthiscond==0) = nan;
                        catch ME
                            if strcmp(ME.identifier,'MATLAB:nomem')
                                for i = 1:size(tracesthiscond,4)
                                    tmp = tracesthiscond(:,:,:,i);
                                    tmp(tmp==0) = nan;
                                    tracesthiscond(:,:,:,i) = tmp;
                                end
                            else
                                disp(ME)
                                keyboard
                            end
                        end
                        
                        
                        for i = 1:size(tracesthiscond,4)
                            % Apply brainmask
                            tmp = tracesthiscond(:,:,:,i);
                            tmp(~repmat(brainmask,[1,1,size(tmp,3)]))= nan;
                            tracesthiscond(:,:,:,i) = imgaussfilt(tmp,smoothfact); %Gaussian filter
                            
                        end
                        
                        
                        for i = 1:100:size(tracesthiscond,1) %Can't do the whole dataset in one go, internal conversion to double?
                            % Extract data
                            tmp = tracesthiscond(i:i+99,:,:,:);
                            
                            %dFF
                            if baselinemethod==1 % dFF (After drift correction) a trial specific baseline
                                %Slow trent correction
                                tmp = tmp ./ permute(repmat(BASELINEMAT(i:i+99,:,trialidx),[1,1,1,size(tmp,3)]),[1,2,4,3]);
                                base = single(squeeze(nanmean(tmp(:,:,timeline>=-300 & timeline<0,:),3))); %Baseline
                                base(base==0) = nan; %Remove 0 and make nan; cannot divide by 0
                                tmp = single(tmp(:,:,timeline>=timelim(1) &timeline<=timelim(2),:)); %F
                                tmp = (tmp - permute(repmat(base,[1,1,1,size(tmp,3)]),[1,2,4,3]))./permute(repmat(base,[1,1,1,size(tmp,3)]),[1,2,4,3]); %dF/F
                            elseif baselinemethod==2 % Only drift correction, then use average baseline
                                %Slow trent correction
                                tmp = tmp ./ permute(repmat(BASELINEMAT(i:i+99,:,trialidx),[1,1,1,size(tmp,3)]),[1,2,4,3]);
                                tmp = single(tmp(:,:,timeline>=timelim(1) &timeline<=timelim(2),:)); %F
%                                 tmp = (tmp - repmat(base(i:i+99,:),[1,1,size(tmp,3),size(tmp,4)]))./repmat(base(i:i+99,:),[1,1,size(tmp,3),size(tmp,4)]); %dF/F
                            elseif baselinemethod == 3 %no detrending, only average baseline
                                base = single(squeeze(nanmean(BASELINEMAT(i:i+99,:,:),3))); %Baseline
                                base(base==0) = nan; %Remove 0 and make nan; cannot divide by 0
                                tmp = single(tmp(:,:,timeline>=timelim(1) &timeline<=timelim(2),:)); %F
                                tmp = (tmp - repmat(base,[1,1,size(tmp,3),size(tmp,4)]))./repmat(base,[1,1,size(tmp,3),size(tmp,4)]); %dF/F 
                            end
                            
                            
                            [row,column] = find(seed);
                            trialbytrialavg = cell(size(row,column));
                               for i = 1:size(tmp,4);
                                   for n = 1:length(row)
                                       for m = length(column);
                                           trialbytrialavg(row,column) = nanmean(tmp(row,column,:,i),3);
                                       end
                                   end
                               end
                               
                              
                                           
                                   
                                
                            
                            
                             
                            
                            
                            dFFav{stidx,ridx}(i:i+99,:,:) = nanmean(tmp,4); %average over trials
                            SUMSqr{stidx,ridx}(i:i+99,:,:) = nansum(tmp.^2,4);      %Sum of squared dFF values (for z-score calculations)
                            SumdFF{stidx,ridx}(i:i+99,:,:) = nansum(tmp,4); %Sum of trials
                            nrtPerPix{stidx,ridx}(i:i+99,:,:) = size(tracesthiscond,4) - sum(isnan(tmp),4); %
                            
                            
                            
                            
                            
                            
                            
                           
                        end
                        nrt{stidx,ridx} = size(tracesthiscond,4);
                        
                        % RT
                        meanRT{stidx,ridx} = nanmean([LOG.RT(trialidx)]);
                        stdRT{stidx,ridx} = nanstd([LOG.RT(trialidx)]);
                        
                    end
                end
             
               
            end
        end
    end
end

end
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
     
            
            %% 
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            