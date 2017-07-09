function TrialbyTrialNC(info,miceopt,StorePath,Stim2Check,baselinemethod,trialtypes,takeequalsample)
global UserQuestions

%% Define timebins


nback = 20;
createvideosandfigurespermouse =0;
latencyana = 0;
contrastROIs = 0; %Contrast between conditions to determine ROI?
EvokedActivityROI = 1;
wholebrainana = 1;
newsize = [400 400];
scalefct = 0.5;
cpimrange = [-0.5 0.8];

if strcmp(Stim2Check,'FGTask')
    basel = [-250 -50];
    fgtw = [120 250]; %FOR FG
    VisInit = [50 120];
    bigtw = [200 450];
    TW = {VisInit,fgtw,bigtw}
    % TW = {bigtw}
else
    basel = [-300 -50];
    visint = [50 200];
    vistw = [200 500]; %FOR FG
    Delayearlytw = [600 1350];
    Delaylatetw = [1400 1600];
    Responstw = [1900 2200];
    TW = {basel,visint,vistw,Delayearlytw,Responstw,Delaylatetw}
    %     TW = {vistw};
    
end
micechosen = zeros(1,length(miceopt));


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
            
            %% Extract Alan Brain Model and SEED
            load(fullfile(StorePath,mouse,'brainareamodel.mat'))
            
            evrois = load(fullfile('C:\Users\gillissen\Desktop\EvokedROIMouse',[mouse 'EvokedActivROIs']));
            
            %rightv1 stimulus representation
            rightv1 = zeros(newsize);
            rightv1(sub2ind(size(rightv1),round(rois{1}.xi.*scalefct),round(rois{1}.yi.*scalefct))) =1;
            rightv1 = imfill(rightv1');
            rightv1XY = regionprops(rightv1,'centroid');
            [xgrid, ygrid] = meshgrid(1:newsize(2), 1:newsize(1));
            centerRV1 = ((xgrid-rightv1XY.Centroid(1)).^2 + (ygrid-rightv1XY.Centroid(2)).^2) <= 5.^2;
            
            %leftv1 stimulus representation
            
            leftv1 = zeros(newsize);
            leftv1(sub2ind(size(leftv1),round(rois{2}.xi.*scalefct),round(rois{2}.yi.*scalefct))) =1;
            leftv1 = imfill(leftv1');
            leftv1XY = regionprops(leftv1,'centroid');
            [xgrid, ygrid] = meshgrid(1:newsize(2), 1:newsize(1));
            centerLV1 = ((xgrid-leftv1XY.Centroid(1)).^2 + (ygrid-leftv1XY.Centroid(2)).^2) <= 5.^2;
            
            
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
                
              %% Eye and Motion
              
               %removeidx
                        load(fullfile(StorePath,mouse,[mouse date],[mouse num2str(expnr)],'ThrowAwayIdx.mat'))
                        if strcmp(Stim2Check,'DelayedOriTuningSound')
                            fullfgtr = find(LOG.currentdelay==str2num(trialtypes{id}) & LOG.Gavepassive(LOG.RealTask==1) == 0 & LOG.Ezbox == 0 & LOG.TotalStimDur == 500);
                        else
                            if strcmp(trialtypes{id},'FG')
                                fullfgtr = find(LOG.BGContrast==1 & LOG.Gavepassive==0&LOG.Ezbox==0);
                            else
                                fullfgtr = find(LOG.BGContrast==0 & LOG.Gavepassive==0&LOG.Ezbox==0);
                            end
                        end
                        
                        %Load baselinecorrection
                        load(fullfile(StorePath,mouse,[mouse date],[mouse num2str(expnr)],'BASELINEMAT.mat'))
                        BASELINEMAT = imresize(BASELINEMAT,newsize,'bilinear');
                        %for timeline etc.
                        load(fullfile(StorePath,mouse,[mouse date],[mouse num2str(expnr)],rawdatfiles(strcmp({rawdatfiles(:).name},[mouse num2str(expnr) '_RawData_C' num2str(length(cvec)) '.mat'])).name));
                        clear conddata
                        twidx = find(timeline>=TW{twid}(1)&timeline<=TW{twid}(2));
                        
                        if TW
                        baseidx = find(timeline>=basel(1) & timeline<=basel(2));
                        
                        
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
            
            %% Left vs Right response
            if strcmp(Stim2Check,'DelayedOriTuningSound')
                fullfgtr = find(LOG.currentdelay==trialtype & LOG.Gavepassive(LOG.RealTask==1) == 0 & LOG.Ezbox == 0 & LOG.TotalStimDur == 500);
            else
                if strcmp(trialtype,'FG')
                    fullfgtr = find(LOG.BGContrast==1 & LOG.Gavepassive==0&LOG.Ezbox==0);
                else
                    fullfgtr = find(LOG.BGContrast==0 & LOG.Gavepassive==0&LOG.Ezbox==0);
                end
            end
         
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
                    end
                end
            
            end
        end
    end
    
    SAVE
    
end
end

 

            










