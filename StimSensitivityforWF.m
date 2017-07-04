function StimSensetivityforWF(info,miceopt,StorePath,Stim2Check,baselinemethod,trialtypes,SideOpt,takeequalsample)

global UserQuestions
paths = info.paths;
logs = info.logs;
nrMouse = size(paths,1);

%% USER INPUT
nback = 20;
createvideosandfigurespermouse =0;
latencyana = 0;
contrastROIs = 0; %Contrast between conditions to determine ROI?
EvokedActivityROI = 1;
wholebrainana = 1;
spotlightana = 1;
nrSVM = 12;
if strcmp(Stim2Check,'FGTask')
    basel = [-250 -50];
    fgtw = [120 250]; %FOR FG
    VisInit = [50 120];
    bigtw = [200 450];
    TW = {VisInit,fgtw,bigtw}
    % TW = {bigtw}
else
    basel = [-250 -50];
    vistw = [50 500]; %FOR FG
    Delaytw = [800 1500];
    Responstw = [1900 2200];
    TW = {vistw,Delaytw,Responstw}
    %     TW = {vistw};
    
end
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
    referenceimage = uint8(imread(fullfile(StorePath,mouse,'\RefFilepRF.bmp')));
    xpix = size(referenceimage,1);
    ypix = size(referenceimage,2);
    %Load Alan Brain model
    BrainModel{midx} = load(fullfile(StorePath,mouse,'brainareamodel.mat'))
    
        Hwhole = figure('units','normalized','outerposition',[0 0 1 1],'name',[mouse '_Allsess_' 'stimselectivity' ' wholebrain']);
        Hspot = figure('units','normalized','outerposition',[0 0 1 1],'name',[mouse '_Allsess_'  'stimselectivity' ' Spotlight']);
        for twid = 1:length(TW)
            for id = 1:length(trialtypes)
                sessioncount = 0;
                rightdat = [];
                leftdat = [];
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
                        
                        ConditionNamestmp = ConditionNamestmp(cvec)
                        rawdatfiles = dir(fullfile(StorePath,mouse,[mouse date],[mouse num2str(expnr)],[mouse num2str(expnr) '_RawData*']));
                        
                        %Average over orientations
                        conditionparts = cellfun(@(X) strsplit(X,' '),ConditionNamestmp,'UniformOutput',0);
                        
                        %Find all reactions
                        reaction = cellfun(@(X) X{1},conditionparts,'UniformOutput',0); %Reaction
                        
                        orientation = cellfun(@(X) X{2},conditionparts,'UniformOutput',0); %orientations
                        side = cellfun(@(X) X{4},conditionparts,'UniformOutput',0); %SIdes
                        SideOpt = unique(side,'stable');
                        
                        
                        if ~exist(fullfile(StorePath,mouse,[mouse date],[mouse num2str(expnr)],['Baseline' num2str(baselinemethod) ,'_' trialtypes{id}, '_eqsample' num2str(takeequalsample)],'LEFTVSRIGHT.mat'))
                            disp(fullfile(StorePath,mouse,[mouse date],[mouse num2str(expnr)],['Baseline' num2str(baselinemethod) ,'_' trialtypes{id}, '_eqsample' num2str(takeequalsample)],'LEFTVSRIGHT.mat'))
                            disp('Skipping this session, cause no LEFTvsRIGHT data detected. First run WMspecificAnalysis with this baselinemethod')
                            continue
                        end
                        
                        disp('Loading data...')
                        load(fullfile(StorePath,mouse,[mouse date],[mouse num2str(expnr)],['Baseline' num2str(baselinemethod) ,'_' trialtypes{id}, '_eqsample' num2str(takeequalsample)],'LEFTVSRIGHT.mat'))
                        if ~exist('LEFTVSRIGHT','var')
                            continue
                            disp('Skipping this session, cause no LEFTvsRIGHT data detected. There were no trials of this type in this session')
                        else
                            tmp = ~cellfun(@isempty,LEFTVSRIGHT.nrt);
                            
                            if ~any(tmp(:))
                                disp('Skipping this session, cause no LEFTvsRIGHT data detected. There were no trials of this type in this session')
                                continue
                            end
                        end
                        
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
                        
                        %for timeline etc.
                        load(fullfile(StorePath,mouse,[mouse date],[mouse num2str(expnr)],rawdatfiles(strcmp({rawdatfiles(:).name},[mouse num2str(expnr) '_RawData_C' num2str(length(cvec)) '.mat'])).name));
                        clear conddata
                        twidx = find(timeline>=TW{twid}(1)&timeline<=TW{twid}(2));
                        baseidx = find(timeline>=basel(1) & timeline<=basel(2));
                        
                        % load in rawdata
                        rightidx = find(~cellfun(@isempty,(cellfun(@(X) strfind(X,'right'),ConditionNamestmp,'UniformOutput',0))));
                        if strcmp(Stim2Check,'DelayedOriTuningSound')
                            rightidx(ismember(rightidx,find(~cellfun(@isempty,(cellfun(@(X) strfind(X,'500'),ConditionNamestmp,'UniformOutput',0)))))) = [];
                        end
                        for i = 1:length(rightidx)
                            tmpload = load(fullfile(StorePath,mouse,[mouse date],[mouse num2str(expnr)],rawdatfiles(strcmp({rawdatfiles(:).name},[mouse num2str(expnr) '_RawData_C' num2str(rightidx(i)) '.mat'])).name));
                            try
                                rmtmp = removeidx(1:length(tmpload.ctrials{rightidx(i)}),rightidx(i))';
                                rm2tmp = ~ismember(tmpload.ctrials{rightidx(i)},fullfgtr);
                                rm3tmp = (rmtmp==1 | rm2tmp==1);
                                trialidx = tmpload.ctrials{rightidx(i)};
                                trialidx = trialidx(rm3tmp);
                                
                                rightdattmp =  zeros(800,800,length(twidx),sum(rm3tmp),'single');
                                for j = 1:100:xpix
                                    tmpnw = single(tmpload.conddata(j:j+99,:,:,rm3tmp))./permute(repmat(BASELINEMAT(j:j+99,:,trialidx),[1,1,1,size(tmpload.conddata,3)]),[1,2,4,3]);
                                    rightdattmp(j:j+99,:,:,:) = (tmpnw(:,:,twidx,:)-repmat(nanmean(tmpnw(:,:,baseidx,:),3),[1,1,length(twidx),1]))./repmat(nanmean(tmpnw(:,:,baseidx,:),3),[1,1,length(twidx),1]);
                                end
                                rightdattmp = squeeze(nanmean(rightdattmp,3));
                                rightdat = cat(3,rightdat,rightdattmp);
                            catch ME
                                disp(ME)
                                keyboard
                            end
                            clear tmpload
                        end
                        
                        leftidx = find(~cellfun(@isempty,(cellfun(@(X) strfind(X,'left'),ConditionNamestmp,'UniformOutput',0))));
                        if strcmp(Stim2Check,'DelayedOriTuningSound')
                            leftidx(ismember(leftidx,find(~cellfun(@isempty,(cellfun(@(X) strfind(X,'500'),ConditionNamestmp,'UniformOutput',0)))))) = [];
                        end
                        for i = 1:length(leftidx)
                            tmpload = load(fullfile(StorePath,mouse,[mouse date],[mouse num2str(expnr)],rawdatfiles(strcmp({rawdatfiles(:).name},[mouse num2str(expnr) '_RawData_C' num2str(leftidx(i)) '.mat'])).name));
                            try
                                
                                rmtmp = removeidx(1:length(tmpload.ctrials{leftidx(i)}),leftidx(i))';
                                rm2tmp = ~ismember(tmpload.ctrials{leftidx(i)},fullfgtr);
                                rm3tmp = (rmtmp==1 | rm2tmp==1);
                                trialidx = tmpload.ctrials{leftidx(i)};
                                trialidx = trialidx(rm3tmp);
                                leftdattmp = nan(800,800,length(twidx),sum(rm3tmp),'single');
                                for j = 1:100:xpix
                                    tmpnw =  single(tmpload.conddata(j:j+99,:,:,rm3tmp))./permute(repmat(BASELINEMAT(j:j+99,:,trialidx),[1,1,1,size(tmpload.conddata,3)]),[1,2,4,3]);
                                    leftdattmp(j:j+99,:,:,:) = (tmpnw(:,:,twidx,:)-repmat(nanmean(tmpnw(:,:,baseidx,:),3),[1,1,length(twidx),1]))./repmat(nanmean(tmpnw(:,:,baseidx,:),3),[1,1,length(twidx),1]);
                                end
                                leftdattmp = squeeze(nanmean(leftdattmp,3));
                                leftdat = cat(3,leftdat,leftdattmp);
                            catch ME
                                disp(ME)
                                keyboard
                            end
                            clear tmpload
                            clear leftdattmp
                        end
                    end
                end
                
                
                nleft = size(leftdat,3);
                nright =  size(rightdat,3);
                
                if nleft>nright
                    tn = floor(trainp.*nright); %Number of trial in training set
                    tt = nright-tn; %Number of trials in test set
                else
                    tn = floor(trainp.*nleft); %Number of trial in training set
                    tt = nleft-tn; %Number of trials in test set
                end
                
                %Remove areas
                throwawayareas = find(cellfun(@isempty,BrainModel{midx}.Model.Rnames));
                throwawayareas = [throwawayareas; find(cellfun(@(X) ismember(X,{'OlfactoryBulb','fibrtracts','InfCol','SupColSens'}),BrainModel{midx}.Model.Rnames))];
                
                removepix = zeros(xpix,ypix);
                for areaid = 1:length(throwawayareas)
                    bounds = BrainModel{midx}.Model.Boundaries{throwawayareas(areaid)};
                    
                    for boundid = 1:length(bounds)
                        removepix(poly2mask(bounds{boundid}(:,1),bounds{boundid}(:,2),xpix,ypix)) = 1;
                    end
                end
                removepix = imfill(removepix,'holes');
                removepixvec = reshape(removepix,[xpix*ypix,1]);
                
                    %% Whole-brain approach
                    %WHOLE BRAIN ROC analysis
                
                
                if wholebrainana
                    thistimer = tic;
                    tmpright = reshape(rightdat,[xpix*ypix,nright]);
                    tmpleft = reshape(leftdat,[xpix*ypix,nleft]);
                  
                    tmpleft = tmpleft';
                    tmpright = tmpright';
                    
                    removenanpix = find(squeeze(sum(isnan(tmpleft),1)>0) | squeeze(sum(isnan(tmpright),1)>0) | removepixvec' == 1);
                    tmpleft(:,removenanpix) = [];
                    tmpright(:,removenanpix) = [];
                
                
                    tmpcp = zeros(size(tmphit));
                 
                    %% compute CP for whole brain
                    % how to get back to brainpixelssss???
                   for pixidx = 1:size(tmpleft,1)
                       
                    L1 = size(tmpleft(pixidx),2);
                    L2 = size(tmpright(pixidx),2);
                    labels = [ones(L1,1);zeros(L2,1)];
                    scores = [tmpleft(pixidx,:);tmpright(pixidx,:)];
                    [~,~,~,AUC1] = perfcurve(labels,scores,1);
                    tmpcp(pixidx) = AUC1;
                   end
                   
                   figure(Hwhole)
                   subplot(length(trialtypes),length(TW),(id-1)*length(TW)+twid)
                   aucform = true(xpix*ypix,1);
                   aucform(removenanpix') = 0;
                   newauc = nan(xpix*ypix,1);
                   newauc(aucform) = smooth2a(mean((tmpcp),2),2);% abs
                   h =imagesc(reshape(newauc,xpix,ypix));
                   hold on
                   scatter(BrainModel{midx}.Model.AllX,BrainModel{midx}.Model.AllY,'k.')
                   axis square
                   colorbar
                   set(h,'AlphaData',~isnan(reshape(newbeta,xpix,ypix)));
                   title([num2str(TW{twid}(1)) '-' num2str(TW{twid}(2)) ', ' trialtypes{id}, 'Stimulus Sensitivity per Pixel '])
                    
                    PERF{midx,sessioncount,ridx,id,twid}.CP = newauc;
                    disp(['CP analysis ' num2str(TW{twid}(1)) '-' num2str(TW{twid}(2)) ', ' trialtypes{id} ' took ' num2str(toc(thistimer)./60) ' minutes'])
                   
                    
                    
                  
end
            end
        end
end
end


