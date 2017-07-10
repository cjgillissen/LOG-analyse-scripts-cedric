function CPforWF(info,miceopt,StorePath,Stim2Check,baselinemethod,trialtypes,takeequalsample)

global UserQuestions
paths = info.paths;
logs = info.logs;
nrMouse = length(miceopt);

%% USER INPUT
nback = 20;
createvideosandfigurespermouse =0;
latencyana = 0;
contrastROIs = 0; %Contrast between conditions to determine ROI?
EvokedActivityROI = 1;
wholebrainana = 1;
newsize = [400 400];
scalefct = 0.5;
cpimrange = [0.3 0.7];

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
    TW = {visint,vistw,Delayearlytw,Responstw,Delaylatetw}
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
    xpix = newsize(1);
    ypix = newsize(2);
    %Load Alan Brain model
    BrainModel{midx} = load(fullfile(StorePath,mouse,'brainareamodel.mat'))
    
        for twid = 1:length(TW)
            for id = 1:length(trialtypes)
                sessioncount = 0;
                hitdat = [];
                errordat = [];
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
                        BASELINEMAT = imresize(BASELINEMAT,newsize,'bilinear');
                        %for timeline etc.
                        load(fullfile(StorePath,mouse,[mouse date],[mouse num2str(expnr)],rawdatfiles(strcmp({rawdatfiles(:).name},[mouse num2str(expnr) '_RawData_C' num2str(length(cvec)) '.mat'])).name));
                        clear conddata
                        twidx = find(timeline>=TW{twid}(1)&timeline<=TW{twid}(2));
                        baseidx = find(timeline>=basel(1) & timeline<=basel(2));
                        
                        % load in rawdata
                        hitidx = find(~cellfun(@isempty,(cellfun(@(X) strfind(X,'Hit'),ConditionNamestmp,'UniformOutput',0))));
                        if strcmp(Stim2Check,'DelayedOriTuningSound')
                            hitidx(ismember(hitidx,find(~cellfun(@isempty,(cellfun(@(X) strfind(X,'500'),ConditionNamestmp,'UniformOutput',0)))))) = [];
                        end
                        for i = 1:length(hitidx)
                            tmpload = load(fullfile(StorePath,mouse,[mouse date],[mouse num2str(expnr)],rawdatfiles(strcmp({rawdatfiles(:).name},[mouse num2str(expnr) '_RawData_C' num2str(hitidx(i)) '.mat'])).name));
                            try
                                rmtmp = removeidx(1:length(tmpload.ctrials{hitidx(i)}),hitidx(i))';
                                rm2tmp = ~ismember(tmpload.ctrials{hitidx(i)},fullfgtr);
                                rm3tmp = (rmtmp==1 | rm2tmp==1);
                                trialidx = tmpload.ctrials{hitidx(i)};
                                trialidx = trialidx(rm3tmp);
                                
                                hitdattmp =  zeros(400,400,length(twidx),sum(rm3tmp),'single');
                                resizedconddata = imresize(tmpload.conddata(:,:,:,rm3tmp),newsize,'bilinear');

                                for j = 1:100:newsize(1)
                                    tmpnw = single(resizedconddata(j:j+99,:,:,:))./permute(repmat(BASELINEMAT(j:j+99,:,trialidx),[1,1,1,size(resizedconddata,3)]),[1,2,4,3]);
                                    hitdattmp(j:j+99,:,:,:) = (tmpnw(:,:,twidx,:)-repmat(nanmean(tmpnw(:,:,baseidx,:),3),[1,1,length(twidx),1]))./repmat(nanmean(tmpnw(:,:,baseidx,:),3),[1,1,length(twidx),1]);
                                end
                                hitdattmp = squeeze(nanmean(hitdattmp,3));
                                hitdat = cat(3,hitdat,hitdattmp);
                           
                            catch ME
                                disp(ME)
                                keyboard
                            end
                            clear tmpload
                        end
                        
                        erroridx = find(~cellfun(@isempty,(cellfun(@(X) strfind(X,'Error'),ConditionNamestmp,'UniformOutput',0))));
                        if strcmp(Stim2Check,'DelayedOriTuningSound')
                            erroridx(ismember(erroridx,find(~cellfun(@isempty,(cellfun(@(X) strfind(X,'500'),ConditionNamestmp,'UniformOutput',0)))))) = [];
                        end
                        for i = 1:length(erroridx)
                            tmpload = load(fullfile(StorePath,mouse,[mouse date],[mouse num2str(expnr)],rawdatfiles(strcmp({rawdatfiles(:).name},[mouse num2str(expnr) '_RawData_C' num2str(erroridx(i)) '.mat'])).name));
                            try
                                
                                rmtmp = removeidx(1:length(tmpload.ctrials{erroridx(i)}),erroridx(i))';
                                rm2tmp = ~ismember(tmpload.ctrials{erroridx(i)},fullfgtr);
                                rm3tmp = (rmtmp==1 | rm2tmp==1);
                                trialidx = tmpload.ctrials{erroridx(i)};
                                trialidx = trialidx(rm3tmp);
                                errodattmp = nan(400,400,length(twidx),sum(rm3tmp),'single');
                                
                                %bilineair interpolation is equivalent to
                                %binning 
                                resizedconddata = imresize(tmpload.conddata(:,:,:,rm3tmp),newsize,'bilinear');
                                
                                for j = 1:100:newsize(1)
                                    tmpnw =  single(resizedconddata(j:j+99,:,:,:))./permute(repmat(BASELINEMAT(j:j+99,:,trialidx),[1,1,1,size(resizedconddata,3)]),[1,2,4,3]);
                                    errodattmp(j:j+99,:,:,:) = (tmpnw(:,:,twidx,:)-repmat(nanmean(tmpnw(:,:,baseidx,:),3),[1,1,length(twidx),1]))./repmat(nanmean(tmpnw(:,:,baseidx,:),3),[1,1,length(twidx),1]);
                                end
                                errodattmp = squeeze(nanmean(errodattmp,3));
                                errordat = cat(3,errordat,errodattmp);
                            catch ME
                                disp(ME)
                                keyboard
                            end
                            clear tmpload
                            clear leftdattmp
                        end
                    end
                end
                
                trainp = 0.9;
                nleft = size(errordat,3);
                nright =  size(hitdat,3);
                
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
                        removepix(poly2mask(bounds{boundid}(:,1).*scalefct,bounds{boundid}(:,2).*scalefct,xpix,ypix)) = 1;
                    end
                end
                removepix = imfill(removepix,'holes');
                removepixvec = reshape(removepix,[xpix*ypix,1]);
                if wholebrainana
                    thistimer = tic;
                    tmphit = reshape(hitdat,[xpix*ypix,nright]);
                    tmperror = reshape(errordat,[xpix*ypix,nleft]);
                    tmperror = tmperror';
                    tmphit = tmphit';
                    removenanpix = find(squeeze(sum(isnan(tmperror),1)>0) | squeeze(sum(isnan(tmphit),1)>0) | removepixvec' == 1);
                    tmperror(:,removenanpix) = [];
                    tmphit(:,removenanpix) = [];
                    tmpcp = zeros(size(tmperror,2),1);
                    tmplowerbound = zeros(size(tmperror,2),1);
                    tmpupperbound = zeros(size(tmperror,2),1);
                 
                    %% compute CP for whole brain
                    % how to get back to brainpixelssss???
                   for pixidx = 1:size(tmperror,2)
                       
                    L1 = size(tmperror(:,pixidx),1);
                    L2 = size(tmphit(:,pixidx),1);
                    labels = [ones(L1,1);zeros(L2,1)];
                    scores = [tmperror(:,pixidx);tmphit(:,pixidx)];
                    [~,~,~,AUC1] = perfcurve(labels,scores,1);
%                   [~,~,~,AUC1] = perfcurve(labels,scores,0,'Nboot',500); % error is positive
                    tmpcp(pixidx) = AUC1;
%                     tmplowerbound(pixidx) = AUC1(2);
%                     tmpupperbound(pixidx) = AUC1(3);
                   end
                   
                   Hwhole = figure;
                   aucform = true(xpix*ypix,1);
                   aucform(removenanpix') = 0;
                   newauc = nan(xpix*ypix,1);
                   newauc(aucform) = smooth2a(tmpcp,2,2);
                   h =imagesc(reshape(newauc,xpix,ypix),cpimrange);
                   hold on
                   scatter(BrainModel{midx}.Model.AllX.*scalefct,BrainModel{midx}.Model.AllY.*scalefct,'k.')
                   axis square
                   colormap(ActSupColorMap)
                   colorbar
                   set(h,'AlphaData',~isnan(reshape(newauc,xpix,ypix)));
                   title([num2str(TW{twid}(1)) '-' num2str(TW{twid}(2)) ', ' trialtypes{id}, 'Choice Probabilities per Pixel bilateral' mouse])
                   Perf{midx,id,twid}.CP = newauc;
                   Perf{midx,id,twid}.nrerror = L1;
                   Perf{midx,id,twid}.nrhit = L2;
%                    Perf{midx,id,twid}.lowerbound = tmplowerbound;
%                    Perf{midx,id,twid}.upperbound = tmpupperbound;
                   disp(['CP analysis ' num2str(TW{twid}(1)) '-' num2str(TW{twid}(2)) ', ' trialtypes{id} ' took ' num2str(toc(thistimer)./60) ' minutes'])
                   saveas(Hwhole,fullfile('C:\Users\gillissen\Desktop\Figures CP',['CP HvsE LandR' num2str(TW{twid}(1)) '-' num2str(TW{twid}(2)) trialtypes{id} mouse]))
                
                   % Average cp per area. Also plot confidence bounds. !! 
                   
                end
                
            end
            
        end
          
          save(fullfile('C:\Users\gillissen\Desktop\Figures CP',['Performance CP' strjoin(trialtypes) mouse]),'Perf')

    end
end


 