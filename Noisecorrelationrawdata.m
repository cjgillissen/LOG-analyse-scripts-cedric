function Noisecorrelationrawdata(info,miceopt,StorePath,Stim2Check,baselinemethod,trialtypes,takeequalsample)

paths = info.paths;
logs = info.logs;
nrMouse = size(paths,1);

%% USER INPUT
nback = 20;
createvideosandfigurespermouse =0;
latencyana = 0;
contrastROIs = 0; %Contrast between conditions to determine ROI?
EvokedActivityROI = 1;
newsize = [400 400];
scalefct = 0.5;
cpimrange = [0 0.9];
ReactOptloop = {'Hit','Error'};

normalizemethod = {'Zscore'} % {'substractavg'} {'None'}

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
    try referenceimage = uint8(imread(fullfile(StorePath,mouse,'\RefFile.bmp')));
    catch
        referenceimage = uint8(imread(fullfile(StorePath,mouse,'\RefFilepRF.bmp')));
    end
    
    xpix = newsize(1);
    ypix = newsize(2);
    %Load Alan Brain model
    BrainModel{midx} = load(fullfile(StorePath,mouse,'brainareamodel.mat'));
    
    % Create mask of approximate stimulus representation in v1
    
            evrois = load(fullfile('C:\Users\gillissen\Desktop\EvokedROIMouse',[mouse 'EvokedActivROIs']));
            
            %rightv1 stimulus representation
            rightv1 = zeros(newsize);
            rightv1(sub2ind(size(rightv1),round(evrois.rois{1}.xi.*scalefct),round(evrois.rois{1}.yi.*scalefct))) =1;
            rightv1 = imfill(rightv1');
            rightv1XY = regionprops(rightv1,'centroid');
            [xgrid, ygrid] = meshgrid(1:newsize(2), 1:newsize(1));
            centerRV1 = ((xgrid-rightv1XY.Centroid(1)).^2 + (ygrid-rightv1XY.Centroid(2)).^2) <= 5.^2;
            
            %leftv1 stimulus representation
            leftv1 = zeros(newsize);
            leftv1(sub2ind(size(leftv1),round(evrois.rois{2}.xi.*scalefct),round(evrois.rois{2}.yi.*scalefct))) =1;
            leftv1 = imfill(leftv1');
            leftv1XY = regionprops(leftv1,'centroid');
            [xgrid, ygrid] = meshgrid(1:newsize(2), 1:newsize(1));
            centerLV1 = ((xgrid-leftv1XY.Centroid(1)).^2 + (ygrid-leftv1XY.Centroid(2)).^2) <= 5.^2;
    
    
    
    %% For TW
    
        for twid = 1:length(TW)
            for id = 1:length(trialtypes)
                sessioncount = 0;
             for loopreactionidx = 1:length(ReactOptloop)
                  leftdat = [];
                  rightdat = [];
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
                        ReactionOpt = unique(reaction,'stable');
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
                        load(fullfile(StorePath,mouse,[mouse date],[mouse num2str(expnr)],rawdatfiles(strcmp({rawdatfiles(:).name},[mouse num2str(expnr) '_RawData_C' num2str(length(cvec)) '.mat'])).name),'timeline');
                        twidx = find(timeline>=TW{twid}(1)&timeline<=TW{twid}(2));
                        baseidx = find(timeline>=basel(1) & timeline<=basel(2));
                        
                        %% load in rawdata
                        
                       thistimer = tic;
                        
                        leftidx = find(~cellfun(@isempty,(cellfun(@(X) strfind(X,ReactOptloop{loopreactionidx}),ConditionNamestmp,'UniformOutput',0)))& ~cellfun(@isempty,(cellfun(@(X) strfind(X,'left'),ConditionNamestmp,'UniformOutput',0))));
                        if strcmp(Stim2Check,'DelayedOriTuningSound')
                            leftidx(ismember(leftidx,find(~cellfun(@isempty,(cellfun(@(X) strfind(X,'500'),ConditionNamestmp,'UniformOutput',0)))))) = [];
                        end
                        for i = 1:length(leftidx)
                            tmpload = load(fullfile(StorePath,mouse,[mouse date],[mouse num2str(expnr)],rawdatfiles(strcmp({rawdatfiles(:).name},[mouse num2str(expnr) '_RawData_C' num2str(leftidx(i)) '.mat'])).name));
                            try
                                rmtmp = ~removeidx(1:length(tmpload.ctrials{leftidx(i)}),leftidx(i))';
                                rm2tmp = ismember(tmpload.ctrials{leftidx(i)},fullfgtr);
                                rm3tmp = (rmtmp==1 & rm2tmp==1);
                                if any(rm3tmp==1)
                                trialidx = tmpload.ctrials{leftidx(i)};
                                trialidx = trialidx(rm3tmp);
                                
                                leftdattmp =  nan(400,400,length(twidx),sum(rm3tmp),'single');
                                resizedconddata = imresize(tmpload.conddata(:,:,:,rm3tmp),newsize,'bilinear');

                                if ~any(TW{twid}<0) % if TW is baseline use only baseline division
                                for j = 1:100:newsize(1)
                                    tmpnw = single(resizedconddata(j:j+99,:,:,:))./permute(repmat(BASELINEMAT(j:j+99,:,trialidx),[1,1,1,size(resizedconddata,3)]),[1,2,4,3]);
                                    leftdattmp(j:j+99,:,:,:) = (tmpnw(:,:,twidx,:)-repmat(nanmean(tmpnw(:,:,baseidx,:),3),[1,1,length(twidx),1]))./repmat(nanmean(tmpnw(:,:,baseidx,:),3),[1,1,length(twidx),1]);
                                end
                                else 
                                    
                                for j = 1:100:newsize(1)
                                    tmpnw = single(resizedconddata(j:j+99,:,:,:))./permute(repmat(BASELINEMAT(j:j+99,:,trialidx),[1,1,1,size(resizedconddata,3)]),[1,2,4,3]);
                                    leftdattmp(j:j+99,:,:,:) = tmpnw(:,:,twidx,:);
                                end 
                                end
                                end   
                                leftdattmp = squeeze(nanmean(leftdattmp,3));
                                leftdat = cat(3,leftdat,leftdattmp);
                           
                            catch ME
                                disp(ME)
                                keyboard
                            
                            clear tmpload
                            clear leftdattmp

                            end
                        end
                        
                        rightidx = find(~cellfun(@isempty,(cellfun(@(X) strfind(X,ReactOptloop{loopreactionidx}),ConditionNamestmp,'UniformOutput',0)))& ~cellfun(@isempty,(cellfun(@(X) strfind(X,'right'),ConditionNamestmp,'UniformOutput',0))));
                        if strcmp(Stim2Check,'DelayedOriTuningSound')
                            rightidx(ismember(rightidx,find(~cellfun(@isempty,(cellfun(@(X) strfind(X,'500'),ConditionNamestmp,'UniformOutput',0)))))) = [];
                        end
                        for i = 1:length(rightidx)
                            tmpload = load(fullfile(StorePath,mouse,[mouse date],[mouse num2str(expnr)],rawdatfiles(strcmp({rawdatfiles(:).name},[mouse num2str(expnr) '_RawData_C' num2str(rightidx(i)) '.mat'])).name));
                            try
                                
                                rmtmp = ~removeidx(1:length(tmpload.ctrials{rightidx(i)}),rightidx(i))';
                                rm2tmp = ismember(tmpload.ctrials{rightidx(i)},fullfgtr);
                                rm3tmp = (rmtmp==1 & rm2tmp==1);
                               if any(rm3tmp==1)
                                trialidx = tmpload.ctrials{rightidx(i)};
                                trialidx = trialidx(rm3tmp);
                                rightdattmp = nan(400,400,length(twidx),sum(rm3tmp),'single');
                                resizedconddata = imresize(tmpload.conddata(:,:,:,rm3tmp),newsize,'bilinear');
                                
                                if ~any(TW{twid}<0)
                                for j = 1:100:newsize(1)
                                    tmpnw =  single(resizedconddata(j:j+99,:,:,:))./permute(repmat(BASELINEMAT(j:j+99,:,trialidx),[1,1,1,size(resizedconddata,3)]),[1,2,4,3]);
                                    rightdattmp(j:j+99,:,:,:) = (tmpnw(:,:,twidx,:)-repmat(nanmean(tmpnw(:,:,baseidx,:),3),[1,1,length(twidx),1]))./repmat(nanmean(tmpnw(:,:,baseidx,:),3),[1,1,length(twidx),1]);
                                end
                                else 
                                    for j = 1:100:newsize(1)
                                    tmpnw =  single(resizedconddata(j:j+99,:,:,:))./permute(repmat(BASELINEMAT(j:j+99,:,trialidx),[1,1,1,size(resizedconddata,3)]),[1,2,4,3]);
                                    rightdattmp(j:j+99,:,:,:) = tmpnw(:,:,twidx,:);
                                    end
                                    
                                
                                end
                               end
                               rightdattmp = squeeze(nanmean(rightdattmp,3));
                               rightdat = cat(3,rightdat,rightdattmp);
                               
                            catch ME
                                disp(ME)
                                keyboard
                            end
                            clear tmpload
                            clear rightdattmp
                        end
                        end
                    end
                
                if strcmp(normalizemethod,'Zscore')
                    rightdat = zscore(rightdat,0,3);
                    rightdat(rightdat<-2.5|rightdat>2.5) = nan;
                    leftdat = zscore(leftdat,0,3);
                    leftdat(leftdat<-2.5|leftdat>2.5) = nan;
                elseif strcmp(normalizemethod,'substractavg')
                    meanright = nanmean(rightdat,3);
                    rightdat = bsxfun(@minus,rightdat,repmat(meanright,[1,1,size(rightdat,3)])); %substract only the mean per pixel insrtead of zscoring
                    meanleft = nanmean(leftdat,3);
                    leftdat = bsxfun(@minus,leftdat,repmat(meanleft,[1,1,size(leftdat,3)]));
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
                
                %% compute noise correlations per side. 
                % Reactions are in the parfor loop
                
%                   LEFT V1 seed when STIM is RIGHT
                    leftv1seed = rightdat;
                    seedtmpLv1 = repmat(centerLV1,[1,1,size(rightdat,3)]);
                    leftv1seed(~seedtmpLv1) = nan;
                    leftv1seed = nanmean(reshape(leftv1seed,[size(leftv1seed,1)*size(leftv1seed,2),size(leftv1seed,3)]),1);
                    corrvec = corr(leftv1seed',reshape(rightdat,[size(rightdat,1)*size(rightdat,2),size(rightdat,3)])','rows','pairwise');
                    corrmapRIGHTstimLEFTseed = reshape(corrvec,[size(rightdat,1),size(rightdat,2)]);
                    
%                   RIGHT V1 seed when STIM is LEFT
                    rightv1seed = leftdat;
                    seedtmpRv1 = repmat(centerRV1,[1,1,size(leftdat,3)]);
                    rightv1seed(~seedtmpRv1) = nan;
                    rightv1seed = nanmean(reshape(rightv1seed,[size(rightv1seed,1)*size(rightv1seed,2),size(rightv1seed,3)]),1);
                    corrvec = corr(rightv1seed',reshape(leftdat,[size(leftdat,1)*size(leftdat,2),size(leftdat,3)])','rows','pairwise');
                    corrmapLEFTstimRIGHTseed = reshape(corrvec,[size(leftdat,1),size(leftdat,2)]);
                     
%                   LEFT V1 seed when STIM is LEFT
                    leftv1seed = leftdat;
                    seedtmpLv1 = repmat(centerLV1,[1,1,size(leftdat,3)]);
                    leftv1seed(~seedtmpLv1) = nan;
                    leftv1seed = nanmean(reshape(leftv1seed,[size(leftv1seed,1)*size(leftv1seed,2),size(leftv1seed,3)]),1);
                    corrvec = corr(leftv1seed',reshape(leftdat,[size(leftdat,1)*size(leftdat,2),size(leftdat,3)])','rows','pairwise');
                    corrmapLEFTstimLEFTseed = reshape(corrvec,[size(leftdat,1),size(leftdat,2)]);
                    
% %                   RIGHT V1 seed when STIM is RIGHT
%                     rightv1seed = rightdat;
%                     seedtmpRv1 = repmat(centerRV1,[1,1,size(rightdat,3)]);
%                     rightv1seed(~seedtmpRv1) = nan;
%                     rightv1seed = nanmean(reshape(rightv1seed,[size(rightv1seed,1)*size(rightv1seed,2),size(rightv1seed,3)]),1);
%                     corrvec = corr(rightv1seed',reshape(rightdat,[size(rightdat,1)*size(rightdat,2),size(rightdat,3)]),'rows','pairwise');
%                     corrmapRIGHTstimRIGHTseed = reshape(corrvec,[size(rightdat,1),size(rightdat,2)]);
                    
                    
                    % make plot, average epr area.
                   %RIGHTstimLEFTseed 
                   LEFTV1 = figure;
                   quantvaldffleftv1 = quantile(abs(corrmapRIGHTstimLEFTseed(:)),0.90);
                   links =imagesc(corrmapRIGHTstimLEFTseed,cpimrange);
                   hold on
                   scatter(BrainModel{midx}.Model.AllX.*scalefct,BrainModel{midx}.Model.AllY.*scalefct,'k.')
                   axis square
                   assen = gca;
                   viscircles(assen,round(leftv1XY.Centroid),6,'Color','k','LineStyle','- -')
                   colormap(ActSupColorMap)
                   h = colorbar;
                   ylabel(h,'Pearson correlation coefficient') 
                   set(links,'AlphaData',~isnan(corrmapRIGHTstimLEFTseed));
                   title([num2str(TW{twid}(1)) '-' num2str(TW{twid}(2)) ', ' trialtypes{id}, 'NC Rightstim LEFT seed ' ReactOptloop{loopreactionidx}  mouse ' nrt=' num2str(size(rightdat,3))])
                   NCmat{midx,id,twid,loopreactionidx}.LEFTV1 = corrmapRIGHTstimLEFTseed;
                   NCmat{midx,id,twid,loopreactionidx}.nrtLEFTV1 = size(rightdat,3);
 
                   disp(['NC analysis ' num2str(TW{twid}(1)) '-' num2str(TW{twid}(2)) ', ' trialtypes{id} ' took ' num2str(toc(thistimer)./60) ' minutes'])
                   saveas(LEFTV1,fullfile('C:\Users\gillissen\Desktop\Figures NC',['NC Rightstim LEFT seed' ReactOptloop{loopreactionidx} num2str(TW{twid}(1)) '-' num2str(TW{twid}(2)) trialtypes{id} mouse]))
                   saveas(LEFTV1,fullfile('C:\Users\gillissen\Desktop\Figures NC',['NC Rightstim LEFT seed' ReactOptloop{loopreactionidx} num2str(TW{twid}(1)) '-' num2str(TW{twid}(2)) trialtypes{id} mouse]),'jpg')

                   %LEFTstimRIGHTseed
                   RIGHTV1 = figure;
                   quantvaldffrightv1 = quantile(abs(corrmapLEFTstimRIGHTseed(:)),0.90);
                   rechts =imagesc(corrmapLEFTstimRIGHTseed,cpimrange);
                   hold on
                   scatter(BrainModel{midx}.Model.AllX.*scalefct,BrainModel{midx}.Model.AllY.*scalefct,'k.')
                   axis square
                   assen = gca;
                   viscircles(assen,round(rightv1XY.Centroid),6,'Color','k','LineStyle','- -')
                   colormap(ActSupColorMap)
                   h = colorbar;
                   ylabel(h,'Pearson correlation coefficient')
                   set(rechts,'AlphaData',~isnan(corrmapLEFTstimRIGHTseed));
                   title([num2str(TW{twid}(1)) '-' num2str(TW{twid}(2)) ','' ' trialtypes{id}, 'NC Leftstim RIGHT seed ' ReactOptloop{loopreactionidx}  mouse ' nrt=' num2str(size(leftdat,3))])
                   NCmat{midx,id,twid,loopreactionidx}.RIGHTV1 = corrmapLEFTstimRIGHTseed;
                   NCmat{midx,id,twid,loopreactionidx}.nrtRIGHTV1 = size(leftdat,3);
                  
                   disp(['NC analysis ' num2str(TW{twid}(1)) '-' num2str(TW{twid}(2)) ', ' trialtypes{id} ' took ' num2str(toc(thistimer)./60) ' minutes'])
                   saveas(RIGHTV1,fullfile('C:\Users\gillissen\Desktop\Figures NC',['NC Right V1 seed' ReactOptloop{loopreactionidx} num2str(TW{twid}(1)) '-' num2str(TW{twid}(2)) trialtypes{id} mouse]))
                   saveas(RIGHTV1,fullfile('C:\Users\gillissen\Desktop\Figures NC',['NC Right V1 seed' ReactOptloop{loopreactionidx} num2str(TW{twid}(1)) '-' num2str(TW{twid}(2)) trialtypes{id} mouse]),'jpg')

                  
             end
            end
        end
        
        save(fullfile('C:\Users\gillissen\Desktop\Figures NC','NCmat'),'NCmat');
end

end


% %% Make contrast plots 
% 
% for midx = 1:nrMouse %For this mouse
% 
%     mouse = miceopt{midx};
%     mousecount = mousecount+1
%     try referenceimage = uint8(imread(fullfile(StorePath,mouse,'\RefFile.bmp')));
%     catch
%         referenceimage = uint8(imread(fullfile(StorePath,mouse,'\RefFilepRF.bmp')));
%     end
%     xpix = newsize(1);
%     ypix = newsize(2);
%     %Load Alan Brain model
%     BrainModel{midx} = load(fullfile(StorePath,mouse,'brainareamodel.mat'));
%     % Create mask of approximate stimulus representation in v1
%     evrois = load(fullfile('C:\Users\gillissen\Desktop\EvokedROIMouse',[mouse 'EvokedActivROIs']));
%     %rightv1 stimulus representation
%     rightv1 = zeros(newsize);
%     rightv1(sub2ind(size(rightv1),round(evrois.rois{1}.xi.*scalefct),round(evrois.rois{1}.yi.*scalefct))) =1;
%     rightv1 = imfill(rightv1');
%     rightv1XY = regionprops(rightv1,'centroid');
%     [xgrid, ygrid] = meshgrid(1:newsize(2), 1:newsize(1));
%     centerRV1 = ((xgrid-rightv1XY.Centroid(1)).^2 + (ygrid-rightv1XY.Centroid(2)).^2) <= 5.^2;
%     %leftv1 stimulus representation
%     leftv1 = zeros(newsize);
%     leftv1(sub2ind(size(leftv1),round(evrois.rois{2}.xi.*scalefct),round(evrois.rois{2}.yi.*scalefct))) =1;
%     leftv1 = imfill(leftv1');
%     leftv1XY = regionprops(leftv1,'centroid');
%     [xgrid, ygrid] = meshgrid(1:newsize(2), 1:newsize(1));
%     centerLV1 = ((xgrid-leftv1XY.Centroid(1)).^2 + (ygrid-leftv1XY.Centroid(2)).^2) <= 5.^2;
% 
%         for twid = 1:length(TW)
%             for id = 1:length(trialtypes)
%               for loopreactionidx = 1:length(ReactOptloop)     % not necessary, chooise yourself. 
%               
%                   
% 
% 
% end
% 
%             end
%         end
% end
% end
% 
%           
% %           save(fullfile('C:\Users\gillissen\Desktop\Figures CP',['Performance CP' strjoin(trialtypes) mouse]),'Perf')
% 
% 
%  