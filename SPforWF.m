function SPforWF(info,miceopt,StorePath,Stim2Check,baselinemethod,trialtypes,takeequalsample,resampling)

global UserQuestions
paths = info.paths;
logs = info.logs;
nrMouse = length(miceopt);
removebias = 1;
biaswindow = 15;

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
    BrainModel{midx} = load(fullfile(StorePath,mouse,'brainareamodel.mat'))
    
        for twid = 1:length(TW)
            for id = 1:length(trialtypes)
                sessioncount = 0;
                lefthitdat = [];
                lefterrordat = [];
                righthitdat = [];
                righterrordat = [];
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
                        if removebias
                          [~,biasmat] = extractbiasidx(LOG,ctrials,biaswindow);
                          removeidx(abs(biasmat')>0.3)=1;
                        end
                        
                     
                        % load in rawdata
                        %LEFTHIT
                        lefthitidx = find(~cellfun(@isempty,(cellfun(@(X) strfind(X,'left'),ConditionNamestmp,'UniformOutput',0)))&~cellfun(@isempty,(cellfun(@(X) strfind(X,'Hit'),ConditionNamestmp,'UniformOutput',0))));
                        if strcmp(Stim2Check,'DelayedOriTuningSound')
                            lefthitidx(ismember(lefthitidx,find(~cellfun(@isempty,(cellfun(@(X) strfind(X,'500'),ConditionNamestmp,'UniformOutput',0)))))) = [];
                        end
                        for i = 1:length(lefthitidx)
                            tmpload = load(fullfile(StorePath,mouse,[mouse date],[mouse num2str(expnr)],rawdatfiles(strcmp({rawdatfiles(:).name},[mouse num2str(expnr) '_RawData_C' num2str(lefthitidx(i)) '.mat'])).name));
                            try
                                rmtmp = ~removeidx(1:length(tmpload.ctrials{lefthitidx(i)}),lefthitidx(i))';
                                rm2tmp = ismember(tmpload.ctrials{lefthitidx(i)},fullfgtr);
                                rm3tmp = (rmtmp==1 & rm2tmp==1);
                                if any(rm3tmp==1)
                                trialidx = tmpload.ctrials{lefthitidx(i)};
                                trialidx = trialidx(rm3tmp);
                                
                                leftdattmp =  zeros(400,400,length(twidx),sum(rm3tmp),'single');
                                resizedconddata = imresize(tmpload.conddata(:,:,:,rm3tmp),newsize,'bilinear');

                                if ~any(TW{twid}<0)
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
                             
                                leftdattmp = imgaussfilt(squeeze(nanmean(leftdattmp,3)),3);
                                lefthitdat = cat(3,lefthitdat,leftdattmp);
                                end
                                
                           
                            catch ME
                                disp(ME)
                                keyboard
                            end
                            clear tmpload
                        end
                        
                        % load in rawdata
                        % LEFTERROR
                        lefterroridx = find(~cellfun(@isempty,(cellfun(@(X) strfind(X,'left'),ConditionNamestmp,'UniformOutput',0)))&~cellfun(@isempty,(cellfun(@(X) strfind(X,'Error'),ConditionNamestmp,'UniformOutput',0))));
                        if strcmp(Stim2Check,'DelayedOriTuningSound')
                            lefterroridx(ismember(lefterroridx,find(~cellfun(@isempty,(cellfun(@(X) strfind(X,'500'),ConditionNamestmp,'UniformOutput',0)))))) = [];
                        end
                        for i = 1:length(lefterroridx)
                            tmpload = load(fullfile(StorePath,mouse,[mouse date],[mouse num2str(expnr)],rawdatfiles(strcmp({rawdatfiles(:).name},[mouse num2str(expnr) '_RawData_C' num2str(lefterroridx(i)) '.mat'])).name));
                            try
                                rmtmp = ~removeidx(1:length(tmpload.ctrials{lefterroridx(i)}),lefterroridx(i))';
                                rm2tmp = ismember(tmpload.ctrials{lefterroridx(i)},fullfgtr);
                                rm3tmp = (rmtmp==1 & rm2tmp==1);
                                if any(rm3tmp==1)
                                trialidx = tmpload.ctrials{lefterroridx(i)};
                                trialidx = trialidx(rm3tmp);
                                
                                leftdattmp =  zeros(400,400,length(twidx),sum(rm3tmp),'single');
                                resizedconddata = imresize(tmpload.conddata(:,:,:,rm3tmp),newsize,'bilinear');

                                if ~any(TW{twid}<0)
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
                             
                                leftdattmp = imgaussfilt(squeeze(nanmean(leftdattmp,3)),3);
                                lefterrordat = cat(3,lefterrordat,leftdattmp);
                                end
                                
                           
                            catch ME
                                disp(ME)
                                keyboard
                            end
                            clear tmpload
                        end
                        
                        %RIGHTHIT
                        righthitidx = find(~cellfun(@isempty,(cellfun(@(X) strfind(X,'right'),ConditionNamestmp,'UniformOutput',0)))&~cellfun(@isempty,(cellfun(@(X) strfind(X,'Hit'),ConditionNamestmp,'UniformOutput',0))));
                        if strcmp(Stim2Check,'DelayedOriTuningSound')
                            righthitidx(ismember(righthitidx,find(~cellfun(@isempty,(cellfun(@(X) strfind(X,'500'),ConditionNamestmp,'UniformOutput',0)))))) = [];
                        end
                        for i = 1:length(righthitidx)
                            tmpload = load(fullfile(StorePath,mouse,[mouse date],[mouse num2str(expnr)],rawdatfiles(strcmp({rawdatfiles(:).name},[mouse num2str(expnr) '_RawData_C' num2str(righthitidx(i)) '.mat'])).name));
                            try
                                
                                rmtmp = ~removeidx(1:length(tmpload.ctrials{righthitidx(i)}),righthitidx(i))';
                                rm2tmp = ismember(tmpload.ctrials{righthitidx(i)},fullfgtr);
                                rm3tmp = (rmtmp==1 & rm2tmp==1);
                                if any(rm3tmp==1)
                                trialidx = tmpload.ctrials{righthitidx(i)};
                                trialidx = trialidx(rm3tmp);
                                rigtdattmp = nan(400,400,length(twidx),sum(rm3tmp),'single');
                                
                                %bilineair interpolation is equivalent to
                                %binning 
                                resizedconddata = imresize(tmpload.conddata(:,:,:,rm3tmp),newsize,'bilinear');
                                
                                if ~any(TW{twid}<0)
                                    for j = 1:100:newsize(1)
                                        tmpnw =  single(resizedconddata(j:j+99,:,:,:))./permute(repmat(BASELINEMAT(j:j+99,:,trialidx),[1,1,1,size(resizedconddata,3)]),[1,2,4,3]);
                                        rigtdattmp(j:j+99,:,:,:) = (tmpnw(:,:,twidx,:)-repmat(nanmean(tmpnw(:,:,baseidx,:),3),[1,1,length(twidx),1]))./repmat(nanmean(tmpnw(:,:,baseidx,:),3),[1,1,length(twidx),1]);
                                    end
                                else
                                    for j = 1:100:newsize(1)
                                        tmpnw =  single(resizedconddata(j:j+99,:,:,:))./permute(repmat(BASELINEMAT(j:j+99,:,trialidx),[1,1,1,size(resizedconddata,3)]),[1,2,4,3]);
                                        rigtdattmp(j:j+99,:,:,:) = tmpnw(:,:,twidx,:);
                                    end
                                end
                                
                                rigtdattmp = imgaussfilt(squeeze(nanmean(rigtdattmp,3)),3);
                                righthitdat = cat(3,righthitdat,rigtdattmp);
                                end
                                
                            catch ME
                                disp(ME)
                                keyboard
                            end
                            clear tmpload
                            clear leftdattmp
                        end
                        
                        %RIGHTERROR
                        righterroridx = find(~cellfun(@isempty,(cellfun(@(X) strfind(X,'right'),ConditionNamestmp,'UniformOutput',0)))&~cellfun(@isempty,(cellfun(@(X) strfind(X,'Error'),ConditionNamestmp,'UniformOutput',0))));
                        if strcmp(Stim2Check,'DelayedOriTuningSound')
                            righterroridx(ismember(righterroridx,find(~cellfun(@isempty,(cellfun(@(X) strfind(X,'500'),ConditionNamestmp,'UniformOutput',0)))))) = [];
                        end
                        for i = 1:length(righterroridx)
                            tmpload = load(fullfile(StorePath,mouse,[mouse date],[mouse num2str(expnr)],rawdatfiles(strcmp({rawdatfiles(:).name},[mouse num2str(expnr) '_RawData_C' num2str(righterroridx(i)) '.mat'])).name));
                            try
                                
                                rmtmp = ~removeidx(1:length(tmpload.ctrials{righterroridx(i)}),righterroridx(i))';
                                rm2tmp = ismember(tmpload.ctrials{righterroridx(i)},fullfgtr);
                                rm3tmp = (rmtmp==1 & rm2tmp==1);
                                if any(rm3tmp==1)
                                trialidx = tmpload.ctrials{righterroridx(i)};
                                trialidx = trialidx(rm3tmp);
                                rigtdattmp = nan(400,400,length(twidx),sum(rm3tmp),'single');
                                
                                %bilineair interpolation is equivalent to
                                %binning 
                                resizedconddata = imresize(tmpload.conddata(:,:,:,rm3tmp),newsize,'bilinear');
                                
                                if ~any(TW{twid}<0)
                                    for j = 1:100:newsize(1)
                                        tmpnw =  single(resizedconddata(j:j+99,:,:,:))./permute(repmat(BASELINEMAT(j:j+99,:,trialidx),[1,1,1,size(resizedconddata,3)]),[1,2,4,3]);
                                        rigtdattmp(j:j+99,:,:,:) = (tmpnw(:,:,twidx,:)-repmat(nanmean(tmpnw(:,:,baseidx,:),3),[1,1,length(twidx),1]))./repmat(nanmean(tmpnw(:,:,baseidx,:),3),[1,1,length(twidx),1]);
                                    end
                                else
                                    for j = 1:100:newsize(1)
                                        tmpnw =  single(resizedconddata(j:j+99,:,:,:))./permute(repmat(BASELINEMAT(j:j+99,:,trialidx),[1,1,1,size(resizedconddata,3)]),[1,2,4,3]);
                                        rigtdattmp(j:j+99,:,:,:) = tmpnw(:,:,twidx,:);
                                    end
                                end
                                
                                rigtdattmp = imgaussfilt(squeeze(nanmean(rigtdattmp,3)),3);
                                righterrordat = cat(3,righterrordat,rigtdattmp);
                                end
                                
                            catch ME
                                disp(ME)
                                keyboard
                            end
                            clear tmpload
                            clear leftdattmp
                        end
                    end
                end
                
                %% Resampling of reactions per side
                % resample inside left vs right
                
                if resampling == 1
                    if size(righthitdat,3)>size(righterrordat,3)
                       rightdat = cat(3,righterrordat,righthitdat(:,:,randperm(size(righthitdat,3),size(righterrordat,3))));
                    elseif size(righterrordat,3)>size(righthitdat,3)
                       rightdat = cat(3,righthitdat,righterrordat(:,:,randperm(size(righterrordat,3),size(righthitdat,3))));
                    end

                    if size(lefthitdat,3)>size(lefterrordat,3)
                        leftdat = cat(3,lefterrordat,lefthitdat(:,:,randperm(size(lefthitdat,3),size(lefterrordat,3))));
                    elseif size(lefterrordat,3)>size(lefthitdat,3)
                        leftdat = cat(3,lefthitdat,lefterrordat(:,:,randperm(size(lefterrordat,3),size(lefthitdat,3))));
                    end
                % resample every condition
                elseif resampling == 2
                    
                    minelement = min([size(righthitdat,3) size(righterrordat,3) size(lefthitdat,3) size(lefterrordat,3)]);
                    righterrordat = righterrordat(:,:,randperm(size(righterrordat,3),minelement));
                    righthitdat = righthitdat(:,:,randperm(size(righthitdat,3),minelement));
                    rightdat = cat(3,righterrordat,righthitdat);
                    lefterrordat = lefterrordat(:,:,randperm(size(lefterrordat,3),minelement));
                    lefthitdat = lefthitdat(:,:,randperm(size(lefthitdat,3),minelement));
                    leftdat = cat(3,lefterrordat,lefthitdat);
                
                elseif resampling ==0
                    rightdat = cat(3,righterrordat,righthitdat);
                    leftdat = cat(3,lefterrordat,lefthitdat);
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
                
                %% ROC
                if wholebrainana
                    thistimer = tic;
                    tmpleft = reshape(leftdat,[xpix*ypix,size(leftdat,3)]);
                    tmpright = reshape(rightdat,[xpix*ypix,size(rightdat,3)]);
                    tmpright = tmpright';
                    tmpleft = tmpleft';
                    removenanpix = find(squeeze(sum(isnan(tmpright),1)>0) | squeeze(sum(isnan(tmpleft),1)>0) | removepixvec' == 1);
                    tmpright(:,removenanpix) = [];
                    tmpleft(:,removenanpix) = [];
                    tmpcp = zeros(size(tmpright,2),1);
                    tmplowerbound = zeros(size(tmpright,2),1);
                    tmpupperbound = zeros(size(tmpright,2),1);
                 
                    %% compute CP for whole brain
                     L1 = size(tmpright,1);
                     L2 = size(tmpleft,1);
                     labels = [ones(L1,1);zeros(L2,1)];
                   parfor (pixidx = 1:size(tmpright,2),4)
                    scores = [tmpright(:,pixidx);tmpleft(:,pixidx)];
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
                   quantvallower = quantile(abs(newauc),0.01);
                   quantvalhigher = quantile(abs(newauc),0.99);
                   h =imagesc(reshape(newauc,xpix,ypix),cpimrange);
                   hold on
                   scatter(BrainModel{midx}.Model.AllX.*scalefct,BrainModel{midx}.Model.AllY.*scalefct,'k.')
                   axis square
                   colormap(ActSupColorMap)
                   colorbar
                   set(h,'AlphaData',~isnan(reshape(newauc,xpix,ypix)));
                   title([num2str(TW{twid}(1)) '-' num2str(TW{twid}(2)) ', ' trialtypes{id}, 'Stimulus Probabilities Left vs Right ' mouse])
                   Perf{midx,id,twid}.SP = newauc;
                   Perf{midx,id,twid}.nrright = L1;
                   Perf{midx,id,twid}.nrleft = L2;
%                    Perf{midx,id,twid}.lowerbound = tmplowerbound;
%                    Perf{midx,id,twid}.upperbound = tmpupperbound;
                   disp(['SP analysis ' num2str(TW{twid}(1)) '-' num2str(TW{twid}(2)) ', ' trialtypes{id} ' took ' num2str(toc(thistimer)./60) ' minutes'])
                   saveas(Hwhole,fullfile('C:\Users\gillissen\Desktop\Figures SP',['SP LeftvsRight' num2str(TW{twid}(1)) '-' num2str(TW{twid}(2)) trialtypes{id} mouse]))
                   
                end
                
            end
            
        end
          
          save(fullfile('C:\Users\gillissen\Desktop\Figures SP',['Performance SP' strjoin(trialtypes) mouse]),'Perf')

end

end


 