function ROIselectionFIG(info,miceopt,StorePath,DataDirectory,Stim2Check,baselinemethod,originallim,trialtypes,plotlim,takeequalsample)
global UserQuestions
paths = info.paths;
logs = info.logs;
nrMouse = size(paths,1);

%% USER INPUT
nback = 20;
createvideosandfigurespermouse =1;
latencyana = 0;
contrastROIs = 1; %Contrast between conditions to determine ROI?
EvokedActivityROI = 1;
removeErrors = 0 %If you don't want to include Errors, make it 1
removeHits = 0 %if you don't want to include hits, make it 1
removeMiss = 0 %If you don't want to include misses, make it 1
fgtw = [ 250]; %FOR FG     


%% STOP USER INPUT

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
    referenceimage = uint8(imread(fullfile(StorePath,mouse,'\RefFile.bmp')));
    sessioncount = 0;
    
    %Load Alan Brain model
    BrainModel{midx} = load(fullfile(StorePath,mouse,'brainareamodel.mat'))
    
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
            
            ConditionNamestmp(cvec)
            %% Load LEFTvsRight data and other datafiles
            iddone = zeros(1,length(trialtypes));
            for id = 1:length(trialtypes)
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
                
                rawdatfiles = dir(fullfile(StorePath,mouse,[mouse date],[mouse num2str(expnr)],[mouse num2str(expnr) '_RawData*']));
                iddone(id) = 1;
                
                if isfield(LEFTVSRIGHT,'SUMSqr')
                    %replace empty matrices with zeros
                    [r,c] = find(cell2mat(cellfun(@isempty,LEFTVSRIGHT.SUMSqr,'UniformOutput',0)));
                    [rtempl,ctempl] = find(~cell2mat(cellfun(@isempty,LEFTVSRIGHT.SUMSqr,'UniformOutput',0)),1);
                    if ~isempty(r)
                        for idx = 1:length(r)
                            try
                                LEFTVSRIGHT.SUMSqr{r(idx),c(idx)} = zeros(size(LEFTVSRIGHT.SUMSqr{rtempl,ctempl}));
                            catch ME
                                disp(ME)
                                keyboard
                            end
                            
                        end
                    end
                    %Make 0 of nans, cause summation with 0's is going
                    %wrong.
                    for rr = 1:size(LEFTVSRIGHT.SUMSqr,1)
                        for cc = 1:size(LEFTVSRIGHT.SUMSqr,2)
                            LEFTVSRIGHT.SUMSqr{rr,cc}(isnan(LEFTVSRIGHT.SUMSqr{rr,cc})) = 0;
                        end
                    end
                    if ~isempty(SUMSQALL{id})
                        if size(LEFTVSRIGHT.SUMSqr,2) < size(SUMSQALL{id},2)
                            c = find(~ismember(unique([ReactionOpt{:}]),LEFTVSRIGHT.ReactionOpt));
                            for idx = 1:length(c)
                                for rowidx = 1:size(SUMSQALL{id},1)
                                    LEFTVSRIGHT.SUMSqr{rowidx,c(idx)} = zeros(size(SUMSQALL{id}{rowidx,c(idx)}));
                                end
                            end
                        end
                        try
                            SUMSQALL{id} = cellfun(@(X,Y) X+Y,SUMSQALL{id},LEFTVSRIGHT.SUMSqr(:,~ismember(LEFTVSRIGHT.ReactionOpt,'TooEarly')),'UniformOutput',0);
                        catch ME
                            disp(ME)
                            keyboard
                        end
                        
                    else
                        SUMSQALL{id} = LEFTVSRIGHT.SUMSqr(:,~ismember(LEFTVSRIGHT.ReactionOpt,'TooEarly'));
                    end
                else
                    warning(['No SUM of squared dFF variable, redo FGSpecificAnalysis for this session!'])
                    disp('Skipping session')
                    continue
                end
                % SUMdFF, SUMSQdFF & NRT need to be saved out. throw away
                % LEFTVSRIGHT after that
                if ~isfield(LEFTVSRIGHT,'SumdFF')
                    error(['Redo ConditionAveraginScript for session ' fullfile(StorePath,mouse,[mouse date],[mouse num2str(expnr)])])
                end
                
                %replace empty matrices with zeros
                [r,c] = find(cell2mat(cellfun(@isempty,LEFTVSRIGHT.SumdFF,'UniformOutput',0)));
                [rtempl,ctempl] = find(~cell2mat(cellfun(@isempty,LEFTVSRIGHT.SumdFF,'UniformOutput',0)),1);
                if ~isempty(r)
                    for idx = 1:length(r)
                        LEFTVSRIGHT.SumdFF{r(idx),c(idx)} = zeros(size(LEFTVSRIGHT.SumdFF{rtempl,ctempl}));
                        LEFTVSRIGHT.nrtPerPix{r(idx),c(idx)} = zeros(size(LEFTVSRIGHT.nrtPerPix{rtempl,ctempl}));
                    end
                end
                
                %replace empty with 0
                [r,c] = find(cellfun(@isempty,LEFTVSRIGHT.nrt));
                if ~isempty(r)
                    for idx = 1:length(r)
                        LEFTVSRIGHT.nrt{r(idx),c(idx)}=0;
                    end
                end
                %Make 0 of nans, cause summation with nans is going
                %wrong.
                for rr = 1:size(LEFTVSRIGHT.dFFav,1)
                    for cc = 1:size(LEFTVSRIGHT.dFFav,2)
                        LEFTVSRIGHT.SumdFF{rr,cc}(isnan(LEFTVSRIGHT.SumdFF{rr,cc})) = 0;
                    end
                end
                
                if ~isempty(SumdFF{id})
                    if size(LEFTVSRIGHT.SumdFF,2) < size(SumdFF{id},2)
                        c = find(~ismember(unique([ReactionOpt{:}]),LEFTVSRIGHT.ReactionOpt));
                        for idx = 1:length(c)
                            for rowidx = 1:size(SumdFF{id},1)
                                LEFTVSRIGHT.SumdFF{rowidx,c(idx)} = zeros(size(SumdFF{id}{rowidx,c(idx)}));
                                LEFTVSRIGHT.nrt{rowidx,c(idx)} = 0;
                                LEFTVSRIGHT.nrtPerPix{rowidx,c(idx)} = zeros(size(SumdFF{id}{rowidx,c(idx)}));
                            end
                        end
                    end
                    
                    
                    try
                        SumdFF{id} = cellfun(@(X,Y) X+Y, LEFTVSRIGHT.SumdFF(:,~ismember(LEFTVSRIGHT.ReactionOpt,'TooEarly')),SumdFF{id},'UniformOutput',0);
                        nrtperpix{id} = cellfun(@(X,Y) X+Y, LEFTVSRIGHT.nrtPerPix(:,~ismember(LEFTVSRIGHT.ReactionOpt,'TooEarly')),nrtperpix{id},'UniformOutput',0);
                    catch ME
                        disp(ME)
                        keyboard
                    end
                else
                    SumdFF{id} = LEFTVSRIGHT.SumdFF(:,~ismember(LEFTVSRIGHT.ReactionOpt,'TooEarly'));
                    nrtperpix{id} = LEFTVSRIGHT.nrtPerPix(:,~ismember(LEFTVSRIGHT.ReactionOpt,'TooEarly'));
                end
                
                
                
                if ~isempty(nrt{id})
                    nrt{id} = cellfun(@(X,Y) X+Y,nrt{id},LEFTVSRIGHT.nrt(:,~ismember(LEFTVSRIGHT.ReactionOpt,'TooEarly')),'UniformOutput',0);
                else
                    nrt{id} = LEFTVSRIGHT.nrt(:,~ismember(LEFTVSRIGHT.ReactionOpt,'TooEarly'));
                end
                meanRT{sessioncount,id} = LEFTVSRIGHT.meanRT;
                stdRT{sessioncount,id} = LEFTVSRIGHT.stdRT;
                SideOpt{sessioncount,id} = LEFTVSRIGHT.SideOpt;
                ReactionOpt{sessioncount,id} = LEFTVSRIGHT.ReactionOpt;
                
                clear LEFTVSRIGHT
                load(fullfile(StorePath,mouse,[mouse date],[mouse num2str(expnr)],rawdatfiles(strcmp({rawdatfiles(:).name},[mouse num2str(expnr) '_RawData_C' num2str(length(cvec)) '.mat'])).name));
                clear conddata
                timevectmp{mousecount,sessioncount,id} = timeline(timeline>=originallim(1)&timeline<=originallim(2));
                
            end
        end
    end
    
    %% Differences
    maxnrReactOpt = max(max(cellfun(@length,ReactionOpt)));
    AllReactionOpt = unique([ReactionOpt{:}],'stable');
    maxnrReactOpt =maxnrReactOpt- sum(strcmp(AllReactionOpt,'TooEarly'));
    AllReactionOpt(strcmp(AllReactionOpt,'TooEarly')) = [];
    if removeErrors
        AllReactionOpt(strcmp(AllReactionOpt,'Error')) = [];
    end
    if removeHits
        AllReactionOpt(strcmp(AllReactionOpt,'Hit')) = [];
    end
    if removeMiss
        AllReactionOpt(strcmp(AllReactionOpt,'Miss')) = [];
    end

    maxnrSideOpt = max(max(cellfun(@length,SideOpt)));
    AllSideOpt = unique([SideOpt{:}],'stable');
    fgopt = find(iddone);
    [xpix ypix ~] = size(SumdFF{fgopt(1)}{1,1});
    
    %Check timevecs: if not the same; warning!
    timevec = cat(1,timevectmp{:});
    timevec = unique(timevec,'rows');
    if size(timevec,1) > 1
        error('Multiple versions of timevec...')
    end
    
    for fgid = fgopt
        NewdFF{fgid} = cellfun(@(X,Y) X./Y,SumdFF{fgid},nrtperpix{fgid},'UniformOutput',0);
    end
    
    
    if ~exist(fullfile(StorePath,'Figures',['Baseline' num2str(baselinemethod) '_eqsample' num2str(takeequalsample)]))
        mkdir(fullfile(StorePath,'Figures',['Baseline' num2str(baselinemethod) '_eqsample' num2str(takeequalsample)]))
    end
    
    %% Brainmask
    brainmask = zeros(800,800);
    for i = 1:length(BrainModel{midx}.Model.Regions)
        Borders = BrainModel{midx}.Model.Boundaries{i};
        for j = 1:length(Borders)
            tmp = poly2mask(Borders{j}(:,1),Borders{j}(:,2),800,800);
            tmp = imfill(tmp,'holes');
            brainmask(tmp) = 1;
        end
    end
    brainmask = imfill(brainmask,'holes');
    
    %% latency (Matt 26/7/2016)
    if latencyana
        %LEFT AND RIGHT
        Latencies = cell(max(fgopt),maxnrReactOpt,maxnrSideOpt);
        RSqrdFits = Latencies;
        rscale = 1/4;
        for id = fgopt
            tmpIm = cellfun(@(X) imresize(X,rscale),NewdFF{id},'UniformOutput',0); %Rescale image
            for stidx = 1:maxnrSideOpt
                h = figure('name',[AllSideOpt{stidx}, ' ' trialtypes{id}]);
                for ridx = 1:maxnrReactOpt
                    Latencies{id,ridx,stidx} = nan(size(tmpIm{stidx,ridx},1),size(tmpIm{stidx,ridx},2));
                    RSqrdFits{id,ridx,stidx} = nan(size(tmpIm{stidx,ridx},1),size(tmpIm{stidx,ridx},2));
                    for rowid=1:size(tmpIm{stidx,ridx},1)
                        for colid = 1:size(tmpIm{stidx,ridx},2)
                            ix = find(timevec>=0 & timevec<=500); %only visual
                            bix = find(timevec<0&timevec>=-300);
                            R = squeeze(tmpIm{stidx,ridx}(rowid,colid,ix));
                            Rsd = nanstd(tmpIm{stidx,ridx}(rowid,colid,bix));
                            Rm = nanmean(tmpIm{stidx,ridx}(rowid,colid,bix));
                            if sum(~isnan(R))<length(ix) || nanmax(R)<Rm+2.5*Rsd
                                continue
                            end
                            % figure;plot(T,R)
                            %                         figure,plot(timevec(ix),R); hold on
                            % SR = smooth(R,5);%plot(T(ix),SR(ix))
                            
                            % cutoff = T(ix(find(diff(SR(ix))<0,1,'first')));
                            [~,cutoff] = max(R);
                            cutoff = timevec(ix(cutoff));
                            % h = line([cutoff cutoff],get(gca,'YLim'));
                            % set(h,'Color',[0 0 0])
                            
                            ix = find(timevec>=0 & timevec<=cutoff);
                            
                            %Fit a cumaltive Gauss
                            [y,params] = cumgaussfit(timevec(ix),squeeze(tmpIm{stidx,ridx}(rowid,colid,ix))',[150 20 0.03 0],0);
                            rsqrd = corrcoef(y,squeeze(tmpIm{stidx,ridx}(rowid,colid,ix))').^2;
                            if length(rsqrd)>1
                                RSqrdFits{id,ridx,stidx}(rowid,colid) = rsqrd(2);
                            else
                                RSqrdFits{id,ridx,stidx}(rowid,colid) = rsqrd(1);
                            end
                            % hold on,plot(T(ix),y,'r')
                            
                            %Latency = mean - std
                            Latencies{id,ridx,stidx}(rowid,colid) = params(1)-params(2);
                            
                            %                         if savefigures && ~isnan(rs)
                            %                             %timepoint of max value
                            %                             h = line([cutoff cutoff],get(gca,'YLim'));  set(h,'Color',[0 0 0])
                            %                             % plot cumulative GaussFit
                            %                             hold on,plot(timevec(ix),y,'r')
                            %                             h = line([Lat Lat],get(gca,'YLim'));
                            %                             set(h,'Color',[0 1 0])
                            %                         end
                            
                        end
                    end
                    try
                        qval = quantile(Latencies{id,ridx,stidx}(:),0.9);
                        if ~isnan(qval)
                            subplot(maxnrReactOpt,2,(ridx-1)*2+1); htmp = imagesc(Latencies{id,ridx,stidx},[0 qval]);
                            set(htmp,'Alphadata',~isnan(Latencies{id,ridx,stidx}))
                            hold on
                            plot(round(BrainModel{midx}.Model.AllX*rscale),round(BrainModel{midx}.Model.AllY*rscale),'.k','MarkerSize',5)
                            axis square
                            colorbar
                            if ridx == 1
                                title('Latencies')
                            end
                            ylabel([AllReactionOpt{ridx}])
                        end
                        qval = quantile(RSqrdFits{id,ridx,stidx}(:),0.9);
                        if ~isnan(qval)
                            
                            subplot(maxnrReactOpt,2,(ridx-1)*2+2); htmp = imagesc(RSqrdFits{id,ridx,stidx});
                            set(htmp,'Alphadata',~isnan(Latencies{id,ridx,stidx}))
                            
                            hold on
                            plot(round(BrainModel{midx}.Model.AllX*rscale),round(BrainModel{midx}.Model.AllY*rscale),'.k','MarkerSize',5)
                            
                            axis square
                            colorbar
                            if ridx == 1
                                title('RSqrd')
                            end
                            
                        end
                    end
                end
                saveas(h,fullfile(StorePath,'Figures',['Baseline' num2str(baselinemethod) '_eqsample' num2str(takeequalsample)],[miceopt{midx} '_' AllSideOpt{stidx}, ' ' trialtypes{id} '_Latencies.bmp']))
                saveas(h,fullfile(StorePath,'Figures',['Baseline' num2str(baselinemethod) '_eqsample' num2str(takeequalsample)],[miceopt{midx} '_' AllSideOpt{stidx}, ' ' trialtypes{id} '_Latencies.fig']))
                saveas(h,fullfile(StorePath,'Figures',['Baseline' num2str(baselinemethod) '_eqsample' num2str(takeequalsample)],[miceopt{midx} '_' AllSideOpt{stidx}, ' ' trialtypes{id} '_Latencies.eps']))
                close(h)
            end
        end
        
        % LEFT AND RIGHT AVERAGED
        Latencies = cell(max(fgopt),maxnrReactOpt);
        RSqrdFits = Latencies;
        rscale = 1/4;
        for id = fgopt
            tmpIm = cellfun(@(X) imresize(X,rscale),NewdFF{id},'UniformOutput',0); %Rescale image
            h = figure('name',[trialtypes{id}]);
            for ridx = 1:maxnrReactOpt
                Latencies{id,ridx} = nan(size(tmpIm{1,ridx},1),size(tmpIm{1,ridx},2));
                RSqrdFits{id,ridx} = nan(size(tmpIm{1,ridx},1),size(tmpIm{1,ridx},2));
                for rowid=1:size(tmpIm{1,ridx},1)
                    for colid = 1:size(tmpIm{1,ridx},2)
                        ix = find(timevec>=0 & timevec<=500); %only visual
                        bix = find(timevec<0&timevec>=-300);
                        tmpdat = nanmean(cat(4,tmpIm{1,ridx}(rowid,colid,ix),tmpIm{2,ridx}(rowid,colid,ix)),4);
                        tmpdatbase = nanmean(cat(4,tmpIm{1,ridx}(rowid,colid,bix),tmpIm{2,ridx}(rowid,colid,bix)),4);
                        %                     tmpdat1 = tmpIm{1,ridx}(rowid,colid,ix);
                        %                     tmpdat2 = tmpIm{2,ridx}(rowid,colid,ix);
                        %                     R1 = squeeze(tmpdat1); R2 = squeeze(tmpdat2);
                        %                     Rsd1 = nanstd(R1); Rsd2 = nanstd(R2);
                        %                     Rm1 = nanmean(R1); Rm2 = nanmean(R2);
                        R = squeeze(tmpdat);
                        Rsd = nanstd(tmpdatbase);
                        Rm = nanmean(tmpdatbase);
                        if sum(~isnan(R))<length(ix) || nanmax(R)<Rm+2.5*Rsd % || nanmax(R2)<Rm2+Rsd2
                            continue
                        end
                        % figure;plot(T,R)
                        %                         figure,plot(timevec(ix),R); hold on
                        % SR = smooth(R,5);%plot(T(ix),SR(ix))
                        
                        % cutoff = T(ix(find(diff(SR(ix))<0,1,'first')));
                        [~,cutoff] = max(R);
                        R = R(1:cutoff);
                        cutoff = timevec(ix(cutoff));
                        % h = line([cutoff cutoff],get(gca,'YLim'));
                        % set(h,'Color',[0 0 0])
                        
                        ix = find(timevec>=0 & timevec<=cutoff);
                        
                        %Fit a cumaltive Gauss
                        [y,params] = cumgaussfit(timevec(ix),R',[150 20 0.03 0],0);
                        rsqrd = corrcoef(y,R').^2;
                        if length(rsqrd)>1
                            RSqrdFits{id,ridx}(rowid,colid) = rsqrd(2);
                        else
                            RSqrdFits{id,ridx}(rowid,colid) = rsqrd(1);
                        end
                        % hold on,plot(T(ix),y,'r')
                        
                        %Latency = mean - std
                        Latencies{id,ridx}(rowid,colid) = params(1)-params(2);
                        
                        %                         if savefigures && ~isnan(rs)
                        %                             %timepoint of max value
                        %                             h = line([cutoff cutoff],get(gca,'YLim'));  set(h,'Color',[0 0 0])
                        %                             % plot cumulative GaussFit
                        %                             hold on,plot(timevec(ix),y,'r')
                        %                             h = line([Lat Lat],get(gca,'YLim'));
                        %                             set(h,'Color',[0 1 0])
                        %                         end
                        
                    end
                end
                try
                    subplot(maxnrReactOpt,2,(ridx-1)*2+1); htmp = imagesc(Latencies{id,ridx},[0 350]);
                    colormap(autumn)
                    set(htmp,'Alphadata',~isnan(Latencies{id,ridx}))
                    hold on
                    plot(round(BrainModel{midx}.Model.AllX*rscale),round(BrainModel{midx}.Model.AllY*rscale),'.k','MarkerSize',5)
                    axis square
                    colorbar
                    if ridx == 1
                        title('Latencies')
                    end
                    ylabel([AllReactionOpt{ridx}])
                    
                    subplot(maxnrReactOpt,2,(ridx-1)*2+2); htmp = imagesc(RSqrdFits{id,ridx},[0 1]);
                    set(htmp,'Alphadata',~isnan(Latencies{id,ridx}))
                    
                    hold on
                    plot(round(BrainModel{midx}.Model.AllX*rscale),round(BrainModel{midx}.Model.AllY*rscale),'.k','MarkerSize',5)
                    
                    axis square
                    colorbar
                    if ridx == 1
                        title('RSqrd')
                    end
                    
                catch ME
                    disp(ME)
                    keyboard
                end
            end
            saveas(h,fullfile(StorePath,'Figures',['Baseline' num2str(baselinemethod) '_eqsample' num2str(takeequalsample)],[miceopt{midx} ' ' trialtypes{id} '_Latencies.bmp']))
            saveas(h,fullfile(StorePath,'Figures',['Baseline' num2str(baselinemethod) '_eqsample' num2str(takeequalsample)],[miceopt{midx} ' ' trialtypes{id} '_Latencies.fig']))
            saveas(h,fullfile(StorePath,'Figures',['Baseline' num2str(baselinemethod) '_eqsample' num2str(takeequalsample)],[miceopt{midx} ' ' trialtypes{id} '_Latencies.eps']))
            close(h)
            
        end
        
        % LEFT - RIGHT AVERAGED
        Latencies = cell(max(fgopt),maxnrReactOpt);
        RSqrdFits = Latencies;
        rscale = 1/4;
        for id = fgopt
            tmpIm = cellfun(@(X) imresize(X,rscale),NewdFF{id},'UniformOutput',0); %Rescale image
            h = figure('name',[trialtypes{id}]);
            for ridx = 1:maxnrReactOpt
                Latencies{id,ridx} = nan(size(tmpIm{1,ridx},1),size(tmpIm{1,ridx},2));
                RSqrdFits{id,ridx} = nan(size(tmpIm{1,ridx},1),size(tmpIm{1,ridx},2));
                for rowid=1:size(tmpIm{1,ridx},1)
                    for colid = 1:size(tmpIm{1,ridx},2)
                        ix = find(timevec>=0 & timevec<=500); %only visual
                        bix = find(timevec<0&timevec>=-300);
                        tmpdat = tmpIm{1,ridx}(rowid,colid,ix)-tmpIm{2,ridx}(rowid,colid,ix);
                        tmpdatbase = tmpIm{1,ridx}(rowid,colid,bix)-tmpIm{2,ridx}(rowid,colid,bix);
                        %                     tmpdat1 = tmpIm{1,ridx}(rowid,colid,ix);
                        %                     tmpdat2 = tmpIm{2,ridx}(rowid,colid,ix);
                        %                     R1 = squeeze(tmpdat1); R2 = squeeze(tmpdat2);
                        %                     Rsd1 = nanstd(R1); Rsd2 = nanstd(R2);
                        %                     Rm1 = nanmean(R1); Rm2 = nanmean(R2);
                        R = squeeze(tmpdat);
                        Rsd = nanstd(tmpdatbase);
                        Rm = nanmean(tmpdatbase);
                        if sum(~isnan(R))<length(ix) || nanmax(R)<Rm+2.5*Rsd % || nanmax(R2)<Rm2+Rsd2
                            continue
                        end
                        % figure;plot(T,R)
                        %                         figure,plot(timevec(ix),R); hold on
                        % SR = smooth(R,5);%plot(T(ix),SR(ix))
                        
                        % cutoff = T(ix(find(diff(SR(ix))<0,1,'first')));
                        [~,cutoff] = max(R);
                        R = R(1:cutoff);
                        cutoff = timevec(ix(cutoff));
                        % h = line([cutoff cutoff],get(gca,'YLim'));
                        % set(h,'Color',[0 0 0])
                        
                        ix = find(timevec>=0 & timevec<=cutoff);
                        
                        %Fit a cumaltive Gauss
                        [y,params] = cumgaussfit(timevec(ix),R',[150 20 0.03 0],0);
                        rsqrd = corrcoef(y,R').^2;
                        if length(rsqrd)>1
                            RSqrdFits{id,ridx}(rowid,colid) = rsqrd(2);
                        else
                            RSqrdFits{id,ridx}(rowid,colid) = rsqrd(1);
                        end
                        % hold on,plot(T(ix),y,'r')
                        
                        %Latency = mean - std
                        Latencies{id,ridx}(rowid,colid) = params(1)-params(2);
                        
                        %                         if savefigures && ~isnan(rs)
                        %                             %timepoint of max value
                        %                             h = line([cutoff cutoff],get(gca,'YLim'));  set(h,'Color',[0 0 0])
                        %                             % plot cumulative GaussFit
                        %                             hold on,plot(timevec(ix),y,'r')
                        %                             h = line([Lat Lat],get(gca,'YLim'));
                        %                             set(h,'Color',[0 1 0])
                        %                         end
                        
                    end
                end
                try
                    subplot(maxnrReactOpt,2,(ridx-1)*2+1); htmp = imagesc(Latencies{id,ridx},[0 350]);
                    colormap(autumn)
                    set(htmp,'Alphadata',~isnan(Latencies{id,ridx}))
                    hold on
                    plot(round(BrainModel{midx}.Model.AllX*rscale),round(BrainModel{midx}.Model.AllY*rscale),'.k','MarkerSize',5)
                    axis square
                    colorbar
                    if ridx == 1
                        title('Latencies')
                    end
                    ylabel([AllReactionOpt{ridx}])
                    
                    subplot(maxnrReactOpt,2,(ridx-1)*2+2); htmp = imagesc(RSqrdFits{id,ridx},[0 1]);
                    set(htmp,'Alphadata',~isnan(Latencies{id,ridx}))
                    
                    hold on
                    plot(round(BrainModel{midx}.Model.AllX*rscale),round(BrainModel{midx}.Model.AllY*rscale),'.k','MarkerSize',5)
                    
                    axis square
                    colorbar
                    if ridx == 1
                        title('RSqrd')
                    end
                    
                catch ME
                    disp(ME)
                    keyboard
                end
            end
            saveas(h,fullfile(StorePath,'Figures',['Baseline' num2str(baselinemethod) '_eqsample' num2str(takeequalsample)],[miceopt{midx} ' ' trialtypes{id} '_LEFT-Right_Latencies.bmp']))
            saveas(h,fullfile(StorePath,'Figures',['Baseline' num2str(baselinemethod) '_eqsample' num2str(takeequalsample)],[miceopt{midx} ' ' trialtypes{id} '_LEFT-Right_Latencies.fig']))
            saveas(h,fullfile(StorePath,'Figures',['Baseline' num2str(baselinemethod) '_eqsample' num2str(takeequalsample)],[miceopt{midx} ' ' trialtypes{id} '_LEFT-Right_Latencies.eps']))
            close(h)
            
        end
        
        LatencyAna{midx}.Latencies = Latencies;
        LatencyAna{midx}.RSqrdFits = RSqrdFits;
        
        %Right - left
        Latencies = cell(max(fgopt),maxnrReactOpt);
        RSqrdFits = Latencies;
        rscale = 1/4;
        for id = fgopt
            tmpIm = cellfun(@(X) imresize(X,rscale),NewdFF{id},'UniformOutput',0); %Rescale image
            h = figure('name',[trialtypes{id}]);
            for ridx = 1:maxnrReactOpt
                Latencies{id,ridx} = nan(size(tmpIm{1,ridx},1),size(tmpIm{1,ridx},2));
                RSqrdFits{id,ridx} = nan(size(tmpIm{1,ridx},1),size(tmpIm{1,ridx},2));
                for rowid=1:size(tmpIm{1,ridx},1)
                    for colid = 1:size(tmpIm{1,ridx},2)
                        ix = find(timevec>=0 & timevec<=500); %only visual
                        bix = find(timevec<0&timevec>=-300);
                        tmpdat = tmpIm{2,ridx}(rowid,colid,ix)-tmpIm{1,ridx}(rowid,colid,ix);
                        tmpdatbase = tmpIm{2,ridx}(rowid,colid,bix)-tmpIm{1,ridx}(rowid,colid,bix);
                        %                     tmpdat1 = tmpIm{1,ridx}(rowid,colid,ix);
                        %                     tmpdat2 = tmpIm{2,ridx}(rowid,colid,ix);
                        %                     R1 = squeeze(tmpdat1); R2 = squeeze(tmpdat2);
                        %                     Rsd1 = nanstd(R1); Rsd2 = nanstd(R2);
                        %                     Rm1 = nanmean(R1); Rm2 = nanmean(R2);
                        R = squeeze(tmpdat);
                        Rsd = nanstd(tmpdatbase);
                        Rm = nanmean(tmpdatbase);
                        if sum(~isnan(R))<length(ix) || nanmax(R)<Rm+2.5*Rsd % || nanmax(R2)<Rm2+Rsd2
                            continue
                        end
                        % figure;plot(T,R)
                        %                         figure,plot(timevec(ix),R); hold on
                        % SR = smooth(R,5);%plot(T(ix),SR(ix))
                        
                        % cutoff = T(ix(find(diff(SR(ix))<0,1,'first')));
                        [~,cutoff] = max(R);
                        R = R(1:cutoff);
                        cutoff = timevec(ix(cutoff));
                        % h = line([cutoff cutoff],get(gca,'YLim'));
                        % set(h,'Color',[0 0 0])
                        
                        ix = find(timevec>=0 & timevec<=cutoff);
                        
                        %Fit a cumaltive Gauss
                        [y,params] = cumgaussfit(timevec(ix),R',[150 20 0.03 0],0);
                        rsqrd = corrcoef(y,R').^2;
                        if length(rsqrd)>1
                            RSqrdFits{id,ridx}(rowid,colid) = rsqrd(2);
                        else
                            RSqrdFits{id,ridx}(rowid,colid) = rsqrd(1);
                        end
                        % hold on,plot(T(ix),y,'r')
                        
                        %Latency = mean - std
                        Latencies{id,ridx}(rowid,colid) = params(1)-params(2);
                        
                        %                         if savefigures && ~isnan(rs)
                        %                             %timepoint of max value
                        %                             h = line([cutoff cutoff],get(gca,'YLim'));  set(h,'Color',[0 0 0])
                        %                             % plot cumulative GaussFit
                        %                             hold on,plot(timevec(ix),y,'r')
                        %                             h = line([Lat Lat],get(gca,'YLim'));
                        %                             set(h,'Color',[0 1 0])
                        %                         end
                        
                    end
                end
                try
                    subplot(maxnrReactOpt,2,(ridx-1)*2+1); htmp = imagesc(Latencies{id,ridx},[0 350]);
                    colormap(autumn)
                    set(htmp,'Alphadata',~isnan(Latencies{id,ridx}))
                    hold on
                    plot(round(BrainModel{midx}.Model.AllX*rscale),round(BrainModel{midx}.Model.AllY*rscale),'.k','MarkerSize',5)
                    axis square
                    colorbar
                    if ridx == 1
                        title('Latencies')
                    end
                    ylabel([AllReactionOpt{ridx}])
                    
                    subplot(maxnrReactOpt,2,(ridx-1)*2+2); htmp = imagesc(RSqrdFits{id,ridx},[0 1]);
                    set(htmp,'Alphadata',~isnan(Latencies{id,ridx}))
                    
                    hold on
                    plot(round(BrainModel{midx}.Model.AllX*rscale),round(BrainModel{midx}.Model.AllY*rscale),'.k','MarkerSize',5)
                    
                    axis square
                    colorbar
                    if ridx == 1
                        title('RSqrd')
                    end
                    
                catch ME
                    disp(ME)
                    keyboard
                end
            end
            saveas(h,fullfile(StorePath,'Figures',['Baseline' num2str(baselinemethod) '_eqsample' num2str(takeequalsample)],[miceopt{midx} ' ' trialtypes{id} '_Right-Left_Latencies.bmp']))
            saveas(h,fullfile(StorePath,'Figures',['Baseline' num2str(baselinemethod) '_eqsample' num2str(takeequalsample)],[miceopt{midx} ' ' trialtypes{id} '_Right-Left_Latencies.fig']))
            saveas(h,fullfile(StorePath,'Figures',['Baseline' num2str(baselinemethod) '_eqsample' num2str(takeequalsample)],[miceopt{midx} ' ' trialtypes{id} '_Right-Left_Latencies.eps']))
            close(h)
            
        end
        
    end
    
    %%  T-score per reaction/side combination, for different timewindows
    try
        if ~strcmp(Stim2Check,'FGTask')
            BaseTW = [-250 0];
            VisTW = [50 120];
            DelayTWearly = [unique(LOG.Stimdur)+100 unique(LOG.Stimdur)+100+max(unique(LOG.currentdelay))/2];
            DelayTWlate = [unique(LOG.Stimdur)+100+max(unique(LOG.currentdelay))/2 unique(LOG.Stimdur)+max(unique(LOG.currentdelay))];
            twname = {'BaseTW','VisTW','DelayTWearly','DelayTWlate'};%'VTotal',,'DTotal'
        else
            BaseTW = [-250 0];
            VisInitTW = [0 120];
            VisTW = [0 500];
            FGModTW = [fgtw];
            twname = {'BaseTW','VisInitTW','VisTW','FGModTW'};%'VTotal',,'DTotal'
        end
        TScore = cell(length(fgopt),maxnrReactOpt,length(twname)-1);
        PScore = TScore;
        TTHRESHOLD = TScore;
        for fgid = fgopt
            FF = figure;
            rowid = 1;
            for ridx = 1:maxnrReactOpt
                try
                    tmp = nansum(cat(4,SUMSQALL{fgid}{:,ridx}),4);
                catch ME
                    disp(ME)
                    keyboard
                end
                if isempty(tmp)
                    continue
                end
                % Apply brainmask
                tmp(~repmat(brainmask,[1,1,size(tmp,3)])) = nan;
                
                %Average dFF
                tmp2 = nanmean(cat(4,NewdFF{fgid}{:,ridx}),4); %Average dFF all sessions, trials concatenated
                tmp2(~repmat(brainmask,[1,1,size(tmp2,3)])) = nan;
                
                %NRT
                ntmp = nansum([nrt{fgid}{:,ridx}]);
                
                for twid = 1:length(twname)-1
                    bsidx = timevec>=eval([twname{1} '(1)']) & timevec<eval([twname{1} '(2)']);
                    ttidx = timevec>=eval([twname{twid+1} '(1)']) & timevec<=eval([twname{twid+1} '(2)']);
                    
                    part1 = nanmean(tmp2(:,:,ttidx),3) - nanmean(tmp2(:,:,bsidx),3);
                    part2a = (1/(ntmp-1)).*nanmean(tmp(:,:,ttidx),3);
                    part2b = (ntmp/(ntmp-1)).*nanmean(tmp2(:,:,ttidx),3).^2;
                    part2 = part2a - part2b;
                    %If <0 make nan
                    part2(part2<0) =nan;
                    part2 = sqrt(part2);
                    TTHRESHOLD{fgid,ridx,twid} = ntmp-1;
                    TScore{fgid,ridx,twid} = (part1./part2)*sqrt(ntmp);
                    PScore{fgid,ridx,twid} = tcdf(TScore{fgid,ridx,twid},TTHRESHOLD{fgid,ridx,twid});
                    
                    subplot(maxnrReactOpt,length(twname)-1,(rowid-1)*(length(twname)-1)+twid)
                    h = imagesc(PScore{fgid,ridx,twid},[0 1]);
                    set(h,'alphadata',~isnan(PScore{fgid,ridx,twid}))
                    colormap(PSCOREMAP)
                    axis square
                    
                    if rowid == 1
                        title(twname{twid+1})
                    end
                    if twid ==1
                        ylabel([AllReactionOpt{ridx}])
                    end
                    
                end
                
                
                oldsize = get(gca,'Position');
                colorbar
                set(gca,'Position',oldsize)
                rowid = rowid+1;
            end
            saveas(FF,fullfile(StorePath,'Figures',['Baseline' num2str(baselinemethod) '_eqsample' num2str(takeequalsample)],[miceopt{midx} '_' trialtypes{fgid} '_PSCOREMAPS.bmp']))
            saveas(FF,fullfile(StorePath,'Figures',['Baseline' num2str(baselinemethod) '_eqsample' num2str(takeequalsample)],[miceopt{midx} '_' trialtypes{fgid} '_PSCOREMAPS.fig']))
            saveas(FF,fullfile(StorePath,'Figures',['Baseline' num2str(baselinemethod) '_eqsample' num2str(takeequalsample)],[miceopt{midx} '_' trialtypes{fgid} '_PSCOREMAPS.eps']))
        end
        
    catch ME
        disp(ME)
        keyboard
    end
    
    %% Create ROIs from T-score maps:
    if ~UserQuestions && EvokedActivityROI
        if strcmp(Stim2Check,'FGTask')
            fgid = find(strcmp(trialtypes,'GREY'));
            ridx = find(strcmp(AllReactionOpt,'Miss'));
        else
            if any(fgopt == 2)
                fgid = 2;
            else
                fgid = 1;
            end
            ridx = 1:length(AllReactionOpt);
        end
        
        rois = {};
        roicount = 1;
        for tid = 1:length(twname)-1
            FF = figure; subplot(1,2,1);
            signTWandClus = false(size(PScore{fgid,ridx(1),tid}));
            try
                ntmp = sum([nrt{fgid}{:,ridx}]);
                krit = (0.05./2); %2-sided
                lowertr = krit./ntmp;
                highertr = 1-(krit./ntmp);
                tmphigh = nanmax(cat(4,PScore{fgid,:,:}),[],4);
                tmplow = nanmin(cat(4,PScore{fgid,:,:}),[],4);
                signTWandClus(tmplow<lowertr|tmphigh>highertr) = 1;
            catch ME
                disp(ME)
                keyboard
            end
            h = imagesc(nanmean(signTWandClus,3));
            set(h,'alphadata',nanmean(signTWandClus,3)~=0)
            hold on
            hh = plot(BrainModel{midx}.Model.AllX,BrainModel{midx}.Model.AllY,'k.','MarkerSize',0.5);
            title(['Significant T-scores, ' twname{tid+1}])
            
            for tp =1:size(signTWandClus,3)
                
                tmp1 =  bwmorph(signTWandClus(:,:,tp),'tophat');
                tmp1 = imgaussfilt(single(imcomplement(tmp1)),5);
                tmp2 = bwmorph(signTWandClus(:,:,tp),'open',inf);
                tmp3 = bwmorph(tmp1==1&tmp2==1,'spur',inf);
                signTWandClus(:,:,tp) = tmp3;
                
            end
            subplot(1,2,2);
            h = imagesc(nanmean(signTWandClus,3));
            set(h,'alphadata',nanmean(signTWandClus,3)~=0)
            hold on
            hh = plot(BrainModel{midx}.Model.AllX,BrainModel{midx}.Model.AllY,'k.','MarkerSize',0.5);
            title('Removed noise')
            
            signTWandClus =nanmean(signTWandClus,3);
            %fit ROIS around areas, also defined by alan brain areas
            regio2take = find(~strcmp(BrainModel{midx}.Model.Rnames,''));
            for roiidx = 1:length(regio2take)
                
                Borders = BrainModel{midx}.Model.Boundaries{regio2take(roiidx)};
                count = 1;
                for roi2dx = 1:length(Borders)
                    mask = poly2mask(Borders{roi2dx}(:,1),Borders{roi2dx}(:,2),xpix,ypix);
                    %Shrink to not have border effects
                    mask = bwmorph(mask,'shrink',1);
                    singTWtmp = signTWandClus;
                    singTWtmp(~mask) = 0;
                    
                    singTWtmp = roicolor(singTWtmp,0.2,1);
                    singTWtmp = bwmorph(singTWtmp,'fill',inf);
                    boundartmps =  bwboundaries(singTWtmp,'noholes');
                    boundartmps(cellfun(@(X) size(X,1)<100,boundartmps))=[];
                    
                    if ~isempty(boundartmps)
                        for ii = 1:length(boundartmps)
                            roi.xi = boundartmps{ii}(:,2);
                            roi.yi = boundartmps{ii}(:,1);
                            rois{roicount} =roi;
                            hold on
                            plot(rois{roicount}.xi,rois{roicount}.yi,'r-','LineWidth',1.5)
                            ht(i) = text(mean(roi.xi(:)), mean(roi.yi(:)),num2str(roicount),'color',[0 0 0]);
                            
                            roicount = roicount+1;
                            
                        end
                    end
                end
            end
            saveas(FF,fullfile(StorePath,'Figures',['Baseline' num2str(baselinemethod) '_eqsample' num2str(takeequalsample)],[mouse 'Timewindow '  twname{tid+1} '_' trialtypes{fgid} '_AutomaticROImapEvokedActiv.bmp']))
            saveas(FF,fullfile(StorePath,'Figures',['Baseline' num2str(baselinemethod) '_eqsample' num2str(takeequalsample)],[mouse 'Timewindow '  twname{tid+1} '_' trialtypes{fgid} '_AutomaticROImapEvokedActiv.fig']))
        end
        
        %COMBINE OVERLAPPING ROIS
        FF = figure;
        imagesc(referenceimage)
        hold on
        newrois = {};
        if ~exist('newroiscount','var')
            newroiscount = 1;
        end
        hh = plot(BrainModel{midx}.Model.AllX,BrainModel{midx}.Model.AllY,'k.','MarkerSize',0.5);
        hold on
        regio2take = find(~strcmp(BrainModel{midx}.Model.Rnames,''));
        for roiidx = 1:length(regio2take)
            Borders = BrainModel{midx}.Model.Boundaries{regio2take(roiidx)};
            for roi2dx = 1:length(Borders)
                mask = poly2mask(Borders{roi2dx}(:,1),Borders{roi2dx}(:,2),xpix,ypix);
                roisidx =find(cell2mat(cellfun(@(X) sum(mask(poly2mask(X.xi,X.yi,xpix,ypix)==1))>0,rois,'UniformOutput',0))==1); %find which ROIs fall into this region and combine them
                if isempty(roisidx)
                    continue
                end
                tmpmask = zeros(xpix,ypix);
                for ii = 1:length(roisidx)
                    tmpmask(poly2mask(rois{roisidx(ii)}.xi,rois{roisidx(ii)}.yi,xpix,ypix))=1;
                end
                
                %Try to make it one ROI
                tmpmask = bwmorph(tmpmask,'fill',inf);
                tmpmask = bwmorph(tmpmask,'bridge',inf);
                boundartmps =  bwboundaries(tmpmask,'noholes');
                boundartmps(cellfun(@(X) size(X,1)<100,boundartmps))=[];
                
                if ~isempty(boundartmps)
                    for ii = 1:length(boundartmps)
                        roi.xi = boundartmps{ii}(:,2);
                        roi.yi = boundartmps{ii}(:,1);
                        newrois{newroiscount} =roi;
                        hold on
                        plot(newrois{newroiscount}.xi,newrois{newroiscount}.yi,'r-','LineWidth',1.5)
                        ht(i) = text(mean(roi.xi(:)), mean(roi.yi(:)),num2str(newroiscount),'color',[0 0 0]);
                        drawnow
                        newroiscount = newroiscount+1;
                        
                    end
                end
            end
        end
        saveas(FF,fullfile(StorePath,'Figures',['Baseline' num2str(baselinemethod) '_eqsample' num2str(takeequalsample)],[mouse '_AutomaticEvokedActivROImapAll.bmp']))
        saveas(FF,fullfile(StorePath,'Figures',['Baseline' num2str(baselinemethod) '_eqsample' num2str(takeequalsample)],[mouse '_AutomaticEvokedActivROImapAll.fig']))
        save(fullfile(StorePath,'Figures',['Baseline' num2str(baselinemethod) '_eqsample' num2str(takeequalsample)],[mouse 'EvokedActivROIs']),'rois')
        
        clear rois
        evokedactivrois = newrois;
        clear newrois
    end
    
    end
end


