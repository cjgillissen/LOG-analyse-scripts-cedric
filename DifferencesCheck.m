function DifferencesCheck(info,miceopt,StorePath,DataDirectory,Stim2Check,baselinemethod,originallim,trialtypes,plotlim,takeequalsample)
global UserQuestions
paths = info.paths;
logs = info.logs;
nrMouse = size(paths,1);

%% USER INPUT
nback = 20;
createvideosandfigurespermouse =0;
latencyana = 0;
contrastROIs = 1; %Contrast between conditions to determine ROI?
EvokedActivityROI = 0;
removeErrors = 0 %If you don't want to include Errors, make it 1
removeHits = 0 %if you don't want to include hits, make it 1
removeMiss = 0 %If you don't want to include misses, make it 1
fgtw = [80 250]; %FOR FG     


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
            VisTW = [80 150];
            DelayTWearly = [unique(LOG.Stimdur)+100 unique(LOG.Stimdur)+100+max(unique(LOG.currentdelay))/2];
            DelayTWlate = [unique(LOG.Stimdur)+100+max(unique(LOG.currentdelay))/2 unique(LOG.Stimdur)+max(unique(LOG.currentdelay))];
            twname = {'BaseTW','VisTW','DelayTWearly','DelayTWlate'};%'VTotal',,'DTotal'
        else
            BaseTW = [-250 0];
            VisInitTW = [0 120];
            VisTW = [80 150];
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




    %% calculate (T-score of) difference between Left & Right stimulus
    %      calculate s per condition
    clear ntmp
    VariancePerCond = cell(2,maxnrReactOpt,maxnrSideOpt);
    TScore = cell(maxnrReactOpt,2); %reaction, fg/grey
    PScore = TScore;
    for fgid=fgopt
        for ridx = 1:maxnrReactOpt
            for sdidx =1:length(AllSideOpt)
                try
                    tmp = SUMSQALL{fgid}{sdidx,ridx};
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
                tmp2 = NewdFF{fgid}{sdidx,ridx}; %Average dFF all sessions, trials concatenated
                tmp2(~repmat(brainmask,[1,1,size(tmp2,3)])) = nan;
                
                %NRT
                try
                    ntmp{sdidx} = nrt{fgid}{sdidx,ridx};
                catch ME
                    disp(ME)
                    keyboard
                end
                aa= ((1/(sum(ntmp{sdidx})-1)).*tmp);
                bb = ((sum(ntmp{sdidx})./(sum(ntmp{sdidx})-1)).*(tmp2.^2));
                try
                    SSqr{sdidx} =aa-bb;
                    VariancePerCond{fgid,ridx,sdidx} = SSqr{sdidx};
                catch ME
                    disp(ME)
                    keyboard
                end
            end
            
            TTHRESHOLD{ridx,fgid} = sum(cell2mat(ntmp))-2;
            PooledSSqr = squeeze(((ntmp{1} - 1).*SSqr{1}+(ntmp{2}-1).*SSqr{2}))./(sum(cell2mat(ntmp))-2);
            PooledSSqr(PooledSSqr<0) = nan; % Just delete those values!
            TScore{ridx,fgid} = squeeze(NewdFF{fgid}{1,ridx}-NewdFF{fgid}{2,ridx})./(sqrt(PooledSSqr).*(sqrt(1/ntmp{1}+1/ntmp{2})));
            PScore{ridx,fgid} = tcdf(TScore{ridx,fgid},TTHRESHOLD{ridx,fgid});        
        end
    end
    
    %Save in struct
    LvsR{midx}.TScore = TScore;
    LvsR{midx}.PScore = PScore;
    AllMicedFF{midx} = NewdFF;
    AllMiceVariance{midx} = VariancePerCond;
    AllMicenrt{midx} = nrt;
    
    %% Make video's: L,R,L-R HITS, T-score
    if createvideosandfigurespermouse
        for ridx = 1:length(AllReactionOpt)
            myVideo = VideoWriter(fullfile(StorePath,'Figures',['Baseline' num2str(baselinemethod) '_eqsample' num2str(takeequalsample)],[miceopt{midx} AllReactionOpt{ridx}]));
            myVideo.FrameRate = 5;
            open(myVideo)
            
            if length(fgopt) == 1
                tmp = cat(4,NewdFF{fgopt}{:,ridx});%,cat(4,NewdFF{2}{:,ridx}));
            else
                tmp = cat(4,cat(4,NewdFF{1}{:,ridx}),cat(4,NewdFF{2}{:,ridx}));
            end
            tmp = tmp(:,:,timevec>=plotlim(1)&timevec<=plotlim(2),:);
            quantvaldff = quantile(abs(tmp(:)),0.98);
            
            tmp = cat(4,TScore{ridx,:});
            tmp = tmp(:,:,timevec>plotlim(1)&timevec<plotlim(2),:);
            quantvalt = quantile(abs(tmp(:)),0.98);
            for tp = 1:length(timevec)
                if timevec(tp)<plotlim(1) || timevec(tp)>plotlim(2)
                    continue
                end
                fh = figure('units','normalized','outerposition',[0 0 1 1]); %Make frame
                countwind = 1;
                for fgid = fgopt
                    %Both sides
                    for sdidx = find(ismember(AllSideOpt,{'left','right'}))
                        tmp = NewdFF{fgid}{sdidx,ridx}(:,:,tp);
                        if sum(~isnan(tmp(:))) == 0
                            countwind = countwind+1;
                            continue
                        end
                        subplot(2,sum(ismember(AllSideOpt,{'left','right'}))+1,countwind)
                        h = imagesc(smooth2a(tmp,5),[-quantvaldff quantvaldff]);
                        axis square
                        set(h,'Alphadata',brainmask&~isnan(tmp))
                        colormap(ActSupColorMap)
                        
                        hold on
                        plot(BrainModel{midx}.Model.AllX,BrainModel{midx}.Model.AllY,'.y','MarkerSize',0.1)
                        
                        if sdidx == 1
                            ylabel(trialtypes{fgid})
                        elseif sdidx==sum(ismember(AllSideOpt,{'left','right'}))
                            oldsize = get(gca,'Position');
                            colorbar
                            set(gca,'Position',oldsize)
                        end
                        if fgid == 2
                            xlabel(AllSideOpt{sdidx})
                        end
                        countwind = countwind+1;
                    end
                    
                    %Left-Right t-score
                    tmp2 = squeeze(TScore{ridx,fgid}(:,:,tp));
                    if sum(~isnan(tmp2(:))) ==0
                        countwind = countwind+1;
                        continue
                    end
                    subplot(2,sum(ismember(AllSideOpt,{'left','right'}))+1,countwind)
                    try
                        
                        tmp2(PScore{ridx,fgid}(:,:,tp)>0.025&PScore{ridx,fgid}(:,:,tp)<0.9750) = 0;
                        
                        h = imagesc(smooth2a(tmp2,5),[-quantvalt quantvalt]);
                    catch ME
                        disp(ME)
                        keyboard
                    end
                    axis square
                    set(h,'Alphadata',brainmask&~isnan(tmp))
                    colormap(ActSupColorMap)
                    oldsize = get(gca,'Position');
                    colorbar
                    set(gca,'Position',oldsize)
                    
                    hold on
                    plot(BrainModel{midx}.Model.AllX,BrainModel{midx}.Model.AllY,'.y','MarkerSize',0.1)
                    
                    if fgid == fgopt(end)
                        xlabel([AllSideOpt{1} '-' AllSideOpt{2} ' significant t-scores'])
                    end
                    countwind = countwind+1;
                end
                
                suplabel(['t = ' num2str(timevec(tp))],'t');
                
                % Save
                try
                frame = getframe(fh);
                writeVideo(myVideo,frame)
                close(fh)
                catch ME
                   disp(ME)
                    keyboard
                end
                
            end
            close(myVideo)
        end
        
        %% Time-Course Plots (left -- Right -- Left - Right) both FG and Grey
        regio2take = find(~strcmp(BrainModel{midx}.Model.Rnames,''));
        for roiidx = 1:length(regio2take)
            FF = figure('name',[miceopt{midx},BrainModel{midx}.Model.Rnames{regio2take(roiidx)} '_TimeCourses']);
            legendname = {};
            
            Borders = BrainModel{midx}.Model.Boundaries{regio2take(roiidx)};
            mask4plot = zeros(size(referenceimage,1),size(referenceimage,2));
            count = 1;
            for roi2dx = 1:length(Borders)
                mask = poly2mask(Borders{roi2dx}(:,1),Borders{roi2dx}(:,2),xpix,ypix);
                %Shrink to not have border effects
                mask = bwmorph(mask,'shrink',1);
                mask4plot(mask) = count;
                hlegend = [];
                for fgid = fgopt
                    subplot(length(Borders)+1,2+1,(fgid-1)*3+roi2dx)
                    clear h
                    counth = 1;
                    for ridx = 1:length(AllReactionOpt)
                        for stidx = find(ismember(AllSideOpt,{'left','right'}))
                            legendname{(fgid-1)*(length(AllReactionOpt)*length(AllSideOpt))+(ridx-1)*length(AllSideOpt)+stidx} =[trialtypes{fgid} ', ' AllReactionOpt{ridx} ' ' AllSideOpt{stidx} ', n = ' num2str(nrt{fgid}{stidx,ridx})];
                            
                            tmp = squeeze(NewdFF{fgid}{stidx,ridx});
                            tmp(~repmat(mask,[1,1,size(tmp,3)]))=nan;
                            tmp = smooth(nanmean(reshape(tmp,[xpix*ypix,size(tmp,3)]),1),3);
                            
                            if stidx == 1
                                colid = 1;
                            elseif stidx == length(AllSideOpt)
                                colid = 5;
                            else
                                colid = 3;
                            end
                            tmperror = VariancePerCond{fgid,ridx,stidx};
                            tmperror(~repmat(mask,[1,1,size(tmperror,3)]))=nan;
                            tmperror = smooth(nanmean(reshape(tmperror,[xpix*ypix,size(tmperror,3)]),1),3);
                            
                            try
                            H=shadedErrorBar(timevec,tmp,tmperror,{'color',LineMap(:,colid,ridx),'LineWidth',2},1);
                            catch ME
                                disp(ME)
                                keyboard
                            end
                            h(counth) = H.mainLine;
                            %                         h(counth) = plot(timevec,tmp,'color',LineMap(:,colid,ridx),'LineWidth',2);
                            hold on
                            counth = counth+1;
                            
                            
                        end
                    end
                    yvals = get(gca,'ylim');
                    if strcmp(Stim2Check,'DelayedOriTuningSound')
                        try
                            hpatch(1) = patch([0 0 min(unique(LOG.Stimdur)) min(unique(LOG.Stimdur))],[yvals(1) yvals(2) yvals(2) yvals(1)],repmat(1,4,1),'FaceColor',[0 1 0],'FaceAlpha',.2);
                            hpatch(2) = patch([min(unique(LOG.Stimdur)) min(unique(LOG.Stimdur)) min(unique(LOG.Stimdur))+str2num(trialtypes{fgid}) min(unique(LOG.Stimdur))+str2num(trialtypes{fgid})], [yvals(1) yvals(2) yvals(2) yvals(1)],repmat(1,4,1),'FaceColor',[1 0 0],'FaceAlpha',.2);
                        catch ME
                            disp(ME)
                            keyboard
                        end
                    end
                    hlegend = [hlegend h];
                    ylabel(trialtypes{fgid})
                    xlim([plotlim(1) plotlim(2)])
                    title(num2str(count))
                end
                if strcmp(Stim2Check,'DelayedOriTuningSound')
                    hlegend = [hlegend hpatch];
                end
                count = count+1;
            end
            %Plot brain area taken
            subplot(length(Borders)+1,2+1,roi2dx+1)
            htmp = imagesc(mask4plot);
            set(htmp,'AlphaData',~(mask4plot==0))
            
            colorbar
            hold on
            hh = plot(BrainModel{midx}.Model.AllX,BrainModel{midx}.Model.AllY,'k.','MarkerSize',0.2);
            
            %legend
            subtmp = subplot(length(Borders)+1,2+1,length(Borders)+1+roi2dx+1);
            posleg = get(subtmp,'position');
            delete(subtmp)
            
            legend(hlegend,{legendname{~cellfun(@isempty,legendname)}},'Position',posleg)
            
            saveas(FF,fullfile(StorePath,'Figures',['Baseline' num2str(baselinemethod) '_eqsample' num2str(takeequalsample)],[miceopt{midx} BrainModel{midx}.Model.Rnames{regio2take(roiidx)} '_TimeCoursesLvsR.bmp']))
            saveas(FF,fullfile(StorePath,'Figures',['Baseline' num2str(baselinemethod) '_eqsample' num2str(takeequalsample)],[miceopt{midx} BrainModel{midx}.Model.Rnames{regio2take(roiidx)} '_TimeCoursesLvsR.fig']))
            
            close(FF)
        end
    end
    %% T_SCORE Method for ROIS
    if ~UserQuestions && contrastROIs
        if strcmp(Stim2Check,'FGTask')
            fgid = find(strcmp(trialtypes,'GREY'));
            ridx = find(strcmp(AllReactionOpt,'Miss'));
            TW = {fgtw};
        else
            if any(fgopt == 2)
                fgid = 2;
            else
                fgid = 1;
            end
            ridx = 2;
            TW = {[80 500],[800 1500]};
        end
        
        rois = {};
        roicount = 1;
        for tid = 1:length(TW)
            tw = TW{tid};
            FF = figure; subplot(1,2,1);
            signTWandClus = false(size(PScore{ridx,fgid}));
            try
                ntmp = sum([nrt{fgid}{:,ridx}]);
                krit = (0.05./2); %2-sided
                lowertr = krit./ntmp;
                highertr = 1-(krit./ntmp);
                signTWandClus(PScore{ridx,fgid}<lowertr|PScore{ridx,fgid}>highertr) = 1;
                signTWandClus = (signTWandClus(:,:,timevec>=tw(1)&timevec<=tw(2)));
            catch ME
                disp(ME)
                keyboard
            end
            h = imagesc(nanmean(signTWandClus,3));
            set(h,'alphadata',nanmean(signTWandClus,3)~=0)
            hold on
            hh = plot(BrainModel{midx}.Model.AllX,BrainModel{midx}.Model.AllY,'k.','MarkerSize',0.5);
            title(['Significant T-scores left vs right, t = ' num2str(tw(1)) '-' num2str(tw(2))])
            
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
            saveas(FF,fullfile(StorePath,'Figures',['Baseline' num2str(baselinemethod) '_eqsample' num2str(takeequalsample)],[mouse 'Timewindow ' num2str(tw(1)) '-' num2str(tw(2)) '_' trialtypes{fgid} '_AutomaticROImapLeftvsRight.bmp']))
            saveas(FF,fullfile(StorePath,'Figures',['Baseline' num2str(baselinemethod) '_eqsample' num2str(takeequalsample)],[mouse 'Timewindow ' num2str(tw(1)) '-' num2str(tw(2)) '_' trialtypes{fgid} '_AutomaticROImapLeftvsRight.fig']))
            save(fullfile(StorePath,miceopt{midx},'Averages','ROIsLvsR'),'rois')
        end
    else
        %% Make ROIs (and/or add)
        if exist(fullfile(StorePath,miceopt{midx},'Averages','ROIs.mat'))
            load(fullfile(StorePath,miceopt{midx},'Averages','ROIs.mat'))
        else
            rois = {};
        end
        FG_TW = [fgtw];
        TWidx = find(timevec>=FG_TW(1) & timevec<=FG_TW(2));
        
        % Select extra ROIS
        tmp = squeeze(nanmean(nanmean(cat(3,NewdFF{1}{stidx,ridx}(:,:,TWidx)))));
        quantval = quantile(tmp(:),0.99);
        
        figure; hh = imagesc(tmp,[-quantval quantval]);
        set(hh,'Alphadata',~isnan(tmp))
        hold on
        plot(BrainModel{midx}.Model.AllX,BrainModel{midx}.Model.AllY,'k.')
        done = 0;
        while ~done
            %toc
            a = length(rois);
            col = jet(a+1);
            if a > 0
                for i = 1:a
                    roi = rois{i};
                    hold on
                    h(i) = plot(roi.xi, roi.yi, 'color', col(i,:), 'linewidth', 2);
                    hold off
                    ht(i) = text(mean(roi.xi(:)), mean(roi.yi(:)),num2str(i),'color',col(i,:));
                end
            end
            
            button = questdlg('Do you want more or less ROIs?','ROI selection', 'More','Less','Done','More');
            if strcmp(button, 'Done')
                done = 1;
            elseif strcmp(button, 'More')
                [roi.log, roi.xi, roi.yi] = roipoly;
                rois{a+1} = roi;
                
                try
                    delete(h(:))
                    delete(ht(:))
                catch
                end
            elseif strcmp(button,'Less')
                button = inputdlg('Enter ROIs to delete (separate with comma (,)):');
                data = str2num(button{:});
                rois(data) = [];
                
                delete(h(:))
                delete(ht(:))
            end
        end
        if ~exist(fullfile(StorePath,miceopt{midx},'Averages'),'dir')
            mkdir(fullfile(StorePath,miceopt{midx},'Averages'))
        end
        save(fullfile(StorePath,miceopt{midx},'Averages','ROIs'),'rois')
        export_fig(fullfile(StorePath,miceopt{midx},'Averages','RoiMap.png'))
    end
    clear ntmp   
    
    %% Same analysis for Hits-Errors
    % calculate (T-score of) difference between Hits & Error trials
    %      calculate s per condition
    %Only works when there is more than one type of reaction
    if prod(size(SUMSQALL{1})) >= 4 || prod(size(SUMSQALL{2})) >= 4
        
        TScore = cell(3,2); %reaction, fg/grey
        PScore = TScore;
        for fgid=fgopt
            for ridx = 1:maxnrReactOpt
                try
                    tmp = sum(cat(4,SUMSQALL{fgid}{:,ridx}),4);
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
                try
                    ntmp{ridx} = sum(cat(2,nrt{fgid}{:,ridx}));
                catch ME
                    disp(ME)
                    keyboard
                end
                aa= ((1/(sum(ntmp{ridx})-1)).*tmp);
                bb = ((sum(ntmp{ridx})./(sum(ntmp{ridx})-1)).*(tmp2.^2));
                try
                    SSqr{ridx} =aa-bb;
                catch ME
                    disp(ME)
                    keyboard
                end
            end
            
            for ridx = 1:maxnrReactOpt
                TTHRESHOLD{ridx,fgid} = sum(cell2mat(ntmp))-2;
                ridx2 = ridx+1;
                if ridx2 > maxnrReactOpt
                    ridx2 = 1;
                end
                if strcmp(AllReactionOpt{ridx2},'Hit')
                    ridxtmp = ridx2;
                    ridx2 = ridx;
                else
                    ridxtmp = ridx;
                end
                PooledSSqr = squeeze(((ntmp{ridxtmp} - 1).*SSqr{ridxtmp}+(ntmp{ridx2}-1).*SSqr{ridx2}))./(sum(cell2mat(ntmp))-2);
                PooledSSqr(PooledSSqr<0) = nan; % Just delete those values!
                tmp1 = nanmean(cat(4,NewdFF{fgid}{:,ridxtmp}),4); %Average dFF all sessions, trials concatenated
                tmp1(~repmat(brainmask,[1,1,size(tmp1,3)])) = nan;
                
                tmp2 = nanmean(cat(4,NewdFF{fgid}{:,ridx2}),4); %Average dFF all sessions, trials concatenated
                tmp2(~repmat(brainmask,[1,1,size(tmp2,3)])) = nan;
                
                TScore{ridx,fgid} = squeeze(tmp1-tmp2)./(sqrt(PooledSSqr).*(sqrt(1/ntmp{ridxtmp}+1/ntmp{ridx2})));
                PScore{ridx,fgid} = tcdf(TScore{ridx,fgid},TTHRESHOLD{ridx,fgid});
            end
        end
        %Save in struct
        EvsHvsM{midx}.TScore = TScore;
        EvsHvsM{midx}.PScore = PScore;
    else
        EvsHvsM{midx} = 'Not more than 1 reactiontypes available';
    end
    if createvideosandfigurespermouse
        %% Video's
        if prod(size(SUMSQALL{1})) >= 4 || prod(size(SUMSQALL{2})) >= 4
            
            for fgid = fgopt
                
                myVideo = VideoWriter(fullfile(StorePath,'Figures',['Baseline' num2str(baselinemethod) '_eqsample' num2str(takeequalsample)],[miceopt{midx} 'trialtype' trialtypes{fgid}]));
                myVideo.FrameRate = 5;
                open(myVideo)
                
                if length(fgopt) == 1
                    tmp = cat(4,NewdFF{fgopt}{:,ridx});%,cat(4,NewdFF{2}{:,ridx}));
                else
                    tmp = cat(4,cat(4,NewdFF{1}{:,ridx}),cat(4,NewdFF{2}{:,ridx}));
                end
                tmp = tmp(:,:,timevec>=plotlim(1)&timevec<=plotlim(2),:);
                quantvaldff = quantile(abs(tmp(:)),0.98);
                
                tmp = cat(4,TScore{:,:});
                tmp = tmp(:,:,timevec>plotlim(1)&timevec<plotlim(2),:);
                quantvalt = quantile(abs(tmp(:)),0.98);
                
                for tp = 1:length(timevec)
                    if timevec(tp)<plotlim(1) || timevec(tp)>plotlim(2)
                        continue
                    end
                    fh = figure('units','normalized','outerposition',[0 0 1 1]); %Make frame
                    countwind = 1;
                    %Both sides
                    for ridx = 1:maxnrReactOpt
                        ridx2 = ridx+1;
                        if ridx2 > maxnrReactOpt
                            ridx2 = 1;
                        end
                        if strcmp(AllReactionOpt{ridx2},'Hit')
                            ridxtmp = ridx2;
                            ridx2 = ridx;
                        else
                            ridxtmp = ridx;
                        end
                        avdFF = nanmean(cat(4,NewdFF{fgid}{:,ridxtmp}),4);
                        tmp1 = avdFF(:,:,tp);
                        
                        avdFF = nanmean(cat(4,NewdFF{fgid}{:,ridx2}),4);
                        tmp2 = avdFF(:,:,tp);
                        
                        if sum(~isnan(tmp(:))) == 0
                            countwind = countwind+1;
                            continue
                        end
                        subplot(3,3,countwind)
                        h = imagesc(smooth2a(tmp1,5),[-quantvaldff quantvaldff]);
                        axis square
                        set(h,'Alphadata',brainmask&~isnan(tmp1))
                        colormap(ActSupColorMap)
                        
                        hold on
                        plot(BrainModel{midx}.Model.AllX,BrainModel{midx}.Model.AllY,'.y','MarkerSize',0.1)
                        
                        
                        xlabel(AllReactionOpt{ridxtmp})
                        countwind = countwind+1;
                        
                        subplot(3,3,countwind)
                        h = imagesc(smooth2a(tmp2,5),[-quantvaldff quantvaldff]);
                        axis square
                        set(h,'Alphadata',brainmask&~isnan(tmp2))
                        colormap(ActSupColorMap)
                        
                        
                        hold on
                        plot(BrainModel{midx}.Model.AllX,BrainModel{midx}.Model.AllY,'.y','MarkerSize',0.1)
                        xlabel(AllReactionOpt{ridx2})
                        
                        oldsize = get(gca,'Position');
                        colorbar
                        set(gca,'Position',oldsize)
                        
                        countwind = countwind+1;
                        
                        %Left-Right t-score
                        tmp3 = squeeze(TScore{ridx,fgid}(:,:,tp));
                        if sum(~isnan(tmp3(:))) ==0
                            countwind = countwind+1;
                            continue
                        end
                        subplot(3,3,countwind)
                        try
                            
                            tmp3(PScore{ridx,fgid}(:,:,tp)>0.025&PScore{ridx,fgid}(:,:,tp)<0.9750) = 0;
                            
                            h = imagesc(smooth2a(tmp3,5),[-quantvalt quantvalt]);
                        catch ME
                            disp(ME)
                            keyboard
                        end
                        axis square
                        set(h,'Alphadata',brainmask&~isnan(tmp3))
                        colormap(ActSupColorMap)
                        oldsize = get(gca,'Position');
                        colorbar
                        set(gca,'Position',oldsize)
                        
                        hold on
                        plot(BrainModel{midx}.Model.AllX,BrainModel{midx}.Model.AllY,'.y','MarkerSize',0.1)
                        
                        xlabel([AllReactionOpt{ridxtmp} '-' AllReactionOpt{ridx2} ' significant t-scores'])
                        countwind = countwind+1;
                    end
                    
                    suplabel(['t = ' num2str(timevec(tp))],'t');
                    
                    % Save
                    frame = getframe(fh);
                    writeVideo(myVideo,frame)
                    close(fh)
                    
                end
                close(myVideo)
                
            end
        end
    end
    %% ADD ROIS
    if ~UserQuestions && contrastROIs
        if strcmp(Stim2Check,'FGTask')
            fgid = find(strcmp(trialtypes,'GREY'));
            ridx = find(strcmp(AllReactionOpt,'Miss'));
            TW = {fgtw};
        else
            if any(fgopt) == 2
                fgid = 2;
            else
                fgid = 1;
            end
            TW = {[80 500],[800 1500]};
            ridx =1:length(AllReactionOpt);
            
        end
        
        for tid = 1:length(TW)
            tw = TW{tid};
            FF = figure; subplot(1,2,1);
            tmphigh = nanmax(cat(4,PScore{ridx,fgid}),[],4);
            tmplow = nanmin(cat(4,PScore{ridx,fgid}),[],4);
            signTWandClus = false(size(tmphigh));
            
            ntmp = sum([nrt{fgid}{:,ridx}]);
            krit = (0.05./2); %2-sided
            lowertr = krit./ntmp;
            highertr = 1-(krit./ntmp);
            signTWandClus(tmplow<lowertr|tmphigh>highertr) = 1;
            signTWandClus = (signTWandClus(:,:,timevec>=tw(1)&timevec<=tw(2)));

            
            h = imagesc(nanmean(signTWandClus,3));
            set(h,'alphadata',nanmean(signTWandClus,3)~=0)
            hold on
            hh = plot(BrainModel{midx}.Model.AllX,BrainModel{midx}.Model.AllY,'k.','MarkerSize',0.5);
            title(['Significant T-scores Reaction Differences, t = ' num2str(tw(1)) '-' num2str(tw(2))])
            
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
            saveas(FF,fullfile(StorePath,'Figures',['Baseline' num2str(baselinemethod) '_eqsample' num2str(takeequalsample)],[mouse 'Timewindow ' num2str(tw(1)) '-' num2str(tw(2)) '_' trialtypes{fgid} '_AutomaticROImapResponse.bmp']))
            saveas(FF,fullfile(StorePath,'Figures',['Baseline' num2str(baselinemethod) '_eqsample' num2str(takeequalsample)],[mouse 'Timewindow ' num2str(tw(1)) '-' num2str(tw(2)) '_' trialtypes{fgid} '_AutomaticROImapResponse.fig']))
            
        end
        
    end
    
    %% Combine ROIS that overlap:
    if contrastROIs
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
        saveas(FF,fullfile(StorePath,'Figures',['Baseline' num2str(baselinemethod) '_eqsample' num2str(takeequalsample)],[mouse '_AutomaticContrROImapAll.bmp']))
        saveas(FF,fullfile(StorePath,'Figures',['Baseline' num2str(baselinemethod) '_eqsample' num2str(takeequalsample)],[mouse '_AutomaticContrROImapAll.fig']))
        save(fullfile(StorePath,'Figures',['Baseline' num2str(baselinemethod) '_eqsample' num2str(takeequalsample)],[mouse 'ContrRois']),'newrois')
        
        clear rois
        ContrRois = newrois;
        clear newrois
    end
    
    
    if strcmp(Stim2Check,'FGTask')
        fgid = find(strcmp(trialtypes,'FG'));
        TW = {fgtw};
    else
        if any(fgopt == 2)
            fgid = 2;
        else
            fgid = 1;
        end
        TW = {[80 500],[800 1500]};
    end
    FF=figure;
    
    for ridx = 1:size(PScore,1)
        subplot(1,3,ridx)
        rois = {};
        roicount = 1;
        for tid = 1:length(TW)
            tw = TW{tid};
            try
            signTWandClus = false(size(PScore{ridx,fgid}));
            catch ME
                disp(ME)
                keyboard
            end
            try
                ntmp = sum([nrt{fgid}{:,ridx}]);
                krit = (0.05./2); %2-sided
                lowertr = krit./ntmp;
                highertr = 1-(krit./ntmp);
                signTWandClus(PScore{ridx,fgid}<lowertr|PScore{ridx,fgid}>highertr) = 1;
                signTWandClus = (signTWandClus(:,:,timevec>=tw(1)&timevec<=tw(2)));
            catch ME
                disp(ME)
                keyboard
            end
            h = imagesc(nanmean(signTWandClus,3));
            set(h,'alphadata',nanmean(signTWandClus,3)~=0)
            hold on
            hh = plot(BrainModel{midx}.Model.AllX,BrainModel{midx}.Model.AllY,'k.','MarkerSize',0.5);
            title(['Significant T-scores left vs right, t = ' num2str(tw(1)) '-' num2str(tw(2))])
            
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
                            size(rois)
                            hold on
                            plot(rois{roicount}.xi,rois{roicount}.yi,'r-','LineWidth',1.5)
                            ht(i) = text(mean(roi.xi(:)), mean(roi.yi(:)),num2str(roicount),'color',[0 0 0]);
                            save(fullfile(StorePath,'Figures',['Baseline' num2str(baselinemethod) '_eqsample' num2str(takeequalsample)],[mouse 'ContrROIs' num2str(roicount)]),'rois')
                            roicount = roicount+1;
                            
                        end
                    end
                end
            end
        end
        roissave = rois;
        title([AllReactionOpt{ridx}])
    end
    
    saveas(FF,fullfile(StorePath,'Figures',['Baseline' num2str(baselinemethod) '_eqsample' num2str(takeequalsample)],[mouse 'Timewindow ' num2str(tw(1)) '-' num2str(tw(2)) '_LeftvsRightFG.bmp']))
    saveas(FF,fullfile(StorePath,'Figures',['Baseline' num2str(baselinemethod) '_eqsample' num2str(takeequalsample)],[mouse 'Timewindow ' num2str(tw(1)) '-' num2str(tw(2)) '_LeftvsRightFG.fig']))
    save(fullfile(StorePath,'Figures',['Baseline' num2str(baselinemethod) '_eqsample' num2str(takeequalsample)],[mouse 'ContrROIs']),'rois')
    save(fullfile(StorePath,'Figures',['Baseline' num2str(baselinemethod) '_eqsample' num2str(takeequalsample)],[mouse 'ContrROIs']),'roissave')

end


   