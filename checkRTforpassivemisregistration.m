%% LOAD data

clear all
close all

addpath(genpath('C:\Users\gillissen\Documents\GitHub\Mouse'))
cd('C:\Users\gillissen\Desktop\InternshipCédric\TrainingLogs')
%% Select Mouse
TrainingList = dir(fullfile(cd,'*.mat'));
TrainingList_Parts = cellfun(@(X) strsplit(X,'_'), {TrainingList.name},'UniformOutput',false);
AllMice = unique(cellfun(@(X) X(1), TrainingList_Parts))';
Mice=AllMice(listdlg('ListString',AllMice));
RemoveMice = setdiff(AllMice,Mice);
for i = 1:length(RemoveMice)
    TrainingList(find(~cellfun(@isempty,strfind({TrainingList.name},RemoveMice{i}))))=[];
end

%% Select Days
TrainingList_Parts = cellfun(@(X) strsplit(X,'_'), {TrainingList.name},'UniformOutput',false);
AllDays = unique(cellfun(@(X) X(2), TrainingList_Parts));
Days=AllDays(listdlg('ListString',AllDays));
RemoveDays = setdiff(AllDays,Days);
for i = 1:length(RemoveDays)
    TrainingList(find(~cellfun(@isempty,strfind({TrainingList.name},RemoveDays{i}))))=[];
end
TrainingList_Parts = cellfun(@(X) strsplit(X,'_'), {TrainingList.name},'UniformOutput',false);


%% check first reaction times for correctReaction
%if < 500 always then something wrong in passive condition arduino <->
%matlab
RTleftCell = cell(length(Mice),length(Days));
RTrightCell = cell(length(Mice),length(Days));
MINRTLEFT = cell(size(length(Mice),length(Days))); % the smallest RT left
MINRTRIGHT = cell(size(length(Mice),length(Days))); % the smallest RT right
nrPassives = nan(size(length(Mice),length(Days))); % nr of passives given ... Might want to control for end of sessions
nrPassWrong = nan(size(length(Mice),length(Days))); % nr of passives that are wrongly registered, check index number with X
percPassives = nan(size(length(Mice),length(Days))); % percentage of passives given in whole session, so nr passives / nr trials
percWrong = nan(size(length(Mice),length(Days))); % percentage wrongly registered troughout whole session for all conditions
percWrongisPass = nan(size(length(Mice),length(Days))); % of percWrong, what is the percentage passives...
nrincorrect = nan(size(length(Mice),length(Days))); % nr of wrongly registered trials per session per animal


for Mouseidx = 1:length(Mice) % Mouse loop
    
    for Dayidx = 1: length(Days) % Session loop
        Fileidx = find(~cellfun(@isempty,cellfun(@(X) strfind(X,[Mice{Mouseidx} '_' Days{Dayidx}]),{TrainingList.name},'UniformOutput',false)));
        
           if size(Fileidx)==1 % I don't want to bother concatenating sessions sorry
            load(TrainingList(Fileidx).name); % find in traininglist the string element that has the name of the file you want to load
            idxincorrect = []; % initiliaze the index numbers of trials that have an incongruent correctReaction - Reaction
     
    if exist('LOG.correctReaction') % If LOG.correctReaction doesn't exist for file than create     
    
    X = strcmp(LOG.correctReaction,LOG.Reaction);
    idxincorrect = find(~X);
    nrincorrect(Mouseidx,Dayidx) = length(idxincorrect);
    percWrong(Mouseidx,Dayidx) = CheckReactions(TrainingList(Fileidx).name);
    idxpassive = find(LOG.Gavepassive);
    percWrongisPass(Mouseidx,Dayidx) = length(idxincorrect)/(length(find(LOG.Gavepassive(idxincorrect)))+length(idxincorrect));
    nrPassives(Mouseidx,Dayidx) = length(find(LOG.Gavepassive));
    nrPassWrong(Mouseidx,Dayidx) = length(find(LOG.Gavepassive(idxincorrect))); 
    percPassives(Mouseidx,Dayidx) = nrPassives(Mouseidx,Dayidx)/length(LOG.Reaction);

    if length(idxincorrect)>=1 % otherwise you don't know whether 0 because of zero or nan
                percWrongisPass(Mouseidx,Dayidx) = length(idxincorrect)/(length(find(LOG.Gavepassive(idxincorrect)))+length(idxincorrect));
            else
                percWrongisPass(Mouseidx,Dayidx) = nan;
    end

    RTLEFT = nan(size(idxincorrect));
    RTRIGHT = nan(size(idxincorrect));

    for i = 1:length(idxincorrect)

                if ~isempty(LOG.RTleftVec{i})
                    RTLEFT(i) = min(LOG.RTleftVec{i});
                end
                if ~isempty(LOG.RTrightVec{i})
                    RTRIGHT(i) = min(LOG.RTrightVec{i});
                end

     end
            
            MINRTLEFT(Mouseidx,Dayidx) = RTLEFT;
            MINRTRIGHT(Mouseidx,Dayidx) = RTRIGHT;
        else
            percWrong(Mouseidx,Dayidx) = CheckReactions(TrainingList(Fileidx).name);
            load(TrainingList(Fileidx).name)

            X = strcmp(LOG.correctReaction,LOG.Reaction);
            idxincorrect = find(~X);
            nrincorrect(Mouseidx,Dayidx) = length(idxincorrect);
            idxpassive = find(LOG.Gavepassive);
            nrPassives(Mouseidx,Dayidx) = length(find(LOG.Gavepassive));
            nrPassWrong(Mouseidx,Dayidx) = length(find(LOG.Gavepassive(idxincorrect))); 
            percPassives(Mouseidx,Dayidx) = nrPassives(Mouseidx,Dayidx)/length(LOG.Reaction);
            
            if length(idxincorrect)>=1
                percWrongisPass(Mouseidx,Dayidx) = length(idxincorrect)/(length(find(LOG.Gavepassive(idxincorrect)))+length(idxincorrect));
            else
                percWrongisPass(Mouseidx,Dayidx) = nan;
            end

            RTLEFT = nan(size(idxincorrect));
            RTRIGHT = nan(size(idxincorrect));

            for i = 1:length(idxincorrect)

                if ~isempty(LOG.RTleftVec{i})
                    RTLEFT(i) = min(LOG.RTleftVec{i});
                end
                if ~isempty(LOG.RTrightVec{i})
                RTRIGHT(i) = min(LOG.RTrightVec{i});
                end
                
                MINRTLEFT(Mouseidx,Dayidx) = mat2cell(RTLEFT,1);
                MINRTRIGHT(Mouseidx,Dayidx) = mat2cell(RTLEFT,1);

            end
    end

    end
    end
end

%% Create reactiontime bins

edges = 1:50:3000;
rtfrequency = cellfun(@(x) histc(x, edges), MINRTLEFT, 'UniformOutput',false); % now just count the numbers at the specific bin
emptyidx = find(cellfun(@isempty,rtfrequency)); % not used
rtmat = zeros(1,60);

for i = 1:size(rtfrequency,1)
    
   
    for j = 1:size(rtfrequency,2)
        
        if ~isempty(rtfrequency{i,j})
            rtmat = rtmat+rtfrequency{i,j};
        else
    end
end
end



x = plot(rtmat);
xlabel('Earliest RT of a misregistered trial (ms)')
ylabel('Frequency')
lol = [1 10 20 30 40 50 60];
xticks(lol)
xticklabels({'1','451','951','1451','1951','2451','2951'})



figure    
for i = 1:length(Mice)
    hold all
    plot(percPassives(i,:))
    xlabel('Session')
    ylabel('Proportion passives')
end

figure
for i = 1:length(Mice)
    hold all
    plot(percWrongisPass(i,:))
    xlabel('Session')
    ylabel('Proportion wrong registered is passive')
end

figure
for i = 1:length(Mice)
    hold all
    plot(nrPassWrong(i,:))
    xlabel('Session')
    ylabel('# passives wrongly registered')
end

figure
for i = 1:length(Mice)
    hold all
    plot(percWrong(i,:))
    xlabel('Session')
    ylabel('Proportion wrongly registered')
end
        
        
figure
for i = 1:length(Mice)
    hold all
    plot(percWrong(i,:))
    xlabel('Session')
    ylabel('# wrongly registered')
end