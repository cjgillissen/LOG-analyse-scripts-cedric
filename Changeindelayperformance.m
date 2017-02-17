clear all
close all

addpath(genpath('C:\Users\gillissen\Desktop\InternshipCédric\Scripts\Mouse'))
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



%% Check trials which have a different delay period than t(i-1)

% Take trials whose cue delay is different than the previous trial

ChangePerf = NaN(length(Mice),length(Days));
NormalPerf = NaN(length(Mice),length(Days));
LongerPerf = NaN(length(Mice),length(Days));
SmallerPerf = NaN(length(Mice),length(Days));
ZeroCuePerf = NaN(length(Mice),length(Days));

for Mouseidx = 1:length(Mice);
    
    
    
    for Dayidx = 1: length(Days);
        
        RTleft=[];RTright=[]; %Pre allocate Licks
        ResponseCell = {}; %Preallocate for LOG.Reaction concatenation
        Side = {};
        Phase = {};
        Fileidx = find(~cellfun(@isempty,cellfun(@(X) strfind(X,[Mice{Mouseidx} '_' Days{Dayidx}]),{TrainingList.name},'UniformOutput',false)));
        currentdelay = [];
        
        if ~isempty(Fileidx)
            
            if length(Fileidx)>1;   % If there are more than 1 blocks concatinate blocks 
                
                for j = Fileidx
                         percwrong = CheckReactions(TrainingList(j).name);
                         load(TrainingList(j).name)
                         LOG.Reaction = LOG.correctReaction;%Use improved reactions
                         RTleft = [RTleft LOG.RTleftVec];
                         RTright = [RTright LOG.RTrightVec]; 
                         ResponseCell = [ResponseCell LOG.Reaction];
                         Side = [Side LOG.Side];
                         Phase = [Phase LOG.CurrentPhase];
                         currentdelay = [currentdelay LOG.currentdelay];
                         
                         if length(unique(LOG.currentdelay)) >1;
                     
                     diffidx = find(diff(LOG.currentdelay))+1;
                     nrmlidx = find(diff(LOG.currentdelay)==0)+1;
                     longeridx = find(diff(LOG.currentdelay)==1500)+1;
                     smallerdly = find(diff(LOG.currentdelay)==-1500)+1;
                     
                     delayperf = sum(find(strcmp(LOG.Reaction(diffidx),'Hit')))./(sum(find(strcmp(LOG.Reaction(diffidx),'Error')))+sum(find(strcmp(LOG.Reaction(diffidx),'Hit'))));
                     recurrperf = sum(find(strcmp(LOG.Reaction(nrmlidx),'Hit')))./(sum(find(strcmp(LOG.Reaction(nrmlidx),'Error')))+sum(find(strcmp(LOG.Reaction(nrmlidx),'Hit'))));
                     longerdlyperf = sum(find(strcmp(LOG.Reaction(longeridx),'Hit')))./(sum(find(strcmp(LOG.Reaction(longeridx),'Error')))+sum(find(strcmp(LOG.Reaction(longeridx),'Hit'))));
                     smallerdlyperf = sum(find(strcmp(LOG.Reaction(smallerdly),'Hit')))./(sum(find(strcmp(LOG.Reaction(smallerdly),'Error')))+sum(find(strcmp(LOG.Reaction(smallerdly),'Hit'))));
                     zerocuedelayperf = sum(strcmp(LOG.Reaction,'Hit')& LOG.currentdelay ==0)./(sum(strcmp(LOG.Reaction,'Error')&LOG.currentdelay ==0)+sum(strcmp(LOG.Reaction,'Hit')& LOG.currentdelay ==0));
   
                         if isnan(delayperf)~=1;
                             ChangePerf(Mouseidx,Dayidx)=delayperf;
                             NormalPerf(Mouseidx,Dayidx) = recurrperf;
                             LongerPerf(Mouseidx,Dayidx) = longerdlyperf;
                             SmallerPerf(Mouseidx,Dayidx) = smallerdlyperf;
                             ZeroCuePerf(Mouseidx,Dayidx) = zerocuedelayperf;
                         end
                 end
                         
                                 
                                 
                                 
                end
            else
                percwrong = CheckReactions(TrainingList(Fileidx).name);
                load(TrainingList(Fileidx).name)
                LOG.Reaction = LOG.correctReaction;%Use improved reactions
                RTleft = LOG.RTleftVec;
                RTright = LOG.RTrightVec;
                ResponseCell = LOG.Reaction;
                Side = LOG.Side;
                Phase = LOG.CurrentPhase;
                
                 if length(unique(LOG.currentdelay)) >1;
                     
                     diffidx = find(diff(LOG.currentdelay))+1;
                     nrmlidx = find(diff(LOG.currentdelay)==0)+1;
                     longeridx = find(diff(LOG.currentdelay)==1500)+1;
                     smallerdly = find(diff(LOG.currentdelay)==-1500)+1;
                     
                     delayperf = sum(find(strcmp(LOG.Reaction(diffidx),'Hit')))./(sum(find(strcmp(LOG.Reaction(diffidx),'Error')))+sum(find(strcmp(LOG.Reaction(diffidx),'Hit'))));
                     recurrperf = sum(find(strcmp(LOG.Reaction(nrmlidx),'Hit')))./(sum(find(strcmp(LOG.Reaction(nrmlidx),'Error')))+sum(find(strcmp(LOG.Reaction(nrmlidx),'Hit'))));
                     longerdlyperf = sum(find(strcmp(LOG.Reaction(longeridx),'Hit')))./(sum(find(strcmp(LOG.Reaction(longeridx),'Error')))+sum(find(strcmp(LOG.Reaction(longeridx),'Hit'))));
                     smallerdlyperf = sum(find(strcmp(LOG.Reaction(smallerdly),'Hit')))./(sum(find(strcmp(LOG.Reaction(smallerdly),'Error')))+sum(find(strcmp(LOG.Reaction(smallerdly),'Hit'))));
                     zerocuedelayperf = sum(strcmp(LOG.Reaction,'Hit')& LOG.currentdelay ==0)./(sum(strcmp(LOG.Reaction,'Error')&LOG.currentdelay ==0)+sum(strcmp(LOG.Reaction,'Hit')& LOG.currentdelay ==0));
                         if isnan(delayperf)~=1;
                             ChangePerf(Mouseidx,Dayidx)=delayperf;
                             NormalPerf(Mouseidx,Dayidx) = recurrperf;
                             LongerPerf(Mouseidx,Dayidx) = longerdlyperf;
                             SmallerPerf(Mouseidx,Dayidx) = smallerdlyperf;
                             ZeroCuePerf(Mouseidx,Dayidx) = zerocuedelayperf;
                         end
                 end
                
                
                
                
                
            end
        end
    end
end

clrspec = {'b','r','g','p'};
       
     
           
%             superbar(superbar(nanmean(ChangePerf)))
%             bar(mean(ChangePerf(i,:)));
%             bar(mean(NormalPerf(i,:)));
%             bar(mean(LongerPerf(i,:)));
%             bar(mean(SmallerPerf(i,:)));
     
    Alladin = [nanmean(ChangePerf(1,:)) nanmean(NormalPerf(1,:)) nanmean(LongerPerf(1,:)) nanmean(SmallerPerf(1,:)) nanmean(ZeroCuePerf(1,:))];
    Chief = [nanmean(ChangePerf(2,:)) nanmean(NormalPerf(2,:)) nanmean(LongerPerf(2,:)) nanmean(SmallerPerf(2,:)) nanmean(ZeroCuePerf(2,:))];
    Esmeralda = [nanmean(ChangePerf(3,:)) nanmean(NormalPerf(3,:)) nanmean(LongerPerf(3,:)) nanmean(SmallerPerf(3,:)) nanmean(ZeroCuePerf(3,:))];
    Frey = [nanmean(ChangePerf(4,:)) nanmean(NormalPerf(4,:)) nanmean(LongerPerf(4,:)) nanmean(SmallerPerf(4,:)) nanmean(ZeroCuePerf(4,:))];

Scores = [Alladin' Chief' Esmeralda' Frey'];

bar(1:length(Mice),Scores)
    
    
    
