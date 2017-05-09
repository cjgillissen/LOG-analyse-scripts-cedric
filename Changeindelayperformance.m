clear all
close all
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

for Mouseidx = 1:length(Mice)
   for Dayidx = 1: length(Days)
        
        RTleft=[];RTright=[]; %Pre allocate Licks
        ResponseCell = {}; %Preallocate for ResponseCell concatenation
        Side = {};
        Phase = {};
        Fileidx = find(~cellfun(@isempty,cellfun(@(X) strfind(X,[Mice{Mouseidx} '_' Days{Dayidx}]),{TrainingList.name},'UniformOutput',false)));
        currentdelay = [];
        
        if ~isempty(Fileidx)
            
            %% Concatinate blocks
            
            if length(Fileidx)>1   % If there are more than 1 blocks concatinate blocks 
                
                for j = Fileidx
                         percwrong = CheckReactions(TrainingList(j).name);
                         load(TrainingList(j).name)
%                        ResponseCell = LOG.correctReaction;%Use improved reactions
                         RTleft = [RTleft LOG.RTleftVec];
                         RTright = [RTright LOG.RTrightVec]; 
                         ResponseCell = [ResponseCell LOG.correctReaction];
                         Side = [Side LOG.Side];
                         Phase = [Phase LOG.CurrentPhase];
                         currentdelay = [currentdelay LOG.currentdelay];
                         
                end % end of concatination for loop
            else %if not neeeded to concatinate just load the files
                
                percwrong = CheckReactions(TrainingList(Fileidx).name);
                load(TrainingList(Fileidx).name)
                ResponseCell = LOG.correctReaction;%Use improved reactions
                RTleft = LOG.RTleftVec;
                RTright = LOG.RTrightVec;
                Side = LOG.Side;
                Phase = LOG.CurrentPhase;
                currentdelay = LOG.currentdelay;
         
            end % End the loading and initiliazing of variables
            
            if length(currentdelay)~=length(ResponseCell) %In case of passives at end ofsession currentdly and LOG.Reaction 
                %Don't have same length vectors, concatinate to equalize.
                
                ResponseCell = ResponseCell(1:length(currentdelay));
            end

                 if length(unique(currentdelay)) >1 % if there is more than one delay setting used try to find differences in data
                     
                     diffidx = find(diff(currentdelay))+1; % Al the trials of which the previous cue delay was different than currentdelay
                     nrmlidx = find(diff(currentdelay)==0)+1;   %The trials of which the previous trial had similar cue delay
                     longeridx = find(diff(currentdelay)==1500)+1; % when the delay period is longer than expected
                     smallerdly = find(diff(currentdelay)==-1500)+1; % For when the delay is suddenly smaller than expected
                     
       %Calculate performance of trials defined above
                     
% delayperf = sum(find(strcmp(ResponseCell(diffidx),'Hit')))./(sum(find(strcmp(ResponseCell(diffidx),'Error')))+sum(find(strcmp(ResponseCell(diffidx),'Hit'))));
% recurrperf = sum(find(strcmp(ResponseCell(nrmlidx),'Hit')))./(sum(find(strcmp(ResponseCell(nrmlidx),'Error')))+sum(find(strcmp(ResponseCell(nrmlidx),'Hit'))));
% longerdlyperf = sum(find(strcmp(ResponseCell(longeridx),'Hit')))./(sum(find(strcmp(ResponseCell(longeridx),'Error')))+sum(find(strcmp(ResponseCell(longeridx),'Hit'))));
% smallerdlyperf = sum(find(strcmp(ResponseCell(smallerdly),'Hit')))./(sum(find(strcmp(ResponseCell(smallerdly),'Error')))+sum(find(strcmp(ResponseCell(smallerdly),'Hit'))));
% zerocuedelayperf = sum(strcmp(ResponseCell,'Hit')& currentdelay ==0)./(sum(strcmp(ResponseCell,'Error')&currentdelay ==0)+sum(strcmp(ResponseCell,'Hit')& currentdelay ==0));
% maxcuedelayperf = sum(strcmp(ResponseCell,'Hit')&currentdelay ==1500)./(sum(strcmp(ResponseCell,'Error')&currentdelay ==1500)+sum(strcmp(ResponseCell,'Hit')& currentdelay ==1500));
% 
% 
% ChangePerf(Mouseidx,Dayidx)=   delayperf;
% NormalPerf(Mouseidx,Dayidx)=   recurrperf;
% LongerPerf(Mouseidx,Dayidx)=   longerdlyperf;
% SmallerPerf(Mouseidx,Dayidx)=  smallerdlyperf;
% ZeroCuePerf(Mouseidx,Dayidx)=  zerocuedelayperf;

                     
                 end
        end
                        
    end
end

% clrspec = {'b','r','g','p'};
% namespec = {'DiffDelay','SameDelay','LongerDelay','ShorterDelay','ZeroCueDelay','MaxCueDelay'};
% 
% %             superbar(superbar(nanmean(ChangePerf)))
% %             bar(mean(ChangePerf(i,:)));
% %             bar(mean(NormalPerf(i,:)));
% %             bar(mean(LongerPerf(i,:)));
% %             bar(mean(SmallerPerf(i,:)));
%      
%     Alladin = [nanmean(ChangePerf(1,:)) nanmean(NormalPerf(1,:)) nanmean(LongerPerf(1,:)) nanmean(SmallerPerf(1,:)) nanmean(ZeroCuePerf(1,:))];
%     Chief = [nanmean(ChangePerf(2,:)) nanmean(NormalPerf(2,:)) nanmean(LongerPerf(2,:)) nanmean(SmallerPerf(2,:)) nanmean(ZeroCuePerf(2,:))];
%     Esmeralda = [nanmean(ChangePerf(3,:)) nanmean(NormalPerf(3,:)) nanmean(LongerPerf(3,:)) nanmean(SmallerPerf(3,:)) nanmean(ZeroCuePerf(3,:))];
%     Frey = [nanmean(ChangePerf(4,:)) nanmean(NormalPerf(4,:)) nanmean(LongerPerf(4,:)) nanmean(SmallerPerf(4,:)) nanmean(ZeroCuePerf(4,:))];
%     
%                     h = bar([nanmean(ChangePerf) nanmean(NormalPerf) nanmean(LongerPerf) nanmean(SmallerPerf) nanmean(ZeroCuePerf) nanmean(maxcuedelayperf)]);
%                     xlabel('Condition')
%                     ylabel('Performance')
%                     ylim([0 1]);
%                     xticklabels(namespec);
%                     
% %                    
%                     
                    
%% Plotting bar graph

micenames={'Alladin', 'Esmeralda', 'Frey','Chief'};              
Scores = [Alladin' Esmeralda' Frey' Chief'];
b = bar(Scores);
legend(b,micenames)
xticklabels(namespec)
ylim([0 1]);


%% normal performance line graph




































