clear all
close all
cd('C:\Users\gillissen\Desktop\InternshipCédric\TrainingLogs\FGTask')

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
        
%         RTleft=[];RTright=[]; %Pre allocate Licks
        ResponseCell = {}; %Preallocate for ResponseCell concatenation
        Side = {};
%         Phase = {};
        Fileidx = find(~cellfun(@isempty,cellfun(@(X) strfind(X,[Mice{Mouseidx} '_' Days{Dayidx}]),{TrainingList.name},'UniformOutput',false)));
%         currentdelay = [];
        
        if ~isempty(Fileidx)
            
            %% Concatinate blocks
            
            if length(Fileidx)>1   % If there are more than 1 blocks concatinate blocks 
                
                for j = Fileidx
%                          percwrong = CheckReactions(TrainingList(j).name);
                         load(TrainingList(j).name)
%                        ResponseCell = LOG.correctReaction;%Use improved reactions
                         RTleft = [RTleft LOG.RTleftVec];
                         RTright = [RTright LOG.RTrightVec]; 
                         ResponseCell = [ResponseCell LOG.Reaction];
                         Side = [Side LOG.Side];
%                          Phase = [Phase LOG.CurrentPhase];
%                          currentdelay = [currentdelay LOG.currentdelay];
                         
                end % end of concatination for loop
            else %if not neeeded to concatinate just load the files
                
%                 percwrong = CheckReactions(TrainingList(Fileidx).name);
                load(TrainingList(Fileidx).name)
                ResponseCell = LOG.Reaction;%Use improved reactions
                RTleft = LOG.RTleftVec;
                RTright = LOG.RTrightVec;
                Side = LOG.Side;
%                 Phase = LOG.CurrentPhase;
%                 currentdelay = LOG.currentdelay;
         
            end % End the loading and initiliazing of variables
            
            if length(currentdelay)~=length(ResponseCell) %In case of passives at end ofsession currentdly and LOG.Reaction 
                %Don't have same length vectors, concatinate to equalize.
                
                ResponseCell = ResponseCell(1:length(currentdelay));
            end

                 NormalPerf(Mouseidx,Dayidx) = sum(strcmp(ResponseCell,'Hit'))/(sum(strcmp(ResponseCell,'Hit'))+sum(strcmp(ResponseCell,'Error')));
                 
        end
   end
end


%% normal performance line graph

figure
for mouseidx = 1:size(NormalPerf,1)
    ylim([0 1])
    plot(NormalPerf(mouseidx,:),'DisplayName',sprintf('%s',Mice{mouseidx}))
    hold on
    xlabel('Days')
    ylabel('Performance H/E+H')
    legend('-DynamicLegend');
end
hold on
plot(1:size(NormalPerf,2),repmat(0.5,size(NormalPerf,2),1),'k')
plot(1:size(NormalPerf,2),repmat(0.7,size(NormalPerf,2),1),'g')


   
    
    


































