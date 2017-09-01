%% Biasindex calculation of a couple of sessions

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


BiasSession = NaN(length(Mice),length(Days));


for Mouseidx = 1:length(Mice)
   for Dayidx = 1: length(Days)
        
        RTleft=[];RTright=[]; %Pre allocate Licks
        ResponseCell = {}; %Preallocate for ResponseCell concatenation
        Side = {};
        Phase = {};
        Fileidx = find(~cellfun(@isempty,cellfun(@(X) strfind(X,[Mice{Mouseidx} '_' Days{Dayidx}]),{TrainingList.name},'UniformOutput',false)));
        currentdelay = [];
        Biass=[];
        Biasmat = [];
        if ~isempty(Fileidx)
            
            %% Concatinate blocks
            
            if length(Fileidx)>1   % If there are more than 1 blocks concatinate blocks 
                
                for j = Fileidx
                         load(TrainingList(j).name)
                      [Perfmat,Biasmat] = extractbiasidx(LOG,15);
                      Biass = [Biass Biasmat];
                         
                end % end of concatination for loop
                 
                 BiasSession(Mouseidx,Dayidx) = Biass;
                
            else %if not neeeded to concatinate just load the files
                
                load(TrainingList(Fileidx).name)
                
                
         
            end % End the loading and initiliazing of variables
            
%                   
                     
                 end
        end
                        
    end
end

%                     
                    
%% Plotting bar graph

micenames={'Alladin', 'Esmeralda', 'Frey','Chief'};              
Scores = [Alladin' Esmeralda' Frey' Chief'];
b = bar(Scores);
legend(b,micenames)
xticklabels(namespec)
ylim([0 1]);







