load('C:\Users\gillissen\Desktop\InternshipCédric\MainAnaStorage\Frey\Frey20161121\Frey1\Baseline1_0\LEFTVSRIGHT');
load('C:\Users\gillissen\Desktop\InternshipCédric\MainAnaStorage\Frey\brainareamodel');
load('C:\Users\gillissen\Desktop\InternshipCédric\MainAnaStorage\Frey\Frey20161121\Frey1\Frey1_RawData_C1','timeline');
timelim = [-300 2500];
timevec = timeline(timeline>=timelim(1)&timeline<=timelim(2));
%in real version here the task epochs are also defined. 0-500 visual. 500-1450 delay. 1500- onwards response period. 

areas = Model.Regions;

                dFFav = LEFTVSRIGHT.dFFav;
                nrt = LEFTVSRIGHT.nrt;
                meanRT = LEFTVSRIGHT.meanRT;
                stdRT = LEFTVSRIGHT.stdRT;
                SideOpt = LEFTVSRIGHT.SideOpt;
                ReactionOpt = LEFTVSRIGHT.ReactionOpt;
                
                ConditionNames =  LEFTVSRIGHT.ConditionNames;
                
                conditionparts = cellfun(@(X) strsplit(X,' '),ConditionNames,'UniformOutput',0);
                reaction = cellfun(@(X) X{1},conditionparts,'UniformOutput',0); %Reaction
                orientation = cellfun(@(X) X{2},conditionparts,'UniformOutput',0); %orientations
                side = cellfun(@(X) X{4},conditionparts,'UniformOutput',0); %SIdes
                
  
      ConAreaAv = cell(length(SideOpt),length(ReactionOpt));
      corrFCM1 = cell(length(SideOpt),length(ReactionOpt));
     
 %%               
                  
                
  for ridx = 1:length(ReactionOpt) % for each reaction
    for sideidx = 1:length(SideOpt) % for each side, this way you can iterate trough the 2x3 dFFav cell
          if ~isempty(dFFav{sideidx,ridx})
              areaAv = cell(size(areas)); % this will contain the average timeseries of each area of a particular condition
            for areanr = 1:length(areas) % loop trough areas
                tmp = dFFav{sideidx,ridx}; % load condition dataset
                tmp(~repmat(areas{areanr},[1,1,size(tmp,3)]))=nan; % apply brainmask
                areaAv{areanr} = squeeze(nanmean(nanmean(tmp,1),2)); % average pixels of area
            end
             ConAreaAv{sideidx,ridx} = areaAv; % save all average area specific timeseries in the condition cell
          else 
             ConAreaAv{sideidx,ridx} = nan;
          end
         end
        end
  
  % now the average timeseries for each area per condition is calculated and inside ConAreaAv

    for ridx = 1:length(ReactionOpt)
        for sideidx = 1:length(SideOpt)
            tmp = ConAreaAv{sideidx,ridx};
                if iscell(tmp)
                   FC = nan(size(areas)); % initiliaze vector with correlation measures     
                   
                    for areanr = 1:length(areas)
                        timeseries = tmp{areanr,1};
                        if ~isnan(timeseries)
                            FC(areanr) = corr(tmp{40,1},tmp{areanr,1});
                        else
                            FC(areanr) = nan;
                        end
                    end
                        
                else 
                    FC = nan;
                end
                 corrFCM1{sideidx,ridx} = FC;
                
            end
    end
        
    % next up generalize so different sessions can be pooled together.
    % generalize over different mice. 
    % Use number of trials per condition to weigh correlation. 
    
    % plot the seed correlation vector onto the alan brain map.
    % color code for r etc. 
    
    
  
    %% plot correlations 
    
    
          
  