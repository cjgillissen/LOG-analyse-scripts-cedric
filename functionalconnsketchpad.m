load('C:\Users\gillissen\Desktop\InternshipCédric\MainAnaStorage\Frey\Frey20161121\Frey1\Baseline1_1500\LEFTVSRIGHT');
load('C:\Users\gillissen\Desktop\InternshipCédric\MainAnaStorage\Frey\brainareamodel');
load('C:\Users\gillissen\Desktop\InternshipCédric\MainAnaStorage\Frey\Frey20161121\Frey1\Frey1_RawData_C1','timeline');
% load('C:\Users\gillissen\Desktop\InternshipCédric\MainAnaStorage\Frey\Frey20161121\Frey1\Frey1_RawData_C1','timeline');
% timelim = [-300 2500];
% timevec = timeline(timeline>=timelim(1)&timeline<=timelim(2));
%in real version here the task epochs are also defined. 0-500 visual. 500-1450 delay. 1500- onwards response period. 


brainmask = zeros(800,800);

%% Brainmask
            % create cell with the proper area logicals
            areas = Model.Regions;
            brainmask = zeros(800,800);
            for i = 1:length(Model.Regions)
                Borders = Model.Boundaries{i};
                for j = 1:length(Borders)
                    tmp = poly2mask(Borders{j}(:,1),Borders{j}(:,2),800,800);
                    tmp = imfill(tmp,'holes');
                    brainmask(tmp) = 1+i;
                end
            end
            brainmask = imfill(brainmask,'holes');
            
            
    
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
                
                clear LEFTVSRIGHT
  
      ConAreaAv = cell(length(SideOpt),length(ReactionOpt));
      corrFCM1 = cell(length(SideOpt),length(ReactionOpt));
      cohereFCM1 = cell(length(SideOpt),length(ReactionOpt));
      Seedpixelcorr = cell(length(SideOpt),length(ReactionOpt));
     
 %% Apply brain mask and collect average timeseries per cond per area
 % later per subject! and concatenate sessions!
                  
                
  for ridx = 1:length(ReactionOpt) % for each reaction
    for sideidx = 1:length(SideOpt) % for each side, this way you can iterate trough the 2x3 dFFav cell
          if ~isempty(dFFav{sideidx,ridx})
              areaAv = cell(size(areas)); % this will contain the average timeseries of each area of a particular condition
            for areanr = 1:length(areas) % loop trough areas
                tmp = dFFav{sideidx,ridx}; % load condition dataset
                tmp(~repmat(brainmask==areanr,[1,1,size(tmp,3)]))=nan; % apply brainmask
                areaAv{areanr} = squeeze(nanmean(nanmean(tmp,1),2)); % average pixels of area
            end
             ConAreaAv{sideidx,ridx} = areaAv; % save all average area specific timeseries in the condition cell
          else 
             ConAreaAv{sideidx,ridx} = nan;
         end
     end
 end
  
  % now the average timeseries for each area per condition is calculated and inside ConAreaAv
%% Correlation
% per area
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
    
     %% Correlation per pixel
     
     fig = imagesc(brainmask);
     seed = roipoly;
     seedtemp = seed;
     
      for ridx = 1:length(ReactionOpt)
        for sideidx = 1:length(SideOpt)
            tmp = dFFav{sideidx,ridx};
            tmp2 = tmp;
            seedtemp = repmat(seedtemp,[1,1,size(tmp,3)]);
            tmp2(~seedtemp) = nan;
            seedtemp = nanmean(tmp2,1);
            seedtemp = squeeze(nanmean(seedtemp,2));
            
            FC = nan(size(tmp)); % initiliaze vector with correlation measures     
                for pixelx = 1:size(tmp,1)
                    for pixely = 1:size(tmp,2) 
                     
                    Seedpixelcorr{sideidx,ridx}(pixely,pixelx) = corr(seedtemp,squeeze(tmp(pixelx,pixely,:)));
                    end
                end
        end
      end
                      
                        %% Coherence
    
    for ridx = 1:length(ReactionOpt)
        for sideidx = 1:length(SideOpt)
            tmp = ConAreaAv{sideidx,ridx};
                if iscell(tmp)
                   FC = cell(size(areas)); % initiliaze vector with correlation measures     
                    for areanr = 1:length(areas)
                        timeseries = tmp{areanr,1};
                        if ~isnan(timeseries)
                            FC{areanr} = mscohere(tmp{40,1},tmp{areanr,1});
                        else
                            FC{areanr} = nan;
                        end
                    end
                 else 
                    FC = nan;
                end
                 cohereFCM1{sideidx,ridx} = FC;
         end
    end
      
    %% plot correlations 
    
    %first make figure for each condition seperate
    % transpose corr vector. Concatinate per column. M being the areas and
    % N being the correlations per each codintion.
   
   names  = Model.Rnames;
   xticksss = names(9:30);
   onecor = corrFCM1{1,1}; 
   f =  figure;
   bar(9:43,onecor(9:43));
   
   
   %% plot coherence on refbrain
   
   figure
   
   for i = 1:length(cohereFCM1{1,1}{3,1})
       hold on
       plot(1:129,cohereFCM1{1,1}{i,1})
   end
   
   ref = load('refimg');
   colorlim = [0 95];
   ref = zeros(800,800);
   for areanr = 1:length(areas)
       ref(brainmask==areanr) = onecor(areanr);
   end
   imagesc(ref)
   
   %% Plot correlation on refbrain
   
   figure
   
%    for i = 1:length(corrFCM1{1,1}(3,1))
%        hold on
%        plot(1:129,corrFCM1{1,1}(i,1))
%    end
   colorlim = [-1 99];
   
   for ridx = 1:length(ReactionOpt)
       for sideidx = 1:length(SideOpt)
            figure
            ref = zeros(800,800);
            for areanr = 1:length(areas)
            ref(brainmask==areanr) = corrFCM1{sideidx,ridx}(areanr,1);
            end
            
                imagesc(ref)
                title(sprintf('correlation with M1 as seed %s %s',ReactionOpt{ridx},SideOpt{sideidx}))
                
       end
   end
   
           
  
   %% between conditions coherence V1 vs M1
   figure 
   n = 0;
   x=cell(length(SideOpt),length(ReactionOpt)); % doesn't work
   for ridx = 1:length(ReactionOpt)
       for sideidx = 1:length(SideOpt)
           n=n+1;
           hold on
           plot(cohereFCM1{sideidx,ridx}{11,:})
           x{n} = sprintf('%s %s ',ReactionOpt{ridx},SideOpt{sideidx});
       end
   end
           legend(ConditionNames{1:(length(SideOpt)*length(ReactionOpt))});
           
   
         
   
   
   
   
   
   
    
   
    
    
          
  