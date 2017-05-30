outputfile = '\\vc2nin\WBimaging\MainAnalysis\Outputs.mat';
DataDirectory = '\\vc2nin\WBimaging\'
storepath = '\\vc2nin\WBimaging\MainAnalysis\'
load(outputfile)
storefigures = '\\vc2nin\WBimaging\MainAnalysis\Figures\Baseline1_eqsample0\SumFigures\';
%% Amplitude pRF
TW = [150 500];
twidx = find(timeline>=TW(1)&timeline<=TW(2));
for midx = 1:length(miceopt)
    %load pRF
    %Check with mapping if available
    load(fullfile(storepath,miceopt{midx},'brainareamodel.mat'))
    
    pRFfile = fullfile(DataDirectory,'pRF Results',miceopt{midx},'4Mapping','pRFmaps');
    pRFbrain = fullfile(DataDirectory,'pRF Results',miceopt{midx},'4Mapping','RefImg');
    
    %Load up teh reference image from the pRF maps
    load(pRFfile)
    load(pRFbrain)
    
    if ~exist('FGAmplPerRegio','var')
        FGAmplPerRegio = nan(length(Model.Regions),length(miceopt));
    end
    
    mask = abs(out.SIGNMAPt);
    mask = cat(2,zeros([size(mask,1),size(Model.Regions{1},2)-size(mask,2)]),mask);
    for regidx = 1:length(Model.Regions)
        
        mask2 = false(size(mask,1),size(mask,2));
        mask2(Model.Regions{regidx}' == 1 & mask == 1) = true;
        
        if sum(mask2(:)) <= 20
            continue
        end
        figure; subplot(2,3,1); imagesc(Model.Regions{regidx}');
        axis square
        subplot(2,3,2)
        imagesc(mask);         axis square
        subplot(2,3,3);
        imagesc(mask2);
        axis square
        
        title(Model.Rnames{regidx})
        mask2 = repmat(mask2,[1,1,37]);
        %         tmp = abs(AllMicedFF{midx}{1}{1}-AllMicedFF{midx}{1}{2}); %FG-Mod
        %         tmp(~mask2) = nan;
        tmp1=AllMicedFF{midx}{1}{1};
       
        %Half the cortex is where lambda y coord is is: 
        halfidx = round(Model.Lambda(2));
        tmp1(~mask2) = nan;
        tmp2 = AllMicedFF{midx}{1}{2};
        tmp2(~mask2)= nan;
        
        tmp = nan(size(tmp1));
        %Right side of the cortex (L-R)
        tmp(1:halfidx,:) = (tmp1(1:halfidx,:) - tmp2(1:halfidx,:));
        %Left Side of the cortex (R-L)
        tmp(halfidx+1:end,:) = (tmp2(halfidx+1:end,:) - tmp1(halfidx+1:end,:));
        
        %Right side traces
        tmpdFFR = squeeze(nanmean(nanmean(tmp(1:halfidx,:,:),1),2));
        tmpdFF1R = squeeze(nanmean(nanmean(tmp1(1:halfidx,:,:),1),2));
        tmpdFF2R = squeeze(nanmean(nanmean(tmp2(1:halfidx,:,:),1),2));
        subplot(2,3,4)
        plot(timeline(timeline>=-300&timeline<=1500),tmpdFF1R,'g'); hold on; plot(timeline(timeline>=-300&timeline<=1500),tmpdFF2R,'r');
        plot(timeline(timeline>=-300&timeline<=1500),tmpdFFR,'b')
        title('Traces Right Side cortex (L-R)')
        xlim([-200 800])
        ylim([-3*10^-3 17*10^-3])
        
        %Left Side traces
        tmpdFFL = squeeze(nanmean(nanmean(tmp(halfidx:end,:,:),1),2));
        tmpdFF1L = squeeze(nanmean(nanmean(tmp2(halfidx:end,:,:),1),2));
        tmpdFF2L = squeeze(nanmean(nanmean(tmp1(halfidx:end,:,:),1),2));
        subplot(2,3,5)
        plot(timeline(timeline>=-300&timeline<=1500),tmpdFF1L,'g'); hold on; plot(timeline(timeline>=-300&timeline<=1500),tmpdFF2L,'r');
        plot(timeline(timeline>=-300&timeline<=1500),tmpdFFL,'b')
        title('Traces Left Side cortex (R-L)')
        xlim([-200 800])
        ylim([-3*10^-3 17*10^-3])
        
        subtmp = subplot(2,3,6);
        posleg = get(subtmp,'position');
        delete(subtmp)
    
        legend({'Figure','Ground','Modulation'},'Position',posleg)
        FGAmplPerRegio(regidx,midx) = nanmean(nanmean([tmpdFFR(twidx),tmpdFFL(twidx)],2),1);
        
        saveas(gcf,fullfile(storefigures,[miceopt{midx} '_' Model.Rnames{regidx} '_.bmp']))
        saveas(gcf,fullfile(storefigures,[miceopt{midx} '_' Model.Rnames{regidx} '_.fig']))
        saveas(gcf,fullfile(storefigures,[miceopt{midx} '_' Model.Rnames{regidx} '_.eps']))      

    end
    
    %Plotting extra Rois:
    figure; imagesc(brain)
    colormap gray
    hold on
    scatter(Model.AllX,Model.AllY,'y.')
    title(['Press k for enough ROIs or any other key to draw more'])
    roicount = 1;
    stop = 0;
    while ~stop
        BW{roicount} = roipoly;
        roidrawn(roicount) = bwboundaries(BW{roicount});

        scatter(roidrawn{roicount}(:,2),roidrawn{roicount}(:,1),'r.')
        text(nanmean(roidrawn{roicount}(:,2)),nanmean(roidrawn{roicount}(:,1)),num2str(roicount))
        disp('Press keyboard!') 
        waitforbuttonpress
         key = get(gcf,'CurrentCharacter');
         if strcmp(key,'k')
             stop = 1;
         end     
         roicount = roicount+1;
    end
    
      for regidx = 1:length(roidrawn)
        
     
        figure; subplot(2,2,1); imagesc(BW{regidx});
        axis square
      
        
        title(['roi ' num2str(regidx)])
      
        tmp1=AllMicedFF{midx}{1}{1};
        mask2 = repmat(BW{regidx},[1,1,size(tmp1,3)]);
        %Half the cortex is where lambda y coord is is: 
        halfidx = round(Model.Lambda(2));
        tmp1(~mask2) = nan;
        tmp2 = AllMicedFF{midx}{1}{2};
        tmp2(~mask2)= nan;
        
        tmp = nan(size(tmp1));
        %Right side of the cortex (L-R)
        tmp(1:halfidx,:) = (tmp1(1:halfidx,:) - tmp2(1:halfidx,:));
        %Left Side of the cortex (R-L)
        tmp(halfidx+1:end,:) = (tmp2(halfidx+1:end,:) - tmp1(halfidx+1:end,:));
        
        %Right side traces
        tmpdFFR = squeeze(nanmean(nanmean(tmp(1:halfidx,:,:),1),2));
        tmpdFF1R = squeeze(nanmean(nanmean(tmp1(1:halfidx,:,:),1),2));
        tmpdFF2R = squeeze(nanmean(nanmean(tmp2(1:halfidx,:,:),1),2));
        if any(~isnan((unique(tmpdFFR))))
        
        subplot(2,2,3)
        plot(timeline(timeline>=-300&timeline<=1500),tmpdFF1R,'g'); hold on; plot(timeline(timeline>=-300&timeline<=1500),tmpdFF2R,'r');
        plot(timeline(timeline>=-300&timeline<=1500),tmpdFFR,'b')
        title('Traces Right Side cortex (L-R)')
        xlim([-200 800])
        ylim([-3*10^-3 17*10^-3])
        end
        
        %Left Side traces
        tmpdFFL = squeeze(nanmean(nanmean(tmp(halfidx:end,:,:),1),2));
        tmpdFF1L = squeeze(nanmean(nanmean(tmp2(halfidx:end,:,:),1),2));
        tmpdFF2L = squeeze(nanmean(nanmean(tmp1(halfidx:end,:,:),1),2));
        if any(~isnan((unique(tmpdFFL))))
            subplot(2,2,4)
            plot(timeline(timeline>=-300&timeline<=1500),tmpdFF1L,'g'); hold on; plot(timeline(timeline>=-300&timeline<=1500),tmpdFF2L,'r');
            plot(timeline(timeline>=-300&timeline<=1500),tmpdFFL,'b')
            title('Traces Left Side cortex (R-L)')
            xlim([-200 800])
            ylim([-3*10^-3 17*10^-3])
        end
        
        subtmp = subplot(2,2,2);
        posleg = get(subtmp,'position');
        delete(subtmp)
    
        legend({'Figure','Ground','Modulation'},'Position',posleg)
        FGAmplPerRegio(regidx,midx) = nanmean(nanmean([tmpdFFR(twidx),tmpdFFL(twidx)],2),1);
        
        saveas(gcf,fullfile(storefigures,[miceopt{midx} '_ExtraROI' num2str(regidx)  '_.bmp']))
        saveas(gcf,fullfile(storefigures,[miceopt{midx} '_ExtraROI' num2str(regidx) '_.fig']))
        saveas(gcf,fullfile(storefigures,[miceopt{midx} '_ExtraROI' num2str(regidx) '_.eps']))      

    end
    
    
end

FGAmplPerRegio