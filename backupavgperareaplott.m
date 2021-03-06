miceopt = {'Alladin','Chief','Esmeralda','Frey'}%%% %options for mice
basel = [-300 -50];
visint = [50 200];
vistw = [200 500]; %FOR FG
Delayearlytw = [600 1350];
Delaylatetw = [1400 1600];
Responstw = [1900 2200];
TW = {basel,visint,vistw,Delayearlytw,Responstw,Delaylatetw};
StorePath = '\\vcnin\mouse_working_memory\MainAnalysis\';
scalefct = 0.5;
xpix = 400;
ypix = 400;
allmicetmp = cell(length(miceopt),length(TW));
%%
regio2take = [11,12,13,16,18,23,25,26,27,28,29,30,31,32,37,38,39,40];


% rgnames = {'V1','Vpor','Va','M1','M2'};
% regio2take = 10:43;
% regio2take = 11;
avgcp = 1:length(TW);
plotfig =figure('name','Lineplot average cp per area per tw');
imagescfig = figure('name' ,'ImagescFig average cp per area per tw');
SideAreaAVG = cell(length(miceopt),length(TW));
avgtmpperside = nan(length(miceopt),length(TW),length(regio2take));
SideAreaAVG = cellfun(@(X) nan(xpix,ypix),SideAreaAVG,'UniformOutput',0);
BilateralAVG = SideAreaAVG;
avgtmpbilateral = avgtmpperside;

for midx = 1:length(miceopt)
    mouse = miceopt{midx};
    BrainModel{midx} = load(fullfile(StorePath,mouse,'brainareamodel.mat'));
    masktmpperside = zeros(xpix,ypix);
    maskbilateral = zeros(xpix,ypix);
    
    
    for roiidx = 1:length(regio2take)
        disp(['Mouse ' num2str(midx) ', region ' num2str(roiidx) ' of ' num2str(length(regio2take))])
        Borders = BrainModel{midx}.Model.Boundaries{regio2take(roiidx)};
        masktmpperside = zeros(xpix,ypix);
        masktmpbilateral = zeros(xpix,ypix);
        for twidx = 1:length(TW)
            CP = Perf{midx,1,twidx}.SP;
            CP = reshape(CP,xpix,ypix);
        
        for roi2dx = 1:length(Borders)
            mask = poly2mask(Borders{roi2dx}(:,1).*scalefct,Borders{roi2dx}(:,2).*scalefct,xpix,ypix);
            %Shrink to not have border effects
            mask = bwmorph(mask,'shrink',1);
            masktmpperside(mask)=roiidx+roi2dx;
            maskbilateral(mask) = roiidx;
            avgcp(twidx) = nanmean(CP(masktmpperside==roiidx+roi2dx));
            avgtmpperside(midx,twidx,roiidx+roi2dx) = nanmean(CP(masktmpperside==(roiidx+roi2dx)));
            allmicetmp{midx,twidx} = CP;
            SideAreaAVG{midx,twidx}(masktmpperside==(roiidx+roi2dx)) = avgtmpperside(midx,twidx,roiidx+roi2dx);
        end
        
        avgtmpbilateral(midx,twidx) = nanmean(CP(maskbilateral==roiidx));
        BilateralAVG{midx,twidx}(maskbilateral==roiidx) = avgtmpbilateral(midx,twidx,roiidx);
        
        end
        figure(plotfig)
        subplot(3,2,midx)
        plot(1:length(TW),avgcp)
        title([miceopt{midx}])
        xlabel('Timewindow, fix xticks')
        ylabel('AUC')
        hold on
    end
legend(BrainModel{midx}.Model.Rnames{regio2take})

figure(imagescfig);
for twidx = 1:length(TW)
    subplot(length(miceopt),length(TW),(midx-1)*length(TW)+twidx);
    imagesc(SideAreaAVG{midx,twidx},[0.3 0.7])
       
end
end



%% plot bargraph

  for twidx = 1:length(TW)
      




    % plot bargraph with errors
    

        Y = [11 14 13;
            15 12 16];
        E = [ 3  4  2;
            5  2  3];

        % Plot with bar
        subplot(1, 2, 1);
        text(0.5, 0.5, 'not implemented', 'HorizontalAlignment', 'center');
        title('bar');

        % Plot with superbar
        subplot(1, 2, 2);
        superbar(Y, 'E', E);
        title('superbar');

        % Since one does not know the X locations of the individual bars, BAR and
        % ERRORBAR can not be used together to add error bars to grouped bar plots.


        set(gca, 'XTick', []); % place vector
        %         set(gca, 'XTick', [1 2])
        %         set(gca, 'XTickLabel', {'Model1' 'Model2'})

        set(gca, 'YTick', []);

  
    
    
 

%% avg across mice
% ciplot lower upper x
% avgmice = squeeze(nanmean(avgtmp,1));
% stdmice = squeeze(std(avgtmp));
% upperlim = avgmice+stdmice;
% lowerlim = avgmice-stdmice;
% 
% 
% avgmap = cellfun(nanmean,allmicetmp{
% test = nanmean(avgmap,1)

 
avgallmice = squeeze(nanmean(avgtmpperside,1));
stdallmice = std(avgtmpperside,1,1);
collorross = {'g','r','b','y','m'};

figure
for i = 1:size(avgallmice,2)
shadedErrorBar(1:size(avgallmice,1),avgallmice(:,i),std(avgtmpperside(:,:,i),1,1),collorross{i},1)
hold on
end


figure;
hold on
 areaidx = 1
  ciplot(lowerlim(:,areaidx),upperlim(:,areaidx),avgmice(:,areaidx))
  hold on

plot(lowerlim)

% Examples
y=randn(30,80); x=1:size(y,2);
shadedErrorBar(x,mean(y,1),std(y),'g');
shadedErrorBar(x,y,{@median,@std},{'r-o','markerfacecolor','r'});    
shadedErrorBar([],y,{@median,@std},{'r-o','markerfacecolor','r'});    

% Overlay two transparent lines
y=randn(30,80)*10; x=(1:size(y,2))-40;
shadedErrorBar(x,y,{@mean,@std},'-r',1); 
hold on
y=ones(30,1)*x; y=y+0.06*y.^2+randn(size(y))*10;
shadedErrorBar(x,y,{@mean,@std},'-b',1);   
hold off










