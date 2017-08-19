miceopt = {'Alladin','Chief','Esmeralda','Frey'}%%% %options for mice
basel = [-300 -50];
visint = [50 200];
vistw = [200 500]; %FOR FG
Delayearlytw = [600 1350];
Delaylatetw = [1400 1600];
Responstw = [1900 2200];
TW = {basel,visint,vistw,Delayearlytw,Responstw,Delaylatetw};
StorePath = '\\vcnin\mouse_working_memory\MainAnalysis\';
scalefct = 1;
xpix = 800;
ypix = 800;
allmicetmp = cell(length(miceopt),length(TW));
SideOpt = {'Left','Right'};
ReactionOpt = {'Error','Hit'};

%% Colormap

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
%%
% regio2take = [11,12,13,16,18,23,25,26,27,28,29,30,31,32,37,38,39,40];
% regio2take = [11,13,18,23,25,26,27,29,30,31,32,37,38,39];

% rgnames = {'V1','Vpor','Va','M1','M2'};
regio2take = 10:43;
% regio2take = 11;
avgcp = 1:length(TW);
plotfig =figure('name','Lineplot average SP per area per tw');
imagescfig = figure('name' ,'ImagescFig average SP per area per tw');
SideAreaAVG = cell(length(miceopt),length(TW));
avgtmpperside = nan(length(miceopt),length(TW),length(regio2take));
SideAreaAVG = cellfun(@(X) nan(xpix,ypix),SideAreaAVG,'UniformOutput',0);
BilateralAVG = SideAreaAVG;
avgtmpbilateral = avgtmpperside;
minntrials = zeros(1,length(miceopt));

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
       
            CP = Perf{midx,1,twidx}.CP;
            CP = reshape(CP,xpix,ypix);
%             minntrials(midx) = min(Perf{midx,1,twidx}.nrright,Perf{midx,1,twidx}.nrleft);
%             nrleft = Perf{midx,1,twidx}.nrleft;
%             nrright = Perf{midx,1,twidx}.nrright;
            nrhit = Perf{midx,1,twidx}.nrhit;
            nrerror = Perf{midx,1,twidx}.nrerror;
        
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
        
        avgtmpbilateral(midx,twidx,roiidx) = nanmean(CP(maskbilateral==roiidx));
        BilateralAVG{midx,twidx}(maskbilateral==roiidx) = avgtmpbilateral(midx,twidx,roiidx);
        
      
        
        C = linspecer(length(regio2take));
        figure(plotfig)
        subplot(3,2,midx)
        plot(1:length(TW),avgtmpbilateral(midx,:,roiidx),'color',C(roiidx,:),'linewidth',2)
        title([miceopt{midx} ' nr hit ' num2str(nrhit) ', nr nrerror ' num2str(nrerror)]) %num2str(minntrials(midx))
        xlabel('Timewindow, fix xticks')
        ylabel('AUC')
        hold on
        
        
    end
legend(BrainModel{midx}.Model.Rnames{regio2take})

figure(imagescfig);
for twidx = 1:length(TW)
    subplot(length(miceopt),length(TW),(midx-1)*length(TW)+twidx);
    j = imagesc(SideAreaAVG{midx,twidx},[0.1 0.9]);
    colormap(ActSupColorMap)
    set(j,'AlphaData',~isnan(SideAreaAVG{midx,twidx}));
       
end
end


y = squeeze(nanmean(avgtmpbilateral,1));
x = 1:length(TW);
errorrr = squeeze(std(avgtmpbilateral,1));

 figure
for regionidx = 1:length(regio2take)
    errorbar(x,y(:,regionidx),errorrr(:,regionidx),'color',C(regionidx,:),'linewidth',3)
    ylim([0.3,0.7])
    hold on
end
legend(BrainModel{midx}.Model.Rnames{regio2take})


% just take a couple of regions
figure
regions = [11,25,26,38];
reg = find(ismember(regio2take,regions));
for regionidx = reg
    errorbar(x,y(:,regionidx),errorrr(:,regionidx),'color',C(regionidx,:),'linewidth',3)
    ylim([0.3,0.7])
    hold on
end
legend(BrainModel{midx}.Model.Rnames{regions})


%plot bargraph
yy = y(:,ismember(regio2take,regions));
erry = errorrr(:,ismember(regio2take,regions));
h = barwitherr(erry,yy);
xlabel('TimeWindow')
ylabel('AUC')
legend(BrainModel{midx}.Model.Rnames{regions})
ylim([0.3,0.7])
set(gca,'XTickLabel',{num2str(TW{1}),num2str(TW{2}),num2str(TW{3}),num2str(TW{4}),num2str(TW{5}),num2str(TW{6})})
set(h(1),'FaceColor',C(ismember(regio2take,regions(1)),:));
set(h(2),'FaceColor',C(ismember(regio2take,regions(2)),:));
set(h(3),'FaceColor',C(ismember(regio2take,regions(3)),:));
set(h(4),'FaceColor',C(ismember(regio2take,regions(4)),:));

xlim=get(gca,'xlim');
hold on
plot(xlim,[0.5,0.5])
% set(gca,'FaceColor,',C(ismember(regio2take,regions)))



% 
% %% plot bargraph
% 
%       
% 
% 
% 
% 
%     % plot bargraph with errors
%     
% 
%         Y = [11 14 13;
%             15 12 16];
%         E = [ 3  4  2;
%             5  2  3];
% 
%         % Plot with bar
%         subplot(1, 2, 1);
%         text(0.5, 0.5, 'not implemented', 'HorizontalAlignment', 'center');
%         title('bar');
% 
%         % Plot with superbar
%         subplot(1, 2, 2);
%         superbar(Y, 'E', E);
%         title('superbar');
% 
%         % Since one does not know the X locations of the individual bars, BAR and
%         % ERRORBAR can not be used together to add error bars to grouped bar plots.
% 
% 
%         set(gca, 'XTick', []); % place vector
%         %         set(gca, 'XTick', [1 2])
%         %         set(gca, 'XTickLabel', {'Model1' 'Model2'})
% 
%         set(gca, 'YTick', []);
% 
%   
%     
%     
%  
% 
% %% avg across mice
% % ciplot lower upper x
% % avgmice = squeeze(nanmean(avgtmp,1));
% % stdmice = squeeze(std(avgtmp));
% % upperlim = avgmice+stdmice;
% % lowerlim = avgmice-stdmice;
% % 
% % 
% % avgmap = cellfun(nanmean,allmicetmp{
% % test = nanmean(avgmap,1)
% 
%  
% avgallmice = squeeze(nanmean(avgtmpperside,1));
% stdallmice = std(avgtmpperside,1,1);
% collorross = {'g','r','b','y','m'};
% 
% figure
% for i = 1:size(avgallmice,2)
% shadedErrorBar(1:size(avgallmice,1),avgallmice(:,i),std(avgtmpperside(:,:,i),1,1),collorross{i},1)
% hold on
% end
% 
% 
% figure;
% hold on
%  areaidx = 1
%   ciplot(lowerlim(:,areaidx),upperlim(:,areaidx),avgmice(:,areaidx))
%   hold on
% 
% plot(lowerlim)
% 
% % Examples
% y=randn(30,80); x=1:size(y,2);
% shadedErrorBar(x,mean(y,1),std(y),'g');
% shadedErrorBar(x,y,{@median,@std},{'r-o','markerfacecolor','r'});    
% shadedErrorBar([],y,{@median,@std},{'r-o','markerfacecolor','r'});    
% 
% % Overlay two transparent lines
% y=randn(30,80)*10; x=(1:size(y,2))-40;
% shadedErrorBar(x,y,{@mean,@std},'-r',1); 
% hold on
% y=ones(30,1)*x; y=y+0.06*y.^2+randn(size(y))*10;
% shadedErrorBar(x,y,{@mean,@std},'-b',1);   
% hold off
% 
% 








