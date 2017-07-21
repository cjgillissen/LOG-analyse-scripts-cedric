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


% rgnames = {'V1','Vpor','Va','M1','M2'};
% regio2take = 10:43;
regio2take = 11;
avgcp = 1:length(TW);
figure('name','Lineplot average cp per area per tw')
avgtmp = nan(length(miceopt),length(TW),length(rgnames));
for midx = 1:length(miceopt)
    mouse = miceopt{midx};
    BrainModel{midx} = load(fullfile(StorePath,mouse,'brainareamodel.mat'));
    masktmp = zeros(xpix,ypix);
        for roiidx = 1:length(regio2take)
            disp(['Mouse ' num2str(midx) ', region ' num2str(roiidx) ' of ' num2str(length(regio2take))])
            Borders = BrainModel{midx}.Model.Boundaries{regio2take(roiidx)};
            masktmp = zeros(xpix,ypix);
        %         mask = [];r
            for roi2dx = 1:length(Borders)
                mask = poly2mask(Borders{roi2dx}(:,1).*scalefct,Borders{roi2dx}(:,2).*scalefct,xpix,ypix);
                %Shrink to not have border effects
                mask = bwmorph(mask,'shrink',1);
                masktmp(mask)=roiidx;
            end

        for twidx = 1:length(TW)
            CP = Perf{midx,1,twidx}.CP;
            CP = reshape(CP,xpix,ypix);
            avgcp(twidx) = nanmean(CP(masktmp==roiidx));
            avgtmp(midx,twidx,roiidx) = nanmean(CP(masktmp==roiidx));
            allmicetmp{midx,twidx} = CP;
        end
    subplot(3,2,midx)
    plot(1:length(TW),avgcp)
    title([miceopt{midx}])
    xlabel('Timewindow, fix xticks')
    ylabel('AUC')
    hold on
    end
legend(BrainModel{midx}.Model.Rnames{regio2take})
end

%% avg across mice
% ciplot lower upper x
% avgmice = squeeze(nanmean(avgtmp,1));
% stdmice = squeeze(std(avgtmp));
% upperlim = avgmice+stdmice;
% lowerlim = avgmice-stdmice;


avgmap = cellfun(nanmean,allmicetmp{
test = nanmean(avgmap,1)







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










