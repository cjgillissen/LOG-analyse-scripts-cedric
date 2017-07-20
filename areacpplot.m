miceopt = {'Alladin','Chief','Esmeralda','Frey'}%%% %options for mice
basel = [-300 -50];
visint = [50 200];
vistw = [200 500]; %FOR FG
Delayearlytw = [600 1350];
Delaylatetw = [1400 1600];
Responstw = [1900 2200];
TW = {basel,visint,vistw,Delayearlytw,Responstw,Delaylatetw};
StorePath = '\\vcnin\mouse_working_memory\MainAnalysis\';
xpix = 400;
ypix = 400;
scalefct = 0.5;
%%
avgcp = 1:length(TW);
figure('name','Lineplot average cp per area per tw')

for midx = 1:length(miceopt)
    mouse = miceopt{midx};
    BrainModel{midx} = load(fullfile(StorePath,mouse,'brainareamodel.mat'));
        
        rgnames = {'V1','Vpor','Va','M1','M2'};
        regio2take = [11,13,25,38,40];
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
            
        end
        plot(1:length(TW),avgcp)
        hold on
        end
    legend(rgnames)
end
