
%% section 1

cd('\\VCNIN\mouse_working_memory\TrainingLogs')


currentdelay = [];
rewcount = [];
valvedur = [];
guidurleft = [];
guidurright = [];
setup = [];
Side = [];
mice = {};
gavepassive = [];
lefthit = [];
righthit = [];
Reaction = [];
passivefirst = [];


load('Frey_20170314_B1.mat')
CheckReactions('Frey_20170314_B1.mat');

 Side = [Side LOG.Side];
 currentdelay = [currentdelay LOG.currentdelay];
 valvedur = [valvedur LOG.ValveOpenTime];
 setup = [setup repmat(LOG.Setup,[1 size(LOG.ActualRewPeriod,2)])];
 mice = {mice,repmat({LOG.Mouse},[1 size(LOG.ActualRewPeriod,2)])};
 gavepassive = [gavepassive LOG.Gavepassive];
 guidurleft = [guidurleft LOG.rdurleft];
 guidurright = [guidurright LOG.rdurright];
 lefthit = [lefthit strcmp(LOG.Reaction,'Hit')&strcmp(LOG.Side,'left')];
 righthit = [righthit strcmp(LOG.Reaction,'Hit')&strcmp(LOG.Side,'right')];
 Reaction = [Reaction LOG.correctReaction];
 passivefirst = [passivefirst LOG.Passivefirst];
 
 load('Alladin_20170314_B1.mat')
 CheckReactions('Alladin_20170314_B1.mat');

 
Side = [Side LOG.Side];
 currentdelay = [currentdelay LOG.currentdelay];
 valvedur = [valvedur LOG.ValveOpenTime];
 setup = [setup repmat(LOG.Setup,[1 size(LOG.ValveOpenTime,2)])];
 mice = {mice,repmat({LOG.Mouse},[1 size(LOG.ActualRewPeriod,2)])};
 gavepassive = [gavepassive LOG.Gavepassive];
 guidurleft = [guidurleft LOG.rdurleft];
 guidurright = [guidurright LOG.rdurright];
 lefthit = [lefthit strcmp(LOG.Reaction,'Hit')&strcmp(LOG.Side,'left')];
 righthit = [righthit strcmp(LOG.Reaction,'Hit')&strcmp(LOG.Side,'right')];
 Reaction = [Reaction LOG.correctReaction];
  passivefirst = [passivefirst LOG.Passivefirst];

 
 load('Chief_20170314_B1.mat')
 CheckReactions('Chief_20170314_B1.mat');

 Side = [Side LOG.Side];
 currentdelay = [currentdelay LOG.currentdelay];
 valvedur = [valvedur LOG.ValveOpenTime];
 setup = [setup repmat(LOG.Setup,[1 size(LOG.ActualRewPeriod,2)])];
 mice = {mice,repmat({LOG.Mouse},[1 size(LOG.ActualRewPeriod,2)])};
 gavepassive = [gavepassive LOG.Gavepassive];
 guidurleft = [guidurleft LOG.rdurleft];
 guidurright = [guidurright LOG.rdurright];
 lefthit = [lefthit strcmp(LOG.Reaction,'Hit')&strcmp(LOG.Side,'left')];
 righthit = [righthit strcmp(LOG.Reaction,'Hit')&strcmp(LOG.Side,'right')];
 Reaction = [Reaction LOG.correctReaction];
  passivefirst = [passivefirst LOG.Passivefirst];


 
 load('Chief_20170313_B1.mat')
 CheckReactions('Chief_20170313_B1.mat');

 
 Side = [Side LOG.Side];
 currentdelay = [currentdelay LOG.currentdelay];
 valvedur = [valvedur LOG.ValveOpenTime];
 setup = [setup repmat(LOG.Setup,[1 size(LOG.ActualRewPeriod,2)])];
 mice = {mice,repmat({LOG.Mouse},[1 size(LOG.ActualRewPeriod,2)])};
 gavepassive = [gavepassive LOG.Gavepassive];
 guidurleft = [guidurleft LOG.rdurleft];
 guidurright = [guidurright LOG.rdurright];
 lefthit = [lefthit strcmp(LOG.Reaction,'Hit')&strcmp(LOG.Side,'left')];
 righthit = [righthit strcmp(LOG.Reaction,'Hit')&strcmp(LOG.Side,'right')];
 Reaction = [Reaction LOG.correctReaction];
   passivefirst = [passivefirst LOG.Passivefirst];


 
 
 load('Alladin_20170313_B1.mat')
 CheckReactions('Alladin_20170313_B1.mat');

 
 Side = [Side LOG.Side];
 currentdelay = [currentdelay LOG.currentdelay];
 valvedur = [valvedur LOG.ValveOpenTime];
 setup = [setup repmat(LOG.Setup,[1 size(LOG.ActualRewPeriod,2)])];
 mice = {mice,repmat({LOG.Mouse},[1 size(LOG.ActualRewPeriod,2)])};
 gavepassive = [gavepassive LOG.Gavepassive];
 guidurleft = [guidurleft LOG.rdurleft];
 guidurright = [guidurright LOG.rdurright];
 lefthit = [lefthit strcmp(LOG.Reaction,'Hit')&strcmp(LOG.Side,'left')];
 righthit = [righthit strcmp(LOG.Reaction,'Hit')&strcmp(LOG.Side,'right')];
 Reaction = [Reaction LOG.correctReaction];
   passivefirst = [passivefirst LOG.Passivefirst];


 
 load('Frey_20170313_B1.mat')
 CheckReactions('Frey_20170313_B1.mat');

 
 Side = [Side LOG.Side];
 currentdelay = [currentdelay LOG.currentdelay];
 valvedur = [valvedur LOG.ValveOpenTime];
 setup = [setup repmat(LOG.Setup,[1 size(LOG.ActualRewPeriod,2)])];
 mice = {mice,repmat({LOG.Mouse},[1 size(LOG.ActualRewPeriod,2)])};
 gavepassive = [gavepassive LOG.Gavepassive];
 guidurleft = [guidurleft LOG.rdurleft];
 guidurright = [guidurright LOG.rdurright];
 lefthit = [lefthit strcmp(LOG.Reaction,'Hit')&strcmp(LOG.Side,'left')];
 righthit = [righthit strcmp(LOG.Reaction,'Hit')&strcmp(LOG.Side,'right')];
 Reaction = [Reaction LOG.correctReaction];
 passivefirst = [passivefirst LOG.Passivefirst];

 guiL = guidurright;
 guiR = guidurleft;

 guidur = guiR;
 guidur(strcmp(Side,'left'))= guiL(strcmp(Side,'left'));
 

 
 
 %% section 3 
 %compute difference with gui
 
 
 setupopt = [1,3,5];
 sideopt = {'left','right'};
 
     for j = setupopt
       for k = sideopt
         figure;
         histogram(guidur(strcmp(Side,k)&setup==j&~gavepassive&strcmp(Reaction,'Hit')&~passivefirst)-(valvedur(strcmp(Side,k)&setup==j&~gavepassive&strcmp(Reaction,'Hit')&~passivefirst)),0:1:20);
         xlim([0 20])
         title(sprintf('setup %1.0f side %s',j,cell2mat(k)))
         ylabel('Frequency')
         xlabel('Gui value - measured valve opening duration ')
         saveas(gcf,fullfile('C:\Users\gillissen\Desktop\InternshipCédric',sprintf('%s.jpg',get(get(gca,'title'),'string'))))
 
      end
     end
 
 
     %% check passives
     
     
 
 
 
 