
cd('C:\Users\gillissen\Desktop\TrainingLogs')


currentdelay = [];
rewcount = [];
valvedur = [];
setup = [];
Side = [];
mice = {};
gavepassive = [];

load('Frey_20170306_B1.mat')

 Side = [Side LOG.Side];
 currentdelay = [currentdelay LOG.currentdelay];
 valvedur = [valvedur LOG.ValveOpenTime];
 setup = [setup repmat(LOG.Setup,[1 size(LOG.ActualRewPeriod,2)])];
 mice = {mice repmat(LOG.Mouse,[1 size(LOG.ActualRewPeriod,2)])};
 gavepassive = [gavepassive LOG.Gavepassive];
 
 load('Frey_20170303_B1.mat')
 
 
Side = [Side LOG.Side];
 currentdelay = [currentdelay LOG.currentdelay];
 valvedur = [valvedur LOG.ValveOpenTime];
 setup = [setup repmat(LOG.Setup,[1 size(LOG.ValveOpenTime,2)])];
 mice = {mice repmat(LOG.Mouse,[1 size(LOG.ActualRewPeriod,2)])};
 gavepassive = [gavepassive LOG.Gavepassive];

 
 load('Chief_20170306_B1.mat')
 
 Side = [Side LOG.Side];
 currentdelay = [currentdelay LOG.currentdelay];
 valvedur = [valvedur LOG.ValveOpenTime];
 setup = [setup repmat(LOG.Setup,[1 size(LOG.ActualRewPeriod,2)])];
 mice = {mice repmat(LOG.Mouse,[1 size(LOG.ActualRewPeriod,2)])};
 gavepassive = [gavepassive LOG.Gavepassive];

 
 
 load('Chief_20170303_B1.mat')
 
 Side = [Side LOG.Side];
 currentdelay = [currentdelay LOG.currentdelay];
 valvedur = [valvedur LOG.ValveOpenTime];
 setup = [setup repmat(LOG.Setup,[1 size(LOG.ActualRewPeriod,2)])];
 mice = {mice repmat(LOG.Mouse,[1 size(LOG.ActualRewPeriod,2)])};
 gavepassive = [gavepassive LOG.Gavepassive];

 load('Alladin_20170306_B1.mat')
 
 Side = [Side LOG.Side];
 currentdelay = [currentdelay LOG.currentdelay];
 valvedur = [valvedur LOG.ValveOpenTime];
 setup = [setup repmat(LOG.Setup,[1 size(LOG.ActualRewPeriod,2)])];
 mice = {mice repmat(LOG.Mouse,[1 size(LOG.ActualRewPeriod,2)])};
 gavepassive = [gavepassive LOG.Gavepassive];

 
 load('Alladin_20170303_B1.mat')
 
 Side = [Side LOG.Side];
 currentdelay = [currentdelay LOG.currentdelay];
 valvedur = [valvedur LOG.ValveOpenTime];
 setup = [setup repmat(LOG.Setup,[1 size(LOG.ActualRewPeriod,2)])];
 mice = {mice repmat(LOG.Mouse,[1 size(LOG.ActualRewPeriod,2)])};
 gavepassive = [gavepassive LOG.Gavepassive];

 
 load('Frey_20170307_B1.mat')
 
 Side = [Side LOG.Side];
 currentdelay = [currentdelay LOG.currentdelay];
 valvedur = [valvedur LOG.ValveOpenTime];
 setup = [setup repmat(LOG.Setup,[1 size(LOG.ActualRewPeriod,2)])];
 mice = {mice repmat(LOG.Mouse,[1 size(LOG.ActualRewPeriod,2)])};
 gavepassive = [gavepassive LOG.Gavepassive];

 
 
 load('Chief_20170307_B1.mat')
 
 Side = [Side LOG.Side];
 currentdelay = [currentdelay LOG.currentdelay];
 valvedur = [valvedur LOG.ValveOpenTime];
 setup = [setup repmat(LOG.Setup,[1 size(LOG.ActualRewPeriod,2)])];
 mice = {mice repmat(LOG.Mouse,[1 size(LOG.ActualRewPeriod,2)])};
 gavepassive = [gavepassive LOG.Gavepassive];

 
 load('Alladin_20170307_B1.mat')
 
 Side = [Side LOG.Side];
 currentdelay = [currentdelay LOG.currentdelay];
 valvedur = [valvedur LOG.ValveOpenTime];
 setup = [setup repmat(LOG.Setup,[1 size(LOG.ActualRewPeriod,2)])];
 mice = {mice repmat(LOG.Mouse,[1 size(LOG.ActualRewPeriod,2)])};
 gavepassive = [gavepassive LOG.Gavepassive];

 
 
 %% 
 
 
%  values1 = unique(valvedur(setup==1));
%  values1(isnan(values1)) = [];
%  edges1 = values1(1):5:values1(end);
%  
%  figure
%  
%  
%  
%  [S1, Edges1] = histcounts(valvedur(setup==1),30);
%  bar(S1)
%  ylabel('Frequency')
%  xlabel('ms valve is open, setup 1')
%  xticklabels(Edges1(1:5:end))
%  
% 
%  
%  [S3, Edges3] = histcounts(valvedur(setup==3),30);
%  bar(S3)
%  ylabel('Frequency')
%  xlabel('ms valve is open setup 3')
%  xticklabels(Edges3(1:5:end))
%  
%  
%  
%  [S5, Edges5] = histcounts(valvedur(setup==5),30);
%  bar(S5)
%  
%  ylabel('Frequency')
%  xlabel('ms valve is open setup 5')
%  xticklabels(Edges5(1:5:end))
%  
%  
%  %side left setup 1
%  [S1L, Edges1] = histcounts(valvedur(strcmp(Side,'left')&setup==1),30);
%  bar(S1L)
%  title('setup 1 side left')
%  ylabel('Frequency')
%  xlabel('ms valve is open, setup 1 side left')
%  xticklabels(Edges1(1:5:end))
%  
%  %side right setup 1
%  
%  [S1R, Edges1] = histcounts(valvedur(strcmp(Side,'right')&setup==1),30);
%  bar(S1R)
%  title('setup 1 side right')
%  ylabel('Frequency')
%  xlabel('ms valve is open')
%   xticklabels(Edges1(1:5:end))
%   
%   %setup 3 side left
%   
%   [S3L, Edges3] = histcounts(valvedur(strcmp(Side,'left')&setup==3),30);
%  bar(S3L)
%  title('setup 3 side left')
%  ylabel('Frequency')
%  xlabel('ms valve is open ')
%  xticklabels(Edges3(1:5:end))
%  
%  [S3R, Edges3] = histcounts(valvedur(strcmp(Side,'right')&setup==3),30);
%  bar(S3R)
%  title('setup 3 side right')
%  ylabel('Frequency')
%  xlabel('ms valve is open ')
%  xticklabels(Edges3(1:5:end))
%  
%  
%  [S5R, Edges3] = histcounts(valvedur(strcmp(Side,'right')&setup==5),30);
%  bar(S5R)
%  title('setup 5 side right')
%  ylabel('Frequency')
%  xlabel('ms valve is open ')
%  xticklabels(Edges3(1:5:end))
%  
%  [S5L, Edges3] = histcounts(valvedur(strcmp(Side,'left')&setup==5),30);
%  bar(S5L)
%  title('setup 5 side left')
%  ylabel('Frequency')
%  xlabel('ms valve is open ')
%  xticklabels(Edges3(1:5:end))

 
%  
 %% use histogram
 
 figure
 histogram(valvedur(strcmp(Side,'right')&setup==3&~gavepassive),0:5:400);
 xlim([180 400])
 title('setup 3 side right')
 ylabel('Frequency')
 xlabel('ms valve is open ')
 
 figure
 histogram(valvedur(strcmp(Side,'left')&setup==1&~gavepassive),0:5:400);
 xlim([180 400])
 title('setup 1 side left')
 ylabel('Frequency')
 xlabel('ms valve is open')
 
 figure
 histogram(valvedur(strcmp(Side,'right')&setup==1&~gavepassive),0:5:400);
xlim([180 400])
 title('setup 1 side right')
 ylabel('Frequency')
 xlabel('ms valve is open')
 
 figure
 histogram(valvedur(strcmp(Side,'left')&setup==3&~gavepassive),0:5:400);
 xlim([180 400])
 title('setup 3 side left')
 ylabel('Frequency')
 xlabel('ms valve is open')
 
 figure
 histogram(valvedur(strcmp(Side,'left')&setup==5&~gavepassive),0:5:400);
 xlim([180 400])
 title('setup 5 side left')
 ylabel('Frequency')
 xlabel('ms valve is open')
 
 figure
 histogram(valvedur(strcmp(Side,'right')&setup==5&~gavepassive),0:5:400);
 xlim([180 400])
 title('setup 5 side right')
 ylabel('Frequency')
 xlabel('ms valve is open')
 
 
 
 