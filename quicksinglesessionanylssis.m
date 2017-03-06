
cd('C:\Users\gillissen\Desktop\TrainingLogs')


currentdelay = [];
gavepassive = [];
rewcount = [];
valvedur = [];
setup = [];
Side = {};
mice = {};

load('Frey_20170306_B1.mat')

 Side = {Side LOG.Side};
 currentdelay = [currentdelay LOG.currentdelay];
 valvedur = [valvedur LOG.ValveOpenTime];
 setup = [setup repmat(LOG.Setup,[1 size(LOG.ActualRewPeriod,2)])];
 mice = {mice repmat(LOG.Mouse,[1 size(LOG.ActualRewPeriod,2)])};
 
 load('Frey_20170303_B1.mat')
 
 
Side = {Side LOG.Side};
 currentdelay = [currentdelay LOG.currentdelay];
 valvedur = [valvedur LOG.ValveOpenTime];
 setup = [setup repmat(LOG.Setup,[1 size(LOG.ValveOpenTime,2)])];
 mice = {mice repmat(LOG.Mouse,[1 size(LOG.ActualRewPeriod,2)])};
 
 load('Chief_20170306_B1.mat')
 
 Side = {Side LOG.Side};
 currentdelay = [currentdelay LOG.currentdelay];
 valvedur = [valvedur LOG.ValveOpenTime];
 setup = [setup repmat(LOG.Setup,[1 size(LOG.ActualRewPeriod,2)])];
 mice = {mice repmat(LOG.Mouse,[1 size(LOG.ActualRewPeriod,2)])};
 
 
 load('Chief_20170303_B1.mat')
 
 Side = {Side LOG.Side};
 currentdelay = [currentdelay LOG.currentdelay];
 valvedur = [valvedur LOG.ValveOpenTime];
 setup = [setup repmat(LOG.Setup,[1 size(LOG.ActualRewPeriod,2)])];
 mice = {mice repmat(LOG.Mouse,[1 size(LOG.ActualRewPeriod,2)])};
 
 load('Alladin_20170306_B1.mat')
 
 Side = {Side LOG.Side};
 currentdelay = [currentdelay LOG.currentdelay];
 valvedur = [valvedur LOG.ValveOpenTime];
 setup = [setup repmat(LOG.Setup,[1 size(LOG.ActualRewPeriod,2)])];
 mice = {mice repmat(LOG.Mouse,[1 size(LOG.ActualRewPeriod,2)])};
 
 load('Alladin_20170303_B1.mat')
 
 Side = {Side LOG.Side};
 currentdelay = [currentdelay LOG.currentdelay];
 valvedur = [valvedur LOG.ValveOpenTime];
 setup = [setup repmat(LOG.Setup,[1 size(LOG.ActualRewPeriod,2)])];
 mice = {mice repmat(LOG.Mouse,[1 size(LOG.ActualRewPeriod,2)])};
 
 
 values1 = unique(valvedur(setup==1));
 values1(isnan(values1)) = [];
 
 valvehistsetup1 = histc(valvedur(setup==1),values1(~isnan(values)));
 
 