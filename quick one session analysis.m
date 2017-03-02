RTleft=[];RTright=[]; %Pre allocate Licks
        ResponseCell = {}; %Preallocate for ResponseCell concatenation
        Side = {};
        Phase = {};
currentdelay = [];
rewcount = [];
actrewper = [];

load('Frey_20170301_B1')

 RTleft = [RTleft LOG.RTleftVec];
 RTright = [RTright LOG.RTrightVec]; 
 ResponseCell = [ResponseCell LOG.correctReaction];
 Side = [Side LOG.Side];
 Phase = [Phase LOG.CurrentPhase];
 currentdelay = [currentdelay LOG.currentdelay];
 rewcount = LOG.RewardCount;
actrewper = LOG.ActualRewPeriod;
 
 load('Frey_20170301_B2')
 
 RTleft = [RTleft LOG.RTleftVec];
 RTright = [RTright LOG.RTrightVec]; 
 ResponseCell = [ResponseCell LOG.correctReaction];
 Side = [Side LOG.Side];
 Phase = [Phase LOG.CurrentPhase];
 currentdelay = [currentdelay LOG.currentdelay];
 rewcount = [rewcount LOG.RewardCount];
actrewper = [actrewper LOG.ActualRewPeriod];

hits = sum(strcmp(ResponseCell,'Hit'));
errors = sum(strcmp(ResponseCell,'Error'));
