% Specific WM analysis script

Scriptsdir = 'C:\Users\gillissen\Documents\GitHub\Mouse\ImagingAnalysis' %Direction of scripts
storepath = '\\vcnin\mouse_working_memory\MainAnalysis\' %Main storage --> Will contain all processed data & model of the brain
DataDirectory = '\\vcnin\mouse_working_memory\Imaging\'%Raw images.
miceopt = {'Alladin','Chief','Esmeralda','Frey'}%%% %options for mice
Stim2Check = 'DelayedOriTuningSound'%Name of the stimulus as written in the LOG-file
fijiLoc = 'C:\Users\beest\Fiji.app\scripts' %FijiLocation
tempstorage = 'E:\Imaging\RegisteredImages\'%Temporary storage of 
baselinemethod = 1 %1: trial by trial baseline (no detrending necessary per se), 2: filtered baseline (F0), no trial specific baseline (basically F), 3: average baseline (averaged over all trials per condition and then averaged over conditions, after detrending)
MappingDir = '\\vcnin\mouse_working_memory\pRF results\'
takeequalsample = 0;
smoothfact = 2;
RedoAll = 0;
resampling = 2; % 0 no resampling, 1 resample conditions within left vs right, 2 resample all conditions
trialtypes = {'1500'};
%ACHTUNG: you need to install: https://sites.google.com/site/qingzongtseng/template-matching-ij-plugin#downloads

addpath(genpath(Scriptsdir))
global UserQuestions %Defines whether gui's show up to ask what to do with existing datasets. If 0, existing datasets are used and not overwritten
UserQuestions =0;
 
%% Load in data --> 
info = BuildInfo(miceopt,DataDirectory,storepath,Stim2Check); %Get all logs of WM-imaging sessions
save(fullfile('C:\Users\gillissen\Desktop\Figures NC','sessionstruct'),'info')
%% Register data
RegisterSessions(info,storepath,miceopt,fijiLoc,tempstorage,MappingDir,0) %Registeres all images of different sessions of one mouse, including a mapping session, if one available

%% Load Imaging Data
LoadData(info,storepath,miceopt,Stim2Check,3,0) %Loads all data, and saves in cells per condition per mouse, all aligned. Last input: smoothfactor 2d spatial smoothing

%% Load Alan Brain Map onto Mice brains
AllenBrainMappRF(info,storepath,miceopt,DataDirectory,MappingDir,0) %applies a map to the brain of the mouse

%% calculate baseline (for all mice/sessions in info)
BaselineCalculation(info,miceopt,storepath,DataDirectory,Stim2Check,smoothfact,0) %Last factor: trial smoothing for baseline 'fit'

%% Makes Averages + other necessary variables per condition (LEFTVSRIGHT
%structure comes out) %Repeat this script for every 'trialtype' if you have
%different trialtypes
ConditionAveragingScriptWF(info,miceopt,storepath,DataDirectory,Stim2Check,baselinemethod,0,[-300 2500],smoothfact,takeequalsample,0); 
ConditionAveragingScriptWF(info,miceopt,storepath,DataDirectory,Stim2Check,baselinemethod,1500,[-300 2500],smoothfact,takeequalsample,0); 

%% Differences
DifferencesCheck(info,miceopt,storepath,DataDirectory,Stim2Check,baselinemethod,[-300 2500],{'750','1500'},[-300 2500],takeequalsample)
Noisecorrelationrawdata(info,miceopt,storepath,Stim2Check,baselinemethod,trialtypes,takeequalsample)
SPforWF(info,miceopt,storepath,Stim2Check,baselinemethod,trialtypes,takeequalsample,resampling)

%% MVPA 
% MVPAPerSession(info,miceopt,storepath,DataDirectory,Stim2Check,baselinemethod,[-300 1500],{'1500'},[-300 2500],takeequalsample); 
MVPAAllSessions(info,miceopt,storepath,Stim2Check,baselinemethod,{'1500'},{'Error','Hit'},takeequalsample); 

%% Modulation L-R
% FGSpecificAna(info,miceopt,storepath,DataDirectory,Stim2Check,baselinemethod,{'1500'},[-200 1800],smoothfact,0)

%% Effective Pixelsize
% CalcEffPixSize(info,miceopt,StorePath,DataDirectory,Stim2Check,baselinemethod,Scriptsdir)
% %% BaselineDiff  
% See whether there are pre-baseline differences between different stimuli
% BaselineDifference(info,miceopt,storepath,DataDirectory,Stim2Check,baselinemethod) %only works with baselinemethod 2 or 3

% %% Specific Analysis
% WMSpecificAnalysis(info,miceopt,storepath,DataDirectory,Stim2Check,baselinemethod)

%% Z-score analysis
% ZScoreAna(info,miceopt,storepath,DataDirectory,Stim2Check) %It actually is a T-score. but, yeah.. 
%% Motion Separation (ICA)
% MotionSeparator(info,miceopt,storepath,DataDirectory,Stim2Check)

 %% MVPA - Do on server
% MVPA_Decoder(info,miceopt,storepath,DataDirectory,Stim2Check,baselinemethod,Scriptsdir) %First time is same as before. 2nd from when to when prediction should be
%% Functional correlation
% FunctionalCorrelationAnalyses(info,miceopt,storepath,DataDirectory,Stim2Check,baselinemethod,'C:\Users\beest\Documents\MATLAB\fdaMatlab\fdaM');

