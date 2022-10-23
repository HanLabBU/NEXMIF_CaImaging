clear all
close all
clc

% initialize
miceKO = {'2917','4277','605842','605887','605926','606717','611114','611312'};
miceWT = {'2978','2089','605890','611115','612265','612266','612270','612668'}; %'2085',
miceStudy = [miceKO, miceWT];
maskWT = ismember(miceStudy,miceWT);
maskKO = ismember(miceStudy,miceKO);

% pathData = 'U:\eng_research_handata\Rebecca_Mount_1';
pathData = 'U:\eng_research_handata\eng_research_handata2\Rebecca_Mount\Autism';
dataFolder = 'U:\eng_research_handata\Athif Mohamed\nexmif_paper\Data_new2_cleaned_ext2';
savePath = 'U:\eng_research_handata\Athif Mohamed\nexmif_paper\NEXMIF_CaImaging\Figures';


% Limits to data 
badSessions = {'2089_D3','611115_D3','612266_D1','612266_D3','612266_D5'};
badSpeedSessions = {'605890_D1','611115_D1','611115_D3','605926_D1','611114_D1','606717_D5'};
fewMovSessions = {'4277_D5','612270_D1','612270_D5','605890_D3','605890_D5'};
fewMovSessionsM = {'4277','612270','612270','605890','605890'};
fewMovSessionsD = [5,1,5,3,5];

% load('blues and greens.mat')

% load('redpurple.mat')
% colors = rp/255;

% load('bluegreen.mat')
% colors = bg/255;

load('redblue.mat')
colors = rb;

Fs = 20;
dIdx = [1,3,5];
conditionList = {'Ball','Platform','Tonepuff'};
timingList = {'Onset','Offset'};
mType = {'WT','KO'};

% correlation 
splabel = {'Low speed','High speed'};
sh = [0.5 0.8];
corrNames = {'Jaccard Index',...
    'Asymmetric Correlation',...
    'Pearson Correlation'};
corrType = {'Jaccard','Asymmetric', 'Pearson'};

% plotting 
set(0, 'DefaultFigureRenderer', 'painters');