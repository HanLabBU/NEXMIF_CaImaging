clear all
close all
clc

addpath(genpath('U:\eng_research_handata\Athif Mohamed\nexmif_paper\Utils\'))  % Add utilities
init % Initialize data directories and genotypes
close all

mouseIDList = {'605842','612265'}; %WT, KO
dayList = [1,1];
% roi idx to plot:

for pp = 1:2
    mouseID = mouseIDList{pp};
    day = dayList(pp);
    
    mPath = fullfile(dataFolder,sprintf('fullData_%s_D%i_%s',mouseID,day,'ball'));
    
    load(mPath) % load fullData
    cellList = fullData.goodIdx;
    if pp == 1  % 605842
        cellList = [44,66,69,111,141,153,157,163,182,188,207,247,267,269,350];
        color = colors(4,:);
    elseif pp == 2 % 612665
        cellList = [37,38,47,48,63,73,83,93,100,118,131,150,160,173,202];
        color = colors(2,:);
    end
    % cellList = cellList(randperm(numel(cellList),min(numel(cellList),15)));
    
    %% Plot raster with events detected
    
    [fig1] = rasterPlot_ext8(fullData,cellList,color); % plots selected traces, 15
    
    set(gcf,'Color','none')
    set(gca,'Units','inches','InnerPosition',[.3 .3 2 2],'TickDir','out','Color','none','Box','off')
    
    saveas(fig1,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\b_WTKO_traces\',sprintf('traces-%s-%d.png',mouseID,day)))
    saveas(fig1,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\b_WTKO_traces\',sprintf('traces-%s-%d.fig',mouseID,day)))
    saveas(fig1,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\b_WTKO_traces\',sprintf('traces-%s-%d.epsc',mouseID,day)))
end