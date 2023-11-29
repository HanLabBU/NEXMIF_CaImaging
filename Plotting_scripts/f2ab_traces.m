clear all
close all
clc

addpath(genpath('U:\eng_research_handata\Athif Mohamed\nexmif_paper\Utils\'))  % Add utilities
init % Initialize data directories and genotypes

mouseIDList = {'605842','612265'}; % 2.15.23: using 612265 day 1 (WT) and 605842 day 1 (KO) (same as first version)
dayList = [1,1];
% roi idx to plot:

for pp = 1:2
    mouseID = mouseIDList{pp};
    day = dayList(pp);
    
    mPath = fullfile(dataFolder,sprintf('fullData_%s_D%i_%s',mouseID,day,'ball'));
    
    load(mPath) % load fullData
    cellList = fullData.goodIdx;
    if pp == 1  % 605842
%         cellList = [44,66,69,111,141,153,157,163,182,188,207,247,267,269,350]; % included in first submission. event rate is too low compared to 2978
        cellList = [38,44,59,63,66,88,92,97,153,157,182,207,247,267,373];
        color = colors(4,:);
    elseif pp == 2  % 612265            % first version: (too many events, alternate to first version) 612665 d1 cellList = [41,48,100,118,160,170,173,179,194,213,220,225,248,298,300]; (too many events) 2978 d3, cellList = [25,28,35,37,61,70,75,76,95,112,121,124,132,148,157]; 
        cellList = [37,38,47,48,63,73,83,93,100,118,131,150,160,173,202];
        color = colors(2,:);
    end
    % cellList = cellList(randperm(numel(cellList),min(numel(cellList),15)));
    
    %% Plot raster with events detected
    
    [fig1] = rasterPlot_ext9(fullData,cellList,color); % plots selected traces, 15
   
        set(gcf,'Units','inches','Position',[0,0,8.5,11])
        set(gcf,'units','normalized')
        set(gca,'Position',[0.05,0.1,0.85,0.25])
    % plot movement bouts in gray
    yy = ylim;
    for gg = 1:numel(fullData.movBoutStart)
        shadeXX = (fullData.movBoutStart(gg):fullData.movBoutFinish(gg))/20;
        shadeYY = repmat(yy(2),size(shadeXX));
        area(shadeXX,shadeYY,'FaceColor',[0 0 0],'FaceAlpha',0.25,'EdgeAlpha',0)           
        hold on     
    end
  ylim([1,24])
%     set(gcf,'Color','none')
%     set(gca,'Units','inches','InnerPosition',[.3 .3 2 2],'TickDir','out','Color','none','Box','off')
    saveas(fig1,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\b_WTKO_traces\',sprintf('traces_918-%s-%d.pdf',mouseID,day)))

%     saveas(fig1,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\b_WTKO_traces\',sprintf('traces-%s-%d.png',mouseID,day)))
%     saveas(fig1,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\b_WTKO_traces\',sprintf('traces-%s-%d.fig',mouseID,day)))
%     saveas(fig1,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\b_WTKO_traces\',sprintf('traces-%s-%d.epsc',mouseID,day)))
end