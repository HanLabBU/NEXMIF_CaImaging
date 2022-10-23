% raster plot of traces
% 20 frames/sec

% ext 1 - RAM 9.28.21
% plots only first 100 s of data

function[f1,f2,f3] = rasterPlot_paper_zoomin(fullData) 
close all 

mouseID = '605842';
day = 1;

% %tmax = size(fullData.roi_list_minusBG_new(1).trace,2);
tmax = size(fullData.trace,2);     
tMaxPlot = 200*20;          % 200 s of data
highlights = [40,44,55,66,69,85,111,122,141,153,154,157,163,171,182,188,207,247,267,269];     % 605842 d1 ball

f1 = figure(1)
hold on
% plot traces for all ROIs and collect numbers and locations for all traces

tracesToPlot = zeros(length(highlights), tmax);
for i = 1:length(highlights)
    tracesToPlot(i,:) = fullData.roi_list_minusBG_new(highlights(i)).trace(1:tmax);
end

for tr = 1:size(tracesToPlot,1)
    meanTrace = mean(tracesToPlot(tr,:));        %normalize all traces prior to plotting
    trace = tracesToPlot(tr,:) - meanTrace;
    plot([1:tMaxPlot], trace(1:tMaxPlot) + (tr-1)*2, 'Color', [0,.8706,0],'LineWidth',0.5)
end 

set(gca,'Units','inches','InnerPosition',[0.3 0.3 2 2],'TickDir','out','Color','none','Box','off');           % left bottom (from bottom left corner of fig, width height 
hold off

% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\fig1_traces\',sprintf('traces-zoom-%s-%d.png',mouseID,day)))
% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\fig1_traces\',sprintf('traces-zoom-%s-%d.fig',mouseID,day)))
% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\fig1_traces\',sprintf('traces-zoom-%s-%d.epsc',mouseID,day)))

f2 = figure(2)
plot([1:tMaxPlot],fullData.speed(1:tMaxPlot))
set(gca,'Units','inches','InnerPosition',[0.3 0.3 2 1],'TickDir','out','Color','none','Box','off');

% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\fig1_traces\',sprintf('speed-zoom-%s-%d.png',mouseID,day)))
% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\fig1_traces\',sprintf('speed-zoom-%s-%d.fig',mouseID,day)))
% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\fig1_traces\',sprintf('speed-zoom-%s-%d.epsc',mouseID,day)))

f3 = figure(3)
hold on

for tr = 1:size(tracesToPlot,1)
    meanTrace = mean(tracesToPlot(tr,:));        %normalize all traces prior to plotting
    trace = tracesToPlot(tr,:) - meanTrace;
    plot([tmax-tMaxPlot+1:tmax], trace(tmax-tMaxPlot+1:tmax) + (tr-1)*2, 'Color', [0,.8706,0],'LineWidth',0.5)
end 

set(gca,'Units','inches','InnerPosition',[0.3 0.3 2 2],'TickDir','out','Color','none','Box','off');           % left bottom (from bottom left corner of fig, width height 
hold off

% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\fig1_traces\',sprintf('traces-zoom-end-%s-%d.png',mouseID,day)))
% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\fig1_traces\',sprintf('traces-zoom-end-%s-%d.fig',mouseID,day)))
% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\fig1_traces\',sprintf('traces-zoom-end-%s-%d.epsc',mouseID,day)))