clear all
close all
clc

addpath(genpath('U:\eng_research_handata\Athif Mohamed\nexmif_paper\NEXMIF_CaImaging\Utils'))  % Add utilities
init % Initialize data directories and genotypes
close all

mouseID = '605842';
day = 1;
% roi idx to plot:
highlightForPaper = [40,44,55,66,69,85,111,122,141,153,154,157,163,171,182,188,207,247,267,269];  % 605842 d1 

mPath = fullfile(dataFolder,sprintf('fullData_%s_D%i_%s',mouseID,day,'ball'));

roiPath = fullfile('U:\eng_research_handata\eng_research_handata2\Rebecca_Mount\Autism\',...
    sprintf('%s\\Day_%i\\%s\\roi_%s_D%i_%s',mouseID,day,'ball',mouseID,day,'ball')); % load roi

maxminPath = fullfile('U:\eng_research_handata\eng_research_handata2\Rebecca_Mount\Autism\',...
    sprintf('%s\\Day_%i\\%s\\maxmin_%s_D%i_%s',mouseID,day,'ball',mouseID,day,'ball')); % load maxmin

% load(roiPath) % load roiList
load(mPath) % load fullData
% load(maxminPath) % load maxmin
% 
% roi_overlay_paper(roiList, roiList(1,highlightForPaper), maxmin)  %plot the roi
% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\example_traces\',sprintf('roi-%s-%d.png',mouseID,day)))
% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\example_traces\',sprintf('roi-%s-%d.fig',mouseID,day)))

nCells = numel(fullData.goodIdx);
traceLen = numel(fullData.roi_list_minusBG_new(highlightForPaper(1)).trace);
Fs = 20;

%% Plot heat map of all cells with highlighted cells at the bottom

% make a mat of traces -- highlights first, in order 
traceMat = zeros(nCells,traceLen);
for c = 1:length(highlightForPaper)
    traceMat(c,:) = fullData.roi_list_minusBG_new(highlightForPaper(c)).trace;
end

% find non-highlights
in = 1;
nonHighlight = zeros(1,nCells-length(highlightForPaper));
for i = 1:nCells
    if fullData.goodIdx(i) ~= highlightForPaper
        nonHighlight(in) = fullData.goodIdx(i);
        in = in+1;
    end
end

%% add non-highlights to mat
for c = length(highlightForPaper)+1:nCells
    idx = c-length(highlightForPaper);
    traceMat(c,:) = fullData.roi_list_minusBG_new(nonHighlight(idx)).trace;
end

% flip mat so highlights are at the bottom and same order as green traces
traceMat = flip(traceMat,1);
% subtract mean from each trace
avgs = mean(traceMat,2);
traceMat = traceMat - avgs;
% % plot them
figure()
imagesc(traceMat)
set(gca,'Units','inches','InnerPosition',[0.3 0.3 2 3],'XTick',[],'Color','none','Box','off');
colormap('jet')
caxis([0 1])
colorbar('Box','off','TickDirection','out','TickLength',0.03,'LineWidth',2)

saveas(gcf,fullfile(savePath,sprintf('fig1D1_full_heatmap-%s-%d.png',mouseID,day)))
% saveas(gcf,fullfile(savepath,sprintf('full_heatmap-%s-%d.fig',mouseID,day)))
% saveas(gcf,fullfile(savePath,sprintf('full_heatmap-%s-%d.epsc',mouseID,day)))


% %% Plot raster with speed -- use rasterPlot_paper_ext1 in Rebecca scripts
% % master instead 9.28.21
% 
% [~,fig2] = rasterPlot_ext3(fullData,highlightForPaper);
% set(gca,'Units','inches','InnerPosition',[0.3 0.3 2 1],'TickDir','out','LineWidth',2,'Color','none','Box','off');
% % saveas(fig1,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\example_traces\',sprintf('traces-%s-%d.png',mouseID,day)))
% % saveas(fig1,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\example_traces\',sprintf('traces-%s-%d.fig',mouseID,day)))
% 
% saveas(fig2,fullfile(savePath,sprintf('speed-%s-%d.png',mouseID,day)))
% % saveas(fig2,fullfile(savePath,sprintf('speed-%s-%d.fig',mouseID,day)))
% % saveas(fig2,fullfile(savePath,sprintf('speed-%s-%d.epsc',mouseID,day)))


%% Plot raster with speed
[f1,f2,f3] = rasterPlot_paper_zoomin(fullData);

saveas(f1,fullfile(savePath,sprintf('fig1E1_raster-%s-%d.png',mouseID,day)))
saveas(f3,fullfile(savePath,sprintf('fig1E2_raster-%s-%d.png',mouseID,day)))
saveas(f2,fullfile(savePath,sprintf('fig1D2_speed-%s-%d.png',mouseID,day)))

