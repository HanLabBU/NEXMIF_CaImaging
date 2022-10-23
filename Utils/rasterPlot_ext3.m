% RAM, last updated 5.20.21
% raster plot of traces
% 20 frames/sec

% Ext1 - Modified by Athif 5.21.21
% takes in list of good Traces
% Plots groups in multiple figures

% Ext3 - Modified by Athif 6.16.21
% Plots speed

function [fig1,fig2] = rasterPlot_ext3(fullData,goodTraces,tracePerPlot)

tmax = size(fullData.trace,2);
if nargin == 2
    tracePerPlot = numel(goodTraces);
end
nPlots = ceil(numel(goodTraces)/tracePerPlot);

for pp = 1:nPlots
    yMarks = [];
    yLabels = {};
    traceIdxs = tracePerPlot*(pp-1)+1:min(tracePerPlot*(pp),numel(goodTraces));
    
    fig1 = figure;
    set(fig1,'units','normalized','outerposition',[0 0 1 1])

    for ii = 1:numel(traceIdxs)
        idx = goodTraces(traceIdxs(ii));
        trace = fullData.roi_list_minusBG_new(idx).trace;
        plot([1:tmax]/20,trace+(ii)*3);
        
        hold on
        
        eeList  = fullData.roi_list_minusBG_new(idx).event_idx(:,1);
        yMarks(ii) = (ii)*3;
        yLabels{ii} = num2str(idx);
        hold on
    end
    hold on
    
    
    yticks(yMarks)
    yticklabels(yLabels)
    xlabel('Time')
    ylabel('ROI')
end

fig2 = figure;
% set(fig2,'units','normalized','outerposition',[0 0 1 0.2])

speed = fullData.speed;
% speed_norm = bsxfun(@rdivide,bsxfun(@minus, speed, min(speed)),max(speed)-min(speed)); % normalize speed 
plot([1:tmax]/20,speed);

% xlabel('Time')
% ylabel('Speed (cm/s)')