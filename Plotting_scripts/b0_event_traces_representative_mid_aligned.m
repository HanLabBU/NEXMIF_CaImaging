%% Description

% This code runs through all mice - ball data and generates a
% representative event shape 

%% Initialize

addpath(genpath('U:\eng_research_handata\Athif Mohamed\nexmif_paper\Utils\'))  % Add utilities
init % Initialize data directories and genotypes
close all

%% Begin code

miceBad = [];

preIdx = round(4*Fs);
postIdx = round(10*Fs);

blockTraceMeanWTMice = zeros(1, preIdx+postIdx+1);
blockTraceMeanKOMice = zeros(1, preIdx+postIdx+1);

figWT = figure;
figKO = figure;

WT_cnt = 0;
KO_cnt = 0;

for m = 1:numel(miceStudy)
    for d = 1:3
        for c = 1%:3 %:3 %ball
            % Load files
            mPath = fullfile(dataFolder,sprintf('fullData_%s_D%i_%s',miceStudy{m},dIdx(d),lower(conditionList{c})));
            
            if ismember(sprintf('%s_D%i',miceStudy{m},dIdx(d)),badSessions)
                continue
            end
            
            % handle missing files
            try
                load(mPath)
            catch
                miceBad = [miceBad,[miceStudy(m);dIdx(d);{lower(conditionList{c})}]];
            end
            
            % display
            {miceStudy{m}, num2str(dIdx(d)), conditionList{c}}
            
            
            % Fluoroscence trace
            tracesXX = [fullData.roi_list_minusBG_new.trace];
            traces = reshape(tracesXX',[],size(fullData.roi_list_minusBG_new,2));
            
            % Disregard noisy traces (manually marked)
            goodTraces = fullData.goodIdx;
            goodTraces = setdiff(goodTraces,fullData.empty_idx_minusBG_new);
            traces = traces(:,goodTraces);
            n = numel(goodTraces);
            traceLen = size(traces,1);
            
            % Get event indices and start aligning them
            event_idx_all = fullData.roi_list_minusBG_new(:,goodTraces);
            
            blockTraceSumSession = zeros(1, preIdx+postIdx+1);
            for cc = 1:n
                trace = traces(:,cc);
                %                 event_idx = event_idx_all(cc).event_idx(:,2);
                event_idx = round(0.5*(event_idx_all(cc).event_idx(:,1)+ event_idx_all(cc).event_idx(:,2)));
                eventN = numel(event_idx);
                
                blockTraceSum = zeros(1, preIdx+postIdx+1);
                eN = 0;
                for ee = 1:eventN
                    blockIdx = max(1,event_idx(ee)-preIdx):min(traceLen,event_idx(ee)+postIdx);
                    blockTrace = trace(blockIdx)';
                    if numel(blockTrace) == preIdx+postIdx+1
                        eN = eN+1;
                        blockTraceSum = blockTraceSum +blockTrace;
                    end
                end
                blockTraceMean = blockTraceSum./eN;  % Average event shape for a trace
                
                blockTraceMean = (blockTraceMean-min(blockTraceMean))/(max(blockTraceMean)-min(blockTraceMean));
                blockTraceSumSession = blockTraceSumSession +blockTraceMean;
            end
            
            blockTraceMeanSession = blockTraceSumSession./n;
            blockTraceMeanSession = blockTraceMeanSession-min(blockTraceMeanSession);

            if maskWT(m)
                blockTraceMeanWTMice = blockTraceMeanWTMice + blockTraceMeanSession;
                WT_cnt = WT_cnt+1;
                figure(figWT);
                p1 = plot([1:numel(blockTraceMeanSession)]/20*1000,blockTraceMeanSession,'color',[colors(1,:), 0.4],'LineWidth',2);
                hold on
            end
            if ~maskWT(m)
                blockTraceMeanKOMice = blockTraceMeanKOMice + blockTraceMeanSession;
                KO_cnt = KO_cnt+1;
                figure(figKO);
                p2 = plot([1:numel(blockTraceMeanSession)]/20*1000,blockTraceMeanSession,'color',[colors(3,:), 0.4],'LineWidth',2);
                hold on
            end
        end
    end
end

blockTraceMeanWTMice = blockTraceMeanWTMice/WT_cnt;
blockTraceMeanKOMice = blockTraceMeanKOMice/KO_cnt;

%%

figure(figWT);
p3 = plot([1:numel(blockTraceMeanWTMice)]/20*1000,blockTraceMeanWTMice,'color',colors(2,:),'LineWidth',2);
set(gcf,'Units','inches','InnerPosition',[1, 3, 2, 2])
set(gca,'TickDir', 'out','TickLength',[0.03, 0.025], 'Color','none')
box off
%title('Average calcium event WT')
%legend([p1,p3],{'Session Average','WT Average'})
saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\b_WTKO_traces\WT_mid.png'))
saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\b_WTKO_traces\WT_mid.fig'))
saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\b_WTKO_traces\WT_mid.epsc'))
figWT(end)

figure(figKO);
p4 = plot([1:numel(blockTraceMeanKOMice)]/20*1000,blockTraceMeanKOMice,'color',colors(4,:),'LineWidth',2);
set(gcf,'Units','inches','InnerPosition',[1, 3, 2, 2])
set(gca,'TickDir', 'out','TickLength',[0.03, 0.025], 'Color','none')
box off
% title('Average calcium event KO')
% legend([p2,p4],{'Session Average','KO Average'})
saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\b_WTKO_traces\KO_mid.png'))
saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\b_WTKO_traces\KO_mid.fig'))
saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\b_WTKO_traces\KO_mid.epsc'))
figKO(end)
