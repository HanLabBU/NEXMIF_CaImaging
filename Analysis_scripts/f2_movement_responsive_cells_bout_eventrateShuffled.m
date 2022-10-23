% Based on

% Movement_responsive_cells
% before after T test stats

%  based on ext 6
%  use event rates within movement bouts vs not T test
%% code

% close all
addpath('J:\nexmif_paper\Utils')  % Add utilities
init

% cd(pathData)

k = 0;
responsive_cells = struct();
nIter = 1000;

%%
miceBad = [];
for m = 1:numel(miceStudy)
    for d = 1:3
        for c = 1 % ball only
            % Load files
            
            mPath = fullfile(dataFolder,sprintf('fullData_%s_D%i_%s',miceStudy{m},dIdx(d),lower(conditionList{c})));
             
                if ismember(sprintf('%s_D%i',miceStudy{m},dIdx(d)),badSessions)||ismember(sprintf('%s_D%i',miceStudy{m},dIdx(d)),badSpeedSessions)||ismember(sprintf('%s_D%i',miceStudy{m},dIdx(d)),fewMovSessions)
                    continue
                end
            
            try
                load(mPath)
            
            
            % display
            {miceStudy{m}, num2str(dIdx(d)), conditionList{c}}
            
            %% Onsets and offsets
            idx_onset = fullData.movBoutStart;
            time_onset = idx_onset/Fs;
            
            idx_offset = fullData.movBoutFinish;
            time_offset = idx_offset/Fs;
            
            
            %% Trace
            
            % collect traces
            preIdx = 5*Fs;
            traces = fullData.roi_list_minusBG_new;
            nTraces = size(traces,2);
            traceLen = size(traces(1).trace,2);
            
            binTraces = binarizeTrace(traces);
            onsetTraces = binarizeTrace_raster(traces);
            goodTraces = fullData.goodIdx;
            goodTraces = setdiff(goodTraces,fullData.empty_idx_minusBG_new);
            
            %             % plot and see
            %             rasterPlot_ext3(fullData,goodTraces,25)
            %
            for cc = 1:nTraces %[30,44,46,52,58]49:55
                
                if ismember(cc,goodTraces)
                    k = k+1;
                    
                    
                    nBout = numel(fullData.movBoutStart);
                    nRest = numel(fullData.restBoutStart);
                    %% a - Onset raster - Event rate during movement bout vs outside of movement bout
%                     [EventRateRest,EventRateMov] = tTestBoutEventActrate(binTraces,cc,fullData.movBoutStart,fullData.movBoutFinish,fullData.restBoutStart,fullData.restBoutFinish,Fs);
% 
%                     [h,p] = ttest2(EventRateMov,EventRateRest);
%                     [h_pos,p_pos] = ttest2(EventRateMov,EventRateRest,'Tail','right'); % Mov - rest is greater than 0
%                     [h_neg,p_neg] = ttest2(EventRateMov,EventRateRest,'Tail','left'); % Mov - rest is less than 0
                    
                    
                    %% a - bout - full rising within the window as a response
                    [BoutFullRisingShuffled,BoutFullRisingTrue, mov_r,res_r]= shuffledBoutFullRising(binTraces,cc,traceLen,nIter,fullData.movBoutStart,fullData.movBoutFinish,fullData.restBoutStart,fullData.restBoutFinish,Fs);
%                     [BoutFullRisingShuffled,BoutFullRisingTrue, mov_r,res_r]= shuffledBoutFullRising_ca(binTraces,cc,traceLen,nIter,fullData.movBoutIdx,fullData.restBoutIdx,Fs);
%                     prc975Bout  = prctile(BoutFullRisingShuffled,97.5);
%                     prc025Bout  = prctile(BoutFullRisingShuffled,2.5); 
                    prc975Bout  = prctile(BoutFullRisingShuffled,97.5);
                    prc025Bout  = prctile(BoutFullRisingShuffled,2.5);

                    
                    %% Assign basic details
                    responsive_cells(k).animal = miceStudy{m};
                    responsive_cells(k).day = dIdx(d);
                    responsive_cells(k).cell = cc;
                    
                    if maskWT(m)
                        responsive_cells(k).isWT = 1;   %1 WT
                    else
                        responsive_cells(k).isWT = 0;
                    end
                    
                    
                    responsive_cells(k).nBouts = nBout;
                    responsive_cells(k).metric_bout = BoutFullRisingTrue;
                    responsive_cells(k).metric_bout_prc975 = prc975Bout;
                    responsive_cells(k).responsive_bout = BoutFullRisingTrue > prc975Bout;
                    responsive_cells(k).mov_r = mov_r;
                    responsive_cells(k).res_r = res_r;
%                     responsive_cells(k).nRest = nRest;
%                     responsive_cells(k).metric_rest = RestFullRisingTrue;
                    responsive_cells(k).metric_rest_prc025 = prc025Bout;
                    responsive_cells(k).responsive_rest = BoutFullRisingTrue < prc025Bout;
                    
                end
            end
            catch
                miceBad = [miceBad,[miceStudy(m);dIdx(d);{lower(conditionList{c})}]];
            end
        end
    end
end


save('J:\nexmif_paper\code_ball\stats\movement_responsive_cells\responsive_cells_bout_shuffle_08_22','responsive_cells')

%% Plot

load('D:\nexmif_paper\code_ball\stats\movement_responsive_cells\responsive_cells_bout_shuffle','responsive_cells')


% Remove last 1 mouse 688
% responsive_cells_table = responsive_cells_table(1:6007,:);
responsive_cells_table = struct2table(responsive_cells);
miceAll = responsive_cells_table.animal;
miceAllUnique = unique(miceAll);

% Boxplot for animals

featureNames = {'High movement responsive cells','Low movement responsive cells'};
isWT_pop = logical(responsive_cells_table.isWT);
binEdge = 0:0.2:1;

cols = [8,12];  
nEvCols = [5,5]; % number of events
% cols = [8] ; % h_pos
% cols = [9] ; % h_neg



for colIdx = 2% 1:2
    
    featureNames{colIdx}
    nEventsAll = responsive_cells_table{:,nEvCols(colIdx)};
    IdxEv = nEventsAll> 5;
    
    
    col = cols(colIdx);
    
    % Histogram plots of cells care vs not for sessions
    careListPercent = [];
    noCareListPercent = [];
    isWTList = [];
    nList = [];
    nameList = {};
    careList = [];
    noCareList = [];
    for m = 1:numel(miceAllUnique)
        idxM = strcmp(miceAll,miceAllUnique(m));
        isWT = responsive_cells_table{find(idxM,1,'first'),4};
        
        
        idx = idxM & IdxEv;
        if sum(idx) >= 50
            do_cells_care = [responsive_cells_table{idx,col}];
            
            % miceAllUnique(m),d
            % sum(cells_care), numel(cells_care)
            
            if ~isempty(do_cells_care)
                
                care = sum(do_cells_care)/numel(do_cells_care)*100;
                noCare = sum(~ do_cells_care)/numel(do_cells_care)*100;
                
                nList = [nList,numel(do_cells_care)];
                careListPercent = [careListPercent,care];
                isWTList = [isWTList,isWT];
                noCareListPercent = [noCareListPercent,noCare];
                nameList = [nameList, sprintf('%s',miceAllUnique{m})];
                
                careList = [careList,sum(do_cells_care)];
                noCareList = [noCareList,sum(~do_cells_care)];
            end
        end
    end
    
    data = careListPercent;
    group = ~isWTList;
    
    figure
    set(gcf,'units','normalized','outerposition',[0 0 0.25 1])
    
    p1 = boxplot(data,group,'Colors','k','Widths',0.3);
    
    hold on
    yWT = data(logical(isWTList));
    xWT = 1*ones(size(yWT));
    p2 = scatter (xWT,yWT,'filled','MarkerFaceColor',colors(1,:));
    hold on
    yKO = data(~logical(isWTList));
    xKO = 2*ones(size(yKO));
    p3 = scatter (xKO,yKO,'filled','MarkerFaceColor',colors(3,:));
    
    % Stats
    [h_t,p_t] = ttest2(yKO,yWT)
%     meanWT = mean(yWT)
%     meanKO = mean(yKO)
    
    title(featureNames{colIdx})
    set(gca,'xtick',[1 2],'xticklabel',{'WT','KO'})
    xlim([0.5 2.5])
    ylabel('Percentage of movement responsive cells')
    box off
    %     saveas(gcf,fullfile('D:\Autism\Event analysis\Results\Mvmt_Aligned_Activity_Stats_ext4\',sprintf('Session_Care_hist_%s_%s_%i.png',featureNames{f},miceAllUnique{m},d)))
    %     ylim([0 up+0.1])
    
    %     legend([p2,p3,p6,p7,p4,p5],{'WT Sessions','KO Sessions','WT Session Mean','KO Session Mean','WT Population','KO Population'})
    %     saveas(gcf,fullfile('D:\nexmif_paper\plots\movement_responsive_cells\',sprintf('%s_animal_level.png',featureNames{colIdx})))
    saveas(gcf,fullfile('D:\nexmif_paper\code_ball\plots\movement_responsive_cells',sprintf('shuffled %s.png',featureNames{colIdx})))
    saveas(gcf,fullfile('D:\nexmif_paper\code_ball\plots\movement_responsive_cells',sprintf('shuffled %s.fig',featureNames{colIdx})))
   %     updateF fisher mat
    fisherMat = zeros(2);
    fisherMat(1,1) = sum(careList(logical(isWTList)));
    fisherMat(1,2) = sum(noCareList(logical(isWTList)));
    fisherMat(2,1) = sum(careList(logical(~isWTList)));
    fisherMat(2,2) = sum(noCareList(logical(~isWTList)));
    
    
    fisherMat
    [h_f,p_f,stats_f] = fishertest(fisherMat,'Alpha',0.05)%,'Tail','right'
    [h_f_right,p_f_right,stats_f_right] = fishertest(fisherMat,'Tail','right','Alpha',0.05)%
    
    n_WT = fisherMat(1,1)+fisherMat(1,2);
    pr_WT = fisherMat(1,1)/n_WT*100;
    SEP_WT = sqrt(pr_WT*(100-pr_WT)/n_WT);
    n_KO = fisherMat(2,1)+fisherMat(2,2);
    pr_KO = fisherMat(2,1)/n_KO*100;
    SEP_KO = sqrt(pr_KO*(100-pr_KO)/n_KO);
    
    figure
    errorbar(-1,pr_WT,SEP_WT,'o','MarkerFaceColor',colors(1,:),'MarkerEdgeColor',colors(1,:),'Color',colors(1,:));
    hold on
    errorbar(1,pr_KO,SEP_KO,'o','MarkerFaceColor',colors(3,:),'MarkerEdgeColor',colors(3,:),'color',colors(3,:));
    xticks(gca,[-1 1])
    xticklabels(gca, {'WT','KO'})
  
    title('Proportion of resposive cells')
    xlim([-2,2])
    yy = ylim;
    
    text(0,mean(yy),sprintf('Fisher sig: p = %d',p_f),'HorizontalAlignment','center')
    saveas(gcf,fullfile('D:\nexmif_paper\code_ball\plots\movement_responsive_cells',sprintf('shuffled %s fisher.png',featureNames{colIdx})))
    saveas(gcf,fullfile('D:\nexmif_paper\code_ball\plots\movement_responsive_cells',sprintf('shuffled %s fisher.fig',featureNames{colIdx})))
%     PercWT = fisherMat(1,1)/(sum(fisherMat(1,1)+fisherMat(1,2)))
%     PercKO = fisherMat(2,1)/(sum(fisherMat(2,1)+fisherMat(2,2)))
end


%%
% idWT = find(isWTList);
% idKO = setdiff(1:numel(isWTList),idWT);
% 
% figure
% set(gcf,'units','normalized','outerposition',[0 0 0.5 0.9])
% 
% bWTAllList = nList;
% bWTAllList(idKO) = nan;
% bWTAll = bar(categorical(nameList),bWTAllList,'FaceColor','flat');
% bWTAll.CData = repmat(colors(2,:),numel(bWTAllList),1);
% hold on
% 
% bWTCareList = careList;
% bWTCareList(idKO) = nan;
% bWTCare = bar(categorical(nameList),bWTCareList,'FaceColor','flat');
% bWTCare.CData = repmat(colors(1,:),numel(bWTCareList),1);
% hold on
% 
% bKOAllList = nList;
% bKOAllList(idWT) = nan;
% bKOAll = bar(categorical(nameList),bKOAllList,'FaceColor','flat');
% bKOAll.CData = repmat(colors(4,:),numel(bKOAllList),1);
% hold on
% 
% bKOCareList = careList;
% bKOCareList(idWT) = nan;
% bKOCare = bar(categorical(nameList),bKOCareList,'FaceColor','flat');
% bKOCare.CData = repmat(colors(3,:),numel(bKOCareList),1);
% hold on
% 
% text([1:numel(nList)],nList+5,num2cell(nList),'HorizontalAlignment','center')
% hold on
% text([1:numel(nList)],careList+5,num2cell(careList),'HorizontalAlignment','center')
% legend([bWTAll,bWTCare,bKOAll,bKOCare],{'Total # Cells WT','# Responsive Cells WT','Total # Cells KO','# Responsive Cells KO'})
% 
% title('Number of cells that care for movement onsets')
% 
% %%
% 
% nOnsetList = [];
% nameList = {};
% isWTList = [];
% for m = 1:numel(miceAllUnique)
%     idxM = strcmp(miceAll,miceAllUnique(m));
%     for d = [1,3,5]
%         idxD = responsive_cells_table.day == d;
%         idx = idxD & idxM;
%         idxFirst = find(idx,1,'first');
%         nOnset = responsive_cells_table{idxFirst,5};
%         isWT = responsive_cells_table{idxFirst,4};
%         if isempty(nOnset)
%             nOnsetList = [nOnsetList,nan];
%             isWTList = [isWTList,nan];
%         else
%             nOnsetList = [nOnsetList,nOnset];
%             isWTList = [isWTList,isWT];
%         end
%         nameList = [nameList,sprintf('%s-day-%i',miceAllUnique{m},d)];
%     end
% end
% figure
% 
% idWT = find(isWTList);
% idKO = setdiff(1:numel(isWTList),idWT);
% %
% % b1 = bar(categorical(nameList(idWT)),nOnsetList(idWT),'FaceColor','flat');
% % b1.CData = repmat(colors(1,:),numel(idWT),1);
% % hold on
% % b2 = bar(categorical(nameList(idKO)),nOnsetList(idKO),'FaceColor','flat');
% % b2.CData = repmat(colors(3,:),numel(idWT),1);
% % legend('WT','KO')
% 
% b = bar(categorical(nameList),nOnsetList,'FaceColor','flat');
% idWT = find(isWTList);
% idKO = setdiff(1:numel(isWTList),idWT);
% b.CData(idWT,:) = repmat(colors(1,:),numel(idWT),1);
% b.CData(idKO,:) = repmat(colors(3,:),numel(idKO),1);
% 
% 
% %%
% % If the returned test decision h = 1 indicates that fishertest does REJECTS
% % the null hypothesis of no nonrandom association between the categorical variables at the 1% significance level.
% % Since this is a right-tailed hypothesis test,
% % the conclusion is that WT cells DO have greater odds of being responsive to movement onset than
% % KO