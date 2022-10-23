%% Description

% Plots of cells correlated to each other marked by their physical location
% only positive correlation plotted
%% Initialize

addpath('J:\nexmif_paper\Utils')  % Add utilities
init % Initialize data directories and genotypes
% close all

%% Begin code
t_start = tic;
% load correlation stats
% load(fullfile('D:\nexmif_paper\code_ball\stats\correlation','corrStatsTableSpeedAsym'),'corrStatsTable')
load(fullfile('J:\nexmif_paper\code_ball\stats\correlation','corrStatsTableSpeedPearson_08_22'),'corrStatsTable')
fieldNames = corrStatsTable.Properties.VariableNames;
corrStatsTable([6,13,18:19,29,31:35,40:43,45],:) = [];


load('J:\nexmif_paper\code_ball\stats\movement_responsive_cells\responsive_cells_bout_shuffle_08_22','responsive_cells')
responsive_cells_table = struct2table(responsive_cells);



miceBad = {};
listHigh = [];
listLow = [];
listBoth = [];
listHigh_non_sig = [];
listLow_non_sig = [];
listBoth_non_sig = [];
listisKO = [];
listEs = [];

listHighMov = [];
listLowMov = [];
listBothMov = [];
listMov = [];

for k =  1:numel(corrStatsTable.animal)
    mouseName = corrStatsTable{k,1}{:};
    day = corrStatsTable{k,2}{:};
    
    mPath = fullfile(dataFolder,sprintf('fullData_%s_D%i_%s',mouseName,day,'ball'));
    
    % handle missing files
    try
        data_var = load(mPath);
        fullData = data_var.fullData;
    catch
        miceBad = [miceBad,[mouseName,day,'ball']];
    end
    
    % display
    {mouseName, day, 'ball'}
    
    % get cell centroids
    
    goodTraces = fullData.goodIdx;
    goodTraces = setdiff(goodTraces,fullData.empty_idx_minusBG_new);
    centroids = fullData.centroids;
    centroids =  centroids(goodTraces,:);
    
    % get and threshold significant correlations
    
    
    n = numel(goodTraces);
    y_close = double(corrStatsTable{k,10}{:}<20);
    y_close(logical(triu(ones(n)))) = nan;
    
    w_corr_low = corrStatsTable{k,4}{:};
%     w_corr_low(isnan(w_corr_low)) = 0;
    w_corr_low(logical(triu(ones(n)))) = nan;
    w_thresh_low_pos = corrStatsTable{k,5}{:};
    w_thresh_low_neg = corrStatsTable{k,6}{:};
    adj_sig_corr_low_pos = (w_corr_low - w_thresh_low_pos > 1e-5) & (y_close == 0);
    idx_low_pos = find(adj_sig_corr_low_pos);
    adj_sig_corr_low_neg = w_corr_low - w_thresh_low_neg < 1e-5 & (y_close == 0);
    idx_low_neg = find(adj_sig_corr_low_neg);
%     adj_sig_corr_low_non_sig = (w_corr_low < w_thresh_low_pos) & (w_corr_low > w_thresh_low_neg) & (y_close == 0);
    adj_sig_corr_low_non_sig = (w_corr_low < w_thresh_low_pos) & (y_close == 0);
    
    w_corr_high = corrStatsTable{k,7}{:};
%     w_corr_high(isnan(w_corr_high)) = 0;
    w_corr_high(logical(triu(ones(n)))) = nan;
    w_thresh_high_pos = corrStatsTable{k,8}{:};
    w_thresh_high_neg = corrStatsTable{k,9}{:};
    adj_sig_corr_high_pos = (w_corr_high - w_thresh_high_pos > 1e-5) & (y_close == 0);
    idx_high_pos = find(adj_sig_corr_high_pos);
    adj_sig_corr_high_neg = (w_corr_high - w_thresh_high_neg < 1e-5) & (y_close == 0);
    idx_high_neg = find(adj_sig_corr_high_neg);
%     adj_sig_corr_high_non_sig = (w_corr_high < w_thresh_high_pos) & (w_corr_high > w_thresh_high_neg) & (y_close == 0);
     adj_sig_corr_high_non_sig = (w_corr_high < w_thresh_high_pos)  & (y_close == 0);
   
    idx_both_pos = intersect(idx_low_pos,idx_high_pos);
    
%     high_and_neg = adj_sig_corr_high_pos & (w_corr_high < 0);
%     if any(high_and_neg(:))
%         x = 1;
%         break
%     end
    %% Mov responsive
    isM = strcmp({responsive_cells.animal},mouseName);
    isD = [responsive_cells.day] == day;
    isR = [responsive_cells.responsive_bout];
    cell_list = [responsive_cells.cell];
    isR(isnan(isR)) = 0;
    [~,res_cell_list_idx] = ismember(cell_list(isM & isD & isR),goodTraces);
    
    % Both cells have to be mov responsive
    res_cell_mat_both = zeros(n);
    p = nchoosek(res_cell_list_idx,2);
    res_cell_mat_both(sub2ind([n,n],p(:,1),p(:,2))) = 1;
    res_cell_mat_both(sub2ind([n,n],p(:,2),p(:,1))) = 1;
    
    %% Stats
    
    
    listHigh = [listHigh,sum(adj_sig_corr_high_pos(:))];
    listLow = [listLow,sum(adj_sig_corr_low_pos(:))];
    listBoth = [listBoth,sum(adj_sig_corr_low_pos(:) & adj_sig_corr_high_pos(:))];
    
    listHigh_non_sig = [listHigh_non_sig,sum(adj_sig_corr_high_non_sig(:))];
    listLow_non_sig = [listLow_non_sig,sum(adj_sig_corr_low_non_sig(:))];
    listBoth_non_sig = [listBoth_non_sig,sum(adj_sig_corr_low_non_sig(:) & adj_sig_corr_high_non_sig(:))];
    
    listEs = [listEs, n*(n-1)/2 - sum(y_close(:) == 1)];
    listisKO = [listisKO,ismember(mouseName,miceKO)];
    
    listHighMov = [listHighMov,sum(sum(adj_sig_corr_high_pos & res_cell_mat_both))];
    listLowMov = [listLowMov,sum(sum(adj_sig_corr_low_pos & res_cell_mat_both))];
    listBothMov = [listBothMov,sum(sum(adj_sig_corr_low_pos & adj_sig_corr_high_pos & res_cell_mat_both))];
    listMov = [listMov,sum(sum(res_cell_mat_both))/2];
    
    
    
    
    %         %% plot
    %
    %         thresh = 0;
    %
    %         figure
    %         % set(gcf,'windowstyle','docked')
    %         set(gcf,'units','normalized','outerposition',[0 0 1 1])
    %
    %         for ii = 1:numel(idx_low_pos)
    %             [c1,c2] = ind2sub(size(adj_sig_corr_low_pos),idx_low_pos(ii));
    %             w =  w_corr_low(c1,c2);
    %             if w > thresh
    %                 hold on
    %                 l1 = line([centroids(c1,1),centroids(c2,1)],[centroids(c1,2),centroids(c2,2)],'lineWidth',w*10,'color',[0, 0, 1, 0.25]);
    %             end
    %
    %         end
    %
    %         for ii = 1:numel(idx_high_pos)
    %             [c1,c2] = ind2sub(size(adj_sig_corr_high_pos),idx_high_pos(ii));
    %             w =  w_corr_high(c1,c2);
    %             if w > thresh
    %                 hold on
    %                 l2 = line([centroids(c1,1),centroids(c2,1)],[centroids(c1,2),centroids(c2,2)],'lineWidth',w*10,'color',[1, 0, 0, 0.25]);
    %             end
    %         end
    %
    %         for ii = 1:numel(idx_both_pos)
    %             [c1,c2] = ind2sub(size(adj_sig_corr_high_pos),idx_both_pos(ii));
    %             w =  0.5*(w_corr_high(c1,c2)+w_corr_low(c1,c2));
    %             if w > thresh
    %                 hold on
    %                 l3 = line([centroids(c1,1),centroids(c2,1)],[centroids(c1,2),centroids(c2,2)],'lineWidth',w*10,'color',[0, 0, 0, 0.5]);
    %             end
    %         end
    %
    %         s1 = scatter(centroids(:,1),centroids(:,2),'filled','k');
    %         axis square
    %         set(gca,'XColor', 'none','YColor','none')
    %         try
    %         legend([l1,l2,l3,s1],{'Sig. correlated during rest','Sig. correlated during movement','Sig. correlated during both','centroid'})
    %         catch
    %         end
    %         title(sprintf('%s - Day %i - %s - Threshold %.2f',mouseName, day, mType{ismember(mouseName,miceKO)+1},thresh))
    %         saveas(gcf,fullfile('D:\nexmif_paper\code_ball\plots\correlation_analysis\Identity',sprintf('%s - Day %i - %s.fig',mouseName, day, mType{ismember(mouseName,miceKO)+1})));
    %         saveas(gcf,fullfile('D:\nexmif_paper\code_ball\plots\correlation_analysis\Identity',sprintf('%s - Day %i - %s.png',mouseName, day, mType{ismember(mouseName,miceKO)+1})));
    %
    %         pause()
    %         close all
end
%%
listisKO = logical(listisKO);
save(fullfile('J:\nexmif_paper\code_ball\stats\correlation\','identity_stats_pearson_08_22'),'listisKO','listLow','listHigh','listBoth','listLowMov','listHighMov','listBothMov','listEs','listMov','listLow_non_sig','listHigh_non_sig','listBoth_non_sig')

% save(fullfile('J:\nexmif_paper\code_ball\stats\correlation\','corrStatsTableSpeedPearson_08_22'),'corrStatsTable')


ttNames = {'Both','WT','KO'};
for tt = 1:3
    listAll = listLow + listHigh - listBoth;
    switch tt
        case 1
            pLow_1 = listLow./listAll*100;   % _1 - As a fraction of all significant correlations
            pHigh_1 = listHigh./listAll*100;
            pBoth_1 = listBoth./listAll*100;
            
            pLow_2 = listLow./listEs*100;   % _2 - As a fraction of all possible edges
            pHigh_2 = listHigh./listEs*100;
            pBoth_2 = listBoth./listEs*100;
            
            pLow_3 =  listLowMov./listMov*100;   % _2 - As a fraction of all possible movement responsive cell pairs
            pHigh_3 = listHighMov./listMov*100;
            pBoth_3 = listBothMov./listMov*100;
        case 2
            pLow_1 = listLow(~listisKO)./listAll(~listisKO)*100;
            pHigh_1 = listHigh(~listisKO)./listAll(~listisKO)*100;
            pBoth_1 = listBoth(~listisKO)./listAll(~listisKO)*100;
            pLow_2 = listLow(~listisKO)./listEs(~listisKO)*100;
            pHigh_2 = listHigh(~listisKO)./listEs(~listisKO)*100;
            pBoth_2 = listBoth(~listisKO)./listEs(~listisKO)*100;
            pLow_3 =  listLowMov(~listisKO)./listMov(~listisKO)*100;
            pHigh_3 = listHighMov(~listisKO)./listMov(~listisKO)*100;
            pBoth_3 = listBothMov(~listisKO)./listMov(~listisKO)*100;
        case 3
            pLow_1 = listLow(listisKO)./listAll(listisKO)*100;
            pHigh_1 = listHigh(listisKO)./listAll(listisKO)*100;
            pBoth_1 = listBoth(listisKO)./listAll(listisKO)*100;
            pLow_2 = listLow(listisKO)./listEs(listisKO)*100;
            pHigh_2 = listHigh(listisKO)./listEs(listisKO)*100;
            pBoth_2 = listBoth(listisKO)./listEs(listisKO)*100;
            pLow_3 =  listLowMov(listisKO)./listMov(listisKO)*100;
            pHigh_3 = listHighMov(listisKO)./listMov(listisKO)*100;
            pBoth_3 = listBothMov(listisKO)./listMov(listisKO)*100;
    end
    %     subplot(1,3,tt)
    figure
    boxplot([pHigh_1';pLow_1';pBoth_1'],[ones(size(pHigh_1'));2*ones(size(pLow_1'));3*ones(size(pBoth_1'))])
    hold on
    s1 = scatter(1,pHigh_1,'filled');
    cl = reshape([s1(:).CData],3,[])';
    hold on
    s2 = scatter(2,pLow_1,[],cl,'filled');
    hold on
    s3 = scatter(3,pBoth_1,[],cl,'filled');
    
    line([1,2,3],[pHigh_1',pLow_1',pBoth_1'])
    xlim([0.5,3.5])
    xticks([1,2,3])
    xticklabels({'During movement','During rest','Shared'})
    [p,tbl,stats] = anova1([pHigh_1',pLow_1',pBoth_1'],[],'off');
    if p<0.05
    [c,m] = multcompare(stats,'Display','off');
    title(sprintf('Distribution of sig correlations - %s \n P-anova: %.3f \n P-mov-rest: %.3f, P-mov-shared: %.3f,  P-rest-shared: %.3f,' ,ttNames{tt},p, c(1,6),c(2,6),c(3,6)))
    else
        title(sprintf('Distribution of sig correlations - %s \n P-anova: %.3f' ,ttNames{tt},p))
    end
    saveas(gcf,fullfile('D:\nexmif_paper\code_ball\plots\correlation_analysis\Identity',sprintf('Distribution of sig correlations - %s.fig',ttNames{tt})));
    saveas(gcf,fullfile('D:\nexmif_paper\code_ball\plots\correlation_analysis\Identity',sprintf('Distribution of sig correlations - %s.png',ttNames{tt})));
    
    
    figure
    boxplot([pHigh_2';pLow_2';pBoth_2'],[ones(size(pHigh_2'));2*ones(size(pLow_2'));3*ones(size(pBoth_2'))])
    hold on
    s1 = scatter(1,pHigh_2,'filled');
    cl = reshape([s1(:).CData],3,[])';
    hold on
    s2 = scatter(2,pLow_2,[],cl,'filled');
    hold on
    s3 = scatter(3,pBoth_2,[],cl,'filled');
    
    line([1,2,3],[pHigh_2',pLow_2',pBoth_2'])
    xlim([0.5,3.5])
    xticks([1,2,3])
    xticklabels({'During movement','During rest','Shared'})
    [p,tbl,stats] = anova1([pHigh_2',pLow_2',pBoth_2'],[],'off');
    if p<0.05
    [c,m] = multcompare(stats,'Display','off');
    title(sprintf('Perc of sig corr out of all edges - %s \n P-anova: %.3f \n P-mov-rest: %.3f, P-mov-shared: %.3f,  P-rest-shared: %.3f,' ,ttNames{tt},p, c(1,6),c(2,6),c(3,6)))
    else
        title(sprintf('Perc of sig corr out of all edges - %s \n P-anova: %.3f' ,ttNames{tt},p))
    end
    saveas(gcf,fullfile('D:\nexmif_paper\code_ball\plots\correlation_analysis\Identity',sprintf('Perc of sig corr out of all edges - %s.fig',ttNames{tt})));
    saveas(gcf,fullfile('D:\nexmif_paper\code_ball\plots\correlation_analysis\Identity',sprintf('Perc of sig corr out of all edges - %s.png',ttNames{tt})));
    
    figure
    boxplot([pHigh_3';pLow_3';pBoth_3'],[ones(size(pHigh_3'));2*ones(size(pLow_3'));3*ones(size(pBoth_3'))])
    hold on
    s1 = scatter(1,pHigh_3,'filled');
    cl = reshape([s1(:).CData],3,[])';
    hold on
    s2 = scatter(2,pLow_3,[],cl,'filled');
    hold on
    s3 = scatter(3,pBoth_3,[],cl,'filled');
    
    line([1,2,3],[pHigh_3',pLow_3',pBoth_3'])
    xlim([0.5,3.5])
    xticks([1,2,3])
    xticklabels({'During movement','During rest','Shared'})
    [p,tbl,stats] = anova1([pHigh_3',pLow_3',pBoth_3'],[],'off');
    if p<0.05
    [c,m] = multcompare(stats,'Display','off');
    title(sprintf('Perc of sig corr out of mov. responsive edges - %s \n P-anova: %.3f \n P-mov-rest: %.3f, P-mov-shared: %.3f,  P-rest-shared: %.3f,' ,ttNames{tt},p, c(1,6),c(2,6),c(3,6)))
    else
        title(sprintf('Perc of sig corr out of mov. responsive edges - %s \n P-anova: %.3f' ,ttNames{tt},p))
    end
    saveas(gcf,fullfile('D:\nexmif_paper\code_ball\plots\correlation_analysis\Identity',sprintf('Perc of sig corr out of mov. responsive edges - %s.fig',ttNames{tt})));
    saveas(gcf,fullfile('D:\nexmif_paper\code_ball\plots\correlation_analysis\Identity',sprintf('Perc of sig corr out of mov. responsive edges - %s.png',ttNames{tt})));
    
end
% saveas(gcf,fullfile('D:\nexmif_paper\code_ball\plots\correlation_analysis\Identity','fractions.fig'));
% [p,tab,stats] = anova1([listHigh;listLow;listBoth]');
% multcompare(stats);