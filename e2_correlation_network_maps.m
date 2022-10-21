%% Description

% Plots of cells correlated to each other marked by their physical location
% only positive correlation plotted
%% Initialize

addpath(genpath('U:\eng_research_handata\Athif Mohamed\nexmif_paper\Utils\'))  % Add utilities
init % Initialize data directories and genotypes
close all

%% Begin code
t_start = tic;
% load correlation stats
% load(fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\stats\correlation_ball','corrStatsTableSpeedAsym'),'corrStatsTable')
load(fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\stats\correlation_ball','corrStatsTableSpeedPearson_08_22'),'corrStatsTable')
fieldNames = corrStatsTable.Properties.VariableNames;
corrStatsTable([6,13,18:19,29,31:35,40:43,45],:) = [];


load('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\stats\movement_responsive_cells\responsive_cells_bout_shuffle_08_22','responsive_cells')
responsive_cells_table = struct2table(responsive_cells);


% load('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\stats\correlation_ball\identity_stats')
load('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\stats\correlation_ball\identity_stats_pearson_08_22')
miceBad = {};

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
    
    
    %% plot
    
    thresh = 0;
    % colors
    
    if ismember(mouseName,miceKO)
        col1 = colors(3,:);
        col2 = colors(4,:);
    else
        col1 = colors(1,:);
        col2 = colors(2,:);
    end
    
    figure
    % set(gcf,'windowstyle','docked')
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    
%     for ii = 1:numel(idx_low_pos)
%         [c1,c2] = ind2sub(size(adj_sig_corr_low_pos),idx_low_pos(ii));
%         w =  w_corr_low(c1,c2);
%         if w > thresh
%             hold on
%             l1 = line([centroids(c1,1),centroids(c2,1)],[centroids(c1,2),centroids(c2,2)],'lineWidth',1,'color',[col1]);
%         end
%         
%     end
%     
    for ii = 1:numel(idx_high_pos)
        [c1,c2] = ind2sub(size(adj_sig_corr_high_pos),idx_high_pos(ii));
        w =  w_corr_high(c1,c2);
        if w > thresh
            hold on
            l2 = line([centroids(c1,1),centroids(c2,1)],[centroids(c1,2),centroids(c2,2)],'lineWidth',1,'color',[col2]);
        end
    end
    
    for ii = 1:numel(idx_both_pos)
        [c1,c2] = ind2sub(size(adj_sig_corr_high_pos),idx_both_pos(ii));
        w =  0.5*(w_corr_high(c1,c2)+w_corr_low(c1,c2));
        if w > thresh
            hold on
            l3 = line([centroids(c1,1),centroids(c2,1)],[centroids(c1,2),centroids(c2,2)],'lineWidth',1,'color',[0.3, 0.3, 0.3]);
        end
    end
    
    s1 = scatter(centroids(:,1),centroids(:,2),24,'k','filled');
    axis square
    set(gca,'XColor', 'none','YColor','none')
    try
        legend([l1,l2,s1],{'Sig. correlated during running','Sig. correlated during both','centroid'})
    catch
    end
    title(sprintf('%s - Day %i - %s',mouseName, day, mType{ismember(mouseName,miceKO)+1}))
    saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\e2_correlation_network_maps',sprintf('Split high - %s - Day %i - %s.fig',mouseName, day, mType{ismember(mouseName,miceKO)+1})));
    saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\e2_correlation_network_maps',sprintf('Split high - %s - Day %i - %s.png',mouseName, day, mType{ismember(mouseName,miceKO)+1})));

    
      figure
    % set(gcf,'windowstyle','docked')
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    
    for ii = 1:numel(idx_low_pos)
        [c1,c2] = ind2sub(size(adj_sig_corr_low_pos),idx_low_pos(ii));
        w =  w_corr_low(c1,c2);
        if w > thresh
            hold on
            l1 = line([centroids(c1,1),centroids(c2,1)],[centroids(c1,2),centroids(c2,2)],'lineWidth',1,'color',[col1]);
        end
        
    end
    
%     for ii = 1:numel(idx_high_pos)
%         [c1,c2] = ind2sub(size(adj_sig_corr_high_pos),idx_high_pos(ii));
%         w =  w_corr_high(c1,c2);
%         if w > thresh
%             hold on
%             l2 = line([centroids(c1,1),centroids(c2,1)],[centroids(c1,2),centroids(c2,2)],'lineWidth',1,'color',[col2]);
%         end
%     end
    
    for ii = 1:numel(idx_both_pos)
        [c1,c2] = ind2sub(size(adj_sig_corr_high_pos),idx_both_pos(ii));
        w =  0.5*(w_corr_high(c1,c2)+w_corr_low(c1,c2));
        if w > thresh
            hold on
            l3 = line([centroids(c1,1),centroids(c2,1)],[centroids(c1,2),centroids(c2,2)],'lineWidth',1,'color',[0.3, 0.3, 0.3]);
        end
    end
    
    s1 = scatter(centroids(:,1),centroids(:,2),24,'k','filled');
    axis square
    set(gca,'XColor', 'none','YColor','none')
    try
        legend([l1,l3,s1],{'Sig. correlated during rest','Sig. correlated during both','centroid'})
    catch
    end
    title(sprintf('%s - Day %i - %s',mouseName, day, mType{ismember(mouseName,miceKO)+1}))
    saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\e2_correlation_network_maps',sprintf('Split low - %s - Day %i - %s.fig',mouseName, day, mType{ismember(mouseName,miceKO)+1})));
    saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\e2_correlation_network_maps',sprintf('Split low - %s - Day %i - %s.png',mouseName, day, mType{ismember(mouseName,miceKO)+1})));

    close all
end
