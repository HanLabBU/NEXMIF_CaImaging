%% Description

% xxxxxx
%% Initialize

addpath(genpath('U:\eng_research_handata\Athif Mohamed\nexmif_paper\Utils\'))  % Add utilities
init % Initialize data directories and genotypes
close all

%% Begin code
t_start = tic;
% load correlation stats
load(fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\stats\correlation_platform','corrStatsTableAsym'),'corrStatsTable')
fieldNames = corrStatsTable.Properties.VariableNames;


miceBad = {};
listCorr = [];
listisKO = [];
listEs = [];
for k =  1:numel(corrStatsTable.animal)
    mouseName = corrStatsTable{k,1}{:};
    day = corrStatsTable{k,2}{:};
    
    mPath = fullfile(dataFolder,sprintf('fullData_%s_D%i_%s',mouseName,day,'platform'));
    
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
    y_close = double(corrStatsTable{k,7}{:}<20);
    y_close(logical(triu(ones(n)))) = nan;
    
    w_corr = corrStatsTable{k,4}{:};
    w_corr(isnan(w_corr)) = 0;
    w_corr(logical(triu(ones(n)))) = nan;
    w_thresh_pos = corrStatsTable{k,5}{:};
    w_thresh_neg = corrStatsTable{k,6}{:};
    adj_sig_corr_pos = (w_corr > w_thresh_pos) & (y_close == 0);
    idx_pos = find(adj_sig_corr_pos);
    adj_sig_corr_neg = w_corr < w_thresh_neg & (y_close == 0);
    idx_neg = find(adj_sig_corr_neg);
    adj_sig_corr_non_sig = (w_corr < w_thresh_pos) & (w_corr > w_thresh_neg) & (y_close == 0);
%     adj_sig_corr_non_sig = (w_corr < w_thresh_pos) & (y_close == 0);
    
    %%
    
    listCorr = [listCorr,numel(idx_pos)];
    listEs = [listEs, n*(n-1)/2 - sum(y_close(:) == 1)];
    listisKO = [listisKO,ismember(mouseName,miceKO)];
    %% plot 
   
end
listisKO = logical(listisKO);
save('identity_stats_platform','listisKO','listCorr','listEs')
%%

animal = table2cell(corrStatsTable(:,1));
day1 = table2cell(corrStatsTable(:,2));
labels = [animal';day1'];
labels = sprintf('%s-%d ',labels{:});
labels = strsplit(labels,' ');
labels = categorical(labels(1:45));
figure
bar(labels,listCorr)
title('Number of sig.corr cell pairs in each session')
figure
bar(labels,listCorr./listEs)
title('Percentage of sig.corr cell pairs in each session')
