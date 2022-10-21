%% Description

% What happens to significantly correlated pairs when states change % speed low -> speed  high -- corr
% decreases and vice versa
%% Initialize

addpath(genpath('U:\eng_research_handata\Athif Mohamed\nexmif_paper\Utils'))  % Add utilities
init % Initialize data directories and genotypes
% close all

%% collect all values 

t_start = tic;
% load correlation stats
load(fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\stats\correlation_ball','corrStatsTableSpeedPearson_08_22'),'corrStatsTable')
fieldNames = corrStatsTable.Properties.VariableNames;
corrStatsTable([6,13,18:19,29,31:35,40:43,45],:) = [];

miceBad = {};
lists(1).Gtype = 'WT';
lists(1).W1 = []; lists(1).W2 = []; lists(1).W3 = [];lists(1).W4 = [];lists(1).W5 = [];lists(1).W6 = [];lists(1).frac = [];lists(1).n_cells_both = [];lists(1).n_cells_either = [];

lists(2).Gtype = 'KO';
lists(2).W1 = []; lists(2).W2 = []; lists(2).W3 = [];lists(2).W4 = [];lists(2).W5 = [];lists(2).W6 = [];lists(2).frac = [];lists(2).n_cells_both = [];lists(2).n_cells_either = [];

lists_mean = lists;
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
    
    w_corr_high = corrStatsTable{k,7}{:};
    %     w_corr_high(isnan(w_corr_high)) = 0;
    w_corr_high(logical(triu(ones(n)))) = nan;
    w_thresh_high_pos = corrStatsTable{k,8}{:};
    w_thresh_high_neg = corrStatsTable{k,9}{:};
    adj_sig_corr_high_pos = (w_corr_high - w_thresh_high_pos > 1e-5) & (y_close == 0);
    idx_high_pos = find(adj_sig_corr_high_pos);
    idx_both_pos = intersect(idx_low_pos,idx_high_pos);
    idx_either_pos = union(idx_low_pos,idx_high_pos);
    % collect 
    frac_low_pos = numel(idx_low_pos)/(n*(n-1))*100;
    frac_high_pos = numel(idx_high_pos)/(n*(n-1))*100;
    frac_both_pos = numel(idx_both_pos)/numel(idx_either_pos)*100;
    
    if ismember(mouseName,miceWT)
        % take sig low speed correlation cells in a column, get the correlation of the
        % same cell in high
        lists(1).W1 = [lists(1).W1 ;w_corr_low(idx_low_pos)];
        lists(1).W2 = [lists(1).W2;w_corr_high(idx_low_pos)];
        
        % take sig high speed correlation cells in a column, get the correlation of the
        % same cell in low
        
        lists(1).W3 = [lists(1).W3;w_corr_high(idx_high_pos)];
        lists(1).W4 = [lists(1).W4;w_corr_low(idx_high_pos)];
        
        % Both
        lists(1).W5 = [lists(1).W5;w_corr_low(idx_both_pos)];
        lists(1).W6 = [lists(1).W6;w_corr_high(idx_both_pos)];
        
        lists(1).frac = [lists(1).frac, frac_both_pos];
        lists(1).n_cells_both = [lists(1).n_cells_both, numel(idx_both_pos)];
        lists(1).n_cells_either = [lists(1).n_cells_either, numel(idx_either_pos)];


        % means
        lists_mean(1).W1 = [lists_mean(1).W1; mean(w_corr_low(idx_low_pos))];
        lists_mean(1).W2 = [lists_mean(1).W2; mean(w_corr_high(idx_low_pos))];
        
        lists_mean(1).W3 = [lists_mean(1).W3; mean(w_corr_high(idx_high_pos))];
        lists_mean(1).W4 = [lists_mean(1).W4; mean(w_corr_low(idx_high_pos))];
        
        lists_mean(1).W5 = [lists_mean(1).W5; mean(w_corr_low(idx_both_pos))];
        lists_mean(1).W6 = [lists_mean(1).W6; mean(w_corr_high(idx_both_pos))];
        
        
    else
        
        lists(2).W1 = [lists(2).W1 ;w_corr_low(idx_low_pos)];
        lists(2).W2 = [lists(2).W2;w_corr_high(idx_low_pos)];
        lists(2).W3 = [lists(2).W3;w_corr_high(idx_high_pos)];
        lists(2).W4 = [lists(2).W4;w_corr_low(idx_high_pos)];
        lists(2).W5 = [lists(2).W5;w_corr_low(idx_both_pos)];
        lists(2).W6 = [lists(2).W6;w_corr_high(idx_both_pos)];
        lists(2).frac = [lists(2).frac, frac_both_pos];
        lists(2).n_cells_both = [lists(2).n_cells_both, numel(idx_both_pos)];
        lists(2).n_cells_either = [lists(2).n_cells_either, numel(idx_either_pos)];

        % means
        
        
        lists_mean(2).W1 = [lists_mean(2).W1; mean(w_corr_low(idx_low_pos))];
        lists_mean(2).W2 = [lists_mean(2).W2; mean(w_corr_high(idx_low_pos))];
        lists_mean(2).W3 = [lists_mean(2).W3; mean(w_corr_high(idx_high_pos))];
        lists_mean(2).W4 = [lists_mean(2).W4; mean(w_corr_low(idx_high_pos))];
        lists_mean(2).W5 = [lists_mean(2).W5; mean(w_corr_low(idx_both_pos))];
        lists_mean(2).W6 = [lists_mean(2).W6; mean(w_corr_high(idx_both_pos))];
    end
    
end
%% Overall correlated vs uncorrelated 


%%  a - correlated sig in low - what happens in high
va = struct();
va.WTRestSig = lists(1).W1;
va.WTRun = lists(1).W2;
va.KORestSig = lists(2).W1;
va.KORun = lists(2).W2;
figure
va_p = violinplot(va,[],'ShowMean',true);
va_p(1,1).ViolinColor =  colors(1,:);
va_p(1,2).ViolinColor = colors(2,:);
va_p(1,1).BoxColor =  colors(1,:);
va_p(1,2).BoxColor =  colors(2,:);
va_p(1,1).EdgeColor =  colors(1,:);
va_p(1,2).EdgeColor =  colors(2,:);
va_p(1,3).ViolinColor =  colors(3,:);
va_p(1,4).ViolinColor = colors(4,:);
va_p(1,3).BoxColor =  colors(3,:);
va_p(1,4).BoxColor =  colors(4,:);
va_p(1,3).EdgeColor =  colors(3,:);
va_p(1,4).EdgeColor =  colors(4,:);



[h_WT,p_WT] = ttest(va.WTRestSig,va.WTRun);
[h_KO,p_KO] = ttest(va.KORestSig,va.KORun);

mean_std = [mean(va.WTRestSig,'omitnan'), std(va.WTRestSig,'omitnan');...
            mean(va.WTRun,'omitnan'), std(va.WTRun,'omitnan');...
            mean(va.KORestSig,'omitnan'), std(va.KORestSig,'omitnan');...
            mean(va.KORun,'omitnan'), std(va.KORun,'omitnan')];

yy = ylim;
% text(1.5,0.9*yy(2),sprintf('p = %.3f',p_WT),'HorizontalAlignment','Center')
% text(3.5,0.9*yy(2),sprintf('p = %.3f',p_KO),'HorizontalAlignment','Center')
title({'Correlation of cells that are significantly correlated' 'During rest'})
ylabel('Correlation strength')
set(gcf,'Color','none')
set(gca,'Units','inches','InnerPosition',[.8 .5 3 2],'TickDir','out','TickLength',[0.03, 0.025],'Color','none','Box','off','LineWidth',1)
% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\f_correlation_states','Cellwise significantly correlated during rest.png'));
% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\f_correlation_states','Cellwise significantly correlated during rest.fig'));
% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\f_correlation_states','Cellwise significantly correlated during rest.epsc'));
% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\f_correlation_states','Cellwise significantly correlated during rest.svg'));

%%
% sessionwise

figure
WTRestSig = lists_mean(1).W1;
WTRun = lists_mean(1).W2;
KORestSig = lists_mean(2).W1;
KORun = lists_mean(2).W2;
b2 = boxchart(2*ones(size(WTRun')),...
    WTRun',...
    'MarkerStyle','none',...
    'BoxFaceColor',colors(2,:),...
    'BoxFaceAlpha',0,...
    'WhiskerLineColor', colors(2,:),...
    'LineWidth',2);
hold on
b1 = boxchart(ones(size(WTRestSig')),...
    WTRestSig',...
    'MarkerStyle','none',...
    'BoxFaceColor',colors(1,:),...
    'BoxFaceAlpha',0,...
    'WhiskerLineColor', colors(1,:),...
    'LineWidth',2);
hold on

s2 = scatter(2,WTRun,[],colors(2,:),'filled');
hold on
s1 = scatter(1,WTRestSig,[],colors(1,:),'filled');
hold on
line([1,2],[WTRestSig,WTRun],'color','k')
b4 = boxchart(4*ones(size(KORun')),...
    KORun',...
    'MarkerStyle','none',...
    'BoxFaceColor',colors(4,:),...
    'BoxFaceAlpha',0,...
    'WhiskerLineColor', colors(4,:),...
    'LineWidth',2);
hold on
b3 = boxchart(3*ones(size(KORestSig')),...
    KORestSig',...
    'MarkerStyle','none',...
    'BoxFaceColor',colors(3,:),...
    'BoxFaceAlpha',0,...
    'WhiskerLineColor', colors(3,:),...
    'LineWidth',2);
hold on
s4 = scatter(4,KORun,[],colors(4,:),'filled');
hold on
s3 = scatter(3,KORestSig,[],colors(3,:),'filled');
hold on
line([3,4],[KORestSig,KORun],'color','k')

xlim([0.5,4.5])
% ylim([0,27])
xticks([1,2,3,4])
xticklabels({'Rest-Sig','Run','Rest-Sig','Run'})
title({'Correlation of cells that are significantly correlated' 'During rest'})
ylabel('Mean correlation')


[p,tbl,stats] = anova1([WTRestSig',WTRun',KORestSig',KORun'],...
    [ones(size(WTRestSig')),2*ones(size(WTRun')),3*ones(size(KORestSig')),4*ones(size(KORun'))],'off');
[hWT,pWT] = ttest(WTRestSig,WTRun);
[hKO,pKO] = ttest(KORestSig,KORun);
mean_std = [mean(WTRestSig,'omitnan'), std(WTRestSig,'omitnan');...
            mean(WTRun,'omitnan'), std(WTRun,'omitnan');...
            mean(KORestSig,'omitnan'), std(KORestSig,'omitnan');...
            mean(KORun,'omitnan'), std(KORun,'omitnan')];

% yy = ylim;
% text(1.5,0.95*yy(2),sprintf('p* = %.3f',pWT),'HorizontalAlignment','Center')
% text(3.5,0.95*yy(2),sprintf('p* = %.3f',pKO),'HorizontalAlignment','Center')

% set(gcf,'Color','none')
set(gca,'Units','inches','InnerPosition',[.8 .5 3 2],'TickDir','out','TickLength',[0.03, 0.025],'Color','none','Box','off','LineWidth',1,'YLim',[0 0.4])
saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\f_correlation_states\8.3.22','Sessionwise significantly correlated during rest.png'));
saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\f_correlation_states\8.3.22','Sessionwise significantly correlated during rest.fig'));
saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\f_correlation_states\8.3.22','Sessionwise significantly correlated during rest.epsc'));


%% b - correlated sig in high - what happens in low
vb = struct();
vb.WTRunSig = lists(1).W3;
vb.WTRest = lists(1).W4;
vb.KORunSig = lists(2).W3;
vb.KORest = lists(2).W4;
figure
vb_p = violinplot(vb,[],'ShowMean',true);
vb_p(1,1).ViolinColor =  colors(2,:);
vb_p(1,2).ViolinColor = colors(1,:);
vb_p(1,1).BoxColor =  colors(2,:);
vb_p(1,2).BoxColor =  colors(1,:);
vb_p(1,1).EdgeColor =  colors(2,:);
vb_p(1,2).EdgeColor =  colors(1,:);
vb_p(1,3).ViolinColor =  colors(4,:);
vb_p(1,4).ViolinColor = colors(3,:);
vb_p(1,3).BoxColor =  colors(4,:);
vb_p(1,4).BoxColor =  colors(3,:);
vb_p(1,3).EdgeColor =  colors(4,:);
vb_p(1,4).EdgeColor =  colors(3,:);

[h_WT,p_WT] = ttest(vb.WTRunSig,vb.WTRest);
[h_KO,p_KO] = ttest(vb.KORunSig,vb.KORest);

mean_std = [mean(vb.WTRest,'omitnan'), std(vb.WTRest,'omitnan');...
            mean(vb.WTRunSig,'omitnan'), std(vb.WTRunSig,'omitnan');...
            mean(vb.KORest,'omitnan'), std(vb.KORest,'omitnan');...
            mean(vb.KORunSig,'omitnan'), std(vb.KORunSig,'omitnan')];

% yy = ylim;
% text(1.5,0.9*yy(2),sprintf('p = %.3f',p_WT),'HorizontalAlignment','Center')
% text(3.5,0.9*yy(2),sprintf('p = %.3f',p_KO),'HorizontalAlignment','Center')

title({'Correlation of cells that are significantly correlated' 'During run'})
ylabel('Correlation strength')
set(gcf,'Color','none')
set(gca,'Units','inches','InnerPosition',[.8 .5 3 2],'TickDir','out','TickLength',[0.03, 0.025],'Color','none','Box','off','LineWidth',1)
% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\f_correlation_states','Cellwise significantly correlated during run.png'));
% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\f_correlation_states','Cellwise significantly correlated during run.fig'));
% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\f_correlation_states','Cellwise significantly correlated during run.epsc'));
% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\f_correlation_states','Cellwise significantly correlated during run.svg'));

%% session wise
figure
WTRunSig = lists_mean(1).W3;
WTRest = lists_mean(1).W4;
KORunSig = lists_mean(2).W3;
KORest = lists_mean(2).W4;

b1 = boxchart(1*ones(size(WTRunSig')),...
    WTRunSig',...
    'MarkerStyle','none',...
    'BoxFaceColor',colors(2,:),...
    'BoxFaceAlpha',0,...
    'WhiskerLineColor', colors(2,:),...
    'LineWidth',2);
hold on
b2 = boxchart(2*ones(size(WTRest')),...
    WTRest',...
    'MarkerStyle','none',...
    'BoxFaceColor',colors(1,:),...
    'BoxFaceAlpha',0,...
    'WhiskerLineColor', colors(1,:),...
    'LineWidth',2);
hold on

s1 = scatter(1,WTRunSig,[],colors(2,:),'filled');
hold on
s2 = scatter(2,WTRest,[],colors(1,:),'filled');
hold on
line([1,2],[WTRunSig,WTRest],'color','k')
b3 = boxchart(3*ones(size(KORunSig')),...
    KORunSig',...
    'MarkerStyle','none',...
    'BoxFaceColor',colors(4,:),...
    'BoxFaceAlpha',0,...
    'WhiskerLineColor', colors(4,:),...
    'LineWidth',2);
hold on
b4 = boxchart(4*ones(size(KORest')),...
    KORest',...
    'MarkerStyle','none',...
    'BoxFaceColor',colors(3,:),...
    'BoxFaceAlpha',0,...
    'WhiskerLineColor', colors(3,:),...
    'LineWidth',2);
hold on
s3 = scatter(3,KORunSig,[],colors(4,:),'filled');
hold on
s4 = scatter(4,KORest,[],colors(3,:),'filled');
hold on
line([3,4],[KORunSig,KORest],'color','k')

xlim([0.5,4.5])
% ylim([0,27])
xticks([1,2,3,4])
xticklabels({'Run-Sig','Rest','Run-Sig','Rest'})
ylabel('Mean correlation')
title({'Correlation of cells that are significantly correlated' 'During run'})
% set(gca,'Units','inches','InnerPosition',[.5 .3 4 4],'TickDir','out','TickLength',[0.03, 0.025],'Box','off','LineWidth',2)
[p,tbl,stats] = anova1([WTRest',WTRunSig',KORest',KORunSig'],...
    [ones(size(WTRest')),2*ones(size(WTRunSig')),3*ones(size(KORest')),4*ones(size(KORunSig'))],'off');
[hWT,pWT] = ttest(WTRest,WTRunSig);
[hKO,pKO] = ttest(KORest,KORunSig);
mean_std = [mean(WTRest,'omitnan'), std(WTRest,'omitnan');...
            mean(WTRunSig,'omitnan'), std(WTRunSig,'omitnan');...
            mean(KORest,'omitnan'), std(KORest,'omitnan');...
            mean(KORunSig,'omitnan'), std(KORunSig,'omitnan')];

% yy = ylim;
% text(1.5,0.95*yy(2),sprintf('p* = %.3f',pWT),'HorizontalAlignment','Center')
% text(3.5,0.95*yy(2),sprintf('p* = %.3f',pKO),'HorizontalAlignment','Center')

set(gcf,'Color','none')
set(gca,'Units','inches','InnerPosition',[.8 .5 3 2],'TickDir','out','TickLength',[0.03, 0.025],'Color','none','Box','off','LineWidth',1,'YLim',[0 0.4])
saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\f_correlation_states\8.3.22','Sessionwise significantly correlated during run.png'));
saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\f_correlation_states\8.3.22','Sessionwise significantly correlated during run.fig'));
saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\f_correlation_states\8.3.22','Sessionwise significantly correlated during run.epsc'));

%% c - correlated sig in both - what happens between low and high
vc = struct();
vc.WTRestBothSig = lists(1).W5;
vc.WTRunBothSig = lists(1).W6;
vc.KORestBothSig = lists(2).W5;
vc.KORunBothSig = lists(2).W6;
figure
vc_p = violinplot(vc,[],'ShowMean',true);
vc_p(1,1).ViolinColor =  colors(1,:);
vc_p(1,2).ViolinColor = colors(2,:);
vc_p(1,1).BoxColor =  colors(1,:);
vc_p(1,2).BoxColor =  colors(2,:);
vc_p(1,1).EdgeColor =  colors(1,:);
vc_p(1,2).EdgeColor =  colors(2,:);
vc_p(1,3).ViolinColor =  colors(3,:);
vc_p(1,4).ViolinColor = colors(4,:);
vc_p(1,3).BoxColor =  colors(3,:);
vc_p(1,4).BoxColor =  colors(4,:);
vc_p(1,3).EdgeColor =  colors(3,:);
vc_p(1,4).EdgeColor =  colors(4,:);

title({'Correlation of cells that are significantly correlated' 'During both rest and run'})
ylabel('Correlation Strength')
set(gca,'Units','inches','InnerPosition',[.8 .5 3 2],'TickDir','out','TickLength',[0.03, 0.025],'Color','none','Box','off','LineWidth',1)
% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\f_correlation_states','Cellwise significantly correlated during both run and rest.png'));
% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\f_correlation_states','Cellwise significantly correlated during both run and rest.fig'));
% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\f_correlation_states','Cellwise significantly correlated during both run and rest.epsc'));
% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\f_correlation_states','Cellwise significantly correlated during both run and rest.svg'));

[h_WT,p_WT] = ttest(vc.WTRestBothSig,vc.WTRunBothSig);
[h_KO,p_KO] = ttest(vc.KORestBothSig,vc.KORunBothSig);
mean_std = [mean(vc.WTRestBothSig,'omitnan'), std(vc.WTRestBothSig,'omitnan');...
            mean(vc.WTRunBothSig,'omitnan'), std(vc.WTRunBothSig,'omitnan');...
            mean(vc.KORestBothSig,'omitnan'), std(vc.KORestBothSig,'omitnan');...
            mean(vc.KORunBothSig,'omitnan'), std(vc.KORunBothSig,'omitnan')];


yy = ylim;
text(1.5,0.9*yy(2),sprintf('p = %.3f',p_WT),'HorizontalAlignment','Center')
text(3.5,0.9*yy(2),sprintf('p = %.3f',p_KO),'HorizontalAlignment','Center')

% sesion
figure
WTRestBothSig = lists_mean(1).W5;
WTRunBothSig = lists_mean(1).W6;
KORestBothSig = lists_mean(2).W5;
KORunBothSig = lists_mean(2).W6;
b2 = boxchart(2*ones(size(WTRunBothSig')),...
    WTRunBothSig',...
    'MarkerStyle','none',...
    'BoxFaceColor',colors(2,:),...
    'BoxFaceAlpha',0,...
    'WhiskerLineColor', colors(2,:),...
    'LineWidth',2);
hold on
b1 = boxchart(ones(size(WTRestBothSig')),...
    WTRestBothSig',...
    'MarkerStyle','none',...
    'BoxFaceColor',colors(1,:),...
    'BoxFaceAlpha',0,...
    'WhiskerLineColor', colors(1,:),...
    'LineWidth',2);
hold on

s2 = scatter(2,WTRunBothSig,[],colors(2,:),'filled');
hold on
s1 = scatter(1,WTRestBothSig,[],colors(1,:),'filled');
hold on
line([1,2],[WTRestBothSig,WTRunBothSig],'color','k')
b4 = boxchart(4*ones(size(KORunBothSig')),...
    KORunBothSig',...
    'MarkerStyle','none',...
    'BoxFaceColor',colors(4,:),...
    'BoxFaceAlpha',0,...
    'WhiskerLineColor', colors(4,:),...
    'LineWidth',2);
hold on
b3 = boxchart(3*ones(size(KORestBothSig')),...
    KORestBothSig',...
    'MarkerStyle','none',...
    'BoxFaceColor',colors(3,:),...
    'BoxFaceAlpha',0,...
    'WhiskerLineColor', colors(3,:),...
    'LineWidth',2);
hold on
s4 = scatter(4,KORunBothSig,[],colors(4,:),'filled');
hold on
s3 = scatter(3,KORestBothSig,[],colors(3,:),'filled');
hold on
line([3,4],[KORestBothSig,KORunBothSig],'color','k')

xlim([0.5,4.5])
% ylim([0,27])
xticks([1,2,3,4])
xticklabels({'Rest-Sig','Run-sig','Rest-Sig','Run-sig'})
ylabel('Mean correlation')
title({'Correlation of cells that are significantly correlated' 'During both rest and run'})

% set(gcf,'Color','none')
% set(gca,'Units','inches','InnerPosition',[.5 .3 4 4],'TickDir','out','TickLength',[0.03, 0.025],'Box','off','LineWidth',2)
[p,tbl,stats] = anova1([WTRestBothSig',WTRunBothSig',KORestBothSig',KORunBothSig'],...
    [ones(size(WTRestBothSig')),2*ones(size(WTRunBothSig')),3*ones(size(KORestBothSig')),4*ones(size(KORunBothSig'))],'off');
[hWT,pWT] = ttest(WTRestBothSig,WTRunBothSig);
[hKO,pKO] = ttest(KORestBothSig,KORunBothSig);
mean_std = [mean(WTRestBothSig,'omitnan'), std(WTRestBothSig,'omitnan');...
            mean(WTRunBothSig,'omitnan'), std(WTRunBothSig,'omitnan');...
            mean(KORestBothSig,'omitnan'), std(KORestBothSig,'omitnan');...
            mean(KORunBothSig,'omitnan'), std(KORunBothSig,'omitnan')];



yy = ylim;
% text(1.5,0.95*yy(2),sprintf('p* = %.3f',pWT),'HorizontalAlignment','Center')
% text(3.5,0.95*yy(2),sprintf('p* = %.3f',pKO),'HorizontalAlignment','Center')

 set(gcf,'Color','none')
set(gca,'Units','inches','InnerPosition',[.8 .5 3 2],'TickDir','out','TickLength',[0.03, 0.025],'Color','none','Box','off','LineWidth',1,'YLim',[0 0.4])
saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\f_correlation_states\8.3.22','Sessionwise significantly correlated during both run and rest.png'));
saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\f_correlation_states\8.3.22','Sessionwise significantly correlated during both run and rest.fig'));
saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\f_correlation_states\8.3.22','Sessionwise significantly correlated during both run and rest.epsc'));


%% Fraction of cell pairs correlated in both 


frac_mean_std = [mean(lists(1).frac,'omitnan'), std(lists(1).frac,'omitnan');...
            mean(lists(2).frac,'omitnan'), std(lists(2).frac,'omitnan')];


%% Fisher tests for # cell pairs tha tare correlated in both states vs in either state

% Pairs of cells significantly correlated in...
%_______________________________
%|          |                   |
%|  WT_Both | WT_only_one_state |
%|__________|___________________|
%|          |                   |
%|  KO_Both | KO_only_one_state | 
%|__________|___________________|


fisher_struct.fisherMat = zeros(2);
fisher_struct.fisherMat(1,1) = sum(lists(1).n_cells_both);
fisher_struct.fisherMat(1,2) = sum(lists(1).n_cells_either)-sum(lists(1).n_cells_both);
fisher_struct.fisherMat(2,1) = sum(lists(2).n_cells_both);
fisher_struct.fisherMat(2,2) = sum(lists(2).n_cells_either)-sum(lists(2).n_cells_both);


[fisher_struct.h_f,...
    fisher_struct.p_f,...
    fisher_struct.stats_f] = fishertest(fisher_struct.fisherMat,'Alpha',0.05)%,'Tail','right'

fisher_struct.pr_1 = fisher_struct.fisherMat(1,1)/(fisher_struct.fisherMat(1,1)+ fisher_struct.fisherMat(1,2))*100;
fisher_struct.SEP_1 = 1.96*sqrt(fisher_struct.pr_1*(100-fisher_struct.pr_1)/(fisher_struct.fisherMat(1,1)+ fisher_struct.fisherMat(1,2)));
fisher_struct.pr_2 = fisher_struct.fisherMat(2,1)/(fisher_struct.fisherMat(2,1)+ fisher_struct.fisherMat(2,2))*100;
fisher_struct.SEP_2 = 1.96*sqrt(fisher_struct.pr_2*(100-fisher_struct.pr_2)/(fisher_struct.fisherMat(2,1)+ fisher_struct.fisherMat(2,2)));


figure
b = bar([1,2],[fisher_struct.pr_1,...
    fisher_struct.pr_2]);
b.FaceColor = 'none';
b.EdgeColor = 'flat';
b.LineWidth = 2;
b.CData = colors([2,4],:);
hold on
xticks(gca,[1,2,3,4])
xticklabels(gca, {'WT-Both/WT-Either','KO-Both/KO-Either'})
er = errorbar([1,2],[fisher_struct.pr_1,fisher_struct.pr_2],[fisher_struct.SEP_1,fisher_struct.SEP_2]);
er.Color = [0 0 0];
er.LineStyle = 'none';
er.LineWidth = 2;
xlim([0.5,2.5])
yy = ylim;


text(1.5,0.95*yy(2),sprintf('p* = %.4f',fisher_struct.p_f),'HorizontalAlignment','Center')
% text(3.5,0.95*yy(2),sprintf('p* = %.4f',fisher_struct(4).p_f*4),'HorizontalAlignment','Center')
% text(2,1.15*yy(2),sprintf('p* = %.4f',fisher_struct(1).p_f*4),'HorizontalAlignment','Center')
% text(3,1.15*yy(2),sprintf('p* = %.4f',fisher_struct(2).p_f*4),'HorizontalAlignment','Center')
ylim([0,yy(2)*1.2])
ylabel('Significantly correlated cells (%)')
% set(gcf,'Color','none')
set(gca,'Units','inches','InnerPosition',[.8 .5 3 2],'TickDir','out','TickLength',[0.03, 0.025],'Color','none','Box','off','LineWidth',1)

saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\f_correlation_states\8.3.22',sprintf('Pearson frac cells corr in rest and run Fisher.fig')));
saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\f_correlation_states\8.3.22',sprintf('Pearson frac cells corr in rest and run Fisher.png')));
saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\f_correlation_states\8.3.22',sprintf('Pearson frac cells corr in rest and run Fisher.epsc')));


