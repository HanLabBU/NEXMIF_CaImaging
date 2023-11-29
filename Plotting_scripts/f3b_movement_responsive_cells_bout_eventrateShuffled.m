% Based on

% Movement_responsive_cells
% before after T test stats

%  based on ext 6
%  use event rates within movement bouts vs not T test
%% code
addpath(genpath('U:\eng_research_handata\Athif Mohamed\nexmif_paper\Utils\'))  % Add utilities
init % Initialize data directories and genotypes
close all

%% Plot

load('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\stats\movement_responsive_cells\responsive_cells_bout_shuffle_08_22','responsive_cells')

responsive_cells_table = struct2table(responsive_cells);
miceAll = responsive_cells_table.animal;
miceAllUnique = unique(miceAll);

% Boxplot for animals

featureNames = {'High movement responsive cells','Low movement responsive cells'};
isWT_pop = logical(responsive_cells_table.isWT);
binEdge = 0:0.2:1;

cols = [8,12];
nEvCols = [5,5]; % number of events


for colIdx = 1 % 1 - high movement, 2 - low movement
    
    featureNames{colIdx}
    nEventsAll = responsive_cells_table{:,nEvCols(colIdx)};
    IdxEv = nEventsAll> 5;
    
    
    col = cols(colIdx);
    
    % cells care vs not for sessions
    careListPercent = [];
    noCareListPercent = [];
    isWTList = [];
    nList = [];
    nameList = {};
    careList = [];
    noCareList = [];
    mList ={};
    for m = 1:numel(miceAllUnique)
        idxM = strcmp(miceAll,miceAllUnique(m));
        isWT = responsive_cells_table{find(idxM,1,'first'),4};
        
        idx = idxM & IdxEv;
        if sum(idx) >= 50
            do_cells_care = [responsive_cells_table{idx,col}];
            
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
        
        if ~isempty(idx)
            mList = [mList;miceAllUnique{m}]
        end
    end
    
    
    % update fisher mat
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
    b = bar([-1,1],[pr_WT,pr_KO]);
    b.FaceColor = 'none';
    b.EdgeColor = 'flat';
    b.LineWidth = 2;
    b.CData = colors([2,4],:);
    hold on
    xticks(gca,[-1 1])
    xticklabels(gca, {'WT','KO'})
    er = errorbar([-1,1],[pr_WT,pr_KO],1.96*[SEP_WT,SEP_KO]);
    er.Color = [0 0 0];
    er.LineStyle = 'none';
    er.LineWidth = 2;
    xlim([-2.5,2.5])
    yy = ylim;
    %      text(0,yy(2)*0.99,sprintf('Fisher sig: p = %d',p_f),'HorizontalAlignment','center')
    ylim([0 yy(2)*1.01]);
    ylabel('Modulated cells (%)')
    set(gcf,'Color','none')
    set(gca,'Units','inches','InnerPosition',[.5 .3 2 2],'TickDir','out','TickLength',[0.03, 0.025],'Color','none','Box','off','LineWidth',2)
    
%         saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\d_movement_responsive_cells\8.3.22',sprintf('shuffled %s fisher.png',featureNames{colIdx})))
%         saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\d_movement_responsive_cells\8.3.22',sprintf('shuffled %s fisher.fig',featureNames{colIdx})))
%         saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\d_movement_responsive_cells\8.3.22',sprintf('shuffled %s fisher.epsc',featureNames{colIdx})))
%% plot as stacked bar
    figure
    fullSetWT = [pr_WT;pr_KO];
    fullSetKO = [100-pr_WT;100-pr_KO];
    b = bar([-1,1],[fullSetWT,fullSetKO],"stacked");
    b(1).FaceColor = 'flat';
    b(2).FaceColor = 'flat';
%     b.EdgeColor = 'flat';
%     b.LineWidth = 1;
    b(1).CData = colors([2,4],:);
    b(2).CData = [1 1 1;1 1 1];
    hold on
    xticks(gca,[-1 1])
    xticklabels(gca, {'WT','KO'})
    er = errorbar([-1,1],[pr_WT,pr_KO],1.96*[SEP_WT,SEP_KO],'CapSize',18);
    er.Color = [0 0 0];
    er.LineStyle = 'none';
    er.LineWidth = 0.5;
    xlim([-2.5,2.5])
%     yy = ylim;
%     %      text(0,yy(2)*0.99,sprintf('Fisher sig: p = %d',p_f),'HorizontalAlignment','center')
%     ylim([0 yy(2)*1.01]);
    ylabel('Percentage of Total Cells')
    set(gcf,'Color','none')
    set(gca,'Units','inches','InnerPosition',[.5 .3 2 2],'TickDir','out','TickLength',[0.03, 0.025],'Color','none','Box','off','LineWidth',2)
    
    saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\d_movement_responsive_cells\2.15.23',sprintf('stacked shuffled %s fisher.png',featureNames{colIdx})))
    saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\d_movement_responsive_cells\2.15.23',sprintf('stacked shuffled %s fisher.fig',featureNames{colIdx})))
    saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\d_movement_responsive_cells\2.15.23',sprintf('stacked shuffled %s fisher.epsc',featureNames{colIdx})))
    saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\d_movement_responsive_cells\2.15.23',sprintf('stacked shuffled %s fisher.pdf',featureNames{colIdx})))

%% plot sessionwise
    
    [GWT,MWT,DWT] = findgroups(responsive_cells_table{isWT_pop&IdxEv,1} , responsive_cells_table{isWT_pop&IdxEv,2});
    [GKO,MKO,DKO] = findgroups(responsive_cells_table{~isWT_pop&IdxEv,1} , responsive_cells_table{~isWT_pop&IdxEv,2});
    
    meanWT = splitapply(@count1s,responsive_cells_table{isWT_pop&IdxEv,col},GWT);
    meanKO = splitapply(@count1s,responsive_cells_table{~isWT_pop&IdxEv,col},GKO);
    
    figure
    b1 = boxchart(ones(size(meanWT)),...
        meanWT,...
        'MarkerStyle','none',...
        'BoxFaceColor',colors(2,:),...
        'BoxFaceAlpha',0,...
        'WhiskerLineColor', colors(2,:),...
        'LineWidth',2);
    hold on
    s1 = scatter(1,meanWT,[],colors(2,:),'filled');
    hold on
    b2 = boxchart(2*ones(size(meanKO)),...
        meanKO,...
        'MarkerStyle','none',...
        'BoxFaceColor',colors(4,:),...
        'BoxFaceAlpha',0,...
        'WhiskerLineColor', colors(4,:),...
        'LineWidth',2);
    hold on
    s2 = scatter(2,meanKO,[],colors(4,:),'filled');
    
    
    xlim([0.5,2.5])
    % ylim([0,27])
    xticks([1,2])
    xticklabels({'WT','KO'})
    
    
   
    % set(gcf,'Color','none')
    set(gca,'Units','inches','InnerPosition',[.5 .3 2 2],'TickDir','out','TickLength',[0.03, 0.025],'Color','none','Box','off','LineWidth',2)
    
    
    [h,p] = ttest2(meanWT,meanKO);
%     saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\d_movement_responsive_cells',sprintf('Sessionwise shuffled %s boxplot.png',featureNames{colIdx})))
%     saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\d_movement_responsive_cells',sprintf('Sessionwise shuffled %s boxplot.fig',featureNames{colIdx})))
%     saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\d_movement_responsive_cells',sprintf('Sessionwise shuffled %s boxplot.epsc',featureNames{colIdx})))
end

