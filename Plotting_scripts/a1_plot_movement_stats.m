% Plots of WT vs KO for movement stats

addpath(genpath('U:\eng_research_handata\Athif Mohamed\nexmif_paper\Utils\'))  % Add utilities
init % Initiate variables
close all

% Load data 
load(fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\stats\movement_stats','movStats_08_22'),'movStats');
movStatsTable = struct2table(movStats);
miceAll = movStatsTable.animal;
miceAllUnique = unique(miceAll);
movStatsTableSummary = movStatsTable(1:15,[1,3:end]);

% Summarize animal wise 
for m = 1:numel(miceAllUnique)
    miceIdx = find(ismember(movStatsTable.animal,miceAllUnique(m)));
    movStatsTableSummary{m,3:end} = mean(movStatsTable{miceIdx,4:end},1,'omitnan');
    movStatsTableSummary{m,1} = miceAllUnique(m);
    movStatsTableSummary{m,2} =  movStatsTable.gType(miceIdx(1));
end

fieldNames = movStatsTableSummary.Properties.VariableNames;
featureNames = fieldNames(3:end);

featureNames(1:5) = {'Average speed (cm/s)', 'Movement bout duration(s)','Number of movement bouts','Number of movement onsets','Movement bout speed(cm/s)'};
featureNamesPlot = featureNames;
featureNamesPlot(1:5) = {'Average speed', 'Average movement bout duration','Number of movement bouts','Number of movement onsets','Average speed during movement bout'};

% Plot all 5 features 
stats = struct();
for f = 3%1:5
    figure
    
    data = movStatsTableSummary{:,f+2};
    group = ismember([movStatsTableSummary.animal],miceStudy(maskKO));
    boxplot(data,group,'Colors','k')
    axis square
    
    hold on
    xValue = 0;
    xx =[1,2];
    
    for g = 1:2 % WT, KO
        if g == 1
            maskG = ismember([movStatsTableSummary.animal],miceStudy(maskWT));
        else
            maskG = ismember([movStatsTableSummary.animal],miceStudy(maskKO));
        end
        y = data(maskG);
        x = g*ones(size(y));
        scatter (x,y,'filled','markerFaceColor',colors(2*g,:))
        hold on
    end
    
    % Stats
    maskWTT = ismember([movStatsTableSummary.animal],miceStudy(maskWT));
    maskKOT = ismember([movStatsTableSummary.animal],miceStudy(maskKO));
    yKO = data(maskKOT);
    yWT = data(maskWTT);
    
    yKOIn = yKO(~isoutlier(yKO));
    yWTIn = yWT(~isoutlier(yWT));
    
    [p_rs,h] = ranksum(yKO,yWT); % returns the p-value of a two-sided Wilcoxon rank sum test. ranksum tests the null hypothesis that data in x and y are samples from continuous distributions with equal medians, against the alternative that they are not. The test assumes that the two samples are independent. x and y can have different lengths.
    [h_t,p_t] = ttest2(yKO,yWT);
    
%     up = 1.05*max([yKO',yWT']);
%     hold on
%     text(1.5-0.5,up,sprintf('p-ranksum = %.3f',p_rs),'HorizontalAlignment','center')
%     hold on
%     text(1.5+0.5,up,sprintf('p-ttest = %.3f',p_t),'HorizontalAlignment','center')
    
    ylabel(sprintf('%s', featureNames{f}))
    set(gca,'xtick',xx,'xticklabel',{'WT','KO'},'TickDir', 'out','TickLength',[0.03, 0.025], 'Color','none','LineWidth',1,'YLim',[0 inf])
    box off
    set(gcf,'Units','inches','InnerPosition',[5, 5, 3, 3])
    saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\a1_movement_stats\8.3.22',sprintf('%s.fig',featureNamesPlot{f})))
    saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\a1_movement_stats\8.3.22',sprintf('%s.png',featureNamesPlot{f})))
    saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\a1_movement_stats\8.3.22',sprintf('%s.epsc',featureNamesPlot{f})))

    stats(f).feature = featureNames{f};
    stats(f).h =h_t;
    stats(f).p =p_t;
    stats(f).meanWT = mean(yWT);
    stats(f).meanKO = mean(yKO);
    stats(f).stdWT = std(yWT);
    stats(f).stdKO = std(yKO);
end 
