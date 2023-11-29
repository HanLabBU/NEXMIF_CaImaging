% Plots of WT vs KO for movement stats

addpath(genpath('U:\eng_research_handata\Athif Mohamed\nexmif_paper\Utils\'))  % Add utilities
init % Initiate variables
close all

% Load data 
load(fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\stats\movement_stats','movStats_02_23'),'movStats');
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

featureNames([1:5,9]) = {'Average speed (cm/s)', 'Movement bout duration(s)','Number of movement bouts','Number of movement onsets','Movement bout speed(cm/s)', 'Distance (m)'};
featureNamesPlot = featureNames;
featureNamesPlot([1:5,9]) = {'Average speed', 'Average movement bout duration','Number of movement bouts','Number of movement onsets','Average speed during movement bout','Total distance traveled (m)'};

% Plot all 6features 
stats = struct();
for f = [1:5,9]
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
    
    [p_rs,h_rs] = ranksum(yKO,yWT); % returns the p-value of a two-sided Wilcoxon rank sum test. ranksum tests the null hypothesis that data in x and y are samples from continuous distributions with equal medians, against the alternative that they are not. The test assumes that the two samples are independent. x and y can have different lengths.
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
    saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\Reviews_02_23\plots\a1_movement_stats\',sprintf('%s.fig',featureNamesPlot{f})))
    saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\Reviews_02_23\plots\a1_movement_stats\',sprintf('%s.png',featureNamesPlot{f})))
    saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\Reviews_02_23\plots\a1_movement_stats\',sprintf('%s.epsc',featureNamesPlot{f})))
    

    % Test for normality  - Shapiro - Wilks
    figure;
    featureNames{f}
    names = {'WT';'KO'};
    mat = [1,0;1,1;0,0;0,1];
    isNormal = [];
    for pp = 1:2
        if pp == 1, yData_group = yWT; else; yData_group = yKO; end
        [h_SW, p_SW, W] = swtest(yData_group, 0.05,1);
        subplot(1,2,pp)
        qqplot(yData_group);
        title({featureNames{f}, names{pp}, sprintf('SW test p: %.4f',p_SW)})
        axis square

%         subplot(2,2,pp+2)
%         boxplot(yData_group);
%         title({featureNames{f}, names{pp}, 'Boxplot'})
%         axis square
        isNormal = [isNormal;[~h_SW,p_SW]];
    end

    stats(f).feature = featureNames{f};
    stats(f).hp_t =[h_t,p_t];
    stats(f).hp_rs =[h_rs,p_rs];
    stats(f).mean = [mean(yWT);mean(yKO)];
    stats(f).std = [std(yWT);std(yKO)];
    stats(f).isNormal = isNormal';
end 

save('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\Reviews_02_23\stats\movement_stats.mat','stats')
