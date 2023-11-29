% Compare stats between WT and KO
% For calcium event features for all mice

%%
addpath(genpath('U:\eng_research_handata\Athif Mohamed\nexmif_paper\Utils\'))  % Add utilities
init % Initialize data directories and genotypes
close all

% load events stats new
load(fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\stats\event_features_ball','eventStats_02_23'),'eventStatsAll')
eventStatsTable = struct2table(eventStatsAll);

fieldNames = fieldnames(eventStatsAll);
featureCols = [9,17,18,10,19,20,11,13,14];
featureNames =  {'Rise time','Rise time - Run','Rise time - Rest','Mean FWHM','Mean FWHM - Run','Mean FWHM - Rest','Event rate','Event rate - Run','Event rate - Rest'};
% titleNames = featureNames;
% titleNames(3:end) = {{'Rise time','Ball'},{'Mean FWHM','Ball'},{'Event rate','Ball'},...
%     {'Event activity rate','Ball'},...
%     {'Event rate during high movement','Ball'},{'Event rate during low movement','Ball'},...
%     {'Event activity rate during high movement','Ball'}, {'Event activity rate during low movement','Ball'}};


maskFewM = false(size(eventStatsTable,1),1); % A flag for sessions with less than 10% mov
for ii = 1:numel(fewMovSessionsM)
    maskFewM(ismember([eventStatsTable.animal],fewMovSessionsM(ii))& ismember([eventStatsTable.day],fewMovSessionsD(ii)))= 1;
end

maskBadS = false(size(eventStatsTable,1),1); % A flag for sessions with bad movement trace
for ii = 1:numel(badSpeedSessions)
    ssp = strsplit(badSpeedSessions{ii},'_');
    ssd = str2double(ssp{2}(end));
    maskBadS(ismember([eventStatsTable.animal],ssp(1))& ismember([eventStatsTable.day],ssd))= 1;
end

% xlabels = {'Delta F/F','Inter Event Interval (s)','Rise Time (s)','Event Width - FWHM (s)','Event rate(events/min)',  'Event activity rate', 'Event rate(events/min)',...
%     'Event rate(events/min)','Event activity rate','Event activity rate'};
% Event activity rate = Sum(Ca == 1)/High movement duration*60;

maskG1 = ismember([eventStatsTable.animal],miceStudy(maskWT)); %WT
maskG2 = ismember([eventStatsTable.animal],miceStudy(maskKO)); %KO
maskC = ismember([eventStatsTable.condition],'Ball'); %Ball

%% For each feature - Sessionwise means stats
sessionwise_stats = struct();
ss = 1;
for ii = 1:numel(featureCols)
    f = featureCols(ii);
    if f == 17 || f == 18 || f == 19 || f == 20 || f == 13 || f == 14
        mask1 = maskG1&~maskFewM & ~maskBadS;
        mask2 = maskG2&~maskFewM & ~maskBadS;
        mask3 = ~maskFewM & ~maskBadS;
        [gWT,mWT,dWT] = findgroups(eventStatsTable{mask1,2},eventStatsTable{mask1,3});
        [gKO,mKO,dKO] = findgroups(eventStatsTable{mask2,2},eventStatsTable{mask2,3});
        meanWT = splitapply(@nanmean,eventStatsTable{mask1,f},gWT);
        meanKO = splitapply(@nanmean,eventStatsTable{mask2,f},gKO);
        meanBoth = [];
    else
        [gWT,mWT,dWT] = findgroups(eventStatsTable{maskG1,2} , eventStatsTable{maskG1,3});
        [gKO,mKO,dKO] = findgroups(eventStatsTable{maskG2,2} , eventStatsTable{maskG2,3});
        meanWT = splitapply(@nanmean,eventStatsTable{maskG1,f},gWT);
        meanKO = splitapply(@nanmean,eventStatsTable{maskG2,f},gKO);
        if f == 11
            [gBoth,mBoth,dBoth] = findgroups(eventStatsTable{:,2} , eventStatsTable{:,3});
            meanBoth = splitapply(@nanmean,eventStatsTable{:,f},gBoth);
        else
            meanBoth = [];
        end
    end
    
    [h,p] = ttest2(meanWT,meanKO);
    [p_rs,h_rs] = ranksum(meanWT,meanKO);
    sessionwise_stats(ss).feature = featureNames{ii};
    sessionwise_stats(ss).ttest_h = h;
    sessionwise_stats(ss).ttest_p = p;
    sessionwise_stats(ss).rs_h = h_rs;
    sessionwise_stats(ss).rs_p = p_rs;
    sessionwise_stats(ss).WTmean = mean(meanWT,'omitnan');
    sessionwise_stats(ss).KOmean = mean(meanKO,'omitnan');
    sessionwise_stats(ss).WTstd = std(meanWT,'omitnan');
    sessionwise_stats(ss).KOstd = std(meanKO,'omitnan');
    sessionwise_stats(ss).sessWT = {mWT,dWT};
    sessionwise_stats(ss).sessKO = {mKO,dKO};
    sessionwise_stats(ss).dataKO = meanKO;
    sessionwise_stats(ss).dataWT = meanWT;
    sessionwise_stats(ss).meanBoth = meanBoth;
    ss =ss+1;
end

% Plot figures

for ii =  [2,5,8]
    WTrest = sessionwise_stats(ii+1).dataWT;
    WTrun = sessionwise_stats(ii).dataWT;
    KOrest = sessionwise_stats(ii+1).dataKO;
    KOrun = sessionwise_stats(ii).dataKO;
    titleName = sessionwise_stats(ii-1).feature;
    figure
    b1 = boxchart(ones(size(WTrest)),...
        WTrest,...
        'MarkerStyle','none',...
        'BoxFaceColor',colors(1,:),...
        'BoxFaceAlpha',0,...
        'WhiskerLineColor', colors(1,:),...
        'LineWidth',1);
    hold on
    s1 = scatter(1,WTrest,10,colors(1,:),'filled');
    hold on
    b2 = boxchart(2*ones(size(WTrun)),...
        WTrun,...
        'MarkerStyle','none',...
        'BoxFaceColor',colors(2,:),...
        'BoxFaceAlpha',0,...
        'WhiskerLineColor', colors(2,:),...
        'LineWidth',1);
    hold on
    s2 = scatter(2,WTrun,10,colors(2,:),'filled');
    b3 = boxchart(3*ones(size(KOrest)),...
        KOrest,...
        'MarkerStyle','none',...
        'BoxFaceColor',colors(3,:),...
        'BoxFaceAlpha',0,...
        'WhiskerLineColor', colors(3,:),...
        'LineWidth',1);
    hold on
    s3 = scatter(3,KOrest,10,colors(3,:),'filled');
    hold on
    b4 = boxchart(4*ones(size(KOrun)),...
        KOrun,...
        'MarkerStyle','none',...
        'BoxFaceColor',colors(4,:),...
        'BoxFaceAlpha',0,...
        'WhiskerLineColor', colors(4,:),...
        'LineWidth',1);
    hold on
    s4= scatter(4,KOrun,10,colors(4,:),'filled');

    xlim([0.5,4.5])
    % ylim([0,27])
    xticks([1,2,3,4])
    %     xticklabels({'WTrest','WTrun','KOrest','KOrun'})

    %     title(titleName,'FontSize',11,'FontName','Arial','FontWeight','Bold')
    set(gca,'TickDir', 'out','TickLength',[0.03, 0.025], 'Color','none','LineWidth',1,'YLim',[0 inf])
    box off
    set(gcf,'Units','inches','InnerPosition',[1, 3, 1.5, 1])

    % Save plot
    saveas(gcf,fullfile("U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\Reviews_08_23\Codes\plots\pdfs",sprintf('Fig2 Sessionwise boxplot %s.pdf',titleName)))

    %     saveas(gcf,fullfile("U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\Reviews_02_23\plots\c_event_features_ball",sprintf('Sessionwise boxplot %s.png',titleName)))
    %     saveas(gcf,fullfile("U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\Reviews_02_23\plots\c_event_features_ball",sprintf('Sessionwise boxplot %s.fig',titleName)))
    %     saveas(gcf,fullfile("U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\Reviews_02_23\plots\c_event_features_ball",sprintf('Sessionwise boxplot %s.epsc',titleName)))

    % Do  ANOVA sessionwise
    yData = [WTrest;WTrun;KOrest;KOrun];
    gType = [ones(size(WTrest));ones(size(WTrun));zeros(size(KOrest));zeros(size(KOrun))]; % WT - 1, KO - 0;
    movSp = [zeros(size(WTrest));ones(size(WTrun));zeros(size(KOrest));ones(size(KOrun))]; % Running - 1, Resting - 0;

    % Test for normality  - Shapiro - Wilks
    figure;
    titleName
    names = {'WTrest';'WTrun';'KOrest';'KOrun'};
    mat = [1,0;1,1;0,0;0,1];
    isNormal = [];
    for pp = 1:4
        yData_group = yData(gType == mat(pp,1) & movSp == mat(pp,2));
        [h_SW, p_SW, W] = swtest(yData_group, 0.05,1);
        subplot(2,2,pp)
        qqplot(yData_group);
        title({titleName, names{pp}, sprintf('SW test p: %.4f',p_SW)})
        isNormal = [isNormal;[~h_SW,p_SW]];
    end
    sessionwise_stats(ii).isNormal = isNormal';

    %     saveas(gcf,fullfile("U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\Reviews_02_23\plots\c_event_features_ball",sprintf('Sessionwise QQ plot %s.png',titleName)))


    % Do two way anova
    [sessionwise_stats(ii).ANOVA.p...
        sessionwise_stats(ii).ANOVA.tbl...
        sessionwise_stats(ii).ANOVA.stats...
        sessionwise_stats(ii).ANOVA.terms] = anovan(yData,{gType,movSp},"model","interaction",'display','off');


    % FIT GLM
    sessionwise_stats(ii).GLM = fitglm([gType,movSp],yData,' y ~ x1 + x2 + x1 * x2','Distribution','normal');
    deviance_test = devianceTest(sessionwise_stats(ii).GLM);
    sessionwise_stats(ii).GLM_fit = deviance_test{2,4};

    %     [h,p_AD,adstat,cv] = adtest(yData);
    %     tbl = table(yData,gType,movSp);
    %     tbl.movSp = nominal(tbl.movSp);
    %     tbl.gType = nominal(tbl.gType);
    %     md2 = fitlme(tbl,' yData ~ movSp + (1|gType)');
end

%  save('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\Reviews_02_23\stats\c_event_features_ball.mat','sessionwise_stats')


%% For each feature - Cellwise means

% Get stats
cellwise_stats = struct();
ss = 1;
for ii = 1:numel(featureCols)
    f = featureCols(ii);
    if f == 17 || f == 18 || f == 19 || f == 20 || f == 13 || f == 14
        mask1 = maskG1&~maskFewM & ~maskBadS;
        mask2 = maskG2&~maskFewM & ~maskBadS;
        dataWT = eventStatsTable{mask1,f};
        dataKO = eventStatsTable{mask2,f};
    else
        dataWT = eventStatsTable{maskG1,f};
        dataKO = eventStatsTable{maskG2,f};
    end

    [h,p] = ttest2(dataWT,dataKO);

    cellwise_stats(ss).feature = featureNames{ii};
    cellwise_stats(ss).ttest_h = h;
    cellwise_stats(ss).ttest_p = p;
    cellwise_stats(ss).WTmean = mean(dataWT,'omitnan');
    cellwise_stats(ss).KOmean = mean(dataKO,'omitnan');
    cellwise_stats(ss).WTstd = std(dataWT,'omitnan');
    cellwise_stats(ss).KOstd = std(dataKO,'omitnan');
    cellwise_stats(ss).dataKO = dataKO;
    cellwise_stats(ss).dataWT = dataWT;
    ss =ss+1;
end


% Plot violin plot
for ii = [2,5,8]

    WTrest = cellwise_stats(ii+1).dataWT;
    WTrun = cellwise_stats(ii).dataWT;
    KOrest = cellwise_stats(ii+1).dataKO;
    KOrun = cellwise_stats(ii).dataKO;

    vplot.WTrest = WTrest;
    vplot.WTrun = WTrun;
    vplot.KOrest = KOrest;
    vplot.KOrun = KOrun;
    titleName = cellwise_stats(ii-1).feature;

    figure
    vp = violinplot(vplot);

    vp(1,1).ViolinColor = colors(1,:);
    vp(1,2).ViolinColor = colors(2,:);
    vp(1,3).ViolinColor = colors(3,:);
    vp(1,4).ViolinColor = colors(4,:);
    vp(1,1).BoxColor = colors(1,:);
    vp(1,2).BoxColor = colors(2,:);
    vp(1,3).BoxColor = colors(3,:);
    vp(1,4).BoxColor = colors(4,:);
    vp(1,1).EdgeColor = colors(1,:);
    vp(1,2).EdgeColor = colors(2,:);
    vp(1,3).EdgeColor = colors(3,:);
    vp(1,4).EdgeColor = colors(4,:);

    vp(1,1).ViolinAlpha = 0.1;
    vp(1,2).ViolinAlpha = 0.1;
    vp(1,3).ViolinAlpha = 0.1;
    vp(1,4).ViolinAlpha = 0.1;

    vp(1,1).ShowMean = 1;
    vp(1,1).MeanPlot.Color  = [1 1 1 0.5];
    vp(1,2).ShowMean = 1;
    vp(1,2).MeanPlot.Color  = [1 1 1 0.5];
    vp(1,3).ShowMean = 1;
    vp(1,3).MeanPlot.Color  = [1 1 1 0.5];
    vp(1,4).ShowMean = 1;
    vp(1,4).MeanPlot.Color  = [1 1 1 0.5];

    title(titleName,'FontSize',11,'FontName','Arial','FontWeight','Bold')

    % set(gcf,'Color','none')
    % set(gca,'Units','inches','InnerPosition',[.5 .3 4 2],'TickDir','out','TickLength',[0.03, 0.025],'Color','none','Box','off','LineWidth',1)

    set(gca,'TickDir', 'out','TickLength',[0.03, 0.025], 'Color','none','LineWidth',1,'YLim',[0 inf])
    box off
    set(gcf,'Units','inches','InnerPosition',[1, 3, 4, 2])

    % Save plot
    saveas(gcf,fullfile("U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\Reviews_02_23\plots\c_event_features_ball",sprintf('Cellwise boxplot %s.png',titleName)))
    saveas(gcf,fullfile("U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\Reviews_02_23\plots\c_event_features_ball",sprintf('Cellwise boxplot %s.fig',titleName)))
    saveas(gcf,fullfile("U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\Reviews_02_23\plots\c_event_features_ball",sprintf('Cellwise boxplot %s.epsc',titleName)))
    saveas(gcf,fullfile("U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\Reviews_02_23\plots\c_event_features_ball",sprintf('Cellwise boxplot %s.svg',titleName)))


    % Do  ANOVA cellwise
    yData = [WTrest;WTrun;KOrest;KOrun];
    gType = [ones(size(WTrest));ones(size(WTrun));zeros(size(KOrest));zeros(size(KOrun))]; % WT - 1, KO - 0;
    movSp = [zeros(size(WTrest));ones(size(WTrun));zeros(size(KOrest));ones(size(KOrun))]; % Running - 1, Resting - 0;

    % Test for normality  - Shapiro - Wilks
    figure;
    titleName
    names = {'WTrest';'WTrun';'KOrest';'KOrun'};
    mat = [1,0;1,1;0,0;0,1];
    isNormal = [];
    for pp = 1:4
        yData_group = yData(gType == mat(pp,1) & movSp == mat(pp,2));
        [h_SW, p_SW, W] = swtest(yData_group, 0.05,1);
        subplot(2,2,pp)
        qqplot(yData_group);
        title({titleName, names{pp}, sprintf('SW test p: %.4f',p_SW)})
        isNormal = [isNormal;[~h_SW,p_SW]];
    end
    cellwise_stats(ii).isNormal = isNormal';

    %     saveas(gcf,fullfile("U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\Reviews_02_23\plots\c_event_features_ball",sprintf('Cellwise QQ plot %s.png',titleName)))


    % Do two way anova
    [cellwise_stats(ii).ANOVA.p...
        cellwise_stats(ii).ANOVA.tbl...
        cellwise_stats(ii).ANOVA.stats...
        cellwise_stats(ii).ANOVA.terms] = anovan(yData,{gType,movSp},"model","interaction",'display','off');


    % FIT GLM
    cellwise_stats(ii).GLM = fitglm([gType,movSp],yData,' y ~ x1 + x2 + x1 * x2','Distribution','normal');
    deviance_test = devianceTest(cellwise_stats(ii).GLM);
    cellwise_stats(ii).GLM_fit = deviance_test{2,4};

end

save('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\Reviews_02_23\stats\c_event_features_ball_cellwise.mat','cellwise_stats')

%% Average cellwise mean and sd of event rates of WT and KO

meanWT = mean([WTrest;WTrun])
sdWT = std([WTrest;WTrun])
meanKO = mean([KOrest;KOrun])
sdKO = std([KOrest;KOrun])

%% Get cells in each session 


ss = 1;
f1 = 5;
f2 = 6;

[g,m,d] = findgroups(eventStatsTable{:,2},eventStatsTable{:,3});

allcells = splitapply(@unique,eventStatsTable{:,f1},g);
goodcells = splitapply(@unique,eventStatsTable{:,f2},g);

ncellsTable = table(m,d,allcells,goodcells,'VariableNames',{'Mouse','Day','n_cells', 'n_good_cells'})
