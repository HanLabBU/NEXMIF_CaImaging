%% Description

% This code runs through the list of all correlation values for all
% sessions and plots distribution

% Assymmetric correlation - intersection(A,B)/n(A); % custom (in utils)
% Jaccards Index - intersection(A,B)/union(A,B); % Built in

% Spet 15th 2021 - by Athif Mohamed

%% Initialize

addpath(genpath('U:\eng_research_handata\Athif Mohamed\nexmif_paper\Utils\'))  % Add utilities
init % Initialize data directories and genotypes
close all

%% Begin code

% cd(pathData)
corrTypeId = 3;

% load(fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\stats\correlation_ball','corrStatsTableSpeedAsym'),'corrStatsTable')
load(fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\stats\correlation_ball','corrStatsTableSpeedPearson_08_22'),'corrStatsTable')

fieldNames = corrStatsTable.Properties.VariableNames;
corrStatsTable([6,13,18:19,29,31:35,40:43,45],:) = [];

%% Sessionwise stats
% Get Sessionwise means
corrMat_low = 4; % Low speed correlation all (matrix)
corrThresh_low_pos = 5; % Low speed threshold for positive corr
corrThresh_low_neg = 6; % Low speed threshold for neg corr
corrMat_high = 7;
corrThresh_high_pos = 8;
corrThresh_high_neg = 9;
corrNames = {'Positive significant correlation','Positive random correlation'};

sessionwise_stats = struct;
for typ = 1:2
    yData = [];
    movSp = [];
    gType = [];
    for mov = 1:2
        for rr = 1:size(corrStatsTable,1)
            switch mov
                case 1  % rest
                    corrMat = corrStatsTable{rr,4}{:};
                    thresh = corrStatsTable{rr,5}{:};
                    movSp = [movSp;0];
                case 2  % run
                    corrMat = corrStatsTable{rr,7}{:};
                    thresh = corrStatsTable{rr,8}{:};
                    movSp = [movSp;1];
            end

            n = size(corrMat);
            corrMat(logical(triu(ones(n)))) = nan;
            thresh(logical(triu(ones(n)))) = nan;
            cell_close = double(corrStatsTable{rr,10}{:}<20);
            cell_close(logical(triu(ones(n)))) = nan;

            switch typ
                case 1  % Positive significant correlation
                    yData = [yData;mean(corrMat((corrMat-thresh >1e-5) & (cell_close == 0) & (corrMat >1e-5)))];
                case 2  % Positive random correlations
                    yData = [yData;mean(corrMat((corrMat-thresh<=1e-5) & (cell_close == 0) & (corrMat>1e-5)))];
            end

            gType = [gType;ismember(corrStatsTable{rr,1},miceWT)];
        end
    end

    WTRun = yData(gType == 1 & movSp == 1);
    KORun = yData(gType == 0 & movSp == 1);
    WTRest = yData(gType == 1 & movSp == 0);
    KORest = yData(gType == 0 & movSp == 0);

   figure('Units','Inches','Position',[2 2 6 6]);

    b2 = boxchart(2*ones(size(WTRun')),...
        WTRun',...
        'MarkerStyle','none',...
        'BoxFaceColor',colors(2,:),...
        'BoxFaceAlpha',0,...
        'WhiskerLineColor', colors(2,:),...
        'LineWidth',2);
    hold on
    b1 = boxchart(ones(size(WTRest')),...
        WTRest',...
        'MarkerStyle','none',...
        'BoxFaceColor',colors(1,:),...
        'BoxFaceAlpha',0,...
        'WhiskerLineColor', colors(1,:),...
        'LineWidth',2);
    hold on

    s2 = scatter(2,WTRun,18,colors(2,:),'filled');
    hold on
    s1 = scatter(1,WTRest,18,colors(1,:),'filled');
    hold on
%     line([1,2],[WTRest,WTRun],'color','k')
    b4 = boxchart(4*ones(size(KORun')),...
        KORun',...
        'MarkerStyle','none',...
        'BoxFaceColor',colors(4,:),...
        'BoxFaceAlpha',0,...
        'WhiskerLineColor', colors(4,:),...
        'LineWidth',2);
    hold on
    b3 = boxchart(3*ones(size(KORest')),...
        KORest',...
        'MarkerStyle','none',...
        'BoxFaceColor',colors(3,:),...
        'BoxFaceAlpha',0,...
        'WhiskerLineColor', colors(3,:),...
        'LineWidth',2);
    hold on
    s4 = scatter(4,KORun,18,colors(4,:),'filled');
    hold on
    s3 = scatter(3,KORest,18,colors(3,:),'filled');
    hold on
%     line([3,4],[KORest,KORun],'color','k')

    xlim([0.5,4.5])
    % ylim([0,27])
    xticks([1,2,3,4])
    xticklabels({'Rest','Run','Rest','Run'})
    ylabel('Mean correlation')

%     set(gcf,'Color','none')
    set(gca,'Units','inches','InnerPosition',[.8 .5 1.25 1],'TickDir','out','TickLength',[0.03, 0.025],'Color','none','Box','off','LineWidth',1,'YLim',[0 inf])

%     saveas(gcf,fullfile("U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\Reviews_02_23\plots\e0_correlation\",sprintf("Sessionwise %s.fig",corrNames{typ})))
%     saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\Reviews_02_23\plots\e0_correlation\',sprintf("Sessionwise %s.png",corrNames{typ})))
%     saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\Reviews_02_23\plots\e0_correlation\',sprintf("Sessionwise %s.epsc",corrNames{typ})))
    saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\Reviews_02_23\plots\e0_correlation\',sprintf("Sessionwise %s_small.pdf",corrNames{typ})))

    % Stats
    sessionwise_stats(typ).feature = corrNames{typ};
    sessionwise_stats(typ).means = [mean(WTRest),mean(WTRun),mean(KORest),mean(KORun)]';
    sessionwise_stats(typ).stds = [std(WTRest),std(WTRun),std(KORest),std(KORun)]';

    % Check normality
    figure('Units','Inches','Position',[3 3 5 5]);
    names = {'WTrest';'WTrun';'KOrest';'KOrun'};
    mat = [1,0;1,1;0,0;0,1];
    isNormal = [];
    for pp = 1:4
        yData_group = yData(gType == mat(pp,1) & movSp == mat(pp,2));
        [h_SW, p_SW, W] = swtest(yData_group, 0.05,1);
        subplot(2,2,pp)
        qqplot(yData_group);
        title({names{pp}, sprintf('SW test p: %.4f',p_SW)})
        isNormal = [isNormal;[~h_SW,p_SW]];
        axis square
    end



    sessionwise_stats(typ).isNormal= isNormal';



    % Do two way anova
    [sessionwise_stats(typ).ANOVA.p...
        sessionwise_stats(typ).ANOVA.tbl...
        sessionwise_stats(typ).ANOVA.stats...
        sessionwise_stats(typ).ANOVA.terms] = anovan(yData,{gType,movSp},"model","interaction",'display','off');


    % FIT GLM
    sessionwise_stats(typ).GLM.mdl = fitglm([gType,movSp],yData,' y ~ x1 + x2 + x1 * x2','Distribution','normal');
    deviance_test = devianceTest(sessionwise_stats(typ).GLM.mdl);
    sessionwise_stats(typ).GLM.fit = deviance_test{2,4};
    

    sessionwise_stats(typ).GLM.md2 = fitglm([gType,movSp],yData,' y ~ x1','Distribution','normal');
    deviance_test = devianceTest(sessionwise_stats(typ).GLM.md2);
    sessionwise_stats(typ).GLM.fit2 = deviance_test{2,4};


%    save('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\Reviews_02_23\stats\e0_correlation.mat','sessionwise_stats')



end
