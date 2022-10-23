% Compare stats between WT and KO
% For calcium event features for all mice

%%
addpath(genpath('U:\eng_research_handata\Athif Mohamed\nexmif_paper\Utils\'))  % Add utilities
init % Initialize data directories and genotypes
close all

% load events stats new
load(fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\stats\event_features_ball','eventStats_08_22'),'eventStatsAll')


fieldNames = fieldnames(eventStatsAll);
featureNames = fieldNames(7:end-4);
featureNames(3:end) = {'Rise time','Mean FWHM','Event rate','Event activity rate','Event rate during high movement','Event rate during low movement','Event activity rate during high movement', 'Event activity rate during low movement'};
titleNames = featureNames;
titleNames(3:end) = {{'Rise time','Ball'},{'Mean FWHM','Ball'},{'Event rate','Ball'},...
    {'Event activity rate','Ball'},...
    {'Event rate during high movement','Ball'},{'Event rate during low movement','Ball'},...
    {'Event activity rate during high movement','Ball'}, {'Event activity rate during low movement','Ball'}};
eventStatsTable = struct2table(eventStatsAll);


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

xlabels = {'Delta F/F','Inter Event Interval (s)','Rise Time (s)','Event Width - FWHM (s)','Event rate(events/min)',  'Event activity rate', 'Event rate(events/min)',...
    'Event rate(events/min)','Event activity rate','Event activity rate'};
% Event activity rate = Sum(Ca == 1)/High movement duration*60;

%% For each feature plots violin plot
stats = struct();
sessionwise_stats = struct();
ss = 1;
v_event_rate_struct = struct();
for f = [3,4,8,7,5] % 3:numel(featureNames)
    for c = 1%:3 % ball
        
        maskG1 = ismember([eventStatsTable.animal],miceStudy(maskWT)); %WT
        maskG2 = ismember([eventStatsTable.animal],miceStudy(maskKO)); %KO
        maskC = ismember([eventStatsTable.condition],'Ball'); %Ball
        masksAll = [maskG1,maskG2,maskC];
        
        featureValsAll = eventStatsTable{maskC,f+6};
        featureMaxAll = max(featureValsAll);
        
        legendNames = [];
        y_gen_data = {};
        for subG = 1:2
            switch subG
                case 1
                    nMice = sum(maskWT);
                    miceName = miceStudy(maskWT);
                case 2
                    nMice = sum(maskKO);
                    miceName = miceStudy(maskKO);
            end
            maskG = masksAll(:,subG);
            y_gen = [];
            mList = {};
            dList = {};
            for m = 1:nMice
                maskM = ismember([eventStatsTable.animal], miceName(m));
                %                 maskY =  maskM & maskG;
                if f == 7 || f == 8 || f == 9 || f == 10
                    maskY =  maskM & maskG & ~maskFewM & ~maskBadS;
                else
                    maskY =  maskM & maskG;
                end
                y = eventStatsTable{maskY,f+6};
                
                y(isnan(y)) = [];
                y_gen = [y_gen;y];
                
                maskY(isnan(y)) = 0;
                mAll = eventStatsTable.animal(maskY);
                dAll = eventStatsTable.day(maskY);
                if ~isempty(y)
                    dList = [dList,numel(unique(dAll))];
                    mList = [mList,mAll{1}];
                end
            end
            %             if any(isnan(y_gen))
            %                 xx = 1;
            %                 break
            %             end
            y_gen_data = [y_gen_data,y_gen];
            %% get stats
            
            stats(ss).data(subG).type = mType{subG};
            stats(ss).data(subG).nMice = numel(mList);
            stats(ss).data(subG).nCells = numel(y_gen);
            stats(ss).data(subG).mList = mList;
            stats(ss).data(subG).dList = dList;
            stats(ss).data(subG).mean = mean(y_gen);
            stats(ss).data(subG).std = std(y_gen);
        end
        
        stats(ss).feature = featureNames{f};
        [stats(ss).ttest.h, stats(ss).ttest.p] = ttest2(y_gen_data{1},y_gen_data{2});
        
        %% Plot violin plots
        
        figure;
        
        t = struct();
        t.WT = y_gen_data{1};
        t.KO = y_gen_data{2};
        yylim = [0.9*min([y_gen_data{1};y_gen_data{2}]),1.1*max([y_gen_data{1};y_gen_data{2}])];
        vp = violinplot(t);
        
        if f == 8
            vp(1,1).ViolinColor = colors(1,:);
            vp(1,2).ViolinColor = colors(3,:);
            vp(1,1).BoxColor = colors(1,:);
            vp(1,2).BoxColor = colors(3,:);
            vp(1,1).EdgeColor = colors(1,:);
            vp(1,2).EdgeColor = colors(3,:);
        else
            vp(1,1).ViolinColor = colors(2,:);
            vp(1,2).ViolinColor = colors(4,:);
            vp(1,1).BoxColor = colors(2,:);
            vp(1,2).BoxColor = colors(4,:);
            vp(1,1).EdgeColor = colors(2,:);
            vp(1,2).EdgeColor = colors(4,:);
        end
        vp(1,1).ViolinAlpha = 0.1;
        vp(1,2).ViolinAlpha = 0.1;
        
        
        title(titleNames{f},'FontSize',11,'FontName','Arial','FontWeight','Bold')
        ylabel(xlabels{f})
        set(gcf,'Color','none')
        set(gca,'Units','inches','InnerPosition',[.5 .3 2 2],'TickDir','out','TickLength',[0.03, 0.025],'Color','none','Box','off','LineWidth',1)
        
%         saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\c_event_features_ball',sprintf('violinplot %s.png',featureNames{f})))
%         saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\c_event_features_ball',sprintf('violinplot %s.fig',featureNames{f})))
%         saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\c_event_features_ball',sprintf('violinplot %s.epsc',featureNames{f})))
%         saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\c_event_features_ball',sprintf('violinplot %s.svg',featureNames{f})))
        
        %% collect for event rate violin plots
        if f == 8
            
            v_event_rate_struct.WTRest = y_gen_data{1};
            v_event_rate_struct.KORest = y_gen_data{2};
        elseif f == 7
            
            v_event_rate_struct.WTRun = y_gen_data{1};
            v_event_rate_struct.KORun = y_gen_data{2};
        end
        %% Sessionwise T test
        
        if f == 7 || f == 8 || f == 9 || f == 10 || f == 5
            [GWT,MWT,DWT] = findgroups(eventStatsTable{maskC&maskG1&~maskFewM & ~maskBadS,2},eventStatsTable{maskC&maskG1&~maskFewM & ~maskBadS,3});
            [GKO,MKO,DKO] = findgroups(eventStatsTable{maskC&maskG2&~maskFewM & ~maskBadS,2},eventStatsTable{maskC&maskG2&~maskFewM & ~maskBadS,3});
            meanWT = splitapply(@nanmean,eventStatsTable{maskC&maskG1&~maskFewM & ~maskBadS,f+6},GWT);
            meanKO = splitapply(@nanmean,eventStatsTable{maskC&maskG2&~maskFewM & ~maskBadS,f+6},GKO);
            mWT = splitapply(@unique,eventStatsTable{maskC&maskG1&~maskFewM & ~maskBadS,2},GWT);
            mKO = splitapply(@unique,eventStatsTable{maskC&maskG2&~maskFewM & ~maskBadS,2},GKO);
            dWT = splitapply(@unique,eventStatsTable{maskC&maskG1&~maskFewM & ~maskBadS,3},GWT);
            dKO = splitapply(@unique,eventStatsTable{maskC&maskG2&~maskFewM & ~maskBadS,3},GKO);
        else
            [GWT,MWT,DWT] = findgroups(eventStatsTable{maskC&maskG1,2} , eventStatsTable{maskC&maskG1,3});
            [GKO,MKO,DKO] = findgroups(eventStatsTable{maskC&maskG2,2} , eventStatsTable{maskC&maskG2,3});
            meanWT = splitapply(@nanmean,eventStatsTable{maskC&maskG1,f+6},GWT);
            meanKO = splitapply(@nanmean,eventStatsTable{maskC&maskG2,f+6},GKO);
            mWT = splitapply(@unique,eventStatsTable{maskC&maskG1,2},GWT);
            mKO = splitapply(@unique,eventStatsTable{maskC&maskG2,2},GKO);
            dWT = splitapply(@unique,eventStatsTable{maskC&maskG1,3},GWT);
            dKO = splitapply(@unique,eventStatsTable{maskC&maskG2,3},GKO);
        end
        
        figure
        b1 = boxchart(ones(size(meanWT)),...
            meanWT,...
            'MarkerStyle','none',...
            'BoxFaceColor',colors(2,:),...
            'BoxFaceAlpha',0,...
            'WhiskerLineColor', colors(2,:),...
            'LineWidth',1);
        hold on
        s1 = scatter(1,meanWT,[],colors(2,:),'filled');
        hold on
        b2 = boxchart(2*ones(size(meanKO)),...
            meanKO,...
            'MarkerStyle','none',...
            'BoxFaceColor',colors(4,:),...
            'BoxFaceAlpha',0,...
            'WhiskerLineColor', colors(4,:),...
            'LineWidth',1);
        hold on
        s2 = scatter(2,meanKO,[],colors(4,:),'filled');
        
        
        xlim([0.5,2.5])
        % ylim([0,27])
        xticks([1,2])
        xticklabels({'WT','KO'})
        
        title(titleNames{f},'FontSize',11,'FontName','Arial','FontWeight','Bold')
        ylabel(xlabels{f})
        set(gca,'TickDir', 'out','TickLength',[0.03, 0.025], 'Color','none','LineWidth',1,'YLim',[0 inf])
        box off
        set(gcf,'Units','inches','InnerPosition',[1, 3, 2, 2])

       
        saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\c_event_features_ball\8.3.22',sprintf('Sessionwise boxplot %s.png',featureNames{f})))
        saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\c_event_features_ball\8.3.22',sprintf('Sessionwise boxplot %s.fig',featureNames{f})))
        saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\c_event_features_ball\8.3.22',sprintf('Sessionwise boxplot %s.epsc',featureNames{f})))
         
        [h,p] = ttest2(meanWT,meanKO);
        
        sessionwise_stats(ss).feature = featureNames{f};
        sessionwise_stats(ss).ttest_h = h;
        sessionwise_stats(ss).ttest_p = p;
        sessionwise_stats(ss).WTmean = mean(meanWT,'omitnan');
        sessionwise_stats(ss).KOmean = mean(meanKO,'omitnan');
        sessionwise_stats(ss).WTstd = std(meanWT,'omitnan');
        sessionwise_stats(ss).KOstd = std(meanKO,'omitnan');
        sessionwise_stats(ss).sessWT = {mWT,dWT};
        sessionwise_stats(ss).sessKO = {mKO,dKO};
        sessionwise_stats(ss).dataKO = meanKO;
        sessionwise_stats(ss).dataWT = meanWT;
        ss =ss+1;
    end
end


%% Plot event rate violin plots
v_event_rate_plot.WTRest = v_event_rate_struct.WTRest;
v_event_rate_plot.WTRun = v_event_rate_struct.WTRun;
v_event_rate_plot.KORest = v_event_rate_struct.KORest;
v_event_rate_plot.KORun = v_event_rate_struct.KORun;


meanWTRest = mean(v_event_rate_struct.WTRest);
meanWTRun = mean(v_event_rate_struct.WTRun);
meanKORest = mean(v_event_rate_struct.KORest);
meanKORun = mean(v_event_rate_struct.KORun);

stdWTRest = std(v_event_rate_struct.WTRest);
stdWTRun = std(v_event_rate_struct.WTRun);
stdKORest = std(v_event_rate_struct.KORest);
stdKORun = std(v_event_rate_struct.KORun);


figure
vp = violinplot(v_event_rate_plot);

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

title('Event rate','FontSize',11,'FontName','Arial','FontWeight','Bold')
ylabel('Event rate (ev/min)');

% set(gcf,'Color','none')
set(gca,'Units','inches','InnerPosition',[.5 .3 4 2],'TickDir','out','TickLength',[0.03, 0.025],'Color','none','Box','off','LineWidth',1)


[h_WT,p_WT] = ttest(v_event_rate_plot.WTRest,v_event_rate_plot.WTRun);
[h_KO,p_KO] = ttest(v_event_rate_plot.KORest,v_event_rate_plot.KORun);
[h_Rest,p_Rest] = ttest2(v_event_rate_plot.WTRest,v_event_rate_plot.KORest);
[h_Run,p_Run] = ttest2(v_event_rate_plot.WTRun,v_event_rate_plot.KORun);

% yy  = ylim;
% text(1.5,0.7*yy(2),sprintf('p* = %.4f',p_WT*4),'HorizontalAlignment','Center')
% text(3.5,0.7*yy(2),sprintf('p* = %.4f',p_KO*4),'HorizontalAlignment','Center')
% text(2,0.95*yy(2),sprintf('p* = %.4f',p_Rest*4),'HorizontalAlignment','Center')
% text(3,0.95*yy(2),sprintf('p* = %.4f',p_Run*4),'HorizontalAlignment','Center')

% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\c_event_features_ball',sprintf('Violin plot %s.png','event_Rate_all')))
% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\c_event_features_ball',sprintf('Violin plot %s.fig','event_Rate_all')))
% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\c_event_features_ball',sprintf('Violin plot %s.svg','event_Rate_all')))
% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\c_event_features_ball',sprintf('Violin plot %s.epsc','event_Rate_all')))


%% Calculate number of events total 

% load events stats new
load(fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\stats\event_features_ball','eventStats_List'),'eventStatsList')
eventStatsTable = struct2table(eventStatsList);
maskG1 = ismember([eventStatsTable.animal],miceStudy(maskWT)); %WT
maskG2 = ismember([eventStatsTable.animal],miceStudy(maskKO)); %KO
[GWT,MWT,DWT] = findgroups(eventStatsTable{maskG1,2} , eventStatsTable{maskG1,3});
[GKO,MKO,DKO] = findgroups(eventStatsTable{maskG2,2} , eventStatsTable{maskG2,3});
           
nEvWT = sum(splitapply(@sum,eventStatsTable{maskG1,8},GWT))
nEVKO = sum(splitapply(@sum,eventStatsTable{maskG2,8},GKO))


%% Stats sessionwise boxplots 
box_event_rate.WTRest = sessionwise_stats(3).dataWT;
box_event_rate.WTRun = sessionwise_stats(4).dataWT;
box_event_rate.KORest = sessionwise_stats(3).dataKO;
box_event_rate.KORun = sessionwise_stats(4).dataKO;



[h_WT,p_WT] = ttest(box_event_rate.WTRest,box_event_rate.WTRun);
[h_KO,p_KO] = ttest(box_event_rate.KORest,box_event_rate.KORun);
[h_Rest,p_Rest] = ttest2(box_event_rate.WTRest,box_event_rate.KORest);
[h_Run,p_Run] = ttest2(box_event_rate.WTRun,box_event_rate.KORun);

% yy  = ylim;
% text(1.5,0.7*yy(2),sprintf('p* = %.4f',p_WT*4),'HorizontalAlignment','Center')
% text(3.5,0.7*yy(2),sprintf('p* = %.4f',p_KO*4),'HorizontalAlignment','Center')
% text(2,0.95*yy(2),sprintf('p* = %.4f',p_Rest*4),'HorizontalAlignment','Center')
% text(3,0.95*yy(2),sprintf('p* = %.4f',p_Run*4),'HorizontalAlignment','Center')

% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\c_event_features_ball',sprintf('Violin plot %s.png','event_Rate_all')))
% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\c_event_features_ball',sprintf('Violin plot %s.fig','event_Rate_all')))
% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\c_event_features_ball',sprintf('Violin plot %s.svg','event_Rate_all')))
% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\c_event_features_ball',sprintf('Violin plot %s.epsc','event_Rate_all')))