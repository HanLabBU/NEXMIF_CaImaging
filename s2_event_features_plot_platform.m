% Compare stats between WT and KO
% For calcium event features for all mice

%%
addpath(genpath('U:\eng_research_handata\Athif Mohamed\nexmif_paper\Utils\'))  % Add utilities
init % Initialize data directories and genotypes
close all

% load events stats new
load(fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\stats\event_features_platform','eventStats'),'eventStatsAll')
eventStatsTable = struct2table(eventStatsAll);
fieldNames = fieldnames(eventStatsAll);
featureNames = fieldNames(7:end);
featureNames(3:end) = {'Rise time','Event rate','Mean FWHM','Event activity rate'};
titleNames = featureNames;
titleNames(3:end) = {{'Rise time','Platform'},{'Event rate','Platform'},{'Mean FWHM','Platform'},{'Event activity rate','Platform'}};



maskFewM = false(size(eventStatsTable,1),1); % A flag for sessions with less than 10% mov
for ii = 1:numel(fewMovSessionsM)
    maskFewM(ismember([eventStatsTable.animal],fewMovSessionsM(ii))& ismember([eventStatsTable.day],fewMovSessionsD(ii)))= 1;
end
xlabels = {'dF/F','Inter-event interval (s)','Rise time (s)','FWHM (s)','Event rate (events/min)',  'Event activity rate', 'Event rate (events/min)',...
    'Event rate (events/min)','Event activity rate','Event activity rate'};
% Event activity rate = Sum(Ca == 1)/High movement duration*60;

%% For each feature plots violin plot

for f = [3,4,5,6] % 3:numel(featureNames)
    for c = 1%:3 % ball

        maskG1 = ismember([eventStatsTable.animal],miceStudy(maskWT)); %WT
        maskG2 = ismember([eventStatsTable.animal],miceStudy(maskKO)); %KO
        maskC = ismember([eventStatsTable.condition],'Platform'); %Ball
        masksAll = [maskG1,maskG2,maskC];
        
        featureValsAll = eventStatsTable{maskC,f+6};
        featureMaxAll = max(featureValsAll);
        
        ss = 1;
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
            for m = 1:nMice
                maskM = ismember([eventStatsTable.animal], miceName(m));
             
                    maskY =  maskM & maskG;
                y = eventStatsTable{maskY,f+6};
                y(isnan(y)) = [];
                y_gen = [y_gen;y];
                
            end
            y_gen_data = [y_gen_data,y_gen];
        end
 
        %% Plot violin plots
        
        plot5 = [y_gen_data{1};y_gen_data{2}];
        label5 = [repmat({'WT'},size(y_gen_data{1}));repmat({'KO'},size(y_gen_data{2}))];
        
        figure;
        subplot(1,2,1)

        yylim = [0.9*min([y_gen_data{1};y_gen_data{2}]),1.1*max([y_gen_data{1};y_gen_data{2}])];
        vsWT = violinplot(y_gen_data{1}, repmat({'WT'},size(y_gen_data{1})),'ViolinColor',colors(2,:),'ViolinAlpha',0.1,'BoxColor',colors(2,:),'EdgeColor',colors(2,:));
        xlim([0.5,1.5])
        ylim([0 yylim(2)]);
        ylabel(xlabels{f})
        set(gca,'TickDir','out','Color','none','Box','off','LineWidth',2)
        
        subplot(1,2,2)
        vsKO = violinplot(y_gen_data{2}, repmat({'KO'},size(y_gen_data{2})),'ViolinColor',colors(4,:),'ViolinAlpha',0.1,'BoxColor',colors(4,:),'EdgeColor',colors(4,:));
        xlim([0.5,1.5])
        ylim([0 yylim(2)]);
        sgtitle(titleNames{f},'FontSize',11,'FontName','Helvetica','FontWeight','Bold') 
        ylabel(xlabels{f})
        set(gca,'TickDir','out','Color','none','Box','off','LineWidth',2)
        
%         set(gcf,'Color','none')
        
        
        saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\supplementary\s2',sprintf('violinplot %s.png',featureNames{f})))
        saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\supplementary\s2',sprintf('violinplot %s.fig',featureNames{f})))
        saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\supplementary\s2',sprintf('violinplot %s.svg',featureNames{f})))

    end
end

