% Plots of WT vs KO for movement stats

addpath(genpath('U:\eng_research_handata\Athif Mohamed\nexmif_paper\Utils\'))  % Add utilities
init % Initiate variables
close all

% Load data movement
load(fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\stats\movement_stats','movStats_02_23'),'movStats');
movStatsTable = struct2table(movStats);

% load data age
load(fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\stats\movement_stats','ageStats_09_23'),'ageStatsTable');


% add age-d
for row = 1:size(movStatsTable,1)
    animal_age = cell2mat(cellfun(@str2double,cellstr([ageStatsTable.animal]),'UniformOutput',false));
    [~,idx] = ismember(str2double(movStatsTable{row,"animal"}),animal_age);
    if numel(idx) == 1
        movStatsTable{row,"age_im_d"} = ageStatsTable{idx,"age_im_m"}*28+ movStatsTable{row,"day"}-1;
        movStatsTable{row,"age_virus_d"} = ageStatsTable{idx,"age_virus_w"}*7+ movStatsTable{row,"day"}-1;
    else
        movStatsTable{row,"age_im_d"} = nan;
        movStatsTable{row,"age_virus_d"} = nan;
    end
end


fieldNames = movStatsTable.Properties.VariableNames;
% featureNames = fieldNames(4:end);
featureNames = fieldNames;
featureNames([1:5,9]+3) = {'Average speed (cm/s)', 'Movement bout duration(s)','Number of movement bouts','Number of movement onsets','Movement bout speed(cm/s)', 'Distance (m)'};
featureNamesPlot = featureNames;
featureNamesPlot([1:5,9]+3) = {'Average speed', 'Average movement bout duration','Number of movement bouts','Number of movement onsets','Average speed during movement bout','Total distance traveled'};


marker_list = {'o','+','*','x','square','v','diamond','pentagram'};
WTTable =   movStatsTable(movStatsTable.gType == 1,:);
KOTable =   movStatsTable(movStatsTable.gType == 0,:);

WTUnique = cellfun(@str2double,unique(WTTable.animal));
KOUnique = cellfun(@str2double,unique(KOTable.animal));
figPath = 'U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\Reviews_08_23\Codes\plots\age';
%%
stats = struct();a = 1;
for ff = [1:5,9]+3
    figure('WindowStyle','docked')
    for rr = 1:size(WTTable,1)
        [~,idx] = ismember(str2double(WTTable{rr,"animal"}),WTUnique);
        scatter(WTTable{rr,"age_im_d"},WTTable{rr,ff},18,colors(2,:),marker_list{idx},'LineWidth',1.25);
        hold on
    end
    WTx = WTTable{:,"age_im_d"}; WTy = WTTable{:,ff};
    WTmdl = fitlm(WTx,WTy);
    WTls = feval(WTmdl,WTx);
    hold on
    plot(WTx,WTls,'Color',colors(2,:));

    WTrs = WTmdl.Rsquared.Ordinary;
    text(330,feval(WTmdl,330),sprintf('WT: r^2 = %.2f',WTrs),"Color",colors(2,:));


    legWT = scatter(0,0,[],colors(2,:),'_','Visible','on','LineWidth',2);

    for rr = 1:size(KOTable,1)
        [~,idx] = ismember(str2double(KOTable{rr,"animal"}),KOUnique);
        scatter(KOTable{rr,"age_im_d"},KOTable{rr,ff},18,colors(4,:),marker_list{idx},'LineWidth',1.25);
        hold on
    end
    KOx = KOTable{:,"age_im_d"}; KOy = KOTable{:,ff};
    KOmdl = fitlm(KOx,KOy);
    KOls = feval(KOmdl,KOx);
    hold on
    plot(KOx,KOls,'Color',colors(4,:));
    KOrs = KOmdl.Rsquared.Ordinary;
    text(330,feval(KOmdl,330),sprintf('KO: r^2 = %.2f',KOrs),"Color",colors(4,:));


    legWKO = scatter(0,0,[],colors(4,:),'_','Visible','on','LineWidth',2);
    title(featureNamesPlot{ff})
    xlabel('Age at imaging (days)')
    ylabel(featureNames{ff})
    legend([legWT,legWKO],{'WT','KO'})
    xlim([50,350])
    set(gca,'Units','Inches','InnerPosition',[0.5,0.8,2,1.5],'TickDir','out','TickLength',[0.03, 0.025],'Color','none','Box','off','LineWidth',1)
    % saveas(gcf,fullfile(figPath,sprintf('age_at_imaging_vs_%s.pdf',featureNamesPlot{ff})))
    stats(a).feature = featureNames{ff};
    stats(a).WTr_2 = WTrs;
    stats(a).WTp = WTmdl.Coefficients{2,4};
    stats(a).KOr_2 = KOrs;
    stats(a).KOp = KOmdl.Coefficients{2,4};
  a = a+1;
end

%%
for ff = [1:5,9]+3
    figure('WindowStyle','docked')
    for rr = 1:size(WTTable,1)
        [~,idx] = ismember(str2double(WTTable{rr,"animal"}),WTUnique);
        scatter(WTTable{rr,"age_virus_d"},WTTable{rr,ff},[],colors(2,:),marker_list{idx},'LineWidth',1.25);
        hold on
    end
    legWT = scatter(0,0,[],colors(2,:),'_','Visible','on','LineWidth',2);
    for rr = 1:size(KOTable,1)
        [~,idx] = ismember(str2double(KOTable{rr,"animal"}),KOUnique);
        scatter(KOTable{rr,"age_virus_d"},KOTable{rr,ff},[],colors(4,:),marker_list{idx},'LineWidth',1.25);
        hold on
    end
    legWKO = scatter(0,0,[],colors(4,:),'_','Visible','on','LineWidth',2);
    title(featureNamesPlot{ff})
    xlabel('Age at virus infusion (days)')
    ylabel(featureNames{ff})
    legend([legWT,legWKO],{'WT','KO'})
    xlim([40,250])
    saveas(gcf,fullfile(figPath,sprintf('age_at_virus_infusion_vs_%s.png',featureNamesPlot{ff})))
end
