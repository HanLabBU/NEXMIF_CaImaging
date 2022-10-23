% Plots of WT vs KO for movement stats


% Ext 2
% boxplot to include all animals
% Add stats

% Ext 5
% Fuzzy high vs low speed threshold

addpath('J:\nexmif_paper\Utils')  % Add utilities
init % Initiate variables

%% Add mov bouts and onsets in the fullData
isPlot = 0;
thresh_list = [];
miceBad = [];
for m = 1: numel(miceStudy)
    for d = 1:3
        for c = 1 % ball only
            % Load files
            
            mPath = fullfile(dataFolder,sprintf('fullData_%s_D%i_%s',miceStudy{m},dIdx(d),lower(conditionList{c})));
            
            if ismember(sprintf('%s_D%i',miceStudy{m},dIdx(d)),badSessions)%||ismember(sprintf('%s_D%i',miceStudy{m},dIdx(d)),badSpeedSessions)||ismember(sprintf('%s_D%i',miceStudy{m},dIdx(d)),fewMovSessions)
                continue
            end
            
            % handle missing files
            try
                load(mPath)
            
            
            
            % display
            {miceStudy{m}, num2str(dIdx(d)), conditionList{c}}
            
            
            %% Do fuzzy threshold
            
            speed = fullData.speed;
            
            percent_high = 0.2;
            smooth_filt_width = 30; % 1.5s
            sigmoid_slope = 0.8;
            sigmoid_threshold = 0.1;
            on_dur_threshold = 40;
            off_dur_threshold = 40; % 2s
            
            mouse = miceStudy{m};
            mouseType = mType{maskKO(m)+1};
            day = dIdx(d);
            
            
            if isPlot
                figure
                set(gcf,'units','normalized','outerposition',[0 0 1 1])
                %                 figure
                %                 set(gcf,'WindowStyle','docked')
            end
            speed_smooth = movmean(speed,20);
            
            [start_idx_high,finish_idx_high,start_idx_low,finish_idx_low,movBoutIdx,restBoutIdx] = detectMovementHigh(speed,percent_high,smooth_filt_width,sigmoid_slope,sigmoid_threshold,on_dur_threshold,off_dur_threshold,mouse,mouseType,day,Fs,isPlot);% Detect onst and offset of movement
            if isPlot
                saveas(gcf,fullfile('D:\Autism\Event analysis\Results\SpeedPlots\',sprintf('Speed_%s_D%i_%s.fig',miceStudy{m},dIdx(d),lower(conditionList{c}))))
                saveas(gcf,fullfile('D:\Autism\Event analysis\Results\SpeedPlots\',sprintf('Speed_%s_D%i_%s.png',miceStudy{m},dIdx(d),lower(conditionList{c}))))
            end
            fullData.movBoutStart = start_idx_high;
            fullData.movBoutFinish = finish_idx_high;
            fullData.movBoutIdx = movBoutIdx;
            
            fullData.restBoutStart = start_idx_low;
            fullData.restBoutFinish = finish_idx_low;
            fullData.restBoutIdx = restBoutIdx;
            
            %% Do onset detection
            percent_high = 0.20;
            percent_low = 0.05;
            pre_thresh = 0.75;
            
            upper_thresh =  max(5,min(speed)+ percent_high*(max(speed)-min(speed)));
            lower_thresh =  min(5,min(speed)+ percent_low*(max(speed)-min(speed)));
            
            thresh_list = [thresh_list; [upper_thresh,lower_thresh]];
            onsetIdx = detectMovementOnset_basic(speed,lower_thresh,upper_thresh,pre_thresh,Fs,isPlot);
            if isPlot
                title(sprintf('Speed - %s - %s - Day %i',mouse,mouseType,day))
            end
            
            fullData.onsetIdx = onsetIdx;
            
            if isPlot
                % pause()
                saveas(gcf,fullfile('D:\Autism\Event analysis\Results\SpeedPlots\',sprintf('Onsets_%s_D%i_%s.fig',miceStudy{m},dIdx(d),lower(conditionList{c}))))
                saveas(gcf,fullfile('D:\Autism\Event analysis\Results\SpeedPlots\',sprintf('Onsets_%s_D%i_%s.png',miceStudy{m},dIdx(d),lower(conditionList{c}))))
            end
            save(fullfile(dataFolder,sprintf('fullData_%s_D%i_%s',miceStudy{m},dIdx(d),lower(conditionList{c}))),'fullData');
            catch
                miceBad = [miceBad,[miceStudy(m);dIdx(d);{lower(conditionList{c})}]];
            end
            close all
        end
    end
end

%% Plot distributuiions of threshold

% figure
% subplot(2,1,1)
% histogram(thresh_list(:,1),'BinWidth',1, 'FaceAlpha',0.4)
% subplot(2,1,2)
% histogram(thresh_list(:,2),'BinWidth',1, 'FaceAlpha',0.4)


%% gather movement stats

k = 0;
movStats = struct();
thresh_list = [];
miceBad = [];
for m = 1: numel(miceStudy)
    for d = 1:3
        for c = 1 % ball only
            % Load files
            
            mPath = fullfile(dataFolder,sprintf('fullData_%s_D%i_%s',miceStudy{m},dIdx(d),lower(conditionList{c})));
            
            if ismember(sprintf('%s_D%i',miceStudy{m},dIdx(d)),badSessions)||ismember(sprintf('%s_D%i',miceStudy{m},dIdx(d)),badSpeedSessions)%||ismember(sprintf('%s_D%i',miceStudy{m},dIdx(d)),fewMovSessions)
                continue
            end
            
            % handle missing files
            try
                load(mPath)
          
           
            
            % display
            {miceStudy{m}, num2str(dIdx(d)), conditionList{c}}
            
            durations = (fullData.movBoutFinish - fullData.movBoutStart)/Fs;
            nonMovBoutIdx = fullData.restBoutIdx;
            
            % Assign basic details
            k = k+1;
            
            movStats(k).animal = miceStudy{m};
            movStats(k).day = dIdx(d);
            
            
            if maskWT(m)
                movStats(k).gType = 1;   %1 WT
            else
                movStats(k).gType = 0;
            end
            
            movStats(k).averageSpeed = mean(fullData.speed);
            movStats(k).averageDurations = mean(durations);
            movStats(k).nMovBouts = numel(fullData.movBoutStart);
            movStats(k).nOnsets = numel(fullData.onsetIdx);
            
            movStats(k).movBoutSpeed = mean(fullData.speed(fullData.movBoutIdx));
            movStats(k).nonMovBoutSpeed = mean(fullData.speed(nonMovBoutIdx));
            
            movStats(k).movBoutTimeRatio = sum(fullData.movBoutIdx)/numel(fullData.speed)*100;
            movStats(k).nonMovBoutTimeRatio = sum(nonMovBoutIdx)/numel(fullData.speed)*100;
            catch
                miceBad = [miceBad,[miceStudy(m);dIdx(d);{lower(conditionList{c})}]];
            end
        end
    end
end

save(fullfile('J:\nexmif_paper\code_ball\stats\movement_stats','movStats_8_22'),'movStats');


%% Plots
close all
% durationsAll = [movStats.durations];
% durationsWT = [movStats(find([movStats.gType])).durations];
% durationsKO = [movStats(find(~[movStats.gType])).durations];

movStatsTable = struct2table(movStats);
miceAll = movStatsTable.animal;
miceAllUnique = unique(miceAll);


movStatsTableSummary = movStatsTable(1:14,[1,3:end]);

% Summarize for animals
for m = 1:numel(miceAllUnique)
    miceIdx = find(ismember(movStatsTable.animal,miceAllUnique(m)));
    movStatsTableSummary{m,3:end} = mean(movStatsTable{miceIdx,4:end},1,'omitnan');
    movStatsTableSummary{m,1} = miceAllUnique(m);
    movStatsTableSummary{m,2} =  movStatsTable.gType(miceIdx(1));
end

fieldNames = movStatsTableSummary.Properties.VariableNames;
featureNames = fieldNames(3:end);

featureNames(1:5) = {'Average speed (cm/s)', 'Average movement bout duration(s)','Number of movement bouts','Number of movement onsets','Average speed during movement bout (cm/s)'};
featureNamesPlot = featureNames;
featureNamesPlot(1:5) = {'Average speed', 'Average movement bout duration','Number of movement bouts','Number of movement onsets','Average speed during movement bout'};

for f = 1:5
    figure%('units','normalized','outerposition',[0 0 0.5 0.5])
    
    data = movStatsTableSummary{:,f+2};
    group = ismember([movStatsTableSummary.animal],miceStudy(maskKO));
    
    % subplot(2,1,1)
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
        
        scatter (x,y,'filled','markerFaceColor',colors(2*g-1,:))
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
    [h,p_t] = ttest2(yKO,yWT);
    
    %     [p_rsIn,h] = ranksum(yKOIn,yWTIn); % returns the p-value of a two-sided Wilcoxon rank sum test. ranksum tests the null hypothesis that data in x and y are samples from continuous distributions with equal medians, against the alternative that they are not. The test assumes that the two samples are independent. x and y can have different lengths.
    %     [h,p_tIn] = ttest2(yKOIn,yWTIn);
    
    
    up = 1.05*max([yKO',yWT']);
    %     down = 1.1*max([yKO',yWT']);
    
    hold on
    text(1.5-0.5,up,sprintf('p-ranksum = %.3f',p_rs),'HorizontalAlignment','center')
    hold on
    text(1.5+0.5,up,sprintf('p-ttest = %.3f',p_t),'HorizontalAlignment','center')
    
    %      hold on
    %     text(1.5-0.5,down,sprintf('Inliers p-ranksum = %.4f',p_rsIn),'HorizontalAlignment','center')
    %     hold on
    %     text(1.5+0.5,down,sprintf('Inliers p-ttest = %.4f',p_tIn),'HorizontalAlignment','center')
    
    
    ylabel(sprintf('%s',featureNames{f}))
    set(gca,'xtick',xx,'xticklabel',{'WT','KO'})
    box off
    saveas(gcf,fullfile('D:\nexmif_paper\code_ball\plots\movement_stats',sprintf('%s.fig',featureNamesPlot{f})))
    saveas(gcf,fullfile('D:\nexmif_paper\code_ball\plots\movement_stats',sprintf('%s.png',featureNamesPlot{f})))
end

