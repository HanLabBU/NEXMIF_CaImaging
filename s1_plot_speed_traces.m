%% Initialize

addpath(genpath('U:\eng_research_handata\Athif Mohamed\nexmif_paper\Utils\'))  % Add utilities
init % Initialize data directories and genotypes
close all

for m = 1: numel(miceStudy)
    for d = 1:3
        c = 1;  % ball only
        % Load files
        
        mPath = fullfile(dataFolder,sprintf('fullData_%s_D%i_%s',miceStudy{m},dIdx(d),lower(conditionList{c})));
        
        if ismember(sprintf('%s_D%i',miceStudy{m},dIdx(d)),badSessions)||ismember(sprintf('%s_D%i',miceStudy{m},dIdx(d)),badSpeedSessions)||ismember(sprintf('%s_D%i',miceStudy{m},dIdx(d)),fewMovSessions)
            continue
        end
        
        % handle missing files
        try
            load(mPath)
        catch
            miceBad = [miceBad,[miceStudy(m);dIdx(d);{lower(conditionList{c})}]];
        end
        
        % display
        {miceStudy{m}, num2str(dIdx(d)), conditionList{c}}
        
        % collect traces
        nTraces = size(fullData.roi_list_minusBG_new,2);
        traceLen = size(fullData.roi_list_minusBG_new(1).trace,2);
        
        % color
        
        if ismember(miceStudy(m),miceKO)
            col_l = colors(3,:);
            col_h = colors(4,:);
        else
            col_l = colors(1,:);
            col_h = colors(2,:);
        end
        
        % plot speed
        figure
        set(gcf,'units','normalized','outerposition',[0 0.5 1 0.35])
        
        yy = max(fullData.speed);
        hold on
        for gg = 1:numel(fullData.movBoutStart)
            shadeXX = (fullData.movBoutStart(gg):fullData.movBoutFinish(gg));
            shadeYY = repmat(yy,size(shadeXX));
            area(shadeXX,shadeYY,'FaceColor',col_h,'FaceAlpha',0.75,'EdgeAlpha',0)
            hold on
            
        end
        for gg = 1:numel(fullData.restBoutStart)
            shadeXX = (fullData.restBoutStart(gg):fullData.restBoutFinish(gg));
            shadeYY = repmat(yy,size(shadeXX));
            area(shadeXX,shadeYY,'FaceColor',col_l,'FaceAlpha',0.75,'EdgeAlpha',0)
            hold on
        end
        hold on
        plot(fullData.speed,'k')
        xx = xticks;
        set(gca,'XTickLabel',xx/Fs);
        ylabel('Speed (cm/s)')
        xlabel('Time(s)')
        title(sprintf('%s - %s - Day %i',miceStudy{m},mType{maskKO(m)+1},dIdx(d)))
        xlim([0 600*20])
        saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\supplementary\s1',sprintf('%s - %s - Day %i - responsive cells.png',miceStudy{m},mType{maskKO(m)+1},dIdx(d))))
        saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\supplementary\s1',sprintf('%s - %s - Day %i - responsive cells.fig',miceStudy{m},mType{maskKO(m)+1},dIdx(d))))
        
        close all
        
    end
end




