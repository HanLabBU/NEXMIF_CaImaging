
%% Initialize

addpath(genpath('U:\eng_research_handata\Athif Mohamed\nexmif_paper\Utils\'))  % Add utilities
init % Initialize data directories and genotypes
close all

load('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\stats\movement_responsive_cells\responsive_cells_bout_shuffle','responsive_cells')
responsive_cells = struct2table(responsive_cells);

for m = 9%: numel(miceStudy)             % using 2978 day 3 (WT) and 605842 day 1 (KO)
    for d = 2%1:3
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
        
        binTraces = binarizeTrace(fullData.roi_list_minusBG_new);
        
        
        idx_plotted = [];
        
        % Select movement responsive cells
        
        isM = strcmp([responsive_cells.animal],miceStudy{m});
        isD = [responsive_cells.day] == dIdx(d);
        isR = [responsive_cells.responsive_bout];
        cell_list = [responsive_cells.cell];
        mov_r_list = [responsive_cells.mov_r];
        res_r_list = [responsive_cells.res_r];
        isR(isnan(isR)) = 0;
        res_cell_list_idx = find(isM & isD & isR);
        nonres_cell_list_idx = find(isM & isD & ~isR);
        tracePerPlot = 20;
        
        %% plot responsive
        hh1 = 0.25*numel(res_cell_list_idx)/(numel(res_cell_list_idx)+numel(nonres_cell_list_idx ));
        hh2 = 0.25*numel(nonres_cell_list_idx)/(numel(res_cell_list_idx)+numel(nonres_cell_list_idx ));
        figure %('Units','inches','PaperSize',[8.5 11],'PaperType','usletter')
        set(gcf,'Units','inches','Position',[0,0,8.5,11])
        set(gcf,'units','normalized')
        ax1 = subplot('Position',[0.05, hh2+0.625, 0.4, hh1]);
%      
        imagesc(binTraces(cell_list(res_cell_list_idx),:));

        bluemap = multigradient([1,1,1;colors(2,:)], [1,2]);
        redmap = multigradient([1,1,1;colors(4,:)], [1,2]);
        colormap(ax1,bluemap)
%         colormap(flipud(gray))


        ylabel('Responsive Cells')
        %set(gca,'YTick',1:numel(res_cell_list_idx ),'YTickLabel',cell_list(res_cell_list_idx ));
        xlim([1,600.*Fs])
        xx = xticks;
        set(gca,'XTickLabel',xx/Fs);
        xl = xlim;
        title(sprintf('%s - %s - Day %i - responsive cells',miceStudy{m},mType{maskKO(m)+1},dIdx(d)))
        rightLabel = strsplit(sprintf('%.1f - %.1f_', [mov_r_list(res_cell_list_idx )';res_r_list(res_cell_list_idx )']),'_');
        yt = yticks;
%         text((xl(end)+10)*ones(size(yt)),yt, rightLabel(1:end-1),'FontSize',8);
        
        hold on
        yy = ylim;
        ev_rate_trace = zeros(size(binTraces));
        for gg = 1:numel(fullData.movBoutStart)
            shadeXX = (fullData.movBoutStart(gg):fullData.movBoutFinish(gg));
            shadeYY = repmat(yy(2),size(shadeXX));
            area(shadeXX,shadeYY,'FaceColor',[0 0 0],'FaceAlpha',0.25,'EdgeAlpha',0)
            hold on
            hXX = sum(binTraces(cell_list(res_cell_list_idx),shadeXX),2)/numel(shadeXX)*100 ;
            ev_rate_trace(cell_list(res_cell_list_idx),shadeXX) = repmat(hXX,[1,numel(shadeXX)]);
        end
        
        for gg = 1:numel(fullData.restBoutStart)
            shadeXX = (fullData.restBoutStart(gg):fullData.restBoutFinish(gg));
            shadeYY = repmat(yy(2),size(shadeXX));
%             area(shadeXX,shadeYY,'FaceColor',[1 1 1],'FaceAlpha',0.25,'EdgeAlpha',0)
%             hold on
            hXX = sum(binTraces(cell_list(res_cell_list_idx),shadeXX),2)/numel(shadeXX)*100;
            ev_rate_trace(cell_list(res_cell_list_idx),shadeXX) = repmat(hXX,[1,numel(shadeXX)]);
        end
        set(gca,'TickDir','out','Box','off','LineWidth',2)

        %% plot  non responsive

        ax2 = subplot('Position',[0.05, 0.6, 0.4, hh2]);
        %set(gcf,'units','normalized','outerposition',[0 0 1 1])
        imagesc(binTraces(cell_list(nonres_cell_list_idx),:));
        colormap(ax2,flipud(gray))
        ylabel('Non-responsive Cells')
%         set(gca,'YTick',1:numel(nonres_cell_list_idx ),'YTickLabel',cell_list(nonres_cell_list_idx ));
         xlim([1,600.*Fs])
        xx = xticks;
        set(gca,'XTickLabel',xx/Fs);
        xl = xlim;
        title(sprintf('%s - %s - Day %i - non responsive cells',miceStudy{m},mType{maskKO(m)+1},dIdx(d)))
        rightLabel = strsplit(sprintf('%.1f - %.1f_', [mov_r_list(nonres_cell_list_idx )';res_r_list(nonres_cell_list_idx )']),'_');
        yt = yticks;
%         text((xl(end)+10)*ones(size(yt)),yt, rightLabel(1:end-1),'FontSize',8);
        
        hold on
        yy = ylim;
        ev_rate_trace = zeros(size(binTraces));
        for gg = 1:numel(fullData.movBoutStart)
            shadeXX = (fullData.movBoutStart(gg):fullData.movBoutFinish(gg));
            shadeYY = repmat(yy(2),size(shadeXX));
            area(shadeXX,shadeYY,'FaceColor',[0 0 0],'FaceAlpha',0.25,'EdgeAlpha',0)
            hold on
            hXX = sum(binTraces(cell_list(res_cell_list_idx),shadeXX),2)/numel(shadeXX)*100 ;
            ev_rate_trace(cell_list(res_cell_list_idx),shadeXX) = repmat(hXX,[1,numel(shadeXX)]);
            hXXnon = sum(binTraces(cell_list(nonres_cell_list_idx),shadeXX),2)/numel(shadeXX)*100 ;
            ev_rate_trace(cell_list(nonres_cell_list_idx),shadeXX) = repmat(hXXnon,[1,numel(shadeXX)]);
        end
        
        for gg = 1:numel(fullData.restBoutStart)
            shadeXX = (fullData.restBoutStart(gg):fullData.restBoutFinish(gg));
            shadeYY = repmat(yy(2),size(shadeXX));
%             area(shadeXX,shadeYY,'FaceColor',[1 1 1],'FaceAlpha',0.25,'EdgeAlpha',0)
%             hold on
            hXX = sum(binTraces(cell_list(res_cell_list_idx),shadeXX),2)/numel(shadeXX)*100;
            ev_rate_trace(cell_list(res_cell_list_idx),shadeXX) = repmat(hXX,[1,numel(shadeXX)]);
            hXXnon = sum(binTraces(cell_list(nonres_cell_list_idx),shadeXX),2)/numel(shadeXX)*100 ;
            ev_rate_trace(cell_list(nonres_cell_list_idx),shadeXX) = repmat(hXXnon,[1,numel(shadeXX)]);
        end
        set(gca,'TickDir','out','Box','off','LineWidth',2)
   
   
        % Plot mean event density
    
%         subplot('Position',[0.05, 0.7, 0.45, 0.05])
        subplot('Position',[0.05, 0.525, 0.4, 0.05])
        plot(mean(ev_rate_trace(cell_list(nonres_cell_list_idx),:)),'Color',[0,0,0],'LineWidth',1)      % blue
        hold on 
        plot(mean(ev_rate_trace(cell_list(res_cell_list_idx),:)),'Color',colors(2,:),'LineWidth',1)  % red
        hold on 
%         plot(mean(ev_rate_trace(cell_list([res_cell_list_idx;nonres_cell_list_idx]),:)),'Color',[159 97 168]/255,'LineWidth',1)
%         hold on 
        xlim(xl)
        xx = xticks;
        set(gca,'XTickLabel',xx/Fs);
        ylabel('Averaged dF/F (AU)')
%         legend('Non-responsive cells','Responsive cells','All cells')
        set(gca,'TickDir','out','Box','off','LineWidth',2)
        
          yy = ylim;

         for gg = 1:numel(fullData.movBoutStart)
            shadeXX = (fullData.movBoutStart(gg):fullData.movBoutFinish(gg));
            shadeYY = repmat(yy(2),size(shadeXX));
            area(shadeXX,shadeYY,'FaceColor',[0 0 0],'FaceAlpha',0.15,'EdgeAlpha',0)
            hold on
         end
         ylim(yy)






%         % Plot mean fluoresence
%         fl_traces = {fullData.roi_list_minusBG_new.trace}';
%         fl_traces = cell2mat(fl_traces);
%     
%         subplot('Position',[0.05, 0.2, 0.9, 0.1])
%         plot(mean(fl_traces(cell_list(nonres_cell_list_idx),:)),'Color',[62 111 255]/255,'LineWidth',1)      % blue
%         hold on 
%         plot(mean(fl_traces(cell_list(res_cell_list_idx),:)),'Color',[255 26 26]/255,'LineWidth',1)  % red
%         hold on 
%         plot(mean(fl_traces(cell_list([res_cell_list_idx;nonres_cell_list_idx]),:)),'Color',[159 97 168]/255,'LineWidth',1)
%         hold on 
%         xlim(xl)
%         xx = xticks;
%         set(gca,'XTickLabel',xx/Fs);
%         ylabel('Averaged dF/F (AU)')
%         legend('Non-responsive cells','Responsive cells','All cells')
%         set(gca,'FontSize',12,'TickDir','out','Box','off','LineWidth',2)
        
%         % plot speed
%         subplot('Position',[0.05, 0.05, 0.9, 0.1])
%         plot(fullData.speed)
%         xlim(xl)
%         xx = xticks;
%         set(gca,'XTickLabel',xx/Fs);
%         ylabel('Speed (cm/s)')
%         set(gca,'FontSize',12,'TickDir','out','Box','off','LineWidth',1)
%     
        
%         saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\supplementary\s7',sprintf('%s - %s - Day %i - responsive cells.png',miceStudy{m},mType{maskKO(m)+1},dIdx(d))))
%         saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\supplementary\s7',sprintf('%s - %s - Day %i - responsive cells.fig',miceStudy{m},mType{maskKO(m)+1},dIdx(d))))
%         saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\supplementary\s7',sprintf('%s - %s - Day %i - responsive cells.epsc',miceStudy{m},mType{maskKO(m)+1},dIdx(d))))
        close all

    end
end




