% Gets all events
% Saves new fullData structure with these information
% Generates a event stat structure for each trial
% ball - slow and fast
% platform
% tonepuff first 15 seconds - 2nd 15 seconds

% Extension 1
% For background subtracted

% Extension 2
% Event width analysis

% Extension 3
% For Ball data cleaned

% Extension 4
% Add event rates during slow and fast

% Extension 6
% Update the fuzzy threshold

% ext 7
% new data
%%
addpath(genpath('J:\nexmif_paper\Utils'))  % Add utilities
init

newIdx = 1;
dIdx = [1,3,5];
conditionList = {'Ball','Platform','Tonepuff'};

miceBad = [];
eventStatsAll = struct();
eventStatsList = struct();
k = 0;

for m = 1:numel(miceStudy)
    for d = 1:3
        for c = 1%:2 % Ball
            % Load files

            mPath = fullfile(dataFolder,sprintf('fullData_%s_D%i_%s',miceStudy{m},dIdx(d),lower(conditionList{c})));

            if ismember(sprintf('%s_D%i',miceStudy{m},dIdx(d)),badSessions)
                continue
            end
            % handle missing files
            try
                load(mPath)


                % display
                {miceStudy{m}, num2str(dIdx(d)), conditionList{c}}


                goodTraces = fullData.goodIdx;
                goodTraces = setdiff(goodTraces,fullData.empty_idx_minusBG_new);


                % fuzzy high movement
                movBoutStart = fullData.movBoutStart;
                movBoutFinish = fullData.movBoutFinish;
                restBoutStart = fullData.restBoutStart;
                restBoutFinish = fullData.restBoutFinish;

                % binarised trace
                binTraces = binarizeTrace(fullData.roi_list_minusBG_new);

                for n = 1:size(fullData.roi_list_minusBG_new,2) %randperm(size(fullData.roi_list_minusBG_new,2),3)  %
                    if ismember(n,goodTraces)

                        k = k+1

                        % Assign basic details
                        eventStatsAll(k).condition = conditionList{c};
                        eventStatsAll(k).animal = miceStudy{m};
                        eventStatsAll(k).day = dIdx(d);
                        eventStatsAll(k).cellNumber = n;
                        eventStatsAll(k).nCells = size(fullData.roi_list_minusBG_new,2);
                        eventStatsAll(k).nGoodCells = numel(goodTraces);
                        % Gather event statistics

                        nEv = numel(fullData.roi_list_minusBG_new(n).event_amp);
                        traceLen = numel(fullData.roi_list_minusBG_new(n).trace);
                        EventHeight = mean(fullData.roi_list_minusBG_new(n).event_amp);
                        if nEv >=2
                            evInts = fullData.roi_list_minusBG_new(n).event_time(2:end,1) - fullData.roi_list_minusBG_new(n).event_time(1:end-1,1);
                        else
                            evInts = [];
                        end
                        InterEventInterval = mean(evInts);

                        DeconvRiseTimeList = fullData.roi_list_minusBG_new(n).event_time(:,2) - fullData.roi_list_minusBG_new(n).event_time(:,1);%*1000;
                        DeconvRiseTime = mean(DeconvRiseTimeList);
                        EventRate = nEv/(traceLen/20)*60;
                        % AreaUnderRiseTime =  mean((fullData.roi_list_preproc_minusBG(n).event_time(evIdxs,2) - fullData.roi_list_preproc_minusBG(n).event_time(evIdxs,1)).*...
                        % fullData.roi_list_preproc_minusBG(n).event_amp(evIdxs);
                        fwhmList = eventWidthFWHM_ext2(fullData.roi_list_minusBG_new(n).trace,fullData.roi_list_minusBG_new(n).event_idx,fullData.roi_list_minusBG_new(n).event_fall,Fs);
                        meanFWHM = mean(fwhmList,'omitnan');




                        ca_onsets = fullData.roi_list_minusBG_new(n).event_idx(:,1);
                        ca_binary = binTraces(n,:);
                        EventActRate = sum(ca_binary)/(traceLen/Fs)*60;

                        % event rate during high speed - Fuzzy

                        EventCountHigh = 0;
                        EventActHigh =0;
                        DeconvRiseTimeHighList = [];
                        fwhmHighList = [];
                        for pp = 1: numel(movBoutStart)
                            EventCountHigh = EventCountHigh + sum(ca_onsets>= movBoutStart(pp) & ca_onsets<= movBoutFinish(pp));
                            EventActHigh = EventActHigh + sum(ca_binary(movBoutStart(pp):movBoutFinish(pp)));
                            DeconvRiseTimeHighList = [DeconvRiseTimeHighList; DeconvRiseTimeList(fullData.roi_list_minusBG_new(n).event_idx(:,1)>= movBoutStart(pp)& ...
                                fullData.roi_list_minusBG_new(n).event_idx(:,2)<= movBoutFinish(pp))];
                            fwhmHighList = [fwhmHighList; fwhmList(fullData.roi_list_minusBG_new(n).event_idx(:,1)>= movBoutStart(pp)& ...
                                fullData.roi_list_minusBG_new(n).event_idx(:,2)<= movBoutFinish(pp))];
                        end

                        DurationHigh = sum((movBoutFinish - movBoutStart)/Fs);
                        EventRateHigh = EventCountHigh/DurationHigh*60;
                        EventActRateHigh = EventActHigh/DurationHigh*60;
                        DeconvRiseTimeHigh = mean(DeconvRiseTimeHighList);
                        meanFWHMHigh = mean(fwhmHighList,'omitnan');
                        % event rate during low speed - Fuzzy
                        EventCountLow = 0;
                        EventActLow =0;
                        DeconvRiseTimeLowList = [];
                        fwhmLowList = [];
                        for pp = 1: numel(restBoutStart)
                            EventCountLow = EventCountLow + sum(ca_onsets>= restBoutStart(pp) & ca_onsets<= restBoutFinish(pp));
                            EventActLow = EventActLow + sum(ca_binary(restBoutStart(pp): restBoutFinish(pp)));
                            DeconvRiseTimeLowList = [DeconvRiseTimeLowList; DeconvRiseTimeList(fullData.roi_list_minusBG_new(n).event_idx(:,1)>= restBoutStart(pp)& ...
                                fullData.roi_list_minusBG_new(n).event_idx(:,2)<= restBoutFinish(pp))];
                            fwhmLowList = [fwhmLowList; fwhmList(fullData.roi_list_minusBG_new(n).event_idx(:,1)>= restBoutStart(pp)& ...
                                fullData.roi_list_minusBG_new(n).event_idx(:,2)<= restBoutFinish(pp))];
                        end

                        DurationLow = sum((restBoutFinish - restBoutStart)/Fs);
                        EventRateLow = EventCountLow/DurationLow*60;
                        EventActRateLow = EventActLow/DurationLow*60;
                        DeconvRiseTimeLow = mean(DeconvRiseTimeLowList);
                        meanFWHMLow = mean(fwhmLowList,'omitnan');

                        % Assign event statistics
                        eventStatsAll(k).EventHeight = EventHeight;
                        eventStatsAll(k).InterEventInterval = InterEventInterval;
                        eventStatsAll(k).DeconvRiseTime = DeconvRiseTime;
                        eventStatsAll(k).meanFWHM = meanFWHM;
                        eventStatsAll(k).EventRate = EventRate;
                        eventStatsAll(k).EventActRate = EventActRate;


                        % Assign basic details
                        eventStatsList(k).condition = conditionList{c};
                        eventStatsList(k).animal = miceStudy{m};
                        eventStatsList(k).day = dIdx(d);
                        eventStatsList(k).cellNumber = n;
                        eventStatsList(k).nCells = size(fullData.roi_list_minusBG_new,2);
                        eventStatsList(k).nGoodCells = numel(goodTraces);
                        %  Add Lists
                        eventStatsList(k).FWHM = fwhmList;
                        eventStatsList(k).Risetime = ((fullData.roi_list_minusBG_new(n).event_time(:,2) - fullData.roi_list_minusBG_new(n).event_time(:,1)))';

                        if ismember(sprintf('%s_D%i',miceStudy{m},dIdx(d)),badSpeedSessions)||ismember(sprintf('%s_D%i',miceStudy{m},dIdx(d)),fewMovSessions)
                            eventStatsAll(k).EventRateHigh = nan;
                            eventStatsAll(k).EventRateLow = nan;
                            eventStatsAll(k).EventActRateHigh = nan;
                            eventStatsAll(k).EventActRateLow = nan;
                            eventStatsAll(k).DeconvRiseTimeHigh = nan;
                            eventStatsAll(k).DeconvRiseTimeLow = nan;
                            eventStatsAll(k).meanFWHMHigh = nan;
                            eventStatsAll(k).meanFWHMLow = nan;
                        else
                            eventStatsAll(k).EventRateHigh = EventRateHigh;
                            eventStatsAll(k).EventRateLow = EventRateLow;
                            eventStatsAll(k).EventActRateHigh = EventActRateHigh;
                            eventStatsAll(k).EventActRateLow = EventActRateLow;
                            eventStatsAll(k).DeconvRiseTimeHigh = DeconvRiseTimeHigh;
                            eventStatsAll(k).DeconvRiseTimeLow = DeconvRiseTimeLow;
                            eventStatsAll(k).meanFWHMHigh = meanFWHMHigh;
                            eventStatsAll(k).meanFWHMLow = meanFWHMLow;
                        end

                    end
                end
            catch
                miceBad = [miceBad,[miceStudy(m);dIdx(d);{lower(conditionList{c})}]];
            end
        end
    end
end

save(fullfile('J:\nexmif_paper\code_ball\stats\event_features','eventStats_08_22'),'eventStatsAll')
% save(fullfile('D:\nexmif_paper\code_ball\stats\event_features','eventStats_List'),'eventStatsList')
