%% Description

% This code runs through all mice - ball data and accumulates the pearsons
% c orrelation for delta F/F trace and assymmetric correlation and jaccards
% index for % Binarized trace where 1 == calcium event occurs to peak; 0
% otherwise

% Assymmetric correlation - intersection(A,B)/n(A); % custom (in utils)
% Jaccards Index - intersection(A,B)/union(A,B); % Built in

% use fuzzy speed thresholds new fuzzy parameters

% Jan 18th 2021 - by Athif Mohamed
% new data

% par for loop
%% Initialize

addpath('J:\nexmif_paper\Utils')  % Add utilities
init % Initialize data directories and genotypes
close all

%% Begin code

newIdx = 1;
dIdx = [1,3,5];
conditionList = {'Ball','Platform','Tonepuff'};

t_start = tic;

miceBad = [];
%% create struct
List_animal = cell(1,numel(miceStudy)*3);
List_day = cell(1,numel(miceStudy)*3);
List_condition = cell(1,numel(miceStudy)*3);

List_corrMatPearsonLow = cell(1,numel(miceStudy)*3);
List_corrMatPearsonThreshLow_pos = cell(1,numel(miceStudy)*3);
List_corrMatPearsonThreshLow_neg = cell(1,numel(miceStudy)*3);

List_corrMatPearsonHigh = cell(1,numel(miceStudy)*3);
List_corrMatPearsonThreshHigh_pos = cell(1,numel(miceStudy)*3);
List_corrMatPearsonThreshHigh_neg = cell(1,numel(miceStudy)*3);

parfor k = [46,47,48]% 1:numel(miceStudy)*3
    m = floor((k+2)/3);
    d = mod(k,3);
    if d==0
        d = 3;
    end
    [m,d]
    c = 1;%:3 %:3 %ball
    % Load files
    mPath = fullfile(dataFolder,sprintf('fullData_%s_D%i_%s',miceStudy{m},dIdx(d),lower(conditionList{c})));

    if ismember(sprintf('%s_D%i',miceStudy{m},dIdx(d)),badSessions)||ismember(sprintf('%s_D%i',miceStudy{m},dIdx(d)),badSpeedSessions)||ismember(sprintf('%s_D%i',miceStudy{m},dIdx(d)),fewMovSessions)
        continue
    end
    % handle missing files
    try
        data_var = load(mPath);
        fullData = data_var.fullData;



        % display
        {miceStudy{m}, num2str(dIdx(d)), conditionList{c}}

        goodTraces = fullData.goodIdx;
        goodTraces = setdiff(goodTraces,fullData.empty_idx_minusBG_new);


        %% Get high and low speed times - old
        %                 idx_onset = fullData.idx_onset;
        %                 idx_offset = fullData.idx_offset;
        %                 [paired_onsets,paired_offsets,lone_onsets,lone_offsets] = matchOnsetOffset(idx_onset,idx_offset);
        %                 traceLen = numel(fullData.roi_list_minusBG_new(1).trace);
        %
        %                 speedHighIdx = [];
        %                 for pp = 1:numel(paired_onsets)
        %                     speedHighIdx = [speedHighIdx,idx_onset(paired_onsets(pp)):idx_offset(paired_offsets(pp))];
        %                 end
        %
        %                 speedLowIdx = 1:traceLen;
        %                 speedLowIdx = setdiff(speedLowIdx, speedHighIdx);
        %% Get high and low speed times - fuxxy
        movBoutIdx = fullData.movBoutIdx;
        nonMovBoutIdx = fullData.restBoutIdx;



        %% Correlations for high and low speeds
        % Get binarized trace
        binTracesAll = binarizeTrace(fullData.roi_list_minusBG_new);
        for ss = 1:2  % 1 - low , 2 - high
            switch ss
                case 1
                    idxs =  nonMovBoutIdx;
                case 2
                    idxs =  movBoutIdx;
            end

            %% Binarized trace  - Jaccards Index

            binTraces = binTracesAll(goodTraces,idxs);
            if isempty(binTraces)
                corrMatAsymThresh = [];
                corrMatPearson = [];
            else
                [corrMatPearsonThresh_pos,corrMatPearsonThresh_neg,corrMatPearson] = correlation_circShift_pearson(2000,binTraces);
            end


            %% Add into struct
            List_animal{1,k} = miceStudy{m};
            List_day{1,k} = dIdx(d);
            List_condition{1,k} = conditionList{c};

            switch ss
                case 1
                    List_corrMatPearsonLow{1,k} = corrMatPearson;
                    List_corrMatPearsonThreshLow_pos{1,k} = corrMatPearsonThresh_pos;
                    List_corrMatPearsonThreshLow_neg{1,k} = corrMatPearsonThresh_neg;
                case 2
                    List_corrMatPearsonHigh{1,k} = corrMatPearson;
                    List_corrMatPearsonThreshHigh_pos{1,k} = corrMatPearsonThresh_pos;
                    List_corrMatPearsonThreshHigh_neg{1,k} = corrMatPearsonThresh_neg;
            end

        end
    catch
        miceBad = [miceBad,[miceStudy(m);dIdx(d);{lower(conditionList{c})}]];
    end
end

corrStats = table(List_animal,...
    List_day,...
    List_condition,...
    List_corrMatPearsonLow,...
    List_corrMatPearsonThreshLow_pos,...
    List_corrMatPearsonThreshLow_neg,...
    List_corrMatPearsonHigh ,...
    List_corrMatPearsonThreshHigh_pos,...
    List_corrMatPearsonThreshHigh_neg);


% save(fullfile('D:\Autism\Correlation analysis\Results\CorrstatsSignificant_speed','corrStats'),'corrStats')
% time_elapsed = toc(t_start)
varNames = corrStats.Properties.VariableNames;
corrStatsTable = table();
for tt = 1:numel(varNames)
    corrStatsTable{:,tt} = corrStats{1,tt}';
    corrStatsTable.Properties.VariableNames{tt} = varNames{tt}(6:end);
end
save(fullfile('J:\nexmif_paper\code_ball\stats\correlation\','corrStatsTableSpeedPearson_08_22'),'corrStatsTable')

time_spent = toc(t_start)