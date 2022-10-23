%% Description

%

% Feb 2nd 2021 - by Athif Mohamed
%% Initialize

addpath('D:\nexmif_paper\Utils')  % Add utilities
init % Initialize data directories and genotypes
close all

%% Begin code
% cd(pathData)
t_start = tic;
miceBad = [];

for m = 16%:numel(miceStudy)
    for d = 1:3
        for c = 1 %:3 %ball
            
            mPath = fullfile(dataFolder,sprintf('fullData_%s_D%i_%s',miceStudy{m},dIdx(d),lower(conditionList{c})));
            
            
            % handle missing files
            try
                load(mPath)
           
            
            
            % Load files
%             mPath = fullfile('D:\Autism\Event analysis\Data_new2_cleaned\',sprintf('fullData_%s_D%i_%s',fo(mIdx).name,dIdx(d),lower(conditionList{c})));
                        dPath = fullfile('U:\eng_research_handata\eng_research_handata2\Rebecca_Mount\Autism\',...
                            sprintf('%s\\Day_%i\\%s\\roi_%s_D%i_%s',miceStudy{m},dIdx(d)+1,lower(conditionList{c}),miceStudy{m},dIdx(d)+1,lower(conditionList{c})));
%             dPath = fullfile('U:\eng_research_handata\eng_research_handata2\Rebecca_Mount\Autism\',...
%                 sprintf('%s\\Day_%i\\%s\\roi_%s_D%i_%s',miceStudy{m},dIdx(d),lower(conditionList{c}),miceStudy{m},dIdx(d),lower(conditionList{c})));
            % handle missing files
            %             try

            load(dPath)
            
             
            % display
            {miceStudy{m}, num2str(dIdx(d)), conditionList{c}}
            
            
            %% Get centroid and save
            centroids = [];
            for cc = 1:size(roiList,2)
                idx = roiList(cc).pixel_idx;
                [y,x] = ind2sub([1024,1024],idx);
                centroids = [centroids;[mean(x),mean(y)]];
            end
            
            fullData.centroids = centroids;
            %                 %% test
            %                 load('U:\eng_research_handata\eng_research_handata2\Rebecca_Mount\Autism\2089\Day_1\Ball\maxmin_2089_D1_ball.mat')
            %                 im = mat2gray(maxmin);
            %                 imshow(im)
            %                 hold on
            %                 scatter(centroid(:,1),centroid(:,2));
            %%
            save(mPath,'fullData');
            %                 catch
            %                 miceBad = [miceBad;[m,d]];
            %             end
            catch
                miceBad = [miceBad,[miceStudy(m);dIdx(d);{lower(conditionList{c})}]];
            end
        end
    end
end

% save(fullfile('D:\Autism\Correlation analysis\Results\CorrstatsUnwrapped','corrStats'),'corrStats')
time_elapsed = toc(t_start)