%% TO DO:
%   - Brain plots of stat vals
%   - Bar plot of means across electrodes for each condition
%       - statistics

subj = 'EC131';
rootdir = '/Users/mattleonard/Documents/Research/data';
brain_dir = [rootdir filesep 'MRI']; % [rootdir filesep '..' filesep 'pia/data_store2/imaging/subjects'];

params_sets = [25 2.25 500 ; ...
    25 2.25 250 ; ...
    25 1 250 ; ...
    10 2.25 250];

plot_avg_ERP_flag = 0;              % whether to plot mean ERPs
test_pre_post_dist_byCond_flag = 1; % whether to plot pre/post bar graphs for each condition
    indivmeans_bar_plot_flag = 0;       % whether to plot mean bar plots
    plot_indiv_statval_brain_flag = 1;  % whether to plot statvals on brain
test_pre_post_dist_allConds_flag = 0;   % whether to plot pre/post bar graphs for all conditions
    allmeans_bar_plot_flag = 0;             % whether to plot mean bar plots
    bar_plot_statval_flag = 0;              % whether to plot statval bar plots
    plot_mean_statval_allelecs_flag = 0;    % whether to plot mean statvals across all elecs
ekg_analysis_flag = 0;              % whether to run EKG analysis
    plot_time_domain_figure_flag = 0; % whether to plot time domain results
    plot_freq_domain_figure_flag = 0; % whether to plot freq domain results
    
h5_file_info = [];          % which h5 file(s) to use
recType = 'clinical';       % 'clinical' or 'TDT'
zscore_flag = 1;            % whether to z-score data (will only do once)
find_bad_data_flag = 1;     % whether to reject bad trials/channels
thresh = 10;                % STD threshold for bad trials
nBadTrialThresh = 20;       % nTrials threshold for bad channels
remove_badTrial_flag = 1;   % whether to remove bad trials
twins = [-10 -5 ; 0 5];         % which time windows to use
testStat = 'ranksum';           % which test to use ('ranksum', 'signrank', 'ttest', 'ttest2')
hr_method = 'pan-tompkin';      % which method to use to detect heart rate

%% LOAD AND CONCATENATE DATA

if ~exist('dat','var')
    dat.data = [];
    dat.cond = [];
    for i = 1:size(params_sets,1)
        fprintf('Loading data for %sHz %smA %sus....\n',num2str(params_sets(i,1)),num2str(params_sets(i,2)),num2str(params_sets(i,3)));
        load([rootdir filesep subj filesep 'VNS' filesep subj '_VNS_ERPs_' num2str(params_sets(i,1)) '_' num2str(params_sets(i,2)) '_' num2str(params_sets(i,3)) '.mat']);
        
        % REMOVE UNDESIRED BLOCKS
        if ~isempty(h5_file_info) & ~exist('remove_bad_blocks_flag','var') & length(unique(VNS_dat.file_info)) > 1
            fprintf('Removing undesired blocks\n');
            VNS_dat.raw_erps(:,:,find(~strcmpi(VNS_dat.file_info,h5_file_info))) = [];
            VNS_dat.file_info(find(~strcmpi(VNS_dat.file_info,h5_file_info))) = [];
            remove_bad_blocks_flag = 1;
        end
        
        % Z-SCORE DATA
        if zscore_flag
            if isfield(VNS_dat,'zscore_flag')
                fprintf('Data have already been z-scored.\n');
            else
                VNS_dat = VNS_zscore_data(VNS_dat);
            end
        end
        
        % REJECT BAD DATA
        if find_bad_data_flag
            if isfield(VNS_dat,'badTrials')
                fprintf('Data have already been rejected.\n');
            else
                VNS_dat = VNS_find_bad_data(VNS_dat,thresh,nBadTrialThresh,h5_file_info,remove_badTrial_flag);
            end
        end
        
        dat.data = cat(3,dat.data,VNS_dat.raw_erps);
        dat.cond = cat(1,dat.cond,repmat(i,size(VNS_dat.raw_erps,3),1));
        dat.time_axis = VNS_dat.time_axis;
        dat.anatomy_elecs = VNS_dat.anatomy_elecs;
        dat.sampFreq = VNS_dat.sampFreq;
        dat.EKG_ch(i) = VNS_dat.EKG_ch;
        
        clear VNS_dat;
        fprintf('\n');
    end
end

conds = unique(dat.cond);
for i = 1:size(params_sets,1)
    condNames{i} = [num2str(params_sets(i,1)) 'Hz ' num2str(params_sets(i,2)) 'mA ' num2str(params_sets(i,3)) 'us'];
end
alpha_level = 0.05/size(dat.data,1);

%% PLOT AVERAGE ERPs

if plot_avg_ERP_flag
    figure('Name','Average ERP');
    textprogressbar('Plotting electrode average ERPs: ');
    clear sp
    
    for i = 1:size(dat.data,1)
        textprogressbar((i/size(dat.data,1))*100);
        p = plotGridPosition_new(i,size(dat.data,1),ceil(sqrt(size(dat.data,1))));
        subplot('Position',p);
        
        for c = 1:length(conds)
            trls = find(dat.cond == conds(c));
            
            %     h(i) = shadedErrorBar(taxis,squeeze(mean(VNS_dat.raw_erps(i,:,trls),3)),...
            %         squeeze(ste(VNS_dat.raw_erps(i,:,trls),3)));
            %     h(i).mainLine.Color = 'b';
            %     h(i).patch.FaceColor = 'b';
            %     h(i).patch.FaceAlpha = 0.5;
            plot(dat.time_axis,squeeze(mean(dat.data(i,:,trls),3)));
            hold on;
            
            axis tight;
            sp(i) = gca;
            
            set(gca,'XTickLabel',[],'YTickLabel',[]);
        end
        
        fprintf('\n');
        
    end
    
    yl = cell2mat(get(sp, 'Ylim'));
    ylnew = [min(yl(:,1)) max(yl(:,2))];
    set(sp, 'Ylim', ylnew);
    for i = 1:length(sp)
        set(gcf,'CurrentAxes',sp(i));
        line([0 0],get(gca,'YLim'),'Color','k');
        line(get(gca,'XLim'),[0 0],'Color','k');
        
        text(0,(max(get(gca,'YLim')) - (max(get(gca,'YLim'))*0.1)),num2str(i));
    end
end

%% PLOT PRE/POST MEANS FOR EACH COND SEPARATELY

if test_pre_post_dist_byCond_flag
    
    if indivmeans_bar_plot_flag
        for c = 1:length(conds)
            trls = find(dat.cond == conds(c));
            
            figure('Name',[num2str(params_sets(c,1)) 'Hz ' num2str(params_sets(c,2)) 'mA ' num2str(params_sets(c,3)) 'us']);
            textprogressbar('Plotting electrode pre/post VNS means: ');
            for i = 1:size(dat.data,1)
                textprogressbar((i/size(dat.data,1))*100);
                p = plotGridPosition_new(i,size(dat.data,1),ceil(sqrt(size(dat.data,1))));
                subplot('Position',p);
                
                if ~strcmpi(dat.anatomy_elecs(i),'Left-Cerebral-White-Matter') && ~strcmpi(dat.anatomy_elecs(i),'')
                    
                    testDist1 = squeeze(dat.data(i,find(dat.time_axis==twins(1)):find(dat.time_axis==twins(2)),trls));
                    testDist2 = squeeze(dat.data(i,find(dat.time_axis==twins(3)):find(dat.time_axis==twins(4)),trls));
                    
                    switch testStat
                        case 'ranksum'
                            [pval(i),~,tmpStatVal(i)] = ranksum(testDist1(:),testDist2(:));
                            statVal(i) = tmpStatVal(i).zval;
                            
                            hbar(i) = barwitherr([ste(testDist1(:)),ste(testDist2(:))],...
                                [mean(testDist1(:)),mean(testDist2(:))]);
                        case 'signrank'
                            [pval(i),~,tmpStatVal(i)] = signrank(mean(testDist1),mean(testDist2));
                            statVal(i) = tmpStatVal(i).zval;
                            
                            hbar(i) = barwitherr([ste(mean(testDist1)),ste(mean(testDist2))],...
                                [mean(mean(testDist1)),mean(mean(testDist2))]);
                        case 'ttest'
                            [~,pval(i),~,tmpStatVal(i)] = ttest(mean(testDist1),mean(testDist2));
                            statVal(i) = tmpStatVal(i).tstat;
                            
                            hbar(i) = barwitherr([ste(mean(testDist1)),ste(mean(testDist2))],...
                                [mean(mean(testDist1)),mean(mean(testDist2))]);
                        case 'ttest2'
                            [~,pval(i),~,tmpStatVal(i)] = ttest2(testDist1(:),testDist2(:));
                            statVal(i) = tmpStatVal(i).tstat;
                            
                            hbar(i) = barwitherr([ste(testDist1(:)),ste(testDist2(:))],...
                                [mean(testDist1(:)),mean(testDist2(:))]);
                    end
                    
                    if pval(i) <= alpha_level
                        set(hbar(i),'FaceColor','r');
                    end
                    
                    axis tight;
                    sp(i) = gca;
                    
                else
                    sp(i) = gca;
                    sp(i).YLim = [0 1e-50];
                    pval(i) = NaN;
                    statVal(i) = NaN;
                end
                set(gca,'XTickLabel',[],'YTickLabel',[]);
            end
            fprintf('\n');
            
            textprogressbar('Rescaling axes');
            yl = cell2mat(get(sp, 'Ylim'));
            ylnew = [min(yl(:,1)) max(yl(:,2))];
            set(sp, 'Ylim', ylnew);
            for i = 1:length(sp)
                textprogressbar((i/length(sp))*100);
                set(gcf,'CurrentAxes',sp(i));
                line(get(gca,'XLim'),[0 0],'Color','k');
                
                text(0.5,(max(get(gca,'YLim')) - (max(get(gca,'YLim'))*0.1)),num2str(i));
            end
            fprintf('\n%d/%d sig chans (%2.2g%%)\n',length(find(pval <= alpha_level)),length(pval),...
                (length(find(pval <= alpha_level))/length(pval))*100);
            fprintf('\n');
        end
    end
    
    if plot_indiv_statval_brain_flag
        if ~exist('brain','var')
            fprintf('Loading brain anatomy and elecs....\n');
            brain.lh = load([brain_dir filesep subj filesep 'Meshes' filesep subj '_lh_pial.mat']);
            brain.rh = load([brain_dir filesep subj filesep 'Meshes' filesep subj '_rh_pial.mat']);
            load([brain_dir filesep subj filesep 'elecs' filesep recType '_elecs_all.mat']);
        end

        for c = 1:length(conds)
            trls = find(dat.cond == conds(c));
            
            figure('Name',[num2str(params_sets(c,1)) 'Hz ' num2str(params_sets(c,2)) 'mA ' num2str(params_sets(c,3)) 'us'],'Color','w');
        
            ctmr_gauss_plot(brain.lh.cortex,[0 0 0],0,'lh');
            hold on;
            ctmr_gauss_plot(brain.rh.cortex,[0 0 0],0,'rh');
            alpha(0.25);
            cmap = cbrewer('div','RdBu',101);
            
            textprogressbar('Plotting electrode pre/post VNS means: ');
            for i = 1:size(dat.data,1)
                textprogressbar((i/size(dat.data,1))*100);
                
                if ~strcmpi(dat.anatomy_elecs(i),'Left-Cerebral-White-Matter') && ~strcmpi(dat.anatomy_elecs(i),'')
                    
                    testDist1 = squeeze(dat.data(i,find(dat.time_axis==twins(1)):find(dat.time_axis==twins(2)),trls));
                    testDist2 = squeeze(dat.data(i,find(dat.time_axis==twins(3)):find(dat.time_axis==twins(4)),trls));
                    
                    switch testStat
                        case 'ranksum'
                            [pval(i),~,tmpStatVal(i)] = ranksum(testDist1(:),testDist2(:));
                            statVal(i) = tmpStatVal(i).zval;
                        case 'signrank'
                            [pval(i),~,tmpStatVal(i)] = signrank(mean(testDist1),mean(testDist2));
                            statVal(i) = tmpStatVal(i).zval;
                        case 'ttest'
                            [~,pval(i),~,tmpStatVal(i)] = ttest(mean(testDist1),mean(testDist2));
                            statVal(i) = tmpStatVal(i).tstat;
                        case 'ttest2'
                            [~,pval(i),~,tmpStatVal(i)] = ttest2(testDist1(:),testDist2(:));
                            statVal(i) = tmpStatVal(i).tstat;
                    end                                        
                else
                    pval(i) = NaN;
                    statVal(i) = NaN;
                end
                
                if ~ismember(i,dat.EKG_ch(c))
                    if pval(i) <= alpha_level
                        plotDat(i) = (statVal(i) - min(statVal)) / (max(statVal) - min(statVal));
                        scatter3(elecmatrix(i,1),elecmatrix(i,2),elecmatrix(i,3),...
                            50,...
                            cmap(ceil(plotDat(i)*100)+1,:),'filled');
                        hold on;
                    else
                        scatter3(elecmatrix(i,1),elecmatrix(i,2),elecmatrix(i,3),...
                            50,...
                            'k');
                        hold on;
                    end
                end
            end
            fprintf('statVal [min max]: [%2.4g %2.4g]\n',min(statVal),max(statVal));
            fprintf('\n');
        end
    end
end

%% PLOT PRE/POST MEANS FOR ALL CONDS ON ONE PLOT

if test_pre_post_dist_allConds_flag
    
    figure('Name','pre/post VNS');
    textprogressbar('Plotting electrode pre/post VNS means: ');
    clear hbar statVal pval;
    
    if allmeans_bar_plot_flag
        for i = 1:size(dat.data,1)
            textprogressbar((i/size(dat.data,1))*100);
            p = plotGridPosition_new(i,size(dat.data,1),ceil(sqrt(size(dat.data,1))));
            subplot('Position',p);
            
            for c = 1:length(conds)
                trls = find(dat.cond == conds(c));
                
                if ~strcmpi(dat.anatomy_elecs(i),'Left-Cerebral-White-Matter') && ~strcmpi(dat.anatomy_elecs(i),'')
                    
                    testDist1 = squeeze(dat.data(i,find(dat.time_axis==twins(1)):find(dat.time_axis==twins(2)),trls));
                    testDist2 = squeeze(dat.data(i,find(dat.time_axis==twins(3)):find(dat.time_axis==twins(4)),trls));
                    
                    switch testStat
                        case 'ranksum'
                            [pval(i,c),~,tmpStatVal(i,c)] = ranksum(testDist1(:),testDist2(:));
                            statVal(i,c) = tmpStatVal(i,c).zval;
                        case 'signrank'
                            [pval(i,c),~,tmpStatVal(i,c)] = signrank(mean(testDist1),mean(testDist2));
                            statVal(i,c) = tmpStatVal(i,c).zval;
                        case 'ttest'
                            [~,pval(i,c),~,tmpStatVal(i,c)] = ttest(mean(testDist1),mean(testDist2));
                            statVal(i,c) = tmpStatVal(i,c).tstat;
                        case 'ttest2'
                            [~,pval(i,c),~,tmpStatVal(i,c)] = ttest2(testDist1(:),testDist2(:));
                            statVal(i,c) = tmpStatVal(i,c).tstat;
                    end
                    
                else
                    pval(i,c) = NaN;
                    statVal(i,c) = NaN;
                end
                
                plotDat.mean(i,c,:) = [mean(testDist1(:)) mean(testDist2(:))];
                plotDat.err(i,c,:) = [ste(testDist1(:)) ste(testDist2(:))];
            end
            
            
            if ~strcmpi(dat.anatomy_elecs(i),'Left-Cerebral-White-Matter') && ~strcmpi(dat.anatomy_elecs(i),'')
                
                hbar(i,:,:) = barwitherr(squeeze(plotDat.err(i,:,:)),...
                    squeeze(plotDat.mean(i,:,:)));
                
                axis tight;
                sp(i) = gca;
                
            else
                sp(i) = gca;
                sp(i).YLim = [0 1e-50];
            end
            
            set(gca,'XTickLabel',[],'YTickLabel',[]);
        end
        fprintf('\n');
        
        textprogressbar('Rescaling axes');
        yl = cell2mat(get(sp, 'Ylim'));
        ylnew = [min(yl(:,1)) max(yl(:,2))];
        set(sp, 'Ylim', ylnew);
        for i = 1:length(sp)
            textprogressbar((i/length(sp))*100);
            set(gcf,'CurrentAxes',sp(i));
            line(get(gca,'XLim'),[0 0],'Color','k');
            
            text(0.5,(max(get(gca,'YLim')) - (max(get(gca,'YLim'))*0.1)),num2str(i));
        end
        fprintf('\n');
    end
    
    %% PLOT STAT VAL FOR EACH CONDITION ON EACH ELECTRODE
    
    if bar_plot_statval_flag
        figure('Name','pre vs post test stat');
        textprogressbar('Plotting electrode pre/post VNS stats: ');
        
        clear b statVal pval
        
        for i = 1:size(dat.data,1)
            textprogressbar((i/size(dat.data,1))*100);
            p = plotGridPosition_new(i,size(dat.data,1),ceil(sqrt(size(dat.data,1))));
            subplot('Position',p);
            
            if ~strcmpi(dat.anatomy_elecs(i),'Left-Cerebral-White-Matter') && ~strcmpi(dat.anatomy_elecs(i),'')
                for c = 1:length(conds)
                    trls = find(dat.cond == conds(c));
                    
                    testDist1 = squeeze(dat.data(i,find(dat.time_axis==twins(1)):find(dat.time_axis==twins(2)),trls));
                    testDist2 = squeeze(dat.data(i,find(dat.time_axis==twins(3)):find(dat.time_axis==twins(4)),trls));
                    
                    switch testStat
                        case 'ranksum'
                            [pval(i,c),~,tmpStatVal(i,c)] = ranksum(testDist1(:),testDist2(:));
                            statVal(i,c) = tmpStatVal(i,c).zval;
                        case 'signrank'
                            [pval(i,c),~,tmpStatVal(i,c)] = signrank(mean(testDist1),mean(testDist2));
                            statVal(i,c) = tmpStatVal(i,c).zval;
                        case 'ttest'
                            [~,pval(i,c),~,tmpStatVal(i,c)] = ttest(mean(testDist1),mean(testDist2));
                            statVal(i,c) = tmpStatVal(i,c).tstat;
                        case 'ttest2'
                            [~,pval(i,c),~,tmpStatVal(i,c)] = ttest2(testDist1(:),testDist2(:));
                            statVal(i,c) = tmpStatVal(i,c).tstat;
                    end
                    b(i,c) = bar(c,statVal(i,c)*-1);
                    if pval(i,c) <= alpha_level
                        set(b(i,c),'FaceColor','r');
                    end
                    hold on;
                end
                axis tight;
                sp(i) = gca;
                
            else
                pval(i,:) = NaN;
                statVal(i,:) = NaN;
                sp(i) = gca;
                sp(i).YLim = [0 1e-50];
            end
            set(gca,'XTickLabel',[],'YTickLabel',[]);
        end
        fprintf('\n');
        
        textprogressbar('Rescaling axes');
        yl = cell2mat(get(sp, 'Ylim'));
        ylnew = [min(yl(:,1)) max(yl(:,2))];
        set(sp, 'Ylim', ylnew);
        for i = 1:length(sp)
            textprogressbar((i/length(sp))*100);
            set(gcf,'CurrentAxes',sp(i));
            line(get(gca,'XLim'),[0 0],'Color','k');
            
            text(0.5,(max(get(gca,'YLim')) - (max(get(gca,'YLim'))*0.1)),num2str(i));
        end
        fprintf('\n');
    end
    
    %% PLOT MEAN STATVALS ACROSS ALL ELECS
    
    if plot_mean_statval_allelecs_flag
        figure;
        subplot(2,1,1);
        barwitherr(nanstd(statVal*-1),nanmean(statVal*-1));
        ylabel('Mean z-value (ranksum)');
        set(gca,'XTickLabel',condNames);
        
        subplot(2,1,2);
        barwitherr(nanstd(abs(statVal)),nanmean(abs(statVal)));
        ylabel('Magnitude of z-value (ranksum)');
        
        set(gca,'XTickLabel',condNames);
        xlabel('VNS parameters');
    end
end

%% EKG ANALYSIS

if ekg_analysis_flag
    if plot_freq_domain_figure_flag
        h_freq_fig = figure('Name','EKG Spectra');
    end
    for c = 1:length(conds)
        trls = find(dat.cond == conds(c));
        [qrs_i_raw_pre{c},qrs_i_raw_post{c},aligned_EKG_pre{c},aligned_EKG_post{c}] = VNS_analyze_HRV(dat,...
            trls,...
            params_sets(c,:),...
            c,...
            [min(dat.time_axis) max(dat.time_axis)],...
            hr_method,...
            testStat,...
            plot_time_domain_figure_flag,...
            plot_freq_domain_figure_flag,...
            h_freq_fig);
    end
end