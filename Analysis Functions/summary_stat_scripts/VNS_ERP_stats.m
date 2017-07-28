function [pval,statVal] = VNS_ERP_stats(VNS_dat,subj,brain_dir,varargin)

%% INPUT ARGS

% which electrode montage type
if ~isempty(find(strcmpi(varargin,'recType')));
    recType = varargin{find(strcmpi(varargin,'recType'))+1};
else
    recType = 'clinical';
end

% whether to test distributions of pre vs. post VNS data
if ~isempty(find(strcmpi(varargin,'test_pre_post_dist_flag')));
    test_pre_post_dist_flag = varargin{find(strcmpi(varargin,'test_pre_post_dist_flag'))+1};
else
    test_pre_post_dist_flag = 0;
end

% whether to test specific timepoints of ERP
if ~isempty(find(strcmpi(varargin,'test_ERP_timecourse_flag')));
    test_ERP_timecourse_flag = varargin{find(strcmpi(varargin,'test_ERP_timecourse_flag'))+1};
else
    test_ERP_timecourse_flag = 0;
end

% whether to plot brain
if ~isempty(find(strcmpi(varargin,'plot_brain_flag')));
    plot_brain_flag = varargin{find(strcmpi(varargin,'plot_brain_flag'))+1};
else
    plot_brain_flag = 0;
end

% which test to use ('ranksum', 'signrank', 'ttest', 'ttest2')
if ~isempty(find(strcmpi(varargin,'testStat')));
    testStat = varargin{find(strcmpi(varargin,'testStat'))+1};
else
    testStat = 'ranksum'; % 'ranksum', 'signrank', 'ttest', 'ttest2'
end

% which time windows to use
if ~isempty(find(strcmpi(varargin,'twins')));
    twins = varargin{find(strcmpi(varargin,'twins'))+1};
else
    twins = [-10 -5 ; 0 5];
end

% which alpha level to use for significance
if ~isempty(find(strcmpi(varargin,'alpha_level')));
    alpha_level = varargin{find(strcmpi(varargin,'alpha_level'))+1};
else
    alpha_level = 0.05/size(VNS_dat.raw_erps,1);
end

%% TEST PRE/POST STIM DISTRIBUTIONS

if test_pre_post_dist_flag
    clear pval tmpStatVal statVal
    
    figure;
    textprogressbar(['Running ' testStat ' test ']);
    for i = 1:size(VNS_dat.raw_erps,1)
        textprogressbar((i/size(VNS_dat.raw_erps,1))*100);
        p = plotGridPosition_new(i,size(VNS_dat.raw_erps,1),ceil(sqrt(size(VNS_dat.raw_erps,1))));
        subplot('Position',p);
        
        if ~strcmpi(VNS_dat.anatomy_elecs(i),'Left-Cerebral-White-Matter') && ~strcmpi(VNS_dat.anatomy_elecs(i),'')
            
            % test stats
            testDist1 = squeeze(VNS_dat.raw_erps(i,find(VNS_dat.time_axis==twins(1)):find(VNS_dat.time_axis==twins(2)),:));
            testDist2 = squeeze(VNS_dat.raw_erps(i,find(VNS_dat.time_axis==twins(3)):find(VNS_dat.time_axis==twins(4)),:));
            
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
    fprintf('\n');
    
    %% PLOT ELECS ON BRAIN
    
    if plot_brain_flag
        if ~exist('brain','var')
            fprintf('Loading brain anatomy and elecs....\n');
            brain.lh = load([brain_dir filesep subj filesep 'Meshes' filesep subj '_lh_pial.mat']);
            brain.rh = load([brain_dir filesep subj filesep 'Meshes' filesep subj '_rh_pial.mat']);
            load([brain_dir filesep subj filesep 'elecs' filesep recType '_elecs_all.mat']);
        end
        
        fprintf('Plotting recons and elec data....\n');
        figure('Color','w');;
        ctmr_gauss_plot(brain.lh.cortex,[0 0 0],0,'lh');
        hold on;
        ctmr_gauss_plot(brain.rh.cortex,[0 0 0],0,'rh');
        alpha(0.25);
        
        switch testStat
            case 'ranksum'
                cmap = cbrewer('div','RdBu',101);
            case 'signrank'
                cmap = cbrewer('div','RdBu',101);
            case 'ttest'
                cmap = cbrewer('div','RdBu',101);
            case 'ttest2'
                cmap = cbrewer('div','RdBu',101);
        end
        
        for i = 1:length(statVal)
            if ~ismember(i,VNS_dat.EKG_ch)
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
    end
end

%% TEST ERP TIMECOURSE 

if test_ERP_timecourse_flag
    clear pval tmpStatVal statVal

    figure;
    textprogressbar('Plotting ERPs ');
    for i = 1:size(VNS_dat.raw_erps,1)
        textprogressbar((i/size(VNS_dat.raw_erps,1))*100);
        p = plotGridPosition_new(i,size(VNS_dat.raw_erps,1),ceil(sqrt(size(VNS_dat.raw_erps,1))));
        subplot('Position',p);
        
        if ~strcmpi(VNS_dat.anatomy_elecs(i),'Left-Cerebral-White-Matter') && ~strcmpi(VNS_dat.anatomy_elecs(i),'')
            plot(VNS_dat.time_axis,...
                squeeze(mean(VNS_dat.raw_erps(i,:,:),3)));
            hold on;
            
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
        line([0 0],get(gca,'YLim'),'Color','k');
        
        text(0.5,(max(get(gca,'YLim')) - (max(get(gca,'YLim'))*0.1)),num2str(i));
    end
    fprintf('\n');

    textprogressbar(['Running ' testStat ' test']);
    for i = 1:size(VNS_dat.raw_erps,1)
        textprogressbar((i/size(VNS_dat.raw_erps,1))*100);
        p = plotGridPosition_new(i,size(VNS_dat.raw_erps,1),ceil(sqrt(size(VNS_dat.raw_erps,1))));
        subplot('Position',p);
        
        if ~strcmpi(VNS_dat.anatomy_elecs(i),'Left-Cerebral-White-Matter') && ~strcmpi(VNS_dat.anatomy_elecs(i),'')
            for j = 1:100:size(VNS_dat.raw_erps,2)
                
                % test stats
                baselineDist = squeeze(VNS_dat.raw_erps(i,find(VNS_dat.time_axis==twins(1)):find(VNS_dat.time_axis==twins(2)),:));
                
                switch testStat
                    case 'ranksum'
                        [pval(i,j),~,tmpStatVal(i,j)] = ranksum(testDist1(:),testDist2(:));
                        statVal(i,j) = tmpStatVal(i,j).zval;
                        
                    case 'signrank'
                        [pval(i,j),~,tmpStatVal(i,j)] = signrank(squeeze(VNS_dat.raw_erps(i,j,:)),mean(mean(baselineDist)),'method','approximate');
                        statVal(i,j) = tmpStatVal(i,j).zval;
                        
                    case 'ttest'
                        [~,pval(i,j),~,tmpStatVal(i,j)] = ttest(squeeze(VNS_dat.raw_erps(i,j,:)),mean(mean(baselineDist)));
                        statVal(i,j) = tmpStatVal(i,j).tstat;
                        
                    case 'ttest2'
                        [~,pval(i,j),~,tmpStatVal(i,j)] = ttest2(testDist1(:),testDist2(:));
                        statVal(i,j) = tmpStatVal(i,j).tstat;
                        
                end
                if pval(i,j) <= alpha_level
                    line([VNS_dat.time_axis(j) VNS_dat.time_axis(j)],get(gca,'YLim'),'Color','r');
                end
            end
            
        else
            pval(i) = NaN;
            statVal(i) = NaN;
        end
    end
    fprintf('\n');
end