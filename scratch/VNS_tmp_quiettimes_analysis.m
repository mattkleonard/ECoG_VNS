


%% GET TRIALS AND DOWNSAMPLE IF NECESSARY

if isempty(trls)
    trls = 1:size(VNS_dat.raw_erps,3);
end

if ds_flag
    [p,q] = rat(fsDs/VNS_dat.sampFreq);
    
    VNS_dat.time_axis = VNS_dat.time_axis(1):1/fsDs:VNS_dat.time_axis(end);
    
    for i = 1:size(VNS_dat.raw_erps,1)
        VNS_dat.raw_erps(i,:,:) = resample(squeeze(VNS_dat.raw_erps(i,:,:)),p,q);
    end
end

%% PLOT AVERAGE ERPs

if plot_avg_ERP_flag
    figure('Name','Average ERP');
    textprogressbar('Plotting electrode average ERPs: ');
    clear sp
    
    for i = 1:size(VNS_dat.raw_erps,1)
        textprogressbar((i/size(VNS_dat.raw_erps,1))*100);
        p = plotGridPosition_new(i,size(VNS_dat.raw_erps,1),ceil(sqrt(size(VNS_dat.raw_erps,1))));
        subplot('Position',p);
        
        %     h(i) = shadedErrorBar(taxis,squeeze(mean(VNS_dat.raw_erps(i,:,trls),3)),...
        %         squeeze(ste(VNS_dat.raw_erps(i,:,trls),3)));
        %     h(i).mainLine.Color = 'b';
        %     h(i).patch.FaceColor = 'b';
        %     h(i).patch.FaceAlpha = 0.5;
        plot(VNS_dat.time_axis,squeeze(mean(VNS_dat.raw_erps(i,:,trls),3)),'Color','b');
        hold on;
        
        axis tight;
        sp(i) = gca;
        
        set(gca,'XTickLabel',[],'YTickLabel',[]);
    end
    
    fprintf('\n');

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

%% PLOT SINGLE TRIAL RASTERS

if plot_raster_ERP_flag
    figure('Name','Single Trial ERP');
    textprogressbar('Plotting electrode single trial raster ERPs: ');
    clear sp
    
    for i = 1:size(VNS_dat.raw_erps,1)
        textprogressbar((i/size(VNS_dat.raw_erps,1))*100);
        p = plotGridPosition_new(i,size(VNS_dat.raw_erps,1),ceil(sqrt(size(VNS_dat.raw_erps,1))));
        subplot('Position',p);
        
        imagesc(squeeze(VNS_dat.raw_erps(i,:,trls))');
        hold on;
        
        axis tight;
        sp(i) = gca;
        
        set(gca,'XTickLabel',[],'YTickLabel',[]);
    end
    
    fprintf('\n');

    cl = cell2mat(get(sp,'CLim'));
    cl(VNS_dat.badChans,:) = NaN;
    clnew = [min(cl(:,1)) max(cl(:,2))];
    set(sp, 'CLim', clnew);
    textprogressbar('Rescaling axes: ');
    for i = 1:length(sp)
        textprogressbar((i/length(sp))*100);
        set(gcf,'CurrentAxes',sp(i));
        line([find(VNS_dat.time_axis == 0) find(VNS_dat.time_axis == 0)],get(gca,'YLim'),'Color','k');
        
        text(find(VNS_dat.time_axis == 0),(max(get(gca,'YLim')) - (max(get(gca,'YLim'))*0.1)),num2str(i));
    end
end

%%

if plot_avg_ERP_region_flag
    disp(unique(VNS_dat.anatomy_elecs));
    anat_region = input('Which region do you want to plot (see list above)? ','s');
    figure('Name',['Average ERP: ' anat_region]);
    textprogressbar('Plotting electrode average ERPs: ');
    
    clear sp
    elecs = find(strcmpi(anat_region,VNS_dat.anatomy_elecs));
    for i = 1:length(elecs)
        textprogressbar((i/length(elecs))*100);
        p = plotGridPosition_new(i,length(elecs),ceil(sqrt(length(elecs))));
        subplot('Position',p);
        
            h(i) = shadedErrorBar(VNS_dat.time_axis,squeeze(mean(VNS_dat.raw_erps(elecs(i),:,trls),3)),...
                squeeze(ste(VNS_dat.raw_erps(elecs(i),:,trls),3)));
            h(i).mainLine.Color = 'b';
            h(i).patch.FaceColor = 'b';
            h(i).patch.FaceAlpha = 0.5;
%         plot(VNS_dat.time_axis,squeeze(mean(VNS_dat.raw_erps(elecs(i),:,trls),3)),'Color','b');
        hold on;
        
        axis tight;
        sp(i) = gca;
        
        set(gca,'XTickLabel',[],'YTickLabel',[]);
    end
    
    fprintf('\n');

    yl = cell2mat(get(sp, 'Ylim'));
    ylnew = [min(yl(:,1)) max(yl(:,2))];
    set(sp, 'Ylim', ylnew);
    for i = 1:length(sp)
        set(gcf,'CurrentAxes',sp(i));
        line([0 0],get(gca,'YLim'),'Color','k');
        line(get(gca,'XLim'),[0 0],'Color','k');
        
        text(0,(max(get(gca,'YLim')) - (max(get(gca,'YLim'))*0.1)),num2str(elecs(i)));
    end
end

%% PLOT AVERAGE ERP AND SINGLE TRIAL RASTER FOR ONE CHANNEL

if plot_single_chan_flag
    figure;
    subplot(2,1,1);
    plot(VNS_dat.time_axis,squeeze(mean(VNS_dat.raw_erps(elec,:,trls),3)),'Color','b');
    hold on;
    axis tight;
    line([0 0],get(gca,'YLim'),'Color','k');
    line(get(gca,'XLim'),[0 0],'Color','k');
    
    text(0,(max(get(gca,'YLim')) - (max(get(gca,'YLim'))*0.1)),num2str(elec));
    ylabel('z-score');
    
    subplot(2,1,2);
    imagesc(squeeze(VNS_dat.raw_erps(elec,:,trls))');
    xticklabels = -10:10;
    xticks = linspace(1,size(VNS_dat.raw_erps,2), numel(xticklabels));
    set(gca,'XTick',xticks,'XTickLabel',xticklabels);
    axis xy;
    hold on;
    line([find(VNS_dat.time_axis == 0) find(VNS_dat.time_axis == 0)],get(gca,'YLim'),'Color','k');
    ylabel('trials');
    xlabel('Time (sec)');
end

%% TEST PRE/POST STIM DISTRIBUTIONS
% 
% twins = [-10 0 ; 0 10];
% alpha_level = 0.05/size(VNS_dat.raw_erps,1);
% 
% figure;
% for i = 1:size(VNS_dat.raw_erps,1)
%         p = plotGridPosition_new(i,size(VNS_dat.raw_erps,1),ceil(sqrt(size(VNS_dat.raw_erps,1))));
%         subplot('Position',p);
%         
%         if ~strcmpi(VNS_dat.anatomy_elecs(i),'Left-Cerebral-White-Matter') && ~strcmpi(VNS_dat.anatomy_elecs(i),'')
% 
%         % test stats
%         testDist1 = squeeze(VNS_dat.raw_erps(i,find(VNS_dat.time_axis==twins(1)):find(VNS_dat.time_axis==twins(2)),:));
%         testDist2 = squeeze(VNS_dat.raw_erps(i,find(VNS_dat.time_axis==twins(3)):find(VNS_dat.time_axis==twins(4)),:));
%         pval(i) = ranksum(testDist1(:),testDist2(:));
%         %     pval(i) = signrank(testDist1(:),testDist2(:));
%         %     [~,pval(i)] = ttest2(testDist1(:),testDist2(:));
%         
%         hbar(i) = barwitherr([ste(testDist1(:)),ste(testDist2(:))],...
%             [mean(testDist1(:)),mean(testDist2(:))]);
%         
%         if pval(i) <= alpha_level
%             set(hbar(i),'FaceColor','r');
%         end
%         axis tight;
%         sp(i) = gca;
%         
%         set(gca,'XTickLabel',[],'YTickLabel',[]);
%     end
% end
% yl = cell2mat(get(sp, 'Ylim'));
% ylnew = [min(yl(:,1)) max(yl(:,2))];
% set(sp, 'Ylim', ylnew);
% for i = 1:length(sp)
%     set(gcf,'CurrentAxes',sp(i));
%     line(get(gca,'XLim'),[0 0],'Color','k');
%     
%     text(0.5,(max(get(gca,'YLim')) - (max(get(gca,'YLim'))*0.1)),num2str(i));
% end
% 
% %%
% 
% figure;
% for i = 1:size(VNS_dat.raw_erps,1)
%     p = plotGridPosition_new(i,size(VNS_dat.raw_erps,1),ceil(sqrt(size(VNS_dat.raw_erps,1))));
%     subplot('Position',p);
%     
%     testDist1 = squeeze(VNS_dat.raw_erps(i,find(VNS_dat.time_axis==twins(1)):find(VNS_dat.time_axis==twins(2)),:));
%     testDist2 = squeeze(VNS_dat.raw_erps(i,find(VNS_dat.time_axis==twins(3)):find(VNS_dat.time_axis==twins(4)),:));
% 
%     [pxx,f_axis] = pwelch(mean(testDist1,2),[],[],[],fsDs);
%     plot(f_axis(find(f_axis <= fsDs)),pxx(find(f_axis <= fsDs)));
%     hold on;
%     [pxx,f_axis] = pwelch(mean(testDist2,2),[],[],[],fsDs);
%     plot(f_axis(find(f_axis <= fsDs)),pxx(find(f_axis <= fsDs)));
%     
%     axis tight;
%     
% end