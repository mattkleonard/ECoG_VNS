rootdir = '/Users/mattleonard/Documents/Research/data';
fsDs = 200;

%%
if ~exist('VNS_raw','var')
    load([rootdir '/EC131/VNS/EC131_QuietTime/VNS_notched_data.mat']);
    
    
    trls = 1:100;
    
    [p,q] = rat(fsDs/VNS_raw.sampFreq);
    
    taxis = -10:1/fsDs:10;
    
    for i = 1:size(VNS_raw.raw_erps,1)
        VNS_raw.raw_erps_ds(i,:,:) = resample(squeeze(VNS_raw.raw_erps(i,:,:)),p,q);
    end
end

%%

figure;
for i = 1:size(VNS_raw.raw_erps_ds,1)
    p = plotGridPosition_new(i,size(VNS_raw.raw_erps_ds,1),ceil(sqrt(size(VNS_raw.raw_erps_ds,1))));
    subplot('Position',p);
    
    %     h(i) = shadedErrorBar(taxis,squeeze(mean(VNS_raw.raw_erps_ds(i,:,trls),3)),...
    %         squeeze(ste(VNS_raw.raw_erps_ds(i,:,trls),3)));
    %     h(i).mainLine.Color = 'b';
    %     h(i).patch.FaceColor = 'b';
    %     h(i).patch.FaceAlpha = 0.5;
    plot(taxis,squeeze(mean(VNS_raw.raw_erps_ds(i,:,trls),3)),'Color','b');
    hold on;
    
    axis tight;
    sp(i) = gca;
    
    set(gca,'XTickLabel',[],'YTickLabel',[]);
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

%%

elec = 75;

figure;
subplot(2,1,1);
plot(taxis,squeeze(mean(VNS_raw.raw_erps_ds(elec,:,trls),3)),'Color','b');
hold on;
axis tight;
line([0 0],get(gca,'YLim'),'Color','k');
line(get(gca,'XLim'),[0 0],'Color','k');

text(0,(max(get(gca,'YLim')) - (max(get(gca,'YLim'))*0.1)),num2str(elec));

subplot(2,1,2);
imagesc(squeeze(VNS_raw.raw_erps_ds(elec,:,trls))');
xticklabels = -10:10;
xticks = linspace(1,size(VNS_raw.raw_erps_ds,2), numel(xticklabels));
set(gca,'XTick',xticks,'XTickLabel',xticklabels);
axis xy;
hold on;
line([2000 2000],get(gca,'YLim'),'Color','k');

%% TEST PRE/POST STIM DISTRIBUTIONS

twins = [-10 0 ; 0 10];
alpha_level = 0.05/size(VNS_raw.raw_erps_ds,1);

figure;
for i = 1:size(VNS_raw.raw_erps_ds,1)
        p = plotGridPosition_new(i,size(VNS_raw.raw_erps_ds,1),ceil(sqrt(size(VNS_raw.raw_erps_ds,1))));
        subplot('Position',p);
        
        if ~strcmpi(VNS_raw.anatomy_elecs(i),'Left-Cerebral-White-Matter') && ~strcmpi(VNS_raw.anatomy_elecs(i),'')

        % test stats
        testDist1 = squeeze(VNS_raw.raw_erps_ds(i,find(taxis==twins(1)):find(taxis==twins(2)),:));
        testDist2 = squeeze(VNS_raw.raw_erps_ds(i,find(taxis==twins(3)):find(taxis==twins(4)),:));
        pval(i) = ranksum(testDist1(:),testDist2(:));
        %     pval(i) = signrank(testDist1(:),testDist2(:));
        %     [~,pval(i)] = ttest2(testDist1(:),testDist2(:));
        
        hbar(i) = barwitherr([ste(testDist1(:)),ste(testDist2(:))],...
            [mean(testDist1(:)),mean(testDist2(:))]);
        
        if pval(i) <= alpha_level
            set(hbar(i),'FaceColor','r');
        end
        axis tight;
        sp(i) = gca;
        
        set(gca,'XTickLabel',[],'YTickLabel',[]);
    end
end
yl = cell2mat(get(sp, 'Ylim'));
ylnew = [min(yl(:,1)) max(yl(:,2))];
set(sp, 'Ylim', ylnew);
for i = 1:length(sp)
    set(gcf,'CurrentAxes',sp(i));
    line(get(gca,'XLim'),[0 0],'Color','k');
    
    text(0.5,(max(get(gca,'YLim')) - (max(get(gca,'YLim'))*0.1)),num2str(i));
end

%%

figure;
for i = 1:size(VNS_raw.raw_erps_ds,1)
    p = plotGridPosition_new(i,size(VNS_raw.raw_erps_ds,1),ceil(sqrt(size(VNS_raw.raw_erps_ds,1))));
    subplot('Position',p);
    
    testDist1 = squeeze(VNS_raw.raw_erps_ds(i,find(taxis==twins(1)):find(taxis==twins(2)),:));
    testDist2 = squeeze(VNS_raw.raw_erps_ds(i,find(taxis==twins(3)):find(taxis==twins(4)),:));

    [pxx,f_axis] = pwelch(mean(testDist1,2),[],[],[],fsDs);
    plot(f_axis(find(f_axis <= fsDs)),pxx(find(f_axis <= fsDs)));
    hold on;
    [pxx,f_axis] = pwelch(mean(testDist2,2),[],[],[],fsDs);
    plot(f_axis(find(f_axis <= fsDs)),pxx(find(f_axis <= fsDs)));
    
    axis tight;
    
end