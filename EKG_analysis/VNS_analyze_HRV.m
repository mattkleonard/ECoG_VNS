function [qrs_i_raw_pre,qrs_i_raw_post,aligned_EKG_pre,aligned_EKG_post] = VNS_analyze_HRV(dat,trls,VNS_parms,cond,ERP_times,hr_method,testStat,plot_time_domain_figure_flag,plot_freq_domain_figure_flag,h_freq_fig)

%% TO DO: ALIGN DATA TO R PEAK


% subj = 'EC131';
%
% rootdir = '/Users/mattleonard/Documents/Research/data';
%
% params_sets = [25 2.25 250];
%
% stim_freq = VNS_parms(1);
%
% hr_method = 'pan-tompkin'; % 'pan-tompkin', 'max'

notch_VNS_flag = 0;

%%

clear qrs_amp_raw qrs_i_raw delay rr_int_pre rr_int_post

textprogressbar('Calculating EKG waveform using Pan-Tompkin ');
for i = 1:length(trls)
    textprogressbar((i/length(trls))*100);
    
    if notch_VNS_flag
        data(:,i) = applyLineNoiseNotch_VNS_Harmonics(squeeze(dat.data(dat.EKG_ch,:,trls(i))),dat.sampFreq,'stim_freq',stim_freq);
    else
        data(:,i) = squeeze(dat.data(dat.EKG_ch(cond),:,trls(i)));
    end
    
    preStimDat(i,:) = data(1:find(dat.time_axis == 0),i);
    postStimDat(i,:) = data((find(dat.time_axis == 0)+1):end,i);
    
    switch hr_method
        case 'pan-tompkin'
            [qrs_amp_raw_pre{i},qrs_i_raw_pre{i},delay_pre(i)] = pan_tompkin(preStimDat(i,:),dat.sampFreq,0);
            [qrs_amp_raw_post{i},qrs_i_raw_post{i},delay_post(i)] = pan_tompkin(postStimDat(i,:),dat.sampFreq,0);
        case 'max'
            [qrs_amp_raw_pre{i},qrs_i_raw_pre{i}] = max(preStimDat(i,:));
            [qrs_amp_raw_post{i},qrs_i_raw_post{i}] = max(postStimDat(i,:));
    end
    
    rr_int_pre{i} = diff(qrs_i_raw_pre{i}) / dat.sampFreq;
    rr_int_post{i} = diff(qrs_i_raw_post{i}) / dat.sampFreq;
    %     [qrs_amp_raw{i},qrs_i_raw{i},delay(i)] = pan_tompkin(squeeze(VNS_dat.raw_erps(VNS_dat.EKG_ch,:,i)),VNS_dat.sampFreq,0);
    
    %     rr_int{i,1}
end
fprintf('\n');
aligned_EKG_pre = NaN(size((qrs_i_raw_pre{1}(1)) - (dat.sampFreq/2):(qrs_i_raw_pre{1}(1)) + (dat.sampFreq/2)));
aligned_EKG_post = NaN(size((qrs_i_raw_post{1}(1)) - (dat.sampFreq/2):(qrs_i_raw_post{1}(1)) + (dat.sampFreq/2)));

if plot_time_domain_figure_flag
    figure('Name',[num2str(VNS_parms(1)) 'Hz ' num2str(VNS_parms(2)) 'mA ' num2str(VNS_parms(3)) 'us'],...
        'Position',[11    23   944   930]);
    
    %% ALIGNED TO R-PEAK OF EKG
    
    k = 1;
    for i = 1:length(qrs_i_raw_pre)
        for j = 1:length(qrs_i_raw_pre{i})
            try
                aligned_EKG_pre(k,:) = preStimDat(i,(qrs_i_raw_pre{i}(j)) - (dat.sampFreq/2):(qrs_i_raw_pre{i}(j)) + (dat.sampFreq/2));
                aligned_EKG_post(k,:) = postStimDat(i,(qrs_i_raw_post{i}(j)) - (dat.sampFreq/2):(qrs_i_raw_post{i}(j)) + (dat.sampFreq/2));
            catch
                aligned_EKG_pre(k,:) = NaN;
                aligned_EKG_post(k,:) = NaN;
            end
            k = k + 1;
        end
    end
    
    subplot(5,2,2);
    plot(-0.5:1/dat.sampFreq:0.5,nanmean(aligned_EKG_pre,1),'Color','b');
    hold on;
    plot(-0.5:1/dat.sampFreq:0.5,nanmean(aligned_EKG_post,1),'Color','r');
    line([0 0],get(gca,'YLim'),'Color','k');
    legend({'pre','post'});
    title('R-peak aligned average');
    
    %% HISTOGRAMS
    
    nBins = 50;
    subplot(5,2,1);
    histogram([rr_int_pre{:}],nBins,'EdgeColor','none','FaceColor','b');
    hold on;
    histogram([rr_int_post{:}],nBins,'EdgeColor','none','FaceColor','r');
    xlabel('R-R interval');
    ylabel('count');
    title('R-R interval distributions');
    
    %% SINGLE TRIALS WITH R-PEAKS
    
    subplot(5,2,3);
    imagesc(preStimDat);
    hold on;
    for i = 1:length(qrs_i_raw_pre)
        plot(qrs_i_raw_pre{i},i,'or');
    end
    xlabel('time');
    ylabel('trial');
    title('pre-stim');
    subplot(5,2,4);
    for i = 1:length(qrs_i_raw_pre)
        plot(preStimDat(i,:),'b');
        hold on;
        scatter(qrs_i_raw_pre{i},preStimDat(i,qrs_i_raw_pre{i}),'r');
    end
    xlabel('time');
    ylabel('amplitude');
    title('pre-stim');
    
    subplot(5,2,5);
    imagesc(postStimDat);
    hold on;
    for i = 1:length(qrs_i_raw_post)
        plot(qrs_i_raw_post{i},i,'or');
    end
    xlabel('time');
    ylabel('trial');
    title('post-stim');
    subplot(5,2,6);
    for i = 1:length(qrs_i_raw_post)
        plot(postStimDat(i,:),'b');
        hold on;
        scatter(qrs_i_raw_post{i},postStimDat(i,qrs_i_raw_post{i}),'r');
    end
    xlabel('time');
    ylabel('amplitude');
    title('post-stim');
    
    %% R-R PEAK DIFFERENCE STATISTICS
    
    subplot(5,2,7);
    plot([rr_int_pre{:}],'bo');
    hold on;
    plot([rr_int_post{:}],'ro');
    xlabel('trial');
    ylabel('R-R interval');
    
    subplot(5,2,9);
    plot(movmean([rr_int_pre{:}],5),'b');
    hold on;
    plot(movmean([rr_int_post{:}],5),'r');
    xlabel('trial');
    ylabel('R-R interval');
    title('Moving average window');
    
    subplot(5,2,8);
    barwitherr([std(cellfun(@mean,rr_int_pre)) std(cellfun(@mean,rr_int_post))],...
        [mean(cellfun(@mean,rr_int_pre)) mean(cellfun(@mean,rr_int_post))]);
    switch testStat
        case 'ttest2'
            [h,p,ci,stats] = ttest2(cellfun(@mean,rr_int_pre),cellfun(@mean,rr_int_post));
            title(sprintf('Mean: t(%d)=%2.2g, p=%2.2g',stats.df,stats.tstat,p));
        case 'ranksum'
            [p,~,stats] = ranksum(cellfun(@mean,rr_int_pre),cellfun(@mean,rr_int_post));
            title(sprintf('Mean: z=%2.2g, p=%2.2g',stats.zval,p));
    end
    set(gca,'XTickLabel',{'pre','post'});
    
    subplot(5,2,10);
    barwitherr([std(cellfun(@std,rr_int_pre)) std(cellfun(@std,rr_int_post))],...
        [mean(cellfun(@std,rr_int_pre)) mean(cellfun(@std,rr_int_post))]);
    switch testStat
        case 'ttest2'
            [h,p,ci,stats] = ttest2(cellfun(@std,rr_int_pre),cellfun(@std,rr_int_post));
            title(sprintf('STD: t(%d)=%2.2g, p=%2.2g',stats.df,stats.tstat,p));
        case 'ranksum'
            [p,~,stats] = ranksum(cellfun(@std,rr_int_pre),cellfun(@std,rr_int_post));
            title(sprintf('Mean: z=%2.2g, p=%2.2g',stats.zval,p));
    end
    set(gca,'XTickLabel',{'pre','post'});
    
end

%%

if plot_freq_domain_figure_flag
    lf_range = [0.04 0.15];
    hf_range = [0.15 0.4];
    smth_const = 0.2;
    
    [pxx_pre,f_pre] = plomb([rr_int_pre{:}],1:length([rr_int_pre{:}]));
    [pxx_post,f_post] = plomb([rr_int_post{:}],1:length([rr_int_post{:}]));
    
    figure(h_freq_fig);
    subplot(1,length(unique(dat.cond)),cond);
    
    pxx_pre_line = semilogy(f_pre,smooth(pxx_pre,round(smth_const*length(pxx_pre))),'Color','b');
    hold on;
    pxx_post_line = semilogy(f_post,smooth(pxx_post,round(smth_const*length(pxx_post))),'Color','r');
    set(gca,'XLim',[0 max(hf_range)]);
    
    xlabel('frequency');
    ylabel('power');
    
    patch([lf_range(1) lf_range(1) lf_range(2) lf_range(2)],...
        [sort(get(gca,'YLim')), sort(get(gca,'YLim'),'descend')],...
        [0.8 0.9 0.92], 'FaceAlpha',1,'EdgeColor','none');
    patch([hf_range(1) hf_range(1) hf_range(2) hf_range(2)],...
        [sort(get(gca,'YLim')), sort(get(gca,'YLim'),'descend')],...
        [0.2 0.6 0.5], 'FaceAlpha',1,'EdgeColor','none');
    
    uistack(pxx_pre_line,'top');
    uistack(pxx_post_line,'top');
    legend({'lf','hf','pre','post'});
    
    title([num2str(VNS_parms(1)) 'Hz ' num2str(VNS_parms(2)) 'mA ' num2str(VNS_parms(3)) 'us']);
end
