function [] = analyze_heartRate(VNS_dat)
% This function takes a VNS data structure with the heart rate field
% and uses this and the VNS onset timing information to asses changes in 
% various metrics derived from heart rate measurement.
heart_rate = VNS_dat.heart_rate;
onset_inds = VNS_dat.stim_onsets_inds;
analysis_win = 30;
fs = VNS_dat.sampFreq;

rr_metric = true; % if true use RR, else use BPM
if ~rr_metric
    heart_rate = 60./heart_rate;
end
% clear edge:
onset_inds((onset_inds+analysis_win*fs>length(heart_rate))...
    | (onset_inds-analysis_win*fs < 1)) = [];

%% Get metrics for blocks of VNS ON
hr_std_on = zeros(length(onset_inds),1);
hr_rms_on = zeros(length(onset_inds),1);
hr_mean_on = zeros(length(onset_inds),1);
for k = 1:length(onset_inds)
    heart_blk = heart_rate(onset_inds(k):round(onset_inds(k)+analysis_win*fs));
    %heart_blk(heart_blk < 0.33) = [];
    
    hr_std_on(k) = nanstd(heart_blk);
    hr_rms_on(k) = rms(heart_blk);
    hr_mean_on(k) = nanmean(heart_blk);
    
end

%% Get metrics for blocks of VNS OFF
hr_std_off = zeros(length(onset_inds),1);
hr_rms_off = zeros(length(onset_inds),1);
hr_mean_off = zeros(length(onset_inds),1);
for k = 1:length(onset_inds)
    heart_blk = heart_rate(round(onset_inds(k)-analysis_win*fs):onset_inds(k));
    % heart_blk(heart_blk < 0.33) = [];
    
    hr_std_off(k) = nanstd(heart_blk);
    hr_rms_off(k) = rms(heart_blk);
    hr_mean_off(k) = nanmean(heart_blk);
end

% Compare Data:
% KS_test on RR Interval
figure;
subplot(1,3,1)
histogram(hr_mean_on); hold on; histogram(hr_mean_off)
[h,p] = kstest2(hr_mean_off, hr_mean_on);
legend('On', 'Off');
if rr_metric
    title(['Distribution of Mean R-R Interval; KS pval = ' num2str(p,3)]); xlabel('Mean R-R Interval (s)')
else
    title(['Distribution of Mean Heart Rate; KS pval = ' num2str(p,3)]); xlabel('Mean Heart Rate (BPM)')
end

subplot(1,3,2)
histogram(hr_std_on); hold on; histogram(hr_std_off)
[h,p] = kstest2(hr_std_off, hr_std_on);
legend('On', 'Off');
if rr_metric
    title(['Distribution of STD R-R Interval; KS pval = ' num2str(p,3)]); xlabel('Std R-R Interval (s)')
else
    title(['Distribution of STD Heart Rate; KS pval = ' num2str(p,3)]); xlabel('Std Heart Rate (BPM)')
end

subplot(1,3,3)
histogram(hr_rms_on); hold on; histogram(hr_rms_off)
[h,p] = kstest2(hr_rms_off, hr_rms_on);
legend('On', 'Off');
if rr_metric
    title(['Distribution of RMS R-R Interval; KS pval = ' num2str(p,3)]); xlabel('RMS R-R Interval (s)');
else
    title(['Distribution of RMS Heart Rate; KS pval = ' num2str(p,3)]); xlabel('RMS Heart Rate (BMP)');
end


%% Time Course
figure; plot(hr_std_on); hold on; plot(hr_std_off)
xlabel('Stimulation Block Number'); ylabel('StDev');
legend('On','Off')
title(['Time Series of Heart Rate Variability (RRSTD) ; r = ' num2str(corr(hr_std_on, hr_std_off),3)])


%% Raw Data:
is_stim_on = VNS_dat.is_stim_on;
[h,p] = kstest2(heart_rate(is_stim_on), heart_rate(~is_stim_on));

figure;
histogram(heart_rate(is_stim_on), 'EdgeColor', 'none','Normalization','pdf'); hold on; 
histogram(heart_rate(~is_stim_on),'EdgeColor', 'none','Normalization','pdf')
xlabel('R-R Interval (s)')
ylabel('Normalized Counts')
legend('On', 'Off')
title('Comparison of the Distribution of RR Interval with and without VNS')

std_win = 300;
heart_rate_stdwin = zeros(1,floor(length(heart_rate)/std_win));
is_stim_on_win = false(1,floor(length(heart_rate)/std_win));

for k = 1:floor(length(heart_rate)/std_win)
    heart_rate_stdwin(k) = std(heart_rate((std_win*(k-1)+1):(std_win*k)));
    is_stim_on_win(k) = is_stim_on(k*std_win-floor(std_win/2));
end

figure;
histogram(heart_rate_stdwin(is_stim_on_win),'EdgeColor', 'none','Normalization','pdf'); hold on; 
histogram(heart_rate_stdwin(~is_stim_on_win),'EdgeColor', 'none','Normalization','pdf'); hold on; 
xlabel('HRV (RR-STD)')
ylabel('Normalized Counts')
legend('On', 'Off')
[~,p] = kstest2(heart_rate_stdwin(is_stim_on_win),heart_rate_stdwin(~is_stim_on_win))
title({'Comparison of the Distribution of Heart Rate Variability (RR-STD) with VNS';['KS-Test p-value: p = ' num2str(p)]})



%% Tachogram spectra:
tachogram_x = VNS_dat.tachogram_x;
tachogram_y = VNS_dat.tachogram_y;

% for each block 
lf = [0.04 0.15]; % low frequencies of Tachogram
hf = [0.15 0.4];  % high frequencies of tachogram
max_f = max([hf lf]*1.25);

lh_ratio_on = zeros(length(onset_inds),1); % low to high power ratio, vns on
lh_ratio_off = zeros(length(onset_inds),1); % low:high, vns off


for k = 1:length(onset_inds)
    % ON:
    is_on_blk = (tachogram_x >= (onset_inds(k)/fs)) & (tachogram_x <= (analysis_win + onset_inds(k)/fs));
    
    [pxx_on, f_on] = plomb(tachogram_y(is_on_blk), tachogram_x(is_on_blk),max_f,'normalized');
    is_lf_on = (f_on>= lf(1)) & (f_on<=lf(2));
    is_hf_on = (f_on >= hf(1)) & (f_on<=hf(2));
    
    lh_ratio_on(k) = sum(pxx_on(is_lf_on))/sum(pxx_on(is_hf_on));
    
    % OFF:
    is_off_blk = (tachogram_x <= (onset_inds(k)/fs)) & (tachogram_x >= (-analysis_win + onset_inds(k)/fs));
    
    [pxx_off, f_off] = plomb(tachogram_y(is_off_blk), tachogram_x(is_off_blk),max_f,'normalized');
    
    is_lf_off = (f_off>= lf(1)) & (f_off<=lf(2));
    is_hf_off = (f_off >= hf(1)) & (f_off<=hf(2));
    
    lh_ratio_off(k) = sum(pxx_off(is_lf_off))/sum(pxx_off(is_hf_off));
end

figure;
histogram(lh_ratio_on,linspace(0,8,12)); hold on; histogram(lh_ratio_off,linspace(0,8,12))
legend('On', 'Off'); xlabel('Low to High Power Ratio'); ylabel('Number of Stimulation Blocks')
[~,p] = kstest2(lh_ratio_on, lh_ratio_off);
title({'Distribution of Low Frequency to High Frequency Power in Tachogram Spectrum';...
    ['KS-test p-val = ' num2str(p,3)]})


%% universal comparison of spectrograms:
tachogram_is_on = VNS_dat.is_stim_on(round(tachogram_x*fs));
[pxx_on, f_on] = plomb(tachogram_y(tachogram_is_on), tachogram_x(tachogram_is_on),max_f,'normalized');
[pxx_off, f_off] = plomb(tachogram_y(~tachogram_is_on), tachogram_x(~tachogram_is_on), max_f, 'normalized');
smth_const = 0.002; % portion of total length being smoothed over

figure;
semilogy(f_on, smooth(pxx_on,round(smth_const*length(pxx_on)))) ; hold on; semilogy(f_off, smooth(pxx_off,round(smth_const*length(pxx_off))));
patch([lf(1) lf(1) lf(2) lf(2)], [sort(get(gca,'YLim')), sort(get(gca,'YLim'),'descend')], 'r', 'FaceAlpha',0.25,'EdgeColor','none')
patch([hf(1) hf(1) hf(2) hf(2)], [sort(get(gca,'YLim')), sort(get(gca,'YLim'),'descend')], 'b', 'FaceAlpha',0.25,'EdgeColor','none')
xlim([0 max_f])
legend('On', 'Off','Low Frequency', 'High Frequency')
xlabel('Frequency (Hz)'); ylabel('Normalized Spectral Power Density')
title(['Comparison of Tachogram Spectrum with VNS On vs VNS Off'])



% Extra feature: plot all data:
% % [pxx_all, f_all] = plomb(tachogram_y, tachogram_x, max_f,'normalized');
% % figure;
% % semilogy(f_all, smooth(pxx_all, round(smth_const*length(pxx_all)))); hold on;
% % semilogy(f_on, smooth(pxx_on,round(smth_const*length(pxx_on)))); 
% % semilogy(f_off, smooth(pxx_off,round(smth_const*length(pxx_off))));
% % xlabel('Frequency (Hz)'); ylabel('Normalized Spectral Power Density')
% % legend('All Data', 'VNS ON', 'VNS OFF')
% % title(['Comparison of Tachogram Spectrum with VNS On vs VNS Off'])
% % axis tight;