function [] = HRV_trend_measure(VNS_dat)
% Loads VNS data with RR data added and uses linear regression to 

fs = VNS_dat.sampFreq;
tachogram_y = VNS_dat.tachogram_y;
tachogram_x = VNS_dat.tachogram_x;

is_stim_on = VNS_dat.is_stim_on;
onset_inds = VNS_dat.stim_onsets_inds;
analysis_win = 30;
% clear edge:
onset_inds((onset_inds+analysis_win*fs>length(is_stim_on))...
    | (onset_inds-analysis_win*fs < 1)) = [];


%% Betas:
rr_slope_on = zeros(length(onset_inds),1);
rr_slope_off = zeros(length(onset_inds),1);
for k = 1:length(onset_inds)
    % on
    %is_on_blk = (tachogram_x >= (onset_inds(k)/fs)) & (tachogram_x <= (analysis_win + onset_inds(k)/fs));
    % modify to look at first 1 second:
    block_onset_win = [-3 3];
    is_on_blk = (tachogram_x >= block_onset_win(1)+(onset_inds(k)/fs)) & (tachogram_x <= (block_onset_win(2) + onset_inds(k)/fs));
    y = tachogram_y(is_on_blk)';
    x = [ones(sum(is_on_blk),1) tachogram_x(is_on_blk)'];
    
    beta_on = regress(y,x);
    rr_slope_on(k) = beta_on(2);
    
    % off
    %is_off_blk = (tachogram_x <= (onset_inds(k)/fs)) & (tachogram_x >= (-analysis_win + onset_inds(k)/fs));
    block_offset_win = [-3 3];
    is_off_blk = (tachogram_x >= block_offset_win(1)+34+(onset_inds(k)/fs)) & (tachogram_x <= (block_onset_win(2) + 34+ onset_inds(k)/fs));

    y = tachogram_y(is_off_blk)'; 
    x = [ones(sum(is_off_blk),1) tachogram_x(is_off_blk)'];

    beta_off = regress(y,x);
    rr_slope_off(k) = beta_off(2);
    
end

% % figure; histogram(rr_slope_off, linspace(-0.01,0.01,20)); hold on; histogram(rr_slope_on,linspace(-0.01,0.01,20)); legend('off', 'on')
% % xlabel('Distribution of Heart Rate Trend'); ylabel('Number of Stimulation Cycles')
% % legend({'Off', 'On'}); title('Distribution of RR Interval Beta Weights During VNS Cycle')

% Start/Stop of Block:
figure;
subplot(1,2,1); histogram(rr_slope_off); xlabel('Distribution of RR Trend'); ylabel('Number of Stim Cycles');
[~,p] = ttest(rr_slope_off);
title({'Distribution of Slope of RR-Interval during the 6 second Window About VNS Offset'; ['T-Test \mu \neq 0: p-val = ' num2str(p,3)]});

subplot(1,2,2); histogram(rr_slope_on); xlabel('Distribution of RR Trend'); ylabel('Number of Stim Cycles');
[~,p] = ttest(rr_slope_on);
title({'Distribution of Slope RR-Interval during the 6 second Window About VNS Onset'; ['T-Test \mu \neq 0: p-val = ' num2str(p,3)]});


%% Deltas:
% measure delta in mean and std of HRV between begining and end of stim
% blocks
meas_win_size = 8; % window size to measure mean and std at start and end
rr_delta_on = zeros(length(onset_inds),1);
rr_delta_off = zeros(length(onset_inds),1);

rrstd_delta_on = zeros(length(onset_inds),1);
rrstd_delta_off = zeros(length(onset_inds),1);

for k = 1:length(onset_inds)
    % on
    is_on_blk = (tachogram_x >= (onset_inds(k)/fs)) & (tachogram_x <= (analysis_win + onset_inds(k)/fs));
    y_blk = tachogram_y(is_on_blk);
    
    start_inds = 1:meas_win_size; end_inds = (length(y_blk)+1-meas_win_size):length(y_blk);
    rr_delta_on(k) = mean(y_blk(end_inds)) - mean(y_blk(start_inds));
    rrstd_delta_on(k) = std(y_blk(end_inds)) - std(y_blk(start_inds));
    
    % off
    is_off_blk = (tachogram_x <= (onset_inds(k)/fs)) & (tachogram_x >= (-analysis_win + onset_inds(k)/fs));
    y_blk = tachogram_y(is_off_blk);
    
    start_inds = 1:meas_win_size; end_inds = (length(y_blk)+1-meas_win_size):length(y_blk);
    rr_delta_off(k) = mean(y_blk(end_inds)) - mean(y_blk(start_inds));
    rrstd_delta_off(k) = std(y_blk(end_inds)) - std(y_blk(start_inds));
end

figure;
subplot(1,2,1);
histogram(rr_delta_on); hold on; histogram(rr_delta_off)
legend({'On', 'Off'}); xlabel('Change in RR Interval (s)'); ylabel('Number of Cycles');
[~,p] = kstest2(rr_delta_off, rr_delta_on);
title({'Change in RR Interval During Stimulation Cycle';['KS Test: p-val = ' num2str(p,3)]})
    
subplot(1,2,2);
histogram(rrstd_delta_on); hold on; histogram(rrstd_delta_off); 
legend({'On', 'Off'}); xlabel('Change in HRV (RRstd)'); ylabel('Number of Cycles');
[~,p] = kstest2(rrstd_delta_off, rrstd_delta_on);
title({'Change in RRstd During Stimulation Cycle';['KS Test: p-val = ' num2str(p,3)]})

