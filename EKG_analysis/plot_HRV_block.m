function [] = plot_HRV_block(block_names, ekg_channel_num);
%% Loads a single block or set of blocks of clinical recording
% with VNS stimulation and plots the heart rate variability measuremetns
% 
% Input: 
% block_names: list of data block names
%
% ekg_channel_num: ekg channel number (if 2 elements, second is reference)
% 
% Ex:
% block_names = 'EC131_759989b0-8f7a-457e-ab6d-61d9ccd224bc_080.h5'
% block_names = {'EC131_759989b0-8f7a-457e-ab6d-61d9ccd224bc_079.h5', ...
%                   'EC131_759989b0-8f7a-457e-ab6d-61d9ccd224bc_080.h5'};
% ekg_channel_num = 137;
% ekg_channel_num = 137:138;

fs = 1024; % sample rate

if ~iscell(block_names)
    data = hdf5read(block_names, 'ECoG Array');
else
    data = [];
    for k = 1:length(block_names)
        data_blk =  hdf5read(block_names{k}, 'ECoG Array');
        data = [data data_blk];
    end
end
if length(ekg_channel_num) == 1
    ekg_ch = data(ekg_channel_num,:);
else
    ekg_ch = data(ekg_channel_num(1),:) - data(ekg_channel_num(2),:);
end

is_stim_on = find_vns_stim_on(ekg_ch, fs, fs, 25);
time_axis = linspace(0,length(ekg_ch)/1024, length(ekg_ch));
%ekg_ch(1:6000) = [];

[~,qrs_i_raw,delay]=pan_tompkin(ekg_ch,1024, false);
delay = 4; qrs_i_raw = qrs_i_raw - delay;
rr_interval = diff(qrs_i_raw)/fs;
rr_inds = round((qrs_i_raw(2:end)+qrs_i_raw(1:(end-1)))/2); % mid point between pulses
rr_times = time_axis(rr_inds);
is_on_pulse = is_stim_on(rr_inds);

figure;

%% Plot EKG:
subplot(3,1,1)
plot(time_axis, ekg_ch); hold on; plot(time_axis, 100*is_stim_on,'r')
axis tight
xlabel('Time (s)'); ylabel('EKG Signal')
plot(time_axis(qrs_i_raw), ekg_ch(qrs_i_raw) ,'r.')
title('EKG Channel');
legend('EKG Ch', 'VNS On','R Peaks')

%% Plot RR_Rimes and intervals
subplot(3,1,2)
plot(rr_times(is_on_pulse), rr_interval(is_on_pulse), 'rx', rr_times(~is_on_pulse), rr_interval(~is_on_pulse), 'bx')
axis tight
legend('VNS On', 'VNS Off')
xlabel('Time (s)'); ylabel('RR Interval (s)'); 
title('RR Intervals')


%% Plot HR VAR
subplot(3,1,3)
var_win_size = 5; % num pts before & num pts after
HRV_rrstd = zeros(length(rr_interval)-2*(var_win_size-1),1);
HRV_times = rr_times(var_win_size:(length(rr_times)+1-var_win_size));
is_on_var = is_on_pulse(var_win_size:(length(rr_times)+1-var_win_size));

for k = 1:length(HRV_rrstd)
    HRV_rrstd(k) = std(rr_interval(k:(k+var_win_size)));
end
plot(HRV_times(is_on_var), HRV_rrstd(is_on_var), 'rx', HRV_times(~is_on_var), HRV_rrstd(~is_on_var), 'bx')
axis tight
legend('VNS On', 'VNS Off')
xlabel('Time (s)'); ylabel('RRstd'); 
title(['RR STD over ' num2str(2*var_win_size+1) ' point sliding window'])
