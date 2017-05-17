function [erps,time_axis] = make_vns_erps(data, onset_inds,fs, erp_win)
%% Function to make ERPS time locked to VNS onset:
% (Actually quite general in it's application...)
% Inputs:
% data - N_channels x Timepts data matrix
% onset_inds - n_trialx 1 array of onsets
% fs - sample frequency
% erp_win - time win in seconds about the onset ex ([-10 10]);
%
% Output:
% erps - N_chan x Timepts x n_trials erp matrix
% time_axis - time (in seconds relative to the onset timepts for the time
% dimension of the erp matrix.

pre_event_dur = round(erp_win(1)*fs); % number of points to include before onset
post_event_dur = round(erp_win(2)*fs); % num pts after onset

%% clear trails where the data would be clipped by the beginning/end of the recording
is_OoBounds = (onset_inds+pre_event_dur) < 1 | (onset_inds + post_event_dur) > size(data,2);
onset_inds(is_OoBounds) = [];

% preallocate data:
erps = zeros(size(data,1), post_event_dur - pre_event_dur +1, length(onset_inds));
time_axis = linspace(erp_win(1), erp_win(2), size(erps,2));

for i = 1:length(onset_inds);
    erps(:,:,i) = data(:,(onset_inds(i)+pre_event_dur):(onset_inds(i)+post_event_dur));
end

plot_erp = false;
if plot_erp
    i = 1;
    figure; shadedErrorBar(time_axis, mean(squeeze(erps(i,:,:)),2), nansem(squeeze(erps(i, :, :)), 2));
end


end
