function [tstat, time_axis] = tstat_time_course(ecog, onsets, fs, window_len)
% This function takes ecog data and for each channel, aggregates data
% before and after the onset and finds the Kulbler-Liebleck divergence over
% time as a func
onsets = onsets(2:(end-1)); % Clip end pulses

%% define window:
kl_window_start = -30; % amount of time prior to onsets to start lookiin the window
kl_window_stop = 60;
overlap_frac = 0.1;
time_pts = ceil((1/overlap_frac)*(kl_window_stop-kl_window_start)/window_len); % set time points to occur every half window
time_axis = linspace(kl_window_start, kl_window_stop, time_pts);

%% Define off ecog:
off_start_inds = onsets - round(fs*(window_len - kl_window_start));
off_stop_inds = onsets - round(fs*(-kl_window_start));

ecog_off = [];
for k = 1:length(off_start_inds)
    ecog_off = [ecog_off ecog(:,off_start_inds(k):off_stop_inds(k))];
end
range_off = [min(ecog_off,[],2) max(ecog_off,[],2)];

%% Measure kl Divergence at each time point
tstat = zeros(size(ecog_off,1), length(time_axis));
for t = 1:length(time_axis)
   on_start_inds = onsets + round(fs*time_axis(t));
   on_stop_inds = onsets + round(fs*(time_axis(t)+window_len));
   ecog_on = [];
   for k = 1:length(on_start_inds)
       ecog_on = [ecog_on ecog(:,on_start_inds(k):on_stop_inds(k))];
   end
   range_on = [min(ecog_on,[],2) max(ecog_on,[],2)];
   
   tstat(:,t) = sqrt(size(ecog_on,2))*(mean(ecog_on,2) - mean(ecog_off,2))./sqrt(var(ecog_on,[],2)+var(ecog_off,[],2));
   
 %       num_bins = round(size(ecog_on,2)^(1/3));
 %       edges = linspace(-5,5 ,num_bins);

end
   
   
   
   