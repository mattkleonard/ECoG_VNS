function [kl_div, time_axis] = kl_divergence_time_course(ecog, onsets, fs, window_len)
% This function takes ecog data and for each channel, aggregates data
% before and after the onset and finds the Kulbler-Liebleck divergence over
% time as a func
%onsets = onsets(2:(end-1)); % Clip end pulses

%% define window:
kl_window_start = -30; % amount of time prior to onsets to start lookiin the window
kl_window_stop = 60;
overlap_frac = 0.1;
time_pts = ceil((1/overlap_frac)*(kl_window_stop-kl_window_start)/window_len); % set time points to occur every half window
time_axis = linspace(kl_window_start, kl_window_stop, time_pts);

%% Define off ecog:
%off_start_inds = onsets - round(fs*(window_len - kl_window_start));
%off_stop_inds = onsets - round(fs*(-kl_window_start));
off_start_inds = onsets - round(fs*50);
off_stop_inds = onsets - round(fs*5); 

ecog_off = [];
for k = 1:length(off_start_inds)
    ecog_off = [ecog_off ecog(:,off_start_inds(k):off_stop_inds(k))];
end
range_off = [min(ecog_off,[],2) max(ecog_off,[],2)];

%% Measure kl Divergence at each time point
kl_div = zeros(size(ecog_off,1), length(time_axis));
for t = 1:length(time_axis)
   on_start_inds = onsets + round(fs*time_axis(t)-window_len/2);
   on_stop_inds = onsets + round(fs*(time_axis(t)+window_len/2));
   ecog_on = [];
   for k = 1:length(on_start_inds)
       ecog_on = [ecog_on ecog(:,on_start_inds(k):on_stop_inds(k))];
   end
   
%   if sum(t==[44,35])>0
%        a = 1;
%    end
   range_on = [min(ecog_on,[],2) max(ecog_on,[],2)];
   for k = 1:size(ecog_on,1)
       num_bins = round(size(ecog_on,2)^(1/3));
          %edges = linspace(min([range_on(k,1) range_off(k,1)]), max([range_on(k,2) range_off(k,2)]), num_bins);
       edges = linspace(-5,5 ,num_bins);

%       dist_on = histcounts(ecog_on(k,:),edges)/size(ecog_on,2);
%       dist_off = histcounts(ecog_off(k,:), edges)/size(ecog_off,2);
%       %%%dist_off(dist_off == 0) = 0.1/size(ecog_off,2);
%       null_cals = (dist_on == 0) | (dist_off == 0);
%       probs = dist_on.*log2(dist_on./dist_off);
       
       kl_div(k,t) = t_welch(ecog_on(k,:),ecog_off(k,:));
       %[~,~,ks] = kstest2(ecog_on(k,:),ecog_off(k,:));
       %kl_div(k,t) = ks;
       %kl_div(k,t) = sum(probs(~null_cals));
   end
end
   
   
a = 1;

end
   