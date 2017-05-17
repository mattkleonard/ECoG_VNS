function [kl_div, time_axis] = kl_divergence_time_course_2(ecog, onsets, fs, varargin)
% This function is takes a nchannelxtime matrix of ecog data 
% forms it into erps to use to measure the changes in the distribution of data over time.
% The kl_divergence between the points before onsets and a sliding window
% of times is then measured and compared with the central time.
%
% inputs:
% ecog - n_channel x time_pts matrix of ecog data
% onsets - n_trails array of the onset indecies of VNS
% fs - sample rate of the data
%
% Variable Inputs:
% 
% 1 - window_len - (default: 10) the size of the sliding window used to form the distribution
% 
% 2 - time_course_win (default: [-30 60]) start/stop times of KL time
% course
%
% outputs:
% kl_div: the kl_divergence between various windows of ecog and the
% distribution of ecog for onsets
% time_axis - the corresponding time axis for this data.

%% Variable Inputs:
window_len = 10;
if length(varargin)>0
    if ~isempty(varargin{1})
        window_len = varargin{1};
    end
end
kl_win = [-30 60];
if length(varargin)>1
    if ~isempty(varargin{2})
        kl_win = varargin{2};
    end
end
kl_window_start = kl_win(1); % amount of time prior to onsets to start lookiin the window
kl_window_stop = kl_win(2); % amount of time after onsets to advance the window


[erps, time_erps] = make_vns_erps(ecog, onsets,fs, [kl_window_start-window_len/2, kl_window_stop+window_len/2]);


%% Define off ecog:
ecog_off = erps(:,time_erps<0,:);
dims = size(ecog_off);
ecog_off = reshape(ecog_off, dims(1), prod(dims(2:3)));
range_off = [min(ecog_off,[],2) max(ecog_off,[],2)];


overlap_frac = 0.1;
time_pts = ceil((1/overlap_frac)*(kl_window_stop-kl_window_start)/window_len); % set time points to occur every half window
time_axis = linspace(kl_window_start, kl_window_stop, time_pts);


%% Measure kl Divergence at each time point
kl_div = zeros(size(ecog_off,1), length(time_axis));
for t = 1:length(time_axis)
    in_win = (time_erps >= (time_axis(t)-window_len/2)) & (time_erps <= (time_axis(t)+window_len/2));
    ecog_on = erps(:,in_win,:);
    dims = size(ecog_on);
    ecog_on = reshape(ecog_on, dims(1), prod(dims(2:3)));
    
   range_on = [min(ecog_on,[],2) max(ecog_on,[],2)];
   
   num_bins = round(size(ecog_on,2)^(1/3));
   for k = 1:size(ecog_on,1)   
       %edges = linspace(min([range_on(k,1) range_off(k,1)]), max([range_on(k,2) range_off(k,2)]), num_bins);
       edges = linspace(-5,5 ,num_bins);

       dist_on = histcounts(ecog_on(k,:),edges)/sum(~isnan(ecog_on(1,:)));
       dist_off = histcounts(ecog_off(k,:), edges)/sum(~isnan(ecog_off(1,:)));
       %%%dist_off(dist_off == 0) = 0.1/size(ecog_off,2);
       null_cals = (dist_on == 0) | (dist_off == 0);
       probs = dist_on.*log2(dist_on./dist_off);
       
       %kl_div(k,t) = t_welch(ecog_on(k,:),ecog_off(k,:));
       %[~,~,ks] = kstest2(ecog_on(k,:),ecog_off(k,:));
       %kl_div(k,t) = ks;
       kl_div(k,t) = sum(probs(~null_cals));
       % as I understand it, the KL_divergence measures how surprised by
       % the outcome an observer would be if they expected the observed values 
       % to come from VNS off but they instead came from VNS on, compared to when
       % they expect values to come from VNS off and they do come from VNS
       % off. 
   end
end
   
   
a = 1;

end
   