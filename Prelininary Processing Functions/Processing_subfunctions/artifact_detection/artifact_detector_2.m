function is_bad_timept = artifact_detector_2(VNS_dat)
%% Finds bad time points on the basis of unusually high high gamma activity
% sets a threshold and finds times that are very high and creates a boolean
% around these times. This boolean is then appended onto VNS_dat

%% Threshold High Gamma to find bad trials:
power_thresh = 10; % limit on high gamma power to be used
bad_ch_thresh = 1; % number of channels needed to be declared a bad timepoint

fs = VNS_dat.sampFreq;

env_dat = VNS_dat.ecog_hg_env(VNS_dat.good_channels,:);
%onsets = VNS_dat.stim_onsets_inds;
%onsets = VNS_dat.stim_onsets_inds_full;



high_hg = sum((env_dat > power_thresh),1); % Number of electrodes with high highgamma
is_bad_timept = (high_hg >= bad_ch_thresh);

%% smear bad timepoints:
bad_time_spread = 20; % number of points to spread a bad time point over
is_bad_timept = (smooth(double(is_bad_timept), bad_time_spread) >= (1/bad_time_spread));

%VNS_dat.is_bad_timept = is_bad_timept; % originally appened the data
%directly
end