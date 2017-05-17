function VNS_dat = artifact_detector_3(VNS_dat, varargin)
%% Finds bad time points in the envelope data of the 5 major frequency bands
% This is done by finding times where some number of channels have
% unusually high activity in a given band (as determined by a simple threshold)
% The result will append a set of 5 booleans (1 for each band) listing time
% where the data in a given band was high for each band.
% The function returns the original data structure with the bad time point
% booleans appended
%
% Inputs: 
% VNS_dat - VNS data strucuture that includes envelope measurements
%
% Variable Inputs
% var1 - power_thresh - (default 8) threshold zscore for a time point to be
% considered a bad time point.
%
% var2 - bad_ch_thresh - (default 1) the threshold number of channels that
% exceed the power_thresh for a time point to be considered bad
%
% var3 - bad_time_spread (default 20 (200ms at 100 hz)) the number of
% points centered around a bad time point that are also considered to be
% bad time points
%
% Output:
% VNS_dat - same data structure as input with bad_time_point booleans added
%
%
%% Load Variable Inputs:
% var1
power_thresh = 8; % limit on high gamma power to be used
if length(varargin)>0
    if ~isempty(varargin{1})
        power_thresh = varargin{1};
    end
end
% var2
bad_ch_thresh = 1; % number of channels needed to be declared a bad timepoint
if length(varargin)>1
    if ~isempty(varargin{2})
        bad_ch_thresh = varargin{2};
    end
end

% var3
bad_time_spread = 20; % number of points to spread a bad time point over
if length(varargin)>2
    if ~isempty(varargin{3})
        bad_time_spread = varargin{3};
    end
end


%% Threshold Each Band's ECoG to find bad data points:
bands = {'theta', 'alpha', 'beta', 'lg', 'hg'}; % appelation of various frequency bands

for b = 1:length(bands)
    ecog_dat = getfield(VNS_dat, ['ecog_' bands{b} '_env']);
    env_dat = abs(ecog_dat(VNS_dat.good_channels,:));
    
    exteme_pts = sum((env_dat > power_thresh),1); % Number of electrodes with high highgamma
    is_bad_timept_band = (exteme_pts >= bad_ch_thresh);
    is_bad_timept = (smooth(double(is_bad_timept_band), bad_time_spread) >= (1/bad_time_spread));

    VNS_dat = setfield(VNS_dat,['is_bad_time_' bands{b}], is_bad_timept);

end