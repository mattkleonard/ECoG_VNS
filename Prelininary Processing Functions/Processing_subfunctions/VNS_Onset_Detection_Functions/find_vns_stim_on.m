function is_stim_on = find_vns_stim_on(ekg_data, fs_in, varargin)
%% This function finds the timing of the VNS stimulation. It uses the ekg channel
% as well as data about the frequency of the data, the stimulation and its
% duration and uses this to measure the power of VNS stimulation over time
% and then synchronize this with the best fit model given by the duty cycle
% of the stimulation. The data is also downsampled to match the sample rate
% used by other ecog measurements like high gamma.
%% Inputs:
% ekg_data: 1d array of data from the ekg channel
% fs_in: sample rate of ekg channel
%
%% Variable inputs:
% 1: fs_out (default 100) - the sample rate to down sample VNS onset
% measurement to. Should match sample rate of other ecog measurements
%
% 2: stim_freq (default 25) the frequency of VNS stimulation pulses
%
% 3: dutycycle_params (default [30,2]) the timeing of the profile of the
% theoretical trace of the VNS power, array, first value is the duration of
% peak power, the second value is the duration (in seconds of the
% onset/offset ramping time).

%% Load Variable Inputs:
fs_out = 100;
if length(varargin)>0
    if ~isempty(varargin{1})
        fs_out = varargin{1};
    end
end

stim_freq = 25;
if length(varargin)>1
    if ~isempty(varargin{2})
        stim_freq = varargin{2};
    end
end

dutycycle_params = [30,2];
if length(varargin)>2
    if ~isempty(varargin{3})
        dutycycle_params = varargin{3};
    end
end




%% Find VNS power:
    ecg_data_vns = extract_vns_stim_i(ekg_data,fs_in,stim_freq);
    ecg_data_vns_smth = smooth(ecg_data_vns,2*round(fs_in)); % smooth data to eliminate transient noise artifacts
    % Downsample VNS data
    ecg_data_vns_ds = resample(ecg_data_vns_smth, fs_out, fs_in);
    on_thresh = mean([median(ecg_data_vns_ds(ecg_data_vns_ds>0)), median(ecg_data_vns_ds(ecg_data_vns_ds<0))]); % midpoint between high and low value.
    is_stim_on = (ecg_data_vns_ds > on_thresh);

    %% use convlution with idea stimulation artifact to pinpoint onset:
    ecg_data_vns_ds = resample(ecg_data_vns, fs_out, fs_in);
    
    if length(dutycycle_params) == 1
        vns_pulse_ideal = generate_ideal_vns_pulse(fs_out,dutycycle_params(1));
    elseif length(dutycycle_params) == 2
        vns_pulse_ideal = generate_ideal_vns_pulse(fs_out,dutycycle_params(1),dutycycle_params(2));
    elseif length(dutycycle_params) == 3
        vns_pulse_ideal = generate_ideal_vns_pulse(fs_out,dutycycle_params(1),dutycycle_params(2),dutycycle_params(3));
    else
        vns_pulse_ideal = generate_ideal_vns_pulse(fs_out,dutycycle_params(1),dutycycle_params(2), dutycycle_params(3));
        warning('Too Many inputs on duty cycle - possible error')
    end
    is_stim_on_xcorr = false(size(is_stim_on));
    
    stim_duration_thresh = 27; % minimum time in seconds of stim
    [onsets, offsets] = find_boolean_on(is_stim_on, stim_duration_thresh*fs_out);  
    for n = 1:length(onsets)
        start_pad = 5; % seconds before start on onset to look
        stop_pad = 5; % seconds after start of onset to look
        start_ind = onsets(n)-round(start_pad*fs_out);
        stop_ind = onsets(n) + round((stim_duration_thresh+stop_pad)*fs_out);
        if start_ind < 1
            front_pad = zeros((1-start_ind),1);
            start_ind = 1;
        else
            front_pad = [];
        end
        if stop_ind > length(ecg_data_vns_ds)
            back_pad = zeros((start_ind - length(ecg_data_vns_ds)),1);
            stop_ind = length(ecg_data_vns_ds);
        else
            back_pad = [];
        end
        
        pulse_dat = [front_pad; ecg_data_vns_ds(start_ind:stop_ind); back_pad];
        [r,lag] = xcorr(pulse_dat, vns_pulse_ideal);
        [~,max_ind] = max(r);
        start_ind = onsets(n) + lag(max_ind);
        stop_ind = start_ind + round(dutycycle_params(2)*fs_out);
        if start_ind < 1
            start_ind = 1;
        end
        if stop_ind > length(ecg_data_vns_ds)
            stop_ind = length(ecg_data_vns_ds);
        end
        is_stim_on_xcorr(start_ind:stop_ind) = true;
    end
    
    is_stim_on = is_stim_on_xcorr;
end
 