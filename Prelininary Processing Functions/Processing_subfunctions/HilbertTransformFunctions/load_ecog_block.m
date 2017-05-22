function ecog_dat = load_ecog_block(raw_dat, freq_range, get_env, fs_in, fs_out, dat_length_ds)
%% This function loads the channel by channel ecog data for a block of clinincal
% dat and returns the downsampled processed ecog data
%
% Inputs:
% raw_dat - n_chans x timepts raw ecog data sampled at fs_in
% freq_range - [f_min, f_max] the range of frequencies to be processed
% get_env - boolean: (ON)calculate mean Analytic Amplitude of bands
%                    (OFF) calculate mean Phase Angle of bands
% fs_in - sample rate of input
% fs_out - the downsample rate of output
% dat_length_ds - the length of downsampled data
%
% Output:
% ecog_dat - n_chans x timepts processed data sampled at fs_out

% EX:
%ecog_theta_env = load_ecog_block(raw_dat, [4 8], true, 1024, 100, dat_length_ds(k))

%% INITIALIZE OUTPUT:
ecog_dat = zeros(size(raw_dat,1),dat_length_ds);

%% Anlytic Amplitude
if get_env
    for j = 1:size(raw_dat,1)
        fprintf('Electrode [%d] of [%d]\n',j,size(raw_dat,1));
        tmp.data = raw_dat(j,:);
        tmp.sampFreq = fs_in;
        filteredData = processingHilbertTransform_filterbankGUI_meanout(tmp,fs_in, freq_range);
         ecog_tmp = resample(filteredData.data,fs_out,fs_in); % Resample Data
%         % Sliding z-score
%         window_start = -round(60*fs_out);
%         window_end = -round(30*fs_out); % stop of zscoring window
%         first_pt = -window_start;
%         ecog_dat(j,1:first_pt) = (ecog_tmp(1:first_pt)- mean(ecog_tmp(1:first_pt)))/std(ecog_tmp(1:first_pt));
%         for k = (first_pt+1):length(ecog_tmp);
%             ecog_dat(j,k) = (ecog_tmp(k) - mean(ecog_tmp((k+window_start):(k+window_end))))/std(ecog_tmp((k+window_start):(k+window_end)));
%         end
        ecog_dat(j,:) = (ecog_tmp  - mean(ecog_tmp))/std(ecog_tmp); % z-score data
        ecog_dat(j,1:20) = NaN;
        ecog_dat(j,end-20:end) = NaN;
    end
else
    for j = 1:size(raw_dat,1)
        tmp.data = raw_dat(j,:);
        tmp.sampFreq = fs_in;
        filteredData = processingHilbertTransform_filterbankGUI_mean_phase_out(tmp,fs_in, freq_range);
        ecog_dat(j,:) = angle(resample(exp(sqrt(-1)*filteredData.data),fs_out, fs_in)); % Angular Resample
        ecog_dat(j,1:20) = NaN;
        ecog_dat(j,end-20:end) = NaN;
    end
end