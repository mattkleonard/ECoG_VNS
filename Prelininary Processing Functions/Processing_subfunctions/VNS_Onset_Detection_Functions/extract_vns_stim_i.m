function vns_stim = extract_vns_stim_i(data, varargin)
%% This function extracts the VNS stimulation power.
% This is taken to be the power (as measured by the analytic amplitude) of
% the ekg channel bandpassed filtered for the VNS stimulation frequency at
% its harmonics.
% Assuming that the data occurs effectively as a 25 hz spike with many
% harmonics, takes windowed data around the harmonics and takes the analytic
% ampliftude.
% Inputs:
% data - 1d array of ekg data
%
% Variable Inputs:
% sample rate of data (default 1024)
%
% sample_rate of stimulation (default 25)
fs_dat = 1024;
if length(varargin)>0
    if ~isempty(varargin{1})
        fs_dat = varargin{1};
    end
end
fs_stim = 25;
if length(varargin)>1
    if ~isempty(varargin{2})
        fs_stim = varargin{2};
    end
end

%freq_samples = [25:25:450];
freq_samples = [fs_stim:fs_stim:(fs_dat/2)]; % harmonics

% Cut undesireable freqs:
base_freq = 40; % minimum threshold to observe
wall_freq = 60;
freq_samples(freq_samples < base_freq) = []; % Clear low frequency harmonics (too much noise)
freq_samples(mod(freq_samples,wall_freq)==0) = []; % Clear wall socket harmonics (too much noise)

window_widths = 0.6*ones(size(freq_samples)); %

if mod(length(data),2) == 1; % make data even
    data(end) = [];
    append_tail = true;
else
    append_tail = false;
end

freq_axis = (fs_dat/2)*(0:floor(length(data)/2-1))/floor(length(data)/2-1);
freq_axis = [freq_axis, sort(freq_axis,'descend')];


vns_stim = zeros(length(data),1);
vns_stim_envs = zeros([length(data),length(freq_samples)]);

for k = 1:length(freq_samples)

    %% Filter ECoG Channels
    gaussian_window = @(f) exp(-((f-freq_samples(k)).^2)/(2*window_widths(k).^2));
    dat_fft = fft(data);
    T=length(data);
    h = zeros(1,T);
    if 2*fix(T/2)==T %if T is even
        h([1 T/2+1]) = 1;
        h(2:T/2) = 2;
    else
        h(1) = 1; h(2:(T+1)/2) = 2;
    end
    
    gauss_filter = gaussian_window(freq_axis);
    hilbdata=ifft(dat_fft(end,:).*(gauss_filter.*h),T);
    envData=abs(hilbdata);
    vns_stim_envs(:,k) = envData;
end
time_axis = linspace(0,length(data)/fs_stim,length(data));
% % % dat_filt = abs(dat_filt);


%% Z-score each band and Reduce Dimension

for j = 1:size(vns_stim_envs,2)
    vns_stim_envs(:,j) = (vns_stim_envs(:,j)-mean(vns_stim_envs(:,j)))/(std(vns_stim_envs(:,j)));
end

%pca_mat = pca(vns_stim_envs);
%vns_stim = vns_stim_envs*pca_mat(:,1);

% alt
vns_stim = mean(vns_stim_envs,2);

%% Add back trailing data if necessary (for odd length data)
if append_tail
    vns_stim = [vns_stim; vns_stim(end)];
end

end

