% FineFocus Examination (close look at details)
% VNS_dat = VNS_prelim_script_6(root_dir, 137,[97 34]);
% onsets = VNS_dat.stim_onsets_inds;
% fs = VNS_dat.sampFreq = fs_out;
% ecog_hg_env = VNS_dat.ecog_hg_env;


z_score_range = [-60 -30];

erps_hg = zeros(2,2001,length(onsets));
twin = [-10 10];
figure; hold on;
for k = 1:length(onsets)
    if k<=64
        subplot(1,2,1)
    else
        subplot(1,2,2)
    end
    hold on;
    range_ds = (onsets(k) + round(twin(1)*fs)):(onsets(k) + round(twin(2)*fs));
    range_z = (onsets(k) + round(z_score_range(1)*fs)):(onsets(k) + round(z_score_range(2)*fs));
    time_axis = linspace(twin(1), twin(2), length(range_ds));
    data = gdivide(gsubtract(ecog_hg_env(:,range_ds),mean(ecog_hg_env(:,range_z),2)), std(ecog_hg_env(:,range_z),[],2));
    if mean(data(:))>1.5;
        data = zeros(size(data));
    end
    data_raster = data/2 + k;
    plot(time_axis, data_raster(1,:),'r', time_axis,data_raster(2,:),'b');
    plot(twin, k, 'k')
    
    
   erps_hg(:,:,k) = data;
end
subplot(1,2,1)
legend('Depth', 'STG')
plot([0 0], get(gca,'YLim'),'k')
xlabel('Time (s)')
ylabel('Stimulation')
title('Intensity of High Gamma Activity around VNS Onset')

subplot(1,2,2)
legend('Depth', 'STG')
plot([0 0], get(gca,'YLim'),'k')
xlabel('Time (s)')
ylabel('Stimulation')
title('Intensity of High Gamma Activity around VNS Onset')






twin = [-10 10];
figure; hold on;
onsets_full_ds = round(onsets*fs_in/fs_out);
raw_onsets = [];

erps_raw = zeros(2, 20480,length(onsets_full_ds));
for k = 1:length(onsets_full_ds)
    if k<=64
        subplot(1,2,1)
    else
        subplot(1,2,2)
    end
    hold on;
    range_fs = (onsets_full_ds(k) + round(twin(1)*fs_in)):(onsets_full_ds(k) + round(twin(2)*fs_in));
    time_axis_fs = linspace(twin(1), twin(2), length(range_fs));
    data = all_data(:,range_fs)/2 + k;
    plot(time_axis_fs, data(1,:),'r', time_axis_fs,data(2,:),'b');
    plot(twin, k, 'k')
    raw_onsets = [raw_onsets all_data(:,range_fs)];
    
    erps_raw(:,:,k) = (data(:,1:20480) - k)*2;
end
subplot(1,2,1)
legend('Depth', 'STG')
plot([0 0], get(gca,'YLim'),'k')
xlabel('Time (s)')
ylabel('Stimulation')
title('Intensity of RAW Activity around VNS Onset')

subplot(1,2,2)
legend('Depth', 'STG')
plot([0 0], get(gca,'YLim'),'k')
xlabel('Time (s)')
ylabel('Stimulation')
title('Intensity of RAW Activity around VNS Onset')

%% Plot Distributions:
good_trials = [31:58, 62:79 84:108];
flatten = @(x) x(:);
%hg_before = erps_hg(:,1:1000,good_trials);
%hg_after = erps_hg(:,1001:end,good_trials);
hg_before = erps_hg(:,1:1000,:);
hg_after = erps_hg(:,1001:end,:);
dist_edges = linspace(-5, 5, 80);
figure;
subplot(1,2,1); hold on;
dat1 = flatten(squeeze(hg_before(1,:,:)));
dat1(dat1==0) = [];
dat2 = flatten(squeeze(hg_after(1,:,:)));
dat2(dat2==0) = [];
histogram(dat1,dist_edges,'EdgeColor', 'none');
histogram(dat2,dist_edges,'EdgeColor', 'none');
div1 = kl_div(dat1,dat2,[-5 5])

xlabel('High Gamma Intensity')
ylabel('Frequency')
title({'Distribution of High Gamma Before and After VNS in Hippocampus';'Z-score Pre-Onset, data from Stimulations with Low Noise'})

subplot(1,2,2); hold on;
dat1 = flatten(squeeze(hg_before(2,:,:)));
dat1(dat1==0) = [];
dat2 = flatten(squeeze(hg_after(2,:,:)));
dat2(dat2==0) = [];
histogram(dat1,dist_edges,'EdgeColor', 'none');
histogram(dat2,dist_edges,'EdgeColor', 'none');
div1 = kl_div(dat1,dat2,[-5 5])
xlabel('High Gamma Intensity')
ylabel('Frequency')
title({'Distribution of High Gamma Before and After VNS in STG';'Z-score Pre-Onset, data from Stimulations with Low Noise'})
