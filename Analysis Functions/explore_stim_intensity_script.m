%% High Power (1:15 - 2:35)

timestamps = hdf5read('EC131_b4e25fe4-3fc2-4964-8369-d44b58bdf0d4_280.h5', 'timestamp vector');
tsvec = datevec((timestamps - (7*3600))/86400 + datenum(1970,1,1));

data = hdf5read('EC131_b4e25fe4-3fc2-4964-8369-d44b58bdf0d4_280.h5', 'ECoG Array'); % max power
ekg_ch = data(137,:);
data = data(1:137,:);

is_stim_on = find_vns_stim_on(ekg_ch , 1024, 1024, 25);

data_filt =applyLineNoiseBandpass_VNS_Harmonics(data,1024);

xcorr_on = corrcoef(data_filt(:,is_stim_on)');
ekg_corr_on = xcorr_on(137, 1:136);

xcorr_off = corrcoef(data_filt(:,~is_stim_on)');
ekg_corr_off = xcorr_off(137, 1:136);
figure; plot(ekg_corr_on - ekg_corr_off,'bx')

data_intensity = abs(hilbert(data_filt')');
xcorr_env = corrcoef(data_intensity');
ekg_corr_env = xcorr_env(137, 1:136);
ekg_corr_hp =ekg_corr_env;

%% Min Power (915 - 1035)
timestamps = hdf5read('EC131_b4e25fe4-3fc2-4964-8369-d44b58bdf0d4_220.h5', 'timestamp vector');
tsvec = datevec((timestamps - (7*3600))/86400 + datenum(1970,1,1));

data = hdf5read('EC131_b4e25fe4-3fc2-4964-8369-d44b58bdf0d4_220.h5', 'ECoG Array'); % max power
ekg_ch = data(137,:);
data = data(1:137,:);

is_stim_on = find_vns_stim_on(ekg_ch , 1024, 1024, 25);

data_filt =applyLineNoiseBandpass_VNS_Harmonics(data,1024);

xcorr_on = corrcoef(data_filt(:,is_stim_on)');
ekg_corr_on = xcorr_on(137, 1:136);

xcorr_off = corrcoef(data_filt(:,~is_stim_on)');
ekg_corr_off = xcorr_off(137, 1:136);
figure; plot(ekg_corr_on - ekg_corr_off,'bx')

data_intensity = abs(hilbert(data_filt')');
xcorr_env = corrcoef(data_intensity');
ekg_corr_env = xcorr_env(137, 1:136);

ekg_corr_lp =ekg_corr_env;

figure; plot(1:136, ekg_corr_lp,'bx', 1:136, ekg_corr_hp,'rx')
xlabel('Channel Number'); ylabel('Correlation in VNS Intensity')
title({'Effect of Stimulation Settings on the Correlation Between the Intensity of'; 'VNS Harmonics on Electrodes and EKG Channel'})
legend({'1.0 mA; 250 us PW','2.25 mA; 500 us PW'},'location', 'northwest')