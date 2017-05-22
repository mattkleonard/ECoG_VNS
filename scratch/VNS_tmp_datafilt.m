%% SCRATCH CODE FOR tVNS ANALYSIS

subj = 'EC151';
block = 'B46';

ch = 69;
array_idx = [1 800000];
dsFs = 1024;

cutoff = 200;

rootdir = '/Users/mattleonard/Documents/Research/pia/data_store1/human/prcsd_data';

%% LOAD DATA

if ~exist('dat','var')
    dat = h5read([rootdir '/' subj '/' subj '_' block '/h5/' subj '_' block '.h5'],...
        '/ECoG Array',...
        [ch 1],...
        array_idx);
    fs = h5readatt([rootdir '/' subj '/' subj '_' block '/h5/' subj '_' block '.h5'],...
        '/ECoG Array',...
        'Sampling Rate');
    load([rootdir '/../../../data_store2/imaging/subjects/' subj '/elecs/TDT_elecs_all.mat']);
    
    fprintf('Loading CAR channels....\n');
    wm_ch_idx = find(strcmpi('Left-Cerebral-White-Matter',anatomy(:,4)));
    for i = 1:length(wm_ch_idx)
        fprintf('[%d] of [%d]\n',i,length(wm_ch_idx));
        car_dat(i,:) = h5read([rootdir '/' subj '/' subj '_' block '/h5/' subj '_' block '.h5'],...
            '/ECoG Array',...
            [wm_ch_idx(i),1],...
            array_idx);
    end
end

%% DOWNSAMPLE DATA AND GET CAR DATA

[p,q] = rat(dsFs/fs);
for i = 1:size(dat,1)
    dat_ds(i,:) = resample(dat(i,:),p,q);
end

car_dat_mean = mean(car_dat,1);
car_dat_ds = resample(car_dat_mean,p,q);

figure;
for i = 1:size(dat_ds,1)
    if size(dat_ds,1) > 1
        p = plotGridPosition_new(i,size(dat_ds,1),ceil(sqrt(size(dat_ds,1))));
        subplot('Position',p);
    end

    plot(dat_ds(i,:));
end

figure;
for i = 1:size(dat_ds,1)
    if size(dat_ds,1) > 1
        p = plotGridPosition_new(i,size(dat_ds,1),ceil(sqrt(size(dat_ds,1))));
        subplot('Position',p);
    end
    [Pxx(i,:),F(i,:)] = pwelch(dat_ds(i,:),[],[],[],dsFs);
    plot(F(i,:),Pxx(i,:));
end

%% 60 Hz NOTCH FILTER

for i = 1:size(dat,1)
    dat_ds_notch60(i,:) = apply60HzNotch_filter(dat_ds(i,:),dsFs);
end

figure;
for i = 1:size(dat_ds_notch60,1)
    if size(dat_ds_notch60,1) > 1
        p = plotGridPosition_new(i,size(dat_ds_notch60,1),ceil(sqrt(size(dat_ds_notch60,1))));
        subplot('Position',p);
    end

    plot(dat_ds_notch60(i,:));
end

figure;
for i = 1:size(dat_ds_notch60,1)
    if size(dat_ds_notch60,1) > 1
        p = plotGridPosition_new(i,size(dat_ds_notch60,1),ceil(sqrt(size(dat_ds_notch60,1))));
        subplot('Position',p);
    end
    [Pxx(i,:),F(i,:)] = pwelch(dat_ds_notch60(i,:),[],[],[],dsFs);
    plot(F(i,:),Pxx(i,:));
end

%% COMPARE RAW, NOTCH, AND CAR

nRows = 4;

dat_ds_notch60 = apply60HzNotch_filter(dat_ds,dsFs);
% notchFreq = 60;
% order = 2; % To make a 4th order filter
% display('Notch filtering')
% while notchFreq < cutoff % changed from Fs/2, since bandpass is at 200Hz
%     [b1, a1] = butter(order, [notchFreq-2 notchFreq+2]/ (dsFs/2), 'stop');
%     dat_ds_notch60 = filtfilt(b1,a1,dat_ds);
%     notchFreq = notchFreq + 60; % get harmonics
% end


figure;
subplot(nRows,2,1);
plot(dat_ds);
title('Orig signal');
subplot(nRows,2,3);
plot(dat_ds_notch60);
title('Notch60Hz signal');

[Pxx,F] = pwelch(dat_ds,[],[],[],dsFs);
subplot(nRows,2,2);
plot(F(find(F<cutoff)),Pxx(find(F<cutoff)));
title('Orig spectrum');
[Pxx,F] = pwelch(dat_ds_notch60,[],[],[],dsFs);
subplot(nRows,2,4);
plot(F(find(F<cutoff)),Pxx(find(F<cutoff)));
title('Notch60Hz spectrum');

subplot(nRows,2,5);
plot(dat_ds_notch60 - car_dat_ds);
title('Notch60Hz-CAR signal');
subplot(nRows,2,6);
[Pxx,F] = pwelch(dat_ds_notch60 - car_dat_ds,[],[],[],dsFs);
plot(F(find(F<cutoff)),Pxx(find(F<cutoff)));
title('Notch60Hz-CAR spectrum');

car_dat_ds_notch60 = apply60HzNotch_filter(car_dat_ds,dsFs);
subplot(nRows,2,7);
plot(dat_ds_notch60 - car_dat_ds_notch60);
title('Notch60Hz-CARnotch60Hz signal');
subplot(nRows,2,8);
[Pxx,F] = pwelch(dat_ds_notch60 - car_dat_ds_notch60,[],[],[],dsFs);
plot(F(find(F<cutoff)),Pxx(find(F<cutoff)));
title('Notch60Hz-CARnotch60Hz spectrum');

%% PLOT RAW ERPS
% NEED TO LOAD e.g., EC151_B46_VNS_raw.mat (output of BL's code for raw)
%   Best when run at dbstop point in BL's code (e.g., to get var
%   electrodes)

figure;
for i = 1:length(electrodes)
    p = plotGridPosition_new(i,length(electrodes),ceil(sqrt(length(electrodes))));
    subplot('Position',p);
    
    plot(squeeze(mean(VNS_raw.raw_erps(i,:,:),3)));
    hold on;
    line([round(size(VNS_raw.raw_erps,2)/2) round(size(VNS_raw.raw_erps,2)/2)],get(gca,'YLim'),'Color','k');
end

%% PLOT UN-EPOCHED DATA FOR RAW AND HILBERT DATA FROM BL's CODE
% NEED TO LOAD e.g., EC151_B46_VNS_dat.mat (output of BL's code for hilb)

ch = 1;
tmp = [VNS_dat.ecog_raw(electrodes(ch),:) ;...
    VNS_dat.ecog_theta_env(ch,:) ;...
    VNS_dat.ecog_alpha_env(ch,:) ;...
    VNS_dat.ecog_beta_env(ch,:) ;...
    VNS_dat.ecog_lg_env(ch,:) ;...
    VNS_dat.ecog_hg_env(ch,:)];

flds = {'raw','theta','alpha','beta','lg','hg'};
figure;
for i = 1:size(tmp,1)
    subplot(size(tmp,1),1,i);
    plot(tmp(i,:));
    title(flds{i});
end

%% COMPARE NON-DOWNSAMPLED, NON-NOTCHED DATA IN VNS vs. NON-VNS BLOCKS

nRows = 2;
ch = 69;
fs_orig = fileInfo.Datasets(1).Attributes(2).Value;
cutoff = 200;

nonVNS_data = readhtk('/Users/mattleonard/Documents/Research/pia/data_store1/human/prcsd_data/EC151/EC151_B47/RawHTK/Wav25.htk');

figure;
subplot(nRows,2,1);
plot(data(ch,:));
title('VNS Orig signal');

subplot(nRows,2,2);
[Pxx,F] = pwelch(data(ch,:),[],[],[],fs_orig);
plot(F(find(F<cutoff)),Pxx(find(F<cutoff)));
title('VNS Orig spectrum');

subplot(nRows,2,3);
plot(nonVNS_data(1,:));
title('non-VNS Orig signal');

subplot(nRows,2,4);
[Pxx,F] = pwelch(nonVNS_data(1,:),[],[],[],fs_orig);
plot(F(find(F<cutoff)),Pxx(find(F<cutoff)));
title('non-VNS Orig spectrum');

%% COMPARE VNS ON vs. VNS OFF PERIODS
% NEED TO LOAD e.g., EC151_B46_VNS_dat.mat (output of BL's code for hilb)

nRows = 2;
ch = 69;
trl = 2;
cutoff = 200;

clear trlDat;

for i = 2:length(VNS_dat.stim_onsets_inds)
    trlDat.on(:,:,i) = VNS_dat.ecog_raw(:,(VNS_dat.stim_onsets_inds(i):VNS_dat.stim_onsets_inds(i)+round((28*VNS_dat.sampFreq))));
    trlDat.off(:,:,i) = VNS_dat.ecog_raw(:,(VNS_dat.stim_onsets_inds(i)-round((28*VNS_dat.sampFreq)):VNS_dat.stim_onsets_inds(i)));
end

figure;
subplot(nRows,2,1);
plot(squeeze(trlDat.on(ch,:,trl)));
title('VNS ON Orig signal');

subplot(nRows,2,2);
[Pxx,F] = pwelch(squeeze(trlDat.on(ch,:,trl)),[],[],[],VNS_dat.sampFreq);
plot(F(find(F<cutoff)),Pxx(find(F<cutoff)));
title('VNS ON Orig spectrum');

subplot(nRows,2,3);
plot(squeeze(trlDat.off(ch,:,trl)));
title('VNS OFF Orig signal');

subplot(nRows,2,4);
[Pxx,F] = pwelch(squeeze(trlDat.off(ch,:,trl)),[],[],[],VNS_dat.sampFreq);
plot(F(find(F<cutoff)),Pxx(find(F<cutoff)));
title('VNS OFF Orig spectrum');

%% COMPARE NO CAR vs. CAR

noCAR_data = load('EC151_B46_dat_ds_VNS_60Hz.mat');
CAR_data = load('EC151_B46_dat_ds_VNS_60Hz_45chCAR.mat');

nRows = 2;
ch = 69;
cutoff = 200;

figure;
subplot(nRows,2,1);
plot(noCAR_data.data(ch,:));
title('Notch signal');
subplot(nRows,2,2);
[Pxx,F] = pwelch(noCAR_data.data(ch,:),[],[],[],VNS_raw.sampFreq);
plot(F(find(F<cutoff)),Pxx(find(F<cutoff)));
title('Notch spectrum');

subplot(nRows,2,3);
plot(CAR_data.data(ch,:));
title('Notch+CAR signal');
subplot(nRows,2,4);
[Pxx,F] = pwelch(CAR_data.data(ch,:),[],[],[],VNS_raw.sampFreq);
plot(F(find(F<cutoff)),Pxx(find(F<cutoff)));
title('Notch+CAR spectrum');

%% END tVNS ANALYSIS
%%
%%
%%
%%
%%
%%
%% iVNS

subj = 'EC131';
block = '05104e39-8cff-4ea5-b1d5-6f766084f91c_';
dat_num = '015';

ch = 69;
array_idx = [1 800000];
cutoff = 200;

rootdir = '/Users/mattleonard/Documents/Research/data';

%% LOAD DATA

if ~exist('dat','var')
    dat = h5read([rootdir '/' subj '/VNS/' subj '_QuietTime/h5/' subj '_' block dat_num '.h5'],...
        '/ECoG Array',...
        [1 1],...
        [inf inf]);
    fs = h5readatt([rootdir '/' subj '/VNS/' subj '_QuietTime/h5/' subj '_' block dat_num '.h5'],...
        '/ECoG Array',...
        'Sampling Rate');
    load([rootdir '/../pia/data_store2/imaging/subjects/' subj '/elecs/clinical_elecs_all.mat']);
    
    good_ch = 1:size(anatomy,1);
end

%% PLOT RAW DATA AND PSD

figure;
for i = 1:length(good_ch)
    if size(dat,1) > 1
        p = plotGridPosition_new(i,length(good_ch),ceil(sqrt(length(good_ch))));
        subplot('Position',p);
    end

    plot(dat(i,:));
    set(gca,'XTickLabel',[],'YTickLabel',[]);
    text(0,max(get(gca,'YLim')),num2str(i));
end

figure;
for i = 1:length(good_ch)
    if size(dat,1) > 1
        p = plotGridPosition_new(i,length(good_ch),ceil(sqrt(length(good_ch))));
        subplot('Position',p);
    end

    [Pxx,F] = pwelch(dat(i,:),[],[],[],fs);
    plot(F(find(F<cutoff)),Pxx(find(F<cutoff)));
    text(0,max(get(gca,'YLim')),num2str(i));
end

%% PLOT 60Hz NOTCH DATA AND PSD

for i = 1:length(good_ch)
    fprintf('60Hz notch for elec [%d] of [%d]\n',i,length(good_ch));
    dat_notch60(i,:) = apply60HzNotch_filter(dat(i,:),fs);
end

figure;
for i = 1:length(good_ch)
    if size(dat,1) > 1
        p = plotGridPosition_new(i,length(good_ch),ceil(sqrt(length(good_ch))));
        subplot('Position',p);
    end

    plot(dat_notch60(i,:));
    set(gca,'XTickLabel',[],'YTickLabel',[]);
    text(0,max(get(gca,'YLim')),num2str(i));
end

figure;
for i = 1:length(good_ch)
    if size(dat,1) > 1
        p = plotGridPosition_new(i,length(good_ch),ceil(sqrt(length(good_ch))));
        subplot('Position',p);
    end

    [Pxx,F] = pwelch(dat_notch60(i,:),[],[],[],fs);
    plot(F(find(F<cutoff)),Pxx(find(F<cutoff)));
    text(0,max(get(gca,'YLim')),num2str(i));
end

%% EKG

ch = [137 138];

figure;
plot(dat_notch60(ch(1),:) - dat_notch60(ch(2),:));

figure;
[Pxx,F] = pwelch(dat_notch60(ch(1),:) - dat_notch60(ch(2),:),[],[],[],fs);
plot(F(find(F<cutoff)),Pxx(find(F<cutoff)));

is_stim_on_blk = find_vns_stim_on(dat_notch60(ch(1),:) - dat_notch60(ch(2),:), fs, fs, 25, [26 2 2]);