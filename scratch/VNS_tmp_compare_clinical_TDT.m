ch = 5;
tlim = 10000:100000;
fs = 1024;
cutoff = 200;

dat.day1 = h5read('h5/EC131_05104e39-8cff-4ea5-b1d5-6f766084f91c_007.h5','/ECoG Array');
dat.day2 = h5read('h5/EC131_759989b0-8f7a-457e-ab6d-61d9ccd224bc_015.h5','/ECoG Array');
[dat.tdt1,fs_tdt] = readhtk('/Users/mattleonard/Desktop/EC131_tmp/Wav17.htk');
[dat.tdt2,fs_tdt] = readhtk('/Users/mattleonard/Desktop/EC131_tmp/Wav19.htk');

%%

dat.tdt = dat.tdt1 - dat.tdt2;

[p,q] = rat(fs/fs_tdt);
dat.tdt = resample(dat.tdt,p,q);
%%

figure;
subplot(3,2,1);
plot(dat.day1(ch,tlim));
title('clinical day 1');
subplot(3,2,3);
plot(dat.day2(ch,tlim));
title('clinical day 2');
subplot(3,2,5);
plot(dat.tdt(1,tlim));
title('TDT');

subplot(3,2,2);
[pxx,f_axis] = pwelch(dat.day1(ch,tlim),[],[],[],fs);
plot(f_axis(find(f_axis <= cutoff)),pxx(find(f_axis <= cutoff)));
title('clinical day 1');
ylabel('power');
xlabel('frequency');
subplot(3,2,4);
[pxx,f_axis] = pwelch(dat.day2(ch,tlim),[],[],[],fs);
plot(f_axis(find(f_axis <= cutoff)),pxx(find(f_axis <= cutoff)));
title('clinical day 2');
ylabel('power');
xlabel('frequency');
subplot(3,2,6);
[pxx,f_axis] = pwelch(dat.tdt(1,tlim),[],[],[],fs);
plot(f_axis(find(f_axis <= cutoff)),pxx(find(f_axis <= cutoff)));
title('TDT');
ylabel('power');
xlabel('frequency');

%%

elec = [4 7];

load('/Users/mattleonard/Documents/Research/pia/data_store2/imaging/subjects/EC131/Meshes/EC131_lh_pial.mat');
figure;
ctmr_gauss_plot(cortex,[0 0 0],0,'lh');
hold on;

load('/Users/mattleonard/Documents/Research/pia/data_store2/imaging/subjects/EC131/elecs/clinical_elecs_all.mat');
scatter3(elecmatrix(elec(1),1),elecmatrix(elec(1),2),elecmatrix(elec(1),3),50,'r','filled');

load('/Users/mattleonard/Documents/Research/pia/data_store2/imaging/subjects/EC131/elecs/TDT_elecs_all.mat');
for i = 1:128
    scatter3(elecmatrix(i,1),elecmatrix(i,2),elecmatrix(i,3),50,'b');
end

scatter3(elecmatrix(elec(2),1),elecmatrix(elec(2),2),elecmatrix(elec(2),3),50,'g');

