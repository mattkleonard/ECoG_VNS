function dat_filtered=apply60HzNotch_filter(dat, fs)
% Takes data and filters out harmonics of 60 Hz
if fs<200
    warning('ERROR - function should have been written work for sample frequencies <100Hz')
end
%fprintf('Notch Filters: 60 Hz..')
[b,a]=fir2(1000,[0 59 59.5 60.5 61 (fs/2)]/(fs/2),[1 1 0 0 1 1 ]);
%freqz(b,a,[],400)
dat_filtered=filtfilt(b,a,dat')';

%fprintf('120 Hz..')
[b,a]=fir2(1000,[0 119 119.5 120.5 121 (fs/2)]/(fs/2),[1 1 0 0 1 1 ]);
%freqz(b,a,[],400)
dat_filtered=filtfilt(b,a,dat_filtered')';

%fprintf('180 Hz..')
[b,a]=fir2(1000,[0 179 179.5 180.5 181 (fs/2)]/(fs/2),[1 1 0 0 1 1 ]);
%freqz(b,a,[],400)
dat_filtered=filtfilt(b,a,dat_filtered')';



%fprintf('\nNotch Filters Applied\n')
