function VNS_dat = VNS_find_bad_data(VNS_dat,thresh,nBadTrialThresh,h5_file_info,remove_badTrial_flag)

%% Find bad trials on each channel first

textprogressbar('Finding bad trials ');
for i = 1:size(VNS_dat.raw_erps,1)
    textprogressbar((i/size(VNS_dat.raw_erps,1))*100);
    if ~ismember(i,VNS_dat.EKG_ch)
        for j = 1:size(VNS_dat.raw_erps,3)
            if any(abs(squeeze(VNS_dat.raw_erps(i,:,j))) >= thresh)
                badTrials(i,j) = 1;
            else
                badTrials(i,j) = 0;
            end
        end
    end
end
fprintf('\n');

%% Identify channels with too many bad trials

textprogressbar('Finding bad channels ');
VNS_dat.badChans = [];
for i = 1:size(badTrials,1)
    textprogressbar((i/size(badTrials,1))*100);
    if length(find(badTrials(i,:))) >= nBadTrialThresh
        VNS_dat.badChans = [VNS_dat.badChans i];
    end
end
fprintf('\n');

%% Mark bad trials

textprogressbar('Refining bad trials ');
badTrials(VNS_dat.badChans,:) = 0;
badTrials(VNS_dat.EKG_ch,:) = 0;

VNS_dat.badTrials = [];
for i = 1:size(badTrials,2)
    textprogressbar((i/size(badTrials,2))*100);
    if any(badTrials(:,i))
        VNS_dat.badTrials = [VNS_dat.badTrials i];
    end
end
fprintf('\n');

%% Remove bad trials

if remove_badTrial_flag
    textprogressbar('Removing bad trials ');
    for i = length(VNS_dat.badTrials):-1:1
        textprogressbar((i/length(VNS_dat.badTrials))*100)
        VNS_dat.raw_erps(:,:,VNS_dat.badTrials(i)) = [];
    end
end
fprintf('\n');
