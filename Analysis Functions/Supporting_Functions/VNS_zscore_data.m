function VNS_dat = VNS_zscore_data(VNS_dat)

textprogressbar('Z-scoring data ');

for i = 1:size(VNS_dat.raw_erps,1)
    textprogressbar((i/size(VNS_dat.raw_erps,1))*100);
    if ~ismember(i,VNS_dat.EKG_ch)
        tmpDat = squeeze(VNS_dat.raw_erps(i,:,:));
        datMu = mean(tmpDat(:));
        datStd = std(tmpDat(:));
        
        for j = 1:size(VNS_dat.raw_erps,3)
            VNS_dat.raw_erps(i,:,j) = (VNS_dat.raw_erps(i,:,j) - datMu) / datStd;
        end
        
        VNS_dat.zscore_flag = 1;
    end
end
fprintf('\n');
