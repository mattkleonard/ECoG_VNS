function ecog_filt = applyLineNoiseBandpass_VNS_Harmonics(ecog,sampFreq, varargin)

max_harmonic = min([500, (sampFreq/2 - 1)]);
notchFreq = 25;
if length(varargin)>0
    if ~isempty(varargin{1})
        notchFreq=varargin{1};
    end
end
harmonics = notchFreq:notchFreq:max_harmonic;

%% Generate Notch Labels
fir_string = '[0 ';
intensity_string = '[0 ';
for k = 1:length(harmonics)
    fir_string = [fir_string num2str(harmonics(k)-1) ' ' num2str(harmonics(k)-0.5) ' ' num2str(harmonics(k)+0.5)...
        ' ' num2str(harmonics(k)+1) ' '];
    intensity_string = [intensity_string '0 1 1 0 '];
end
fir_string = [fir_string num2str(sampFreq/2) ']/(' num2str(sampFreq/2) ')'];

%% Convert to numerical arrays
fir_coeffs = str2num(fir_string);
intensity_coeffs = str2num([intensity_string '0]']);

[b,a]=fir2(1000,fir_coeffs,intensity_coeffs);

ecog_filt = filtfilt(b,a,ecog')';


end