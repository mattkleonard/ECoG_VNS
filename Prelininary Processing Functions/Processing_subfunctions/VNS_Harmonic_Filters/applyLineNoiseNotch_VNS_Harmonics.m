function ecog_filt = applyLineNoiseNotch_VNS_Harmonics(ecog,sampFreq,varargin)
%% Remove the stimulation frequency and its harmonics from the nchans x time ecog matrix
%
% Inputs
% ecog - channels x time Ecog data matrix
% sampFreq - the sampling rate of the data
% Variable inputs:
% 1 - notchFrequency the frequency to be notched (default 25hz)

% This function automatically generates a single FIR filter that will include all of the
% harmonics and then applies this to the data.
notchFreq=25;
if length(varargin)>0
    if ~isempty(varargin{1})
        notchFreq = varargin{1};
    end
end

max_harmonic = 200; % set as a constant because we only look at ecog bands from [4, 200] hz
harmonics = notchFreq:notchFreq:max_harmonic;

%% Generate Notch Labels
fir_string = '[0 ';
intensity_string = '[1 ';
for k = 1:length(harmonics)
    fir_string = [fir_string num2str(harmonics(k)-1) ' ' num2str(harmonics(k)-0.5) ' ' num2str(harmonics(k)+0.5)...
        ' ' num2str(harmonics(k)+1) ' '];
    intensity_string = [intensity_string '1 0 0 1 '];
end
fir_string = [fir_string num2str(sampFreq/2) ']/(' num2str(sampFreq/2) ')'];
intensity_string = [intensity_string '1]'];

%% Convert to numerical arrays
fir_coeffs = str2num(fir_string);
intensity_coeffs = str2num(intensity_string);

[b,a]=fir2(1000,fir_coeffs,intensity_coeffs);

ecog_filt = filtfilt(b,a,ecog')';




%% Old Function (massively less efficient)
% % % % notchFreq=25;
% % % % while notchFreq<200 %sampFreq/2
% % % %     %fprintf(['Notch Filters: ' int2str(notchFreq) 'Hz.. \n'])
% % % %     [b,a]=fir2(1000,[0 notchFreq-1 notchFreq-.5 notchFreq+.5 notchFreq+1 sampFreq/2]/(sampFreq/2),[1 1 0 0 1 1 ]);
% % % % 	ecog=filtfilt(b,a,ecog')';
% % % %     notchFreq=notchFreq+25;
% % % % end
% % % % %fprintf('\nNotch Filters Applied\n')
