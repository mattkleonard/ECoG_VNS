function [vns_pulse] = generate_ideal_vns_pulse(sample_rate,varargin)
% generates an theoretical representation of the intensity of the VNS pulse
%
% input - sample rate of the ekg channel
% variable inputs:
% 1 - pulse peak time - the duration (in seconds) of peak VNS power
% 2 - pulse onset/offset time - the duration (seconds) of ramp to peak 
% 
%
% output - vns_pulse (an idealized trace of the power of a single VNS
% stimulation pulse given the characteristic duty cycle and sample rate).

%% Pulse shape parameters
pulse_onset_dur = 2;       % time from start to peak (seconds)
pulse_flat_dur = 30;       % duration of peak stimulation power(seconds)
pulse_offset_dur = 2;      % duration of offset of stimulation (seconds)

if length(varargin)>1
    pulse_onset_dur = varargin{1};
end

if length(varargin)>0
    pulse_flat_dur = varargin{2};
end

if length(varargin)>2
    pulse_offset_dur = varargin{3};
end
pre_pulse_dur = 5;         % duration of silence pre stimulation (seconds)
post_pulse_dur = 5;        % duration of silence post stimulation (seconds)

% make pulse (from -1 to 1 (on/off)
vns_pulse = [-ones(1, round(sample_rate*pre_pulse_dur)) ...
    linspace(-1,1,round(sample_rate*pulse_onset_dur)) ...
     ones(1, round(sample_rate*pulse_flat_dur)) ...
     linspace(1,-1,round(sample_rate*pulse_offset_dur))...
     -ones(1,round(sample_rate*post_pulse_dur))];
 
end