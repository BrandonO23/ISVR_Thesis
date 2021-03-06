clc
clear

SampleRate = 44100;
fstart = 500;
fend = SampleRate./2;
RBW = 50;
Navg = 15;
f = fstart:RBW:fend;

Device = 'USB Sound Device';
Driver = 'CoreAudio';


%       'Device', Device, ...
%       'Driver', Driver, ...


aDR = audioDeviceReader(...
      'SampleRate', SampleRate, ...
      'BitDepth', '16-bit integer', ...
      'Device', Device, ...
      'Driver', Driver, ...
      'ChannelMappingSource', 'Property', ...
      'ChannelMapping', 1);

 aDW = audioDeviceWriter(...
      'SampleRate', SampleRate, ...
      'Device', Device, ...
      'Driver', Driver, ...
      'BitDepth', '16-bit integer', ...
      'ChannelMappingSource', 'Property', ...
      'ChannelMapping', 1);
  
SamplesPerFrame = 2^10;
sineSource  = audioOscillator(...
      'Frequency', fstart, ...
      'SignalType', 'sine', ...
      'SampleRate', SampleRate, ...
      'SamplesPerFrame', SamplesPerFrame);
 


% 
% scope = dsp.SpectrumAnalyzer(...
%       'Method', 'Filter bank', ...
%       'SampleRate', SampleRate, ...
%       'RBWSource', 'Property', 'RBW', RBW, ...
%       'SpectralAverages', Navg, ...
%       'FrequencySpan', 'Start and stop frequencies',...
%       'StartFrequency', 0, ...
%       'StopFrequency', SampleRate/2, ...
%       'ReducePlotRate', false, ...
%       'PlotAsTwoSidedSpectrum', false, ...
%       'FrequencyScale', 'Log', ...
%       'PlotMaxHoldTrace', true, ...
%       'ShowLegend', true, ...
%       'YLimits', [-110 20],...
%       'YLabel', 'Power', ...
%       'Title', 'Audio Device Frequency Response');
  
%%
% tic
% while toc < 5
%     x = sineSource()./2;
%     aDW(x);
%     y = aDR();
% %     scope(y);
% end

count = 1;
readerDrops = 0;
writerDrops = 0;

it = 1;
while true
    if count == Navg
        newFreq = sineSource.Frequency + RBW;
        if newFreq > fend
            break
        end
        sineSource.Frequency = newFreq;
        count = 1;
    end
    x = sineSource()./2;
    writerUnderruns = aDW(x);
    [y,readerOverruns] = aDR();
    tt(it) = thd(y);
    readerDrops = readerDrops + readerOverruns;
    writerDrops = writerDrops + writerUnderruns;
    count = count + 1;
    it = it + 1;
end

release(aDR);
release(aDW);

plot(f,tt);

