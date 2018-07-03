clc
clear

SampleRate = 44100;
Device = 'USB Sound Device';
% Device = 'Built-in Microph';
Driver = 'CoreAudio';

aDR = audioDeviceReader(...
      'SampleRate', SampleRate, ...
      'Device', Device, ...
      'Driver', Driver, ...
      'BitDepth', '16-bit integer', ...
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
      'Frequency', 500, ...
      'SignalType', 'sine', ...
      'SampleRate', SampleRate, ...
      'SamplesPerFrame', SamplesPerFrame);
  

RBW = 50;
Navg = 2;



scope = dsp.SpectrumAnalyzer(...
      'Method', 'Filter bank', ...
      'SampleRate', SampleRate, ...
      'RBWSource', 'Property', 'RBW', RBW, ...
      'SpectralAverages', Navg, ...
      'FrequencySpan', 'Start and stop frequencies',...
      'StartFrequency', 0, ...
      'StopFrequency', SampleRate/2, ...
      'ReducePlotRate', false, ...
      'PlotAsTwoSidedSpectrum', false, ...
      'FrequencyScale', 'Log', ...
      'PlotMaxHoldTrace', true, ...
      'ShowLegend', true, ...
      'YLimits', [-110 20],...
      'YLabel', 'Power', ...
      'Title', 'Audio Device Frequency Response');
  
%%
tic
while toc < 5
    x = sineSource();
    aDW(x);
    y = aDR();
    scope(y);
end


%%
count = 1;
readerDrops = 0;
writerDrops = 0;

while true
    if count == Navg
        scope(y);
        newFreq = sineSource.Frequency + RBW;
        if newFreq > SampleRate/2
%         if newFreq > 12000
            break
        end
        sineSource.Frequency = newFreq;
        count = 1;
    end
    x = sineSource();
    writerUnderruns = aDW(x);
    [y,readerOverruns] = aDR();
    readerDrops = readerDrops + readerOverruns;
    writerDrops = writerDrops + writerUnderruns;
    count = count + 1;
end

release(aDR);
release(aDW);
release(scope);

%%
data = getSpectrumData(scope);
freqVector   = data.FrequencyVector{1};
freqResponse = data.MaxHoldTrace{1};

semilogx(freqVector,20*log10(freqResponse));
xlabel('Frequency (Hz)');
ylabel('Vrms');
title('Audio Device Frequency Response');
grid on
xlim([1000 20000])
ylim([-80 -30])