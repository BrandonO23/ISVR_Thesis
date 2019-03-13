clc 
clear


SampleRate = 44100;
SamplesPerFrame = 2^13;

aDR = audioDeviceReader(...
      'SampleRate', SampleRate, ...
      'BitDepth', '16-bit integer', ...
      'ChannelMappingSource', 'Property', ...
      'SamplesPerFrame', SamplesPerFrame,...
      'ChannelMapping', 1);

 aDW = audioDeviceWriter(...
      'SampleRate', SampleRate, ...
      'BitDepth', '16-bit integer', ...
      'ChannelMappingSource', 'Property', ...
      'ChannelMapping', 1);
  
sineSource  = audioOscillator(...
      'Frequency', 1000, ...
      'SignalType', 'sine', ...
      'SampleRate', SampleRate, ...
      'SamplesPerFrame', SamplesPerFrame);
  
 
len = .2;
fre = 10.^(3:.01:4.3);
harmonicdis = zeros(1,length(fre)-1);

i = 1;
tic
while i <= length(fre)
    s = sineSource();
    aDW(s);
    if toc >= len
        y = aDR();
        harmonicdis(i) = thd(y);
        i = i + 1;
        sineSource.Frequency = fre(i);
        tic
    end
end

plot(fre,harmonicdis);


release(sineSource);
release(aDW);
 
 