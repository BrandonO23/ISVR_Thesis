clc
clear

SampleRate = 44100;
% Device = 'USB Sound Device';
% % Device = 'Built-in Microph';
% Driver = 'CoreAudio';
  
Length = 20;
f0 = 100;
fe = 22000;
t = 0:1/SampleRate:Length;
y = chirp(t,f0,Length,fe,'logarithmic');
sinn = sin(2*pi*1000*t);

%%
Recorder = audiorecorder(SampleRate,16,1,3);
P = audioplayer(y./2,SampleRate,16,2);

disp('Recording Started')
play(P)
recordblocking(Recorder,Length);
disp('Recording finished');


out = getaudiodata(Recorder);

%
figure(1)
spectrogram(out,2^10,[],[],SampleRate,'yaxis');

%%
figure(2)
thd(out,SampleRate)
