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

%%
Recorder = audiorecorder(SampleRate,16,1,3);
P = audioplayer(y./5,SampleRate,16,2);

disp('Recording Started')
play(P)
recordblocking(Recorder,Length);
disp('Recording finished');


out = getaudiodata(Recorder);

%
figure(1)
spectrogram(out,2^10,[],[],SampleRate,'yaxis');

%%
t = 0:1/SampleRate:1;


f = 10.^(3:.006:4);

for i = 1:length(f)
    sinn = sin(2*pi*f(i)*t);
    Recorder = audiorecorder(SampleRate,16,1,3);
    P = audioplayer(sinn./5,SampleRate,16,2);

%     disp('Recording Started')
    play(P)
    recordblocking(Recorder,1);
%   disp('Recording finished');
    out = getaudiodata(Recorder);
    tt(i) = thd(out,SampleRate);
%   pwelch(out,[],[],[],SampleRate)
end

%%
semilogx(f,100*10.^(tt./10))
grid on

%%
    sinn = sin(2*pi*f(i)*t);
    Recorder = audiorecorder(SampleRate,16,1,3);
    P = audioplayer(sinn./2,SampleRate,16,2);

    play(P)
    recordblocking(Recorder,1);
    out = getaudiodata(Recorder);
    thd(outSampleRate);



