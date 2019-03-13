load measured_ir_data.mat

fr = measurementData.MagnitudeResponse(3,1).PowerDb;
f  = measurementData.MagnitudeResponse(3,1).Frequency;

%%
Noct = 12;
Z = iosr.dsp.smoothSpectrum(fr,f,Noct);
%%

subplot(1,2,1)
semilogx(f,Z);
xlabel('Frequency (Hz)');
ylabel('Power (dBm)');
title(' MEMS Frequency Response - 12^{th} octave smoothing');
xlim([1000 20000])
ylim([-40 20])
grid 
subplot(1,2,2)
semilogx(f,fr)
xlabel('Frequency (Hz)');
ylabel('Power (dBm)');
title('MEMS Frequency Response');
xlim([1000 20000])
grid 
ylim([-40 20])