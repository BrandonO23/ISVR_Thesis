load M8_F100_AV25_800_12k

Noct = 3;
Z = iosr.dsp.smoothSpectrum(freqResponse,freqVector,Noct);


semilogx(freqVector,Z);
xlabel('Frequency (Hz)');
ylabel('Power (dBm)');
title('Audio Device Frequency Response');

xlim([1000 20000])
grid 