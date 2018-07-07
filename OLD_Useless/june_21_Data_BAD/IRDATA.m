clc
clear
%load measured_ir_data1
%load measured_ir_data2
load measured_ir_data3

f = measurementData.MagnitudeResponse.Frequency;
y = measurementData.MagnitudeResponse.PowerDb;


Noct = 3;
Z = iosr.dsp.smoothSpectrum(y,f,Noct);

semilogx(f,Z)
grid on
xlim([1000 22050])