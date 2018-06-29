% Radius at 1
% 5 Bright points in Cross shape
% Sphere of 266 equally spaced points
% NR          - No Regularization
% EF          - Broadside Double Layer
% #S          - Numer of sources
% 4S_B#       - Spacing between each layer
% 4S_B1cm_D#  - Spacing between sources

clc
clear 

figure('Name','2 cm spacing','NumberTitle','off')
load 4S_B1cm_D4mm
subplot(1,2,1)
semilogx(freq,AC),title('Acoustic Contrast')
subplot(1,2,2)
semilogx(freq,AE),title('Array Effort')


load 4S_B1cm_D8mm
subplot(1,2,1),hold on
semilogx(freq,AC),title('Acoustic Contrast')
subplot(1,2,2),hold on
semilogx(freq,AE),title('Array Effort')

load 4S_B1cm_D12mm
subplot(1,2,1),hold on
semilogx(freq,AC),title('Acoustic Contrast')
subplot(1,2,2),hold on
semilogx(freq,AE),title('Array Effort')

load 4S_B1cm_D16mm
subplot(1,2,1),hold on
semilogx(freq,AC),title('Acoustic Contrast')
subplot(1,2,2),hold on
semilogx(freq,AE),title('Array Effort')

load 4S_B1cm_D20mm
subplot(1,2,1),hold on
semilogx(freq,AC),title('Acoustic Contrast')
subplot(1,2,2),hold on
semilogx(freq,AE),title('Array Effort')

load 4S_B1cm_D40mm
subplot(1,2,1),hold on
semilogx(freq,AC),title('Acoustic Contrast')
subplot(1,2,2),hold on
semilogx(freq,AE),title('Array Effort')

load 4S_B1cm_D60mm
subplot(1,2,1),hold on
semilogx(freq,AC),title('Acoustic Contrast')
subplot(1,2,2),hold on
semilogx(freq,AE),title('Array Effort')

subplot(1,2,1)
grid
subplot(1,2,2)
grid
legend('4 mm','8 mm','12 mm','16 mm','20 mm','40 mm','60 mm')
