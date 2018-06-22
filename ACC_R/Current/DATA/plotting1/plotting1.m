% Radius at 1
% 5 Bright points in Cross shape
% Sphere of 266 equally spaced points
% NR          - No Regularization
% EF          - End Fire
% NR_EF_#     - Numer of sources
% NR_EF_2_#cm - equal spacing between each

clc
clear 
load freq

figure('Name','2 cm spacing','NumberTitle','off')
load NR_EF_2_2cm
subplot(1,2,1)
semilogx(freq,AC),title('Acoustic Contrast')
subplot(1,2,2)
semilogx(freq,AE),title('Array Effort')


load NR_EF_3_2cm
subplot(1,2,1),hold on
semilogx(freq,AC),title('Acoustic Contrast')
subplot(1,2,2),hold on
semilogx(freq,AE),title('Array Effort')

load NR_EF_4_2cm
subplot(1,2,1),hold on
semilogx(freq,AC),title('Acoustic Contrast')
subplot(1,2,2),hold on
semilogx(freq,AE),title('Array Effort')

subplot(1,2,1)
grid
subplot(1,2,2)
grid
legend('2 Sources','3 Source','4 Sources')


%%
figure('Name','4 cm spacing','NumberTitle','off')
load NR_EF_2_4cm
subplot(1,2,1)
semilogx(freq,AC),title('Acoustic Contrast')
subplot(1,2,2)
semilogx(freq,AE),title('Array Effort')


% figure(2)
load NR_EF_3_4cm
subplot(1,2,1),hold on
semilogx(freq,AC),title('Acoustic Contrast')
subplot(1,2,2),hold on
semilogx(freq,AE),title('Array Effort')

% figure(2)
load NR_EF_4_4cm
subplot(1,2,1),hold on
semilogx(freq,AC),title('Acoustic Contrast')
subplot(1,2,2),hold on
semilogx(freq,AE),title('Array Effort')

subplot(1,2,1)
grid
subplot(1,2,2)
grid
legend('2 Sources','3 Source','4 Sources')

%%
figure('Name','6 cm spacing','NumberTitle','off')
load NR_EF_2_6cm
subplot(1,2,1)
semilogx(freq,AC),title('Acoustic Contrast')
subplot(1,2,2)
semilogx(freq,AE),title('Array Effort')


load NR_EF_3_6cm
subplot(1,2,1),hold on
semilogx(freq,AC),title('Acoustic Contrast')
subplot(1,2,2),hold on
semilogx(freq,AE),title('Array Effort')

load NR_EF_4_6cm
subplot(1,2,1),hold on
semilogx(freq,AC),title('Acoustic Contrast')
subplot(1,2,2),hold on
semilogx(freq,AE),title('Array Effort')

subplot(1,2,1)
grid
subplot(1,2,2)
grid
legend('2 Sources','3 Source','4 Sources')

%%
