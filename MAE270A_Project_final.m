ts=1/50;
fs = 1/ts; % sample rate in hertz
nfft = 250; % sub-record length nfft/fs = 5 seconds
win = hamming(nfft); % use Hamming data tapering window
fnyq = fs/2;
t = [0:nfft-1]*ts;

% load response data
load random_u1.mat
y11 = y1;
y21 = y2;
u11 = u1;
load random_u2.mat
y12 = y1;
y22 = y2;
u22 = u2;
load random_u3.mat
y13 = y1;
y23 = y2;
u33 = u3;

[Su1u1,f] = cpsd(u11,u11,win,[],nfft,fs,'twosided'); % auto-spectral density of u
[Su2u2,~] = cpsd(u22,u22,win,[],nfft,fs,'twosided'); % cross-spectral density of y and u
[Su3u3,~] = cpsd(u33,u33,win,[],nfft,fs,'twosided'); % cross-spectral density of y and u

[Su2u1,~] = cpsd(u22,u11,win,[],nfft,fs,'twosided'); % auto-spectral density of u
[Su3u1,~] = cpsd(u33,u11,win,[],nfft,fs,'twosided'); % cross-spectral density of y and u
[Su3u2,~] = cpsd(u33,u22,win,[],nfft,fs,'twosided'); % cross-spectral density of y and u

[Sy1u1,~] = cpsd(y11,u11,win,[],nfft,fs,'twosided'); % auto-spectral density of u
[Sy2u1,~] = cpsd(y21,u11,win,[],nfft,fs,'twosided'); % auto-spectral density of u

[Sy1u2,~] = cpsd(y12,u22,win,[],nfft,fs,'twosided'); % auto-spectral density of u
[Sy2u2,~] = cpsd(y22,u22,win,[],nfft,fs,'twosided'); % auto-spectral density of u

[Sy1u3,~] = cpsd(y13,u33,win,[],nfft,fs,'twosided'); % auto-spectral density of u
[Sy2u3,~] = cpsd(y23,u33,win,[],nfft,fs,'twosided'); % auto-spectral density of u

H11 = Sy1u1./Su1u1;
H21 = Sy2u1./Su1u1;
H12 = Sy1u2./Su2u2;
H22 = Sy2u2./Su2u2;
H13 = Sy1u3./Su3u3;
H23 = Sy2u3./Su3u3;

h11 = ifft(H11);
h21 = ifft(H21);
h12 = ifft(H12);
h22 = ifft(H22);
h13 = ifft(H13);
h23 = ifft(H23);

h = [h11 h12 h13; h21 h22 h23];

%% Part1_task1
%Q1
figure(1)
subplot(311)
loglog(f,Su1u1)
title('Auto Spectra')
grid on; axis([0.1 fnyq 1e-5 1e-2]); legend('Su_1u_1'); legend('Location', 'southwest'); xlabel(); ylabel();
subplot(312)
loglog(f,Su2u2)
grid on; axis([0.1 fnyq 1e-5 1e-2]); legend('Su_2u_2'); legend('Location', 'southwest'); xlabel(); ylabel();
subplot(313)
loglog(f,Su3u3)
grid on; axis([0.1 fnyq 1e-5 1e-2]); legend('Su_3u_3'); legend('Location', 'southwest'); xlabel(); ylabel();

figure(2)
subplot(311)                             
loglog(f,abs(Su2u1))
title('Cross Spectra')
grid on; axis([0.1 fnyq 1e-5 1e-2]); legend('Su_2u_1'); legend('Location', 'southwest'); xlabel(); ylabel();
subplot(312)
loglog(f,abs(Su3u1))
grid on; axis([0.1 fnyq 1e-5 1e-2]); legend('Su_3u_1'); legend('Location', 'southwest'); xlabel(); ylabel();
subplot(313)
loglog(f,abs(Su3u2))
grid on; axis([0.1 fnyq 1e-5 1e-2]); legend('Su_3u_2'); legend('Location', 'southwest'); xlabel(); ylabel();

%Q2
%Calculating average magnitudes of spectra
avg_mag_Su1u1 = mean(abs(Su1u1));
avg_mag_Su2u2 = mean(abs(Su2u2));
avg_mag_Su3u3 = mean(abs(Su3u3));
avg_mag_Su2u1 = mean(abs(Su2u1));
avg_mag_Su3u1 = mean(abs(Su3u1));
avg_mag_Su3u2 = mean(abs(Su3u2));

% Display average magnitudes
fprintf('Average Magnitude Su1u1: %f\n', avg_mag_Su1u1);
fprintf('Average Magnitude Su2u2: %f\n', avg_mag_Su2u2);
fprintf('Average Magnitude Su3u3: %f\n', avg_mag_Su3u3);
fprintf('Average Magnitude Su2u1: %f\n', avg_mag_Su2u1);
fprintf('Average Magnitude Su3u1: %f\n', avg_mag_Su3u1);
fprintf('Average Magnitude Su3u2: %f\n', avg_mag_Su3u2);

%Comparison
fprintf('\nComparison of Average Magnitudes:\n');
fprintf('Su2u1 is, on average, smaller than Su1u1: %s\n', ...
    string(avg_mag_Su2u1 < avg_mag_Su1u1));
fprintf('Su3u1 is, on average, smaller than Su1u1: %s\n', ...
    string(avg_mag_Su3u1 < avg_mag_Su1u1));
fprintf('Su3u2 is, on average, smaller than Su2u2: %s\n', ...
    string(avg_mag_Su3u2 < avg_mag_Su2u2));


%Q3
Mu_1 = var(u11);
Mu_2 = var(u22);
Mu_3 = var(u33);

MSu_1u_1 = mean(Su1u1)*50;
MSu_2u_2 = mean(Su2u2)*50;
MSu_3u_3 = mean(Su3u3)*50;

%Q4
figure(3)
subplot(211)
loglog(f,abs(H11))
title('Frequency Response Magnitudes for u_1')
grid on; axis([0.1 fnyq 1e-3 1e2]); legend('H11'); legend('Location', 'southwest'); xlabel(); ylabel();
subplot(212)
loglog(f,abs(H21))
grid on; axis([0.1 fnyq 1e-3 1e2]); legend('H21'); legend('Location', 'southwest'); xlabel(); ylabel();

figure(4)
subplot(211)
loglog(f,abs(H12))
title('Frequency Response Magnitudes for u_2')
grid on; axis([0.1 fnyq 1e-3 1e2]); legend('H12'); legend('Location', 'southwest'); xlabel(); ylabel();
subplot(212)
loglog(f,abs(H22))
grid on; axis([0.1 fnyq 1e-3 1e2]); legend('H22'); legend('Location', 'southwest'); xlabel(); ylabel();

figure(5)
subplot(211)
loglog(f,abs(H13))
title('Frequency Response Magnitudes for u_3')
grid on; axis([0.1 fnyq 1e-3 1e2]); legend('H13'); legend('Location', 'southwest'); xlabel(); ylabel();
subplot(212)
loglog(f,abs(H23))
grid on; axis([0.1 fnyq 1e-3 1e2]); legend('H23'); legend('Location', 'southwest'); xlabel(); ylabel();

figure(6)
subplot(211)
semilogx(f,rad2deg(angle(H11)))
title('Phase Plots for u_1')
grid on; axis([0.1 fnyq -200 200]); legend('H11'); legend('Location', 'southwest'); xlabel(); ylabel();
subplot(212)
semilogx(f,rad2deg(angle(H21)))
grid on; axis([0.1 fnyq -200 200]); legend('H21'); legend('Location', 'southwest'); xlabel(); ylabel();

figure(7)
subplot(211)
semilogx(f,rad2deg(angle(H12)))
title('Phase Plots for u_2')
grid on; axis([0.1 fnyq -200 200]); legend('H12'); legend('Location', 'southwest'); xlabel(); ylabel();
subplot(212)
semilogx(f,rad2deg(angle(H22)))
grid on; axis([0.1 fnyq -200 200]); legend('H22'); legend('Location', 'southwest'); xlabel(); ylabel();

figure(8)
subplot(211)
semilogx(f,rad2deg(angle(H13)))
title('Phase Plots for u_3')
grid on; axis([0.1 fnyq -200 200]); legend('H13'); legend('Location', 'southwest'); xlabel(); ylabel();
subplot(212)
semilogx(f,rad2deg(angle(H23)))
grid on; axis([0.1 fnyq -200 200]); legend('H23'); legend('Location', 'southwest'); xlabel(); ylabel();

%% Part1_task2
figure(9)
subplot(211)
plot(t,h11)
title('Pulse Response for u_1')
grid on; axis([0 3 -2 3]); legend('h11'); xlabel(); ylabel();
subplot(212)
plot(t,h21)
grid on; axis([0 3 -2 3]); legend('h21'); xlabel(); ylabel();

figure(10)
subplot(211)
plot(t,h12)
title('Pulse Response for u_2')
grid on; axis([0 3 -2 3]); legend('h12'); xlabel(); ylabel();
subplot(212)
plot(t,h22)
grid on; axis([0 3 -2 3]); legend('h22'); xlabel(); ylabel();

figure(11)
subplot(211)
plot(t,h13)
title('Pulse Response for u_3')
grid on; axis([0 3 -2 3]); legend('h13'); xlabel(); ylabel();
subplot(212)
plot(t,h23)
grid on; axis([0 3 -2 3]); legend('h23'); xlabel(); ylabel();

%% Part2_task1
%Q1
n = 25;
M_n = [];

for i = 1:n
    tmp=[];
    for j = 1:n
        hk = [h11(j+i) h12(j+i) h13(j+i); h21(j+i) h22(j+i) h23(j+i)];
        tmp=[tmp hk];
    end
    M_n = [M_n;tmp];
end

figure(12)
semilogy(svd(M_n), '*')
grid on; axis([0 20 1e-3 1e2]); legend(compose('n= %d',n))

[U,S,V] = svd(M_n);

n_s = 8;

S_n = S(1:n_s, 1:n_s);
U_1 = U(:,1:n_s);
V_1 = V(:,1:n_s);
M_w = U_1*S_n*V_1';
L = U_1;
R = S_n*V_1';

Mbar = [];
for i = 1:n
    tmp=[];
    for j = 1:n
        tmp=[tmp [h11(j+i+1) h12(j+i+1) h13(j+i+1); h21(j+i+1) h22(j+i+1) h23(j+i+1)]];
    end
    Mbar = [Mbar;tmp];
end

[Ubar,Sbar,Vbar] = svd(Mbar);
S_nbar = Sbar(1:n_s, 1:n_s);
U_1bar = Ubar(:,1:n_s);
V_1bar = Vbar(:,1:n_s);
M_wbar = U_1bar*S_nbar*V_1bar';

B = R(:,1:3);
C = L(1:2,:);

A = U_1.'*Mbar*V_1/S_n;
max(abs(eig(A)))

D = [0 0 0; 0 0 0];

%Q2
hI = zeros(2,3*3/ts+1);
hI11(1) = 0;
hI21(1) = 0;
hI12(1) = 0;
hI22(1) = 0;
hI13(1) = 0;
hI23(1) = 0;
for i = 2:3/ts+1
    hk = C*A^(i-1)*B;
    hI(1,3*i-2) = hk (1,1);
    hI(2,3*i-2) = hk(2,1);
    hI(1,3*i-1) = hk(1,2);
    hI(1,3*i-1) = hk(2,2);
    hI(1,3*i) = hk(1,3);
    hI(1,3*i) = hk(2,3);

    hI11(i) = hk(1,1);
    hI21(i) = hk(2,1);
    hI12(i) = hk(1,2);
    hI22(i) = hk(2,2);
    hI13(i) = hk(1,3);
    hI23(i) = hk(2,3);
end

figure(13)
subplot(211)
plot(0:ts:3,hI11)
title('Impulse Response for u_1 (calculated for n_s)')
grid on; axis([0 3 -2 3]); legend('h11'); xlabel(); ylabel();
subplot(212)
plot(0:ts:3,hI21)
grid on; axis([0 3 -2 3]); legend('h21'); xlabel(); ylabel();

figure(14)
subplot(211)
plot(0:ts:3,hI12)
title('Impulse Response for u_2 (calculated for n_s)')
grid on; axis([0 3 -2 3]); legend('h12'); xlabel(); ylabel();
subplot(212)
plot(0:ts:3,hI22)
grid on; axis([0 3 -2 3]); legend('h22'); xlabel(); ylabel();

figure(15)
subplot(211)
plot(0:ts:3,hI13)
title('Impulse Response for u_3 (calculated for n_s)')
grid on; axis([0 3 -2 3]); legend('h13'); xlabel(); ylabel();
subplot(212)
plot(0:ts:3,hI23)
grid on; axis([0 3 -2 3]); legend('h23'); xlabel(); ylabel();

%Q3
w = [0.1:0.1:25];
HI11 = zeros(250,1);
HI21 = zeros(250,1);
HI12 = zeros(250,1);
HI22 = zeros(250,1);
HI13 = zeros(250,1);
HI23 = zeros(250,1);
for i = 1:250
    H = C/(exp(sqrt(-1)*ts*w(i)*(2*pi))*eye(n_s)-A)*B+D;
    HI11(i) = H(1,1);
    HI21(i) = H(2,1);
    HI12(i) = H(1,2);
    HI22(i) = H(2,2);
    HI13(i) = H(1,3);
    HI23(i) = H(2,3);
end

figure(16)
subplot(211)
semilogx(w,rad2deg(angle(HI11)))
title('Phase Plots for u_1 (calculated for n_s)')
grid on; axis([0.1 fnyq -200 200]); legend('H11'); legend('Location', 'southwest');
subplot(212)
semilogx(w,rad2deg(angle(HI21)))
grid on; axis([0.1 fnyq -200 200]); legend('H21'); legend('Location', 'southwest');
% 
figure(17)
subplot(211)
loglog(w,abs(HI11))
title('Frequency Response Magnitudes for u_1 (calculated for n_s)')
grid on; axis([0.1 fnyq 1e-3 1e2]); legend('H11'); legend('Location', 'southwest');
subplot(212)
loglog(w,abs(HI21))
grid on; axis([0.1 fnyq 1e-3 1e2]); legend('H21'); legend('Location', 'southwest');

figure(18)
subplot(211)
semilogx(w,rad2deg(angle(HI12)))
title('Phase Plots for u_2 (calculated for n_s)')
grid on; axis([0.1 fnyq -200 200]); legend('H12'); legend('Location', 'southwest');
subplot(212)
semilogx(w,rad2deg(angle(HI22)))
grid on; axis([0.1 fnyq -200 200]); legend('H22'); legend('Location', 'southwest');

figure(19)
subplot(211)
loglog(w,abs(HI12))
title('Frequency Response Magnitudes for u_2 (calculated for n_s)')
grid on; axis([0.1 fnyq 1e-3 1e2]); legend('H12'); legend('Location', 'southwest');
subplot(212)
loglog(w,abs(HI22))
grid on; axis([0.1 fnyq 1e-3 1e2]); legend('H22'); legend('Location', 'southwest');

figure(20)
subplot(211)
semilogx(w,rad2deg(angle(HI13)))
title('Phase Plots for u_3 (calculated for n_s)')
grid on; axis([0.1 fnyq -200 200]); legend('H13'); legend('Location', 'southwest');
subplot(212)
semilogx(w,rad2deg(angle(HI23)))
% grid on; axis([0.1 fnyq -200 200]); legend('H23'); legend('Location', 'southwest');

figure(21)
subplot(211)
loglog(w,abs(HI13))
title('Frequency Response Magnitudes for u_3 (calculated for n_s)')
grid on; axis([0.1 fnyq 1e-3 1e2]); legend('H13'); legend('Location', 'southwest');
subplot(212)
loglog(w,abs(HI23))
grid on; axis([0.1 fnyq 1e-3 1e2]); legend('H23'); legend('Location', 'southwest');

%% Part2_task2
%Q1&2
figure(22)
plot(real(eig(A)),imag(eig(A)),'x')
hold on
fplot(@(t) sin(t), @(t) cos(t));
title(['Eigenvalues of the system for n_s = ',num2str(n_s)])
grid on; axis square; axis([-1.5 1.5 -1.5 1.5]);

%Q3&4
B11 = B(:,1);
C11 = C(1,:);

B21 = B(:,1);
C21 = C(2,:);

B12 = B(:,2);
C12 = C(1,:);

B22 = B(:,2);
C22 = C(2,:);

B13 = B(:,3);
C13 = C(1,:);

B23 = B(:,3);
C23 = C(2,:);

transmission_zeros11 = tzero(A,B11,-C11,0);
transmission_zeros21 = tzero(A,B21,-C21,0);
transmission_zeros12 = tzero(A,B12,-C12,0);
transmission_zeros22 = tzero(A,B22,-C22,0);
transmission_zeros13 = tzero(A,B13,-C13,0);
transmission_zeros23 = tzero(A,B23,-C23,0);

figure(23)
plot(real(transmission_zeros11), imag(transmission_zeros11),'o')
hold on
plot(real(eig(A)),imag(eig(A)),'x')
fplot(@(t) sin(t), @(t) cos(t));
title(['Transmission zeros for SISO: u_1 and y_1 for n_s = ',num2str(n_s)])
grid on; axis square; axis([-1.5 1.5 -1.5 1.5]);
hold off

figure(24)
plot(real(transmission_zeros21), imag(transmission_zeros21),'o')
hold on
plot(real(eig(A)),imag(eig(A)),'x')
fplot(@(t) sin(t), @(t) cos(t));
title(['Transmission zeros for SISO: u_1 and y_2 for n_s = ',num2str(n_s)])
grid on; axis square; axis([-1.5 1.5 -1.5 1.5]);
hold off

figure(25)
plot(real(transmission_zeros12), imag(transmission_zeros12),'o')
hold on
plot(real(eig(A)),imag(eig(A)),'x')
fplot(@(t) sin(t), @(t) cos(t));
title(['Transmission zeros for SISO: u_2 and y_1 for n_s = ',num2str(n_s)])
grid on; axis square; axis([-1.5 1.5 -1.5 1.5]);
hold off

figure(26)
plot(real(transmission_zeros22), imag(transmission_zeros22),'o')
hold on
plot(real(eig(A)),imag(eig(A)),'x')
fplot(@(t) sin(t), @(t) cos(t));
title(['Transmission zeros for SISO: u_2 and y_2 for n_s = ',num2str(n_s)])
grid on; axis square; axis([-1.5 1.5 -1.5 1.5]);
hold off

figure(27)
plot(real(transmission_zeros13), imag(transmission_zeros13),'o')
hold on
plot(real(eig(A)),imag(eig(A)),'x')
fplot(@(t) sin(t), @(t) cos(t));
title(['Transmission zeros for SISO: u_3 and y_1 for n_s = ',num2str(n_s)])
grid on; axis square; axis([-1.5 1.5 -1.5 1.5]);
hold off

figure(28)
plot(real(transmission_zeros23), imag(transmission_zeros23),'o')
hold on
plot(real(eig(A)),imag(eig(A)),'x')
fplot(@(t) sin(t), @(t) cos(t));
title(['Transmission zeros for SISO: u_3 and y_2 for n_s = ',num2str(n_s)])
grid on; axis square; axis([-1.5 1.5 -1.5 1.5]);
hold off


