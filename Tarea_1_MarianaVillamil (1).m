% Student name:
stname = 'Mariana Carolina Villamil Sastre';
% Student codigo:
codigo = '201413371';
disp(['This is the work of ' stname ' with codigo ' codigo])
disp('')
disp('Homework 1 - DATA ANALYSIS FOR GEOSCIENCES')
%%
%Plotting time series
data=dlmread('water.txt')
time=data(:,1)
precipitation=data(:,2)
airtemp=data(:,3)
watertemp=data(:,4)
salinity=data(:,5)
turbidity=data(:,6)
clorofila=data(:,7)

subplot(3,2,1);
plot(time,precipitation)
title('Plot of Precipitation Vs Time','FontSize', 8)
ylabel('Precipitation (inches)','FontSize', 8)
xlabel('Time (days)','FontSize', 8)


subplot(3,2,2); 
plot(time,airtemp)
title('Plot of Air Temperature Vs Time','FontSize', 8)
ylabel('Air temperature (°C)','FontSize', 8)
xlabel('Time (days)','FontSize', 8)

subplot(3,2,3); 
plot(time,watertemp)
title('Plot of Water Temperature Vs Time','FontSize', 8)
ylabel('Water temperature (°C)','FontSize', 8)
xlabel('Time (days)','FontSize', 8)

subplot(3,2,4); 
plot(time,salinity)
title('Plot of Salinity Vs Time','FontSize', 8)
ylabel('Salinity (per thousands)','FontSize', 8)
xlabel('Time (days)','FontSize', 8)

subplot(3,2,5); 
plot(time,turbidity)
title('Plot of Turbidity Vs Time','FontSize', 8)
ylabel('Turbidity','FontSize', 8)
xlabel('Time (days)','FontSize', 8)

subplot(3,2,6); 
plot(time,clorofila)
title('Plot of Chlorophyll Vs Time','FontSize', 8)
ylabel('Chlorophyll (\mu g/l)','FontSize', 8)
xlabel('Time (days)','FontSize', 8)

%% Fourier transform
%In this step I put the respective value of variables

Ts = 1;                                     % Sampling Interval (days)
Fs = 1/Ts;                                      % Sampling Frequency (samples/day)
Fn = Fs/2;                                      % Nyquist Frequency   
n=length(precipitation);


%In this space I calculate the Fourier Transform of each variable
                        
NFFT1= 2^nextpow2(n);
precipitation1=detrend(precipitation)
airtemp1=detrend(airtemp)
watertemp1=detrend(watertemp)
salinity1=detrend(salinity)
turbidity1=detrend(turbidity)
clorofila1=detrend(clorofila)

fF1 = Fs/2*linspace(0,1,NFFT1/2+1);
SF1 = fft(precipitation1,NFFT1)/length(precipitation1); 
SF2 = fft(airtemp1,NFFT1)/length(precipitation1); 
SF3 = fft(watertemp1,NFFT1)/length(precipitation1); 
SF4 = fft(salinity1,NFFT1)/length(precipitation1); 
SF5 = fft(turbidity1,NFFT1)/length(precipitation1); 
SF6 = fft(clorofila1,NFFT1)/length(precipitation1); 

subplot(3,2,1);
plot(fF1,2*abs(SF1(1:NFFT1/2+1)))
ylim([0 0.08])
title('FFT for Precipitation variable','FontSize', 8)
ylabel('Amplitude','FontSize', 8)
xlabel('Frequency (1/days)','FontSize', 8)


subplot(3,2,2); 
plot(fF1,2*abs(SF2(1:NFFT1/2+1)))
ylim([0 1])
title('FFT for Air temperature variable','FontSize', 8)
ylabel('Amplitude','FontSize', 8)
xlabel('Frequency (1/days)','FontSize', 8)

subplot(3,2,3); 
plot(fF1,2*abs(SF3(1:NFFT1/2+1)))
ylim([0 1])
title('FFT for Water temperature variable','FontSize', 8)
ylabel('Amplitude','FontSize', 8)
xlabel('Frequency (1/days)','FontSize', 8)

subplot(3,2,4); 
plot(fF1,2*abs(SF4(1:NFFT1/2+1)))
ylim([0 1])
title('FFT for Salinity variable','FontSize', 8)
ylabel('Amplitude','FontSize', 8)
xlabel('Frequency (1/days)','FontSize', 8)

subplot(3,2,5); 
plot(fF1,2*abs(SF5(1:NFFT1/2+1)))
ylim([0 1])
title('FFT for Turbidity variable','FontSize', 8)
ylabel('Amplitude','FontSize', 8)
xlabel('Frequency (1/days)','FontSize', 8)

subplot(3,2,6); 
plot(fF1,2*abs(SF6(1:NFFT1/2+1)));
ylim([0 1])
title('FFT for Chlorophyll variable','FontSize', 8)
ylabel('Amplitude','FontSize', 8)
xlabel('Frequency (1/days)','FontSize', 8)

%%
%%Filters
%Filter variables for 1 year
[b,a]=butter(6,((1/365)+0.001)/(Fs/2),'low');
precipitationfilt = filtfilt(b,a,precipitation1);

airtempfilt = filtfilt(b,a,airtemp1);
watertempfilt = filtfilt(b,a,watertemp1);
salinityfilt = filtfilt(b,a,salinity1);
turbidityfilt = filtfilt(b,a,turbidity1);
clorofilafilt = filtfilt(b,a,clorofila1);


subplot(3,2,1);
plot(time,precipitationfilt)
title('Filter of 1 year period of Prep Vs Time','FontSize', 8)
ylabel('Precipitation (inches)','FontSize', 8)
xlabel('Time (days)','FontSize', 8)


subplot(3,2,2); 
plot(time,airtempfilt)
title('Filter of 1 year period of AirTemp Vs Time','FontSize', 8)
ylabel('Air temperature (°C)','FontSize', 8)
xlabel('Time (days)','FontSize', 8)

subplot(3,2,3); 
plot(time,watertempfilt)
title('Filter of 1 year period of WaterTemp Vs Time','FontSize', 8)
ylabel('Water temperature (°C)','FontSize', 8)
xlabel('Time (days)','FontSize', 8)

subplot(3,2,4); 
plot(time,salinityfilt)
title('Filter of 1 year period of Sal Vs Time','FontSize', 8)
ylabel('Salinity (per thousands)','FontSize', 8)
xlabel('Time (days)','FontSize', 8)

subplot(3,2,5); 
plot(time,turbidityfilt)
title('Filter of 1 year period of Turb Vs Time','FontSize', 8)
ylabel('Turbidity','FontSize', 8)
xlabel('Time (days)','FontSize', 8)

subplot(3,2,6); 
plot(time,clorofilafilt)
title('Filter of 1 year period of Chlo Vs Time','FontSize', 8)
ylabel('Chlorophyll (\mu g/l)','FontSize', 8)
xlabel('Time (days)','FontSize', 8)

%Filter variables for 5 days
N = 6; F0 = 0.2; BW = 0.05;
peakspec = fdesign.peak('N,F0,BW',N,F0,BW);
peakfilt = design(peakspec,'SystemObject',true);

precipitationfilt5=peakfilt(precipitation1);

airtempfilt5=peakfilt(airtemp1);
watertempfilt5=peakfilt(watertemp1);
salinityfilt5=peakfilt(salinity1);
turbidityfilt5=peakfilt(turbidity1);
clorofilafilt5=peakfilt(clorofila1);

subplot(3,2,1);
plot(time,precipitationfilt5)
title('Filter of 5 days period of Prep Vs Time','FontSize', 8)
ylabel('Precipitation (inches)','FontSize', 8)
xlabel('Time (days)','FontSize', 8)


subplot(3,2,2); 
plot(time,airtempfilt5)
title('Filter of 5 days period of AirTemp Vs Time','FontSize', 8)
ylabel('Air temperature (°C)','FontSize', 8)
xlabel('Time (days)','FontSize', 8)

subplot(3,2,3); 
plot(time,watertempfilt5)
title('Filter of 5 days period of WaterTemp Vs Time','FontSize', 8)
ylabel('Water temperature (°C)','FontSize', 8)
xlabel('Time (days)','FontSize', 8)

subplot(3,2,4); 
plot(time,salinityfilt5)
title('Filter of 5 days period of Sal Vs Time','FontSize', 8)
ylabel('Salinity (per thousands)','FontSize', 8)
xlabel('Time (days)','FontSize', 8)

subplot(3,2,5); 
plot(time,turbidityfilt5)
title('Filter of 5 days period of Turb Vs Time','FontSize', 8)
ylabel('Turbidity','FontSize', 8)
xlabel('Time (days)','FontSize', 8)

subplot(3,2,6); 
plot(time,clorofilafilt5)
title('Filter of 5 days period of Chlo Vs Time','FontSize', 8)
ylabel('Chlorophyll (\mu g/l)','FontSize', 8)
xlabel('Time (days)','FontSize', 8)

disp('It is possible to see that Air temprature and Chlorophyll variables are anti correlated, also precipitation and salinity and the last one, Water temperature and turbidity ')

%%
%Correlation coefficient
vector=[];
vector=[precipitationfilt airtempfilt watertempfilt salinityfilt turbidityfilt clorofilafilt]

correlation=corrcoef(vector)

vector5=[];
vector5=[precipitationfilt5 airtempfilt5 watertempfilt5 salinityfilt5 turbidityfilt5 clorofilafilt5]

correlation5=corrcoef(vector5)
%%
%Chlorophyll autocorrelation

autocorr(clorofila1)
title('Chlorophyll autocorrelation')
xlabel('lags')
ylabel('Autocorrelation')
%%
%Delay air and temperature

d1 = finddelay(watertempfilt,airtempfilt);
[c,lags] = xcorr(watertempfilt,airtempfilt);
d2 = -(lags(c == max(c)));

%%
%Smoothing time series

suavizado= smoothdata(watertemp1,'gaussian',200);
plot(time,watertemp1,'-o',time,suavizado,'-x');
ylabel('Water temperature (°C)','FontSize', 12);
xlabel('Time (days)','FontSize', 12);
legend('Original Data','Smoothed Data');

disp('This method uses a Gaussian sequence to approximate and smooth the time series, using a window of size 200, this is decision of the coder in how much want to smooth the data, I use 200 because 250 softens the data too much and does not allow to identify any relationship with the original data')

%%
%Time frequency transform
addpath(genpath('tftb-0.2'))
winlength = 61; % Window length in samples
[TFstft,~,Fstft]=tfrstft(turbidity1,1:length(turbidity1),NFFT1,hann(winlength));
figure; imagesc(time,Fstft(1:end/2)*Fs,abs(TFstft(1:end/2,:))); set(gca,'YDir','normal')
colorbar;
xlabel('Time (days)'); ylabel('Frequency (days^{-1})')
title('Turbidity variable analyzed with the STFT')

plot(fF1,2*abs(SF5(1:NFFT1/2+1)))
ylim([0 1])
title('FFT for Turbidity variable','FontSize', 8)
ylabel('Amplitude','FontSize', 8)
xlabel('Frequency (Hz)','FontSize', 8)

%%
%Autoregressive method

p = 6; 
[arcoeffs1,variance1] = arburg(precipitation1,p);
[Pbx1,freqar1] = freqz(variance1,arcoeffs1,(NFFT1/2)+1,Fs);

[arcoeffs2,variance2] = arburg(airtemp1,p);
[Pbx2,freqar2] = freqz(variance2,arcoeffs2,(NFFT1/2)+1,Fs);

[arcoeffs3,variance3] = arburg(watertemp1,p);
[Pbx3,freqar3] = freqz(variance3,arcoeffs3,(NFFT1/2)+1,Fs);

[arcoeffs4,variance4] = arburg(salinity1,p);
[Pbx4,freqar4] = freqz(variance4,arcoeffs4,(NFFT1/2)+1,Fs);

[arcoeffs5,variance5] = arburg(turbidity1,p);
[Pbx5,freqar5] = freqz(variance5,arcoeffs5,(NFFT1/2)+1,Fs);

[arcoeffs6,variance6] = arburg(clorofila1,p);
[Pbx6,freqar6] = freqz(variance6,arcoeffs6,(NFFT1/2)+1,Fs);

subplot(3,2,1);
plot(freqar1,abs(Pbx1),'k')
xlabel('Frequency (days^{-1})','FontSize', 8); ylabel('Amplitude','FontSize', 8)
title('AR method for Precipitation frequencies','FontSize', 6)


subplot(3,2,2); 
plot(freqar2,abs(Pbx2),'k')
xlabel('Frequency (days^{-1})','FontSize', 8); ylabel('Amplitude','FontSize', 8)
title('AR method for Airtemp frequencies','FontSize', 6)


subplot(3,2,3); 
plot(freqar3,abs(Pbx3),'k')
xlabel('Frequency (days^{-1})','FontSize', 8); ylabel('Amplitude','FontSize', 8)
title('AR method for WaterTemp frequencies','FontSize', 6)

subplot(3,2,4); 
plot(freqar4,abs(Pbx4),'k')
xlabel('Frequency (days^{-1})','FontSize', 8); ylabel('Amplitude','FontSize', 8)
title('AR method for Salinity frequencies','FontSize', 6)

subplot(3,2,5); 
plot(freqar5,abs(Pbx5),'k')
xlabel('Frequency (days^{-1})','FontSize', 8); ylabel('Amplitude','FontSize', 8)
title('AR method for Turbidity frequencies','FontSize', 6)

subplot(3,2,6); 
plot(freqar6,abs(Pbx6),'k')
xlabel('Frequency (days^{-1})','FontSize', 8); ylabel('Amplitude','FontSize', 8)
title('AR method for Chlorophyll frequencies','FontSize', 6)


