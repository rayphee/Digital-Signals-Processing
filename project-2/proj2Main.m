%% Chebyshev I

clear

% Filter Specs specified in FIR format
fs = 44100;
passbandEdge = 2500/fs; % normalized for easy input
stopbandEdge = 4000/fs; % normalized for easy input
maxPassbandGain = 40;
minPassbandGain = 37;
maxStopbandGain = -55;

% Filter Specs specified for IIR Filter design
passbandRipple = maxPassbandGain - minPassbandGain;
stopbandRipple = maxPassbandGain - maxStopbandGain;

[order, bleh] = cheb1ord(passbandEdge, stopbandEdge, passbandRipple, stopbandRipple);

order

[b, a] = cheby1(order, passbandRipple, passbandEdge);

genFilter = dfilt.df1(b, a);
c = cost(genFilter); 
numberOfOperations = c.nmult

[magResponseLin, genFreq] = freqz(b, a);
grdResponse = grpdelay(b, a);

normalizedFreq = genFreq/pi;

magResponseLog = 20*log10(abs(magResponseLin));

figure
subplot(2,2,1); plot(normalizedFreq,magResponseLog); title('Chebyshev I Magnitude (in dB)')
ylabel('Magnitude (dB)')
xlabel('Normalized Frequency (\times\pi rad/sample)');
subplot(2,2,2); plot(normalizedFreq,abs(magResponseLin)); title('Chebyshev I Magnitude (in Linear)')
ylabel('Magnitude')
xlabel('Normalized Frequency (\times\pi rad/sample)');
subplot(2,2,3); plot(normalizedFreq,grdResponse); title('Chebyshev I Group Delay')
ylabel('Delay')
xlabel('Normalized Frequency (\times\pi rad/sample)');

figure
subplot(2,1,1); zplane(b, a); title('Chebyshev I Pole Zero Plot');

impulseResponse = filter(b, a, [1, zeros(1,99)]);
subplot(2,1,2); stem(impulseResponse); title('Chebyshev I Impulse Response');

%fvtool(b,a); %verification that filter works

%% Chebyshev II

clear

% Filter Specs specified in FIR format
fs = 44100;
passbandEdge = 2500/fs; % normalized for easy input
stopbandEdge = 4000/fs; % normalized for easy input
maxPassbandGain = 40;
minPassbandGain = 37;
maxStopbandGain = -55;

% Filter Specs specified for IIR Filter design
passbandRipple = maxPassbandGain - minPassbandGain;
stopbandRipple = maxPassbandGain - maxStopbandGain;

[order, bleh] = cheb2ord(passbandEdge, stopbandEdge, passbandRipple, stopbandRipple);

order

[b, a] = cheby2(order, stopbandRipple, stopbandEdge);

genFilter = dfilt.df1(b, a);
c = cost(genFilter); 
numberOfOperations = c.nmult

[magResponseLin, genFreq] = freqz(b, a);
grdResponse = grpdelay(b, a);

normalizedFreq = genFreq/pi;

magResponseLog = 20*log10(abs(magResponseLin));

figure
subplot(2,2,1); plot(normalizedFreq,magResponseLog); title('Chebyshev II Magnitude (in dB)')
ylabel('Magnitude (dB)')
xlabel('Normalized Frequency (\times\pi rad/sample)');
subplot(2,2,2); plot(normalizedFreq,abs(magResponseLin)); title('Chebyshev II Magnitude (in Linear)')
ylabel('Magnitude')
xlabel('Normalized Frequency (\times\pi rad/sample)');
subplot(2,2,3); plot(normalizedFreq,grdResponse); title('Chebyshev II Group Delay')
ylabel('Delay')
xlabel('Normalized Frequency (\times\pi rad/sample)');

figure
subplot(2,1,1); zplane(b, a); title('Chebyshev II Pole Zero Plot');

impulseResponse = filter(b, a, [1, zeros(1,99)]);
subplot(2,1,2); stem(impulseResponse); title('Chebyshev II Impulse Response');

%fvtool(b,a); %verification that filter works

%% Kaiser

clear

% Filter Specs specified in FIR format
fs = 44100;
passbandEdge = 2500/fs; % normalized for easy input
stopbandEdge = 4000/fs; % normalized for easy input
maxPassbandGain = 40;
minPassbandGain = 37;
maxStopbandGain = -55;

% Filter Specs specified for IIR Filter design
passbandRipple = maxPassbandGain - minPassbandGain;
stopbandRipple = maxPassbandGain - maxStopbandGain;

specifiedFrequencyEdges = [2500 4000]
desiredAmplitude = [1 0];
desiredRipple = [(10^(passbandRipple/20)-1)/(10^(passbandRipple/20)+1) 10^(-stopbandRipple/20)];

[order, genFrequencyEdges, beta, magType] = kaiserord(specifiedFrequencyEdges, desiredAmplitude, desiredRipple, fs);

order

b = fir1(order, genFrequencyEdges,magType,kaiser(order+1,beta),'noscale');

[magResponseLin, genFreq] = freqz(b, 1);
grdResponse = grpdelay(b, 1);

genFilter = dfilt.df1(b, 1);
c = cost(genFilter); 
numberOfOperations = c.nmult

normalizedFreq = genFreq/pi;

magResponseLog = 20*log10(abs(magResponseLin));

figure
subplot(2,2,1); plot(normalizedFreq,magResponseLog); title('Kaiser Magnitude (in dB)')
ylabel('Magnitude (dB)')
xlabel('Normalized Frequency (\times\pi rad/sample)');
subplot(2,2,2); plot(normalizedFreq,abs(magResponseLin)); title('Kaiser Magnitude (in Linear)')
ylabel('Magnitude')
xlabel('Normalized Frequency (\times\pi rad/sample)');
subplot(2,2,3); plot(normalizedFreq,grdResponse); title('Kaiser Group Delay')
ylabel('Delay')
xlabel('Normalized Frequency (\times\pi rad/sample)');

figure
subplot(2,1,1); zplane(b, 1); title('Kaiser Pole Zero Plot');

impulseResponse = filter(b, 1, [1, zeros(1,99)]);
subplot(2,1,2); stem(impulseResponse); title('Kaiser Impulse Response');

%fvtool(b,1); %verification that filter works

%% Parks-McClellan

clear

% Filter Specs specified in FIR format
fs = 44100;
passbandEdge = 2500/fs; % normalized for easy input
stopbandEdge = 4000/fs; % normalized for easy input
maxPassbandGain = 40;
minPassbandGain = 37;
maxStopbandGain = -55;

% Filter Specs specified for IIR Filter design
passbandRipple = maxPassbandGain - minPassbandGain;
stopbandRipple = maxPassbandGain - maxStopbandGain;

specifiedFrequencyEdges = [2500 4000]
desiredAmplitude = [1 0];
desiredRipple = [(10^(passbandRipple/20)-1)/(10^(passbandRipple/20)+1) 10^(-stopbandRipple/20)];

[order, genFrequencyEdges, genAmplitudes, weights] = firpmord(specifiedFrequencyEdges, desiredAmplitude, desiredRipple, fs);

order


b = firpm(order, genFrequencyEdges, genAmplitudes, weights);

[magResponseLin, genFreq] = freqz(b, 1);
grdResponse = grpdelay(b, 1);

genFilter = dfilt.df1(b, 1);
c = cost(genFilter); 
numberOfOperations = c.nmult

normalizedFreq = genFreq/pi;

magResponseLog = 20*log10(abs(magResponseLin));

figure
subplot(2,2,1); plot(normalizedFreq,magResponseLog); title('Parks-McClellan Magnitude (in dB)')
ylabel('Magnitude (dB)')
xlabel('Normalized Frequency (\times\pi rad/sample)');
subplot(2,2,2); plot(normalizedFreq,abs(magResponseLin)); title('Parks-McClellan Magnitude (in Linear)')
ylabel('Magnitude')
xlabel('Normalized Frequency (\times\pi rad/sample)');
subplot(2,2,3); plot(normalizedFreq,grdResponse); title('Parks-McClellan Group Delay')
ylabel('Delay')
xlabel('Normalized Frequency (\times\pi rad/sample)');

figure
subplot(2,1,1); zplane(b, 1); title('Parks-McClellan Pole Zero Plot');

impulseResponse = filter(b, 1, [1, zeros(1,99)]);
subplot(2,1,2); stem(impulseResponse); title('Parks-McClellan Impulse Response');

%fvtool(b,1); %verification that filter works

%% Butterworth

clear

% Filter Specs specified in FIR format
fs = 44100;
passbandEdge = 2500/fs; % normalized for easy input
stopbandEdge = 4000/fs; % normalized for easy input
maxPassbandGain = 40;
minPassbandGain = 37;
maxStopbandGain = -55;

% Filter Specs specified for IIR Filter design
passbandRipple = maxPassbandGain - minPassbandGain;
stopbandRipple = maxPassbandGain - maxStopbandGain;

[order, cutoffFrequencies] = buttord(passbandEdge, stopbandEdge, passbandRipple, stopbandRipple); 

order

[b, a] = butter(order, cutoffFrequencies);
DF1 = dfilt.df1(b, a);
c1 = cost(DF1); 
numberOfOperationsDF1 = c1.nmult

DF2 = dfilt.df2(b, a);
c2 = cost(DF2); 
numberOfOperationsDF2 = c2.nmult

[z,p,k] = butter(order, cutoffFrequencies);
sections = zp2sos(z, p);
DF2SOS = dfilt.df2sos(sections);
c3 = cost(DF2SOS); 
numberOfOperationsDF2SOS = c3.nmult

DF2SOST = dfilt.df2tsos(sections);
c4 = cost(DF2SOST); 
numberOfOperationsDF2SOST = c4.nmult

[DF1Mag, genFreqDF1] = freqz(DF1);
normalizedFreq = genFreqDF1/pi;
[DF2Mag, genFreqDF2] = freqz(DF2);
[DF2SOSMag, genFreqDF2SOS] = freqz(DF2SOS);
[DF2SOSTMag, genFreqDF2SOST] = freqz(DF2SOST);

DF1GrD = grpdelay(DF1, genFreqDF1);
DF2GrD = grpdelay(DF2, genFreqDF1);
DF2SOSGrD = grpdelay(DF2SOS, genFreqDF1);
DF2SOSTGrD = grpdelay(DF2SOST, genFreqDF1);

figure
subplot(2,1,1); plot(normalizedFreq,abs(DF1Mag),normalizedFreq,abs(DF2Mag),normalizedFreq,abs(DF2SOSMag),normalizedFreq,abs(DF2SOSTMag)); title('Butterworth Magnitude')
ylabel('Magnitude')
xlabel('Normalized Frequency (\times\pi rad/sample)');
legend('Direct Form I','Direct Form II','Direct Form II SOS','Direct Form II SOS Transpose');
subplot(2,1,2); plot(normalizedFreq,DF1GrD,normalizedFreq,DF2GrD,normalizedFreq,DF2SOSGrD,normalizedFreq,DF2SOSTGrD); title('Butterworth Group Delay')
ylabel('Delay')
xlabel('Normalized Frequency (\times\pi rad/sample)');
legend('Direct Form I','Direct Form II','Direct Form II SOS','Direct Form II SOS Transpose');

[zDF1, pDF1, kDF1] = DF1.zpk;
[zDF2, pDF2, kDF2] = DF2.zpk;
[zDF2SOS, pDF2SOS, kDF2SOS] = DF2SOS.zpk;
[zDF2SOST, pDF2SOST, kDF2SOST] = DF2SOST.zpk;

figure
subplot(2,2,1); zplane(zDF1, pDF1); title('Butterworth Direct Form I');
subplot(2,2,2); zplane(zDF2, pDF2); title('Butterworth Direct Form II');
subplot(2,2,3); zplane(zDF2SOS, pDF2SOS); title('Butterworth Direct Form II SOS');
subplot(2,2,4); zplane(zDF2SOST, pDF2SOST); title('Butterworth Direct Form II SOS Transpose');

%fvtool(DF1); %verification that filter works
%fvtool(DF2); %verification that filter works
%fvtool(DF2SOS); %verification that filter works
%fvtool(DF2SOST); %verification that filter works

%% Elliptical

clear

% Filter Specs specified in FIR format
fs = 44100;
passbandEdge = 2500/fs; % normalized for easy input
stopbandEdge = 4000/fs; % normalized for easy input
maxPassbandGain = 40;
minPassbandGain = 37;
maxStopbandGain = -55;

% Filter Specs specified for IIR Filter design
passbandRipple = maxPassbandGain - minPassbandGain;
stopbandRipple = maxPassbandGain - maxStopbandGain;

[order, cutoffFrequencies] = ellipord(passbandEdge, stopbandEdge, passbandRipple, stopbandRipple); 

order

[b, a] = ellip(order, passbandRipple, stopbandRipple, passbandEdge);
DF1 = dfilt.df1(b, a);
c1 = cost(DF1); 
numberOfOperationsDF1 = c1.nmult

DF2 = dfilt.df2(b, a);
c2 = cost(DF2); 
numberOfOperationsDF2 = c2.nmult

[z,p,k] = ellip(order, passbandRipple, stopbandRipple, cutoffFrequencies);
sections = zp2sos(z, p);
DF2SOS = dfilt.df2sos(sections);
c3 = cost(DF2SOS); 
numberOfOperationsDF2SOS = c3.nmult

DF2SOST = dfilt.df2tsos(sections);
c4 = cost(DF2SOST); 
numberOfOperationsDF2SOST = c4.nmult

[DF1Mag, genFreqDF1] = freqz(DF1);
normalizedFreq = genFreqDF1/pi;
[DF2Mag, genFreqDF2] = freqz(DF2);
[DF2SOSMag, genFreqDF2SOS] = freqz(DF2SOS);
[DF2SOSTMag, genFreqDF2SOST] = freqz(DF2SOST);

DF1GrD = grpdelay(DF1, genFreqDF1);
DF2GrD = grpdelay(DF2, genFreqDF1);
DF2SOSGrD = grpdelay(DF2SOS, genFreqDF1);
DF2SOSTGrD = grpdelay(DF2SOST, genFreqDF1);

figure
subplot(2,1,1); plot(normalizedFreq,abs(DF1Mag),normalizedFreq,abs(DF2Mag),normalizedFreq,abs(DF2SOSMag),normalizedFreq,abs(DF2SOSTMag)); title('Elliptical Magnitude')
ylabel('Magnitude')
xlabel('Normalized Frequency (\times\pi rad/sample)');
legend('Direct Form I','Direct Form II','Direct Form II SOS','Direct Form II SOS Transpose');
subplot(2,1,2); plot(normalizedFreq,DF1GrD,normalizedFreq,DF2GrD,normalizedFreq,DF2SOSGrD,normalizedFreq,DF2SOSTGrD); title('Elliptical Group Delay')
ylabel('Delay')
xlabel('Normalized Frequency (\times\pi rad/sample)');
legend('Direct Form I','Direct Form II','Direct Form II SOS','Direct Form II SOS Transpose');

[zDF1, pDF1, kDF1] = DF1.zpk;
[zDF2, pDF2, kDF2] = DF2.zpk;
[zDF2SOS, pDF2SOS, kDF2SOS] = DF2SOS.zpk;
[zDF2SOST, pDF2SOST, kDF2SOST] = DF2SOST.zpk;

figure
subplot(2,2,1); zplane(zDF1, pDF1); title('Elliptical Direct Form I');
subplot(2,2,2); zplane(zDF2, pDF2); title('Elliptical Direct Form II');
subplot(2,2,3); zplane(zDF2SOS, pDF2SOS); title('Elliptical Direct Form II SOS');
subplot(2,2,4); zplane(zDF2SOST, pDF2SOST); title('Elliptical Direct Form II SOS Transpose');

%fvtool(DF1); %verification that filter works
%fvtool(DF2); %verification that filter works
%fvtool(DF2SOS); %verification that filter works
%fvtool(DF2SOST); %verification that filter works