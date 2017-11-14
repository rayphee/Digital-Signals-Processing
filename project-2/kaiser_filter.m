%% Kaiser

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