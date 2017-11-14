%% Chebyshev I

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