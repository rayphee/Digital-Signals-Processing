%% Elliptical

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