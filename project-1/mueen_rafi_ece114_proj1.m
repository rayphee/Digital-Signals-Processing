%% No multi-staging or polyphasing
[audioFile, readSamplingFreq]= audioread('Wagner.wav');
outputAudio = srconvert(audioFile);
y = srconvert([1 zeros(1,3000)]);
verify(y);