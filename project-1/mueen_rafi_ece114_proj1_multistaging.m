%% Multistaging, no polyphasing
[audioFile, readSamplingFreq]= audioread('Wagner.wav');
outputAudio = srconvert_multi(audioFile);
y = srconvert_multi([1 zeros(1,3000)]);
verify(y);