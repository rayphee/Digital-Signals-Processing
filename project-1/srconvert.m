function outputAudio = srconvert(audioFile)
    [b, a] = ellip(5, 0.05, 65, 1/320);
    upsampledAudio = upsample(audioFile, 320);
    filteredUpsampledAudio = filter(b, a, upsampledAudio);
    outputAudio = downsample(filteredUpsampledAudio, 147);
    soundsc(outputAudio, 24000);
    audiowrite('Output.wav',outputAudio, 24000);
end