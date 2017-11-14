function finalOutputAudio = srconvert_multi(audioFile)
    [b2, a2] = ellip(5, 0.1, 70, 1/2);
    upsampledAudio = upsample(audioFile, 2);
    upsampledAudio = filter(b2, a2, upsampledAudio);
    for i=1:5
        upsampledAudio = upsample(upsampledAudio, 2);
        upsampledAudio = filter(b2, a2, upsampledAudio);
    end
    [b5, a5] = ellip(5, 0.1, 70, 1/5);
    finalUpsampledAudio = upsample(upsampledAudio, 5);
    finalUpsampledAudio = filter(b5, a5, finalUpsampledAudio);
    
    outputAudio = downsample(finalUpsampledAudio, 7);
    outputAudio = downsample(outputAudio, 7);
    finalOutputAudio = downsample(outputAudio, 3);
    %soundsc(finalOutputAudio, 24000);
    audiowrite('Output.wav',finalOutputAudio, 24000);
end