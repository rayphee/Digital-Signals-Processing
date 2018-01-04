function finalSignal = estimatedSignal(shortTimeFT)
    specSamples = size(shortTimeFT, 2);
    inverseMatrix = ifft(shortTimeFT, 1024);
    finalSignal = real(inverseMatrix(1:128,1));
    for sampleNumber = 1:(specSamples-1)
        concatTimeSlice = cat(3, real(inverseMatrix(129:256,sampleNumber)), real(inverseMatrix(1:128,sampleNumber+1)));
        newSampleSlice = mean(concatTimeSlice, 3);
        finalSignal = cat(1, finalSignal, newSampleSlice);
    end
    finalSignal = cat(1, finalSignal, real(inverseMatrix(129:256,specSamples)));
end