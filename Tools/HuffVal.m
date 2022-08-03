function [HUFFVAL]=HuffVal(Encoder,Coes)
BandLength = Encoder.Se-Encoder.Ss;
HUFFVAL = [];
Run_of_Length = 0;
for blockid = 1:length(Coes)
    R = 0;
    for bandid = 1:BandLength
        if Coes(bandid,blockid)
            if Run_of_Length
                EOBRUN = EnsureGategory(Run_of_Length)-1;
            end
            if EOBRUN > 14
                HUFFVAL = [HUFFVAL bitshift(14,4)];
                EOBRUN = EOBRUN-14;
            end
            RS = bitshift(EOBRUN,4);
            HUFFVAL = [HUFFVAL RS];
            Run_of_Length = 0;
            S = EnsureGategory(Coes(bandid,blockid));
            RS = bitshift(R,4)+S;
            HUFFVAL = [HUFFVAL RS];
            R = 0;
        else
            R = R + 1;
        end
    end
    if R == BandLength
        Run_of_Length = Run_of_Length + 1;
    elseif R == 0
        continue
    else
        HUFFVAL = [HUFFVAL 0];
    end
end
end