clear
load BabyECGData;
[a,d] = Haar1D(HR);

function [a,D]=Haar1D(data)
    D = {};
    level = log2(length(data));
    scaling_filter = [1/sqrt(2);1/sqrt(2)];
    wavelet_filter = [1/sqrt(2);-1/sqrt(2)];
    for i = 1:level
        a = conv(data,scaling_filter);
        d = conv(data,wavelet_filter);
        a = a(2:2:end);
        d = d(2:2:end);
        D{i} = d; 
        data = a;
    end
end

function signal = InvHaar1D(a,d,level)
    signal = zeros(length(d{1})*2);
    index = length(signal);
    signallength = length(d{1});
    divisor 
    for i = length(d):(level+1)
        signal()
    end
end