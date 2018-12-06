clear
load BabyECGData;
data = HR;
[a,d] = Haar1D(data);
HRrec3 = ihaart(a,d,3);
HRrec5 = ihaart(a,d,5);
HRrec8 = ihaart(a,d,8);
figure
subplot(2,2,1)
plot(times, HR)
xlabel("time")
ylabel("voltage")
title('Original Signal')
subplot(2,2,2)
plot(times, HRrec3)
xlabel("time")
ylabel("voltage")
title('Level 3 Compression')
subplot(2,2,3)
plot(times, HRrec5)
xlabel("time")
ylabel("voltage")
title('Level 5 Compression')
subplot(2,2,4)
plot(times, HRrec8)
xlabel("time")
ylabel("voltage")
title('Level 8')
rec10 = reconstruction(a,d,10);

%%
figure
subplot(2,2,1)
plot(times, rec10)
xlabel("time")
ylabel("voltage")
title('threshold = 10')
rec50 = reconstruction(a,d,50);
subplot(2,2,2)
plot(times, rec50)
xlabel("time")
ylabel("voltage")
title('threshold = 50')
rec20 = reconstruction(a,d,20);
subplot(2,2,3)
plot(times, rec20)
xlabel("time")
ylabel("voltage")
title('threshold = 20')
rec30 = reconstruction(a,d,30);
subplot(2,2,4)
plot(times, rec30)
xlabel("time")
ylabel("voltage")
title('threshold = 30')
threshold = FindThreshold(0.9756,a,d);

%% Compute PSNR
ratio = (0.01:0.005:1)';
psnr = ComputePSNRs(data,a,d,ratio);
figure
plot(ratio,psnr,'x')
xlabel('compression ratio')
ylabel('PSNR')

%% Test Threshold
threshold = FindThreshold(0.05,a,d);
rec5 = reconstruction(a,d,threshold);
figure
%subplot(2,1,1)
hold on
plot(times, rec5)
xlabel("time")
ylabel("voltage")
title('threshold = 5')
%subplot(2,1,2)
threshold = FindThreshold(0.7,a,d);
rec6 = reconstruction(a,d,threshold);
plot(times, rec6)
xlabel("time")
ylabel("voltage")
title('threshold = 6')

hold off
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
function rec = InvHaar1D(a,d,level)
    rec = [a,d,level];
end
function rec = reconstruction(a,D,threshold)
    count = 0;
    if abs(a) < threshold
        a = 0;
        count = count + 1;
    end
    for i = 1:length(D)
        for j = 1:length(D{i})
            if abs(D{i}(j)) < threshold
                D{i}(j) = 0;
                count = count + 1;
            end
        end
    end
    rec = ihaart(a,D);
    count/(length(D{1})*2)
end
function threshold = FindThreshold(ratio,a,d)
    data = [a];
    for i=1:length(d)
        for j = 1:length(d{i})
            data = [data,d{i}(j)];
        end
    end
    data_sort = sort(abs(data));
    threshold = data_sort(ceil(ratio*length(data_sort)));
end
function PSNRs = ComputePSNRs(data,a,d,ratio)
    PSNRs = [];
    for i = 1:length(ratio)
       threshold = FindThreshold(ratio(i),a,d);
       rec = reconstruction(a,d,threshold);
       MSE = sum((rec - data).^2)./length(data);
       PSNR = 10*log10(1./MSE);
       PSNRs = [PSNRs, PSNR];
    end
end