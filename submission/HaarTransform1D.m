clear
load BabyECGData; % MATLAB R2018
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

%% Reconstruction at different threshold
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

%% Compute and plot PSNR, MSE and SSIM 
ratio = (0.01:0.005:1)';
[psnr, mse, ssim] = ComputePSNRs1D(data,a,d,ratio);

figure
plot(ratio,psnr,'x')
xlabel('compression ratio')
ylabel('PSNR')

figure
plot(ratio,mse,'o')
xlabel('compression ratio')
ylabel('MSE')

figure
plot(ratio,ssim,'-')
xlabel('compression ratio')
ylabel('SSIM')

%% Function to obtain approximation and detail coefficients for a 1D signal
% data: a 1D array signal
% a: approximation coefficient
% D: object that stores details coefficient at all levels
function [a,D]=Haar1D(data)
    D = {}; % D stores the detail coefficients of all levels
    level = log2(length(data)); % The max level that can be obtained
    scaling_filter = [1/sqrt(2);1/sqrt(2)]; % Scaling function
    wavelet_filter = [1/sqrt(2);-1/sqrt(2)]; % Wavelet function
    for i = 1:level % At each level i
        a = conv(data,scaling_filter); % Convolute data with scaling function
        d = conv(data,wavelet_filter); % Convolute data with wavelet function
        a = a(2:2:end); % Downsample to obtain approximation coefficients 
        d = d(2:2:end); % Downsample to obtain detail coefficients
        D{i} = d; % Save the result
        data = a; % Update the data used for the next level
    end
end
function rec = reconstruction1D(a,D,threshold)
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
function threshold = FindThreshold1D(ratio,a,d)
    data = [a];
    for i=1:length(d)
        for j = 1:length(d{i})
            data = [data,d{i}(j)];
        end
    end
    data_sort = sort(abs(data));
    threshold = data_sort(ceil(ratio*length(data_sort)));
end
function [PSNRs, MSEs, SSIMs] = ComputePSNRs1D(data,a,d,ratio)
    PSNRs = zeros(length(ratio),1);
    MSEs = zeros(length(ratio),1);
    SSIMs = zeros(length(ratio),1);
    for i = 1:length(ratio)
       threshold = FindThreshold(ratio(i),a,d);
       rec = reconstruction(a,d,threshold);
       MSE = sum((rec - data).^2)./length(data);
       PSNR = 10*log10(1./MSE);
       SSIM = ssim(rec,data);
       MSEs(i) = MSE;
       PSNRs(i) = PSNR;
       SSIMs(i) = SSIM;
    end
end