image = imread('cheetah.png');
grayimage = rgb2gray(image);
imshow(grayimage)
title('compression ratio = 0%')
[a,h,v,d] = haart2(grayimage);

%% Plot Different Compression Ratios
% [threshold5, data_sort5] = FindThreshold2D(a,h,v,d,0.5);
% [rec5,a_rec5,h_rec5,v_rec5,d_rec5] = reconstruction2D(a,h,v,d,threshold5);
% figure
% subplot(1,2,1)
% imshow(uint8(rec5))
% title('compression ratio = 50%')
% [threshold8, data_sort8] = FindThreshold2D(a,h,v,d,0.8);
% [rec8,a_rec8,h_rec8,v_rec8,d_rec8] = reconstruction2D(a,h,v,d,threshold8);
% subplot(1,2,2)
% imshow(uint8(rec8))
% title('compression ratio = 80%')
% 
% figure
% [threshold9, data_sort9] = FindThreshold2D(a,h,v,d,0.9);
% [rec9,a_rec9,h_rec9,v_rec9,d_rec9] = reconstruction2D(a,h,v,d,threshold9);
% subplot(1,2,1)
% imshow(uint8(rec9))
% title('compression ratio = 90%')
% [threshold98, data_sort98] = FindThreshold2D(a,h,v,d,0.98);
% [rec98,a_rec98,h_rec98,v_rec98,d_rec98] = reconstruction2D(a,h,v,d,threshold98);
% subplot(1,2,2)
% imshow(uint8(rec98))
% title('compression ratio = 98%')
% 
% figure
% [threshold99, data_sort99] = FindThreshold2D(a,h,v,d,0.99);
% [rec99,a_rec99,h_rec99,v_rec99,d_rec99] = reconstruction2D(a,h,v,d,threshold99);
% subplot(1,2,1)
% imshow(uint8(rec99))
% title('compression ratio = 99%')
% [threshold995, data_sort995] = FindThreshold2D(a,h,v,d,0.995);
% [rec995,a_rec995,h_rec995,v_rec995,d_rec995] = reconstruction2D(a,h,v,d,threshold995);
% subplot(1,2,2)
% imshow(uint8(rec995))
% title('compression ratio = 99.5%')
% 
% figure
% [threshold999, data_sort999] = FindThreshold2D(a,h,v,d,0.999);
% [rec999,a_rec999,h_rec999,v_rec999,d_rec999] = reconstruction2D(a,h,v,d,threshold999);
% subplot(1,2,1)
% imshow(uint8(rec999))
% title('compression ratio = 99.9%')
% [threshold9995, data_sort9995] = FindThreshold2D(a,h,v,d,0.9995);
% [rec9995,a_rec9995,h_rec9995,v_rec9995,d_rec9995] = reconstruction2D(a,h,v,d,threshold9995);
% subplot(1,2,2)
% imshow(uint8(rec9995))
% title('compression ratio = 99.95%')
% 
% figure
% [threshold9997, data_sort9997] = FindThreshold2D(a,h,v,d,0.9997);
% [rec9997,a_rec9997,h_rec9997,v_rec9997,d_rec9997] = reconstruction2D(a,h,v,d,threshold9997);
% subplot(1,2,1)
% imshow(uint8(rec9997))
% title('compression ratio = 99.97%')
% [threshold9999, data_sort9999] = FindThreshold2D(a,h,v,d,0.9999);
% [rec9999,a_rec9999,h_rec9999,v_rec9999,d_rec9999] = reconstruction2D(a,h,v,d,threshold9999);
% subplot(1,2,2)
% imshow(uint8(rec9999))
% title('compression ratio = 99.99%')
%% Compute PSNR, MSE and SSIM
ratio = (0.01:0.005:1)';
ratio = [ratio; (0.9950:0.0001:0.9999)'];
[psnr,mse,ssim] = ComputePSNRs2D(double(grayimage),a,h,v,d,ratio);

figure
plot(ratio,psnr,'x')
xlabel('compression ratio')
ylabel('PSNR')

figure
plot(ratio,mse,'x')
xlabel('compression ratio')
ylabel('MSE')

figure
plot(ratio,ssim,'x')
xlabel('compression ratio')
ylabel('SSIM')

% Function to reconstruct a 2D image given a specific threshold
% INPUTS
%   a,h,v,d: detail and approximation coefficient
%   threshold: coefficients smaller than threshold will be set to zero
% OUTPUTS
%   rec: reconstructed signal in 2D
%   a,h,v,d: updated detail and approximation coefficient
function [rec,a,h,v,d] = reconstruction2D(a,h,v,d,threshold)
    count = 0;
    total = 0; % total nonzero elements
    if a ~= 0
        if abs(a) < threshold
            a = 0;
            count = count + 1;
        end
        total = total + 1; 
    end
    for i = 1:length(d)
        for j = 1:(size(d{i},1)*size(d{i},2))
            if d{i}(j) ~= 0
                if abs(d{i}(j)) < threshold 
                    d{i}(j) = 0;
                    count = count + 1;
                end
                total = total +1;
            end
            
            if h{i}(j) ~= 0
                if abs(h{i}(j)) < threshold
                    h{i}(j) = 0;
                    count = count + 1;
                end
                total = total + 1;
            end
            
            if v{i}(j) ~= 0
                if abs(v{i}(j)) < threshold
                    v{i}(j) = 0;
                    count = count + 1;
                end
                total = total + 1;
            end
        end
    end
    rec = ihaart2(a,h,v,d);
    compression_rate = count/total
end
% Function to find a desired threshold given a specific compression ratio
% INPUTS
%   a,h,v,d: detail and approximation coefficient
%   ratio: desired compression ratio
% OUTPUTS
%   threshold: desired threshold value
%   data_sort: 1D array of sorted coefficients
function [threshold,data_sort] = FindThreshold2D(a,h,v,d,ratio)
    n = size(h{1},1)*2;
    data = zeros(n,n);
    data(1) = a;
    index = 2;
    for i=1:length(d)
        for j = 1:(size(d{i},1)*size(d{i},2))
            data(index:index+2) = [h{i}(j), v{i}(j),d{i}(j)];
            index = index + 3;
        end
    end
    data_sort = sort(abs(data(:)));
    start_index = sum(data_sort==0);
    threshold = data_sort(start_index + ceil(ratio*(length(data_sort)-start_index)));
end
% Function to compute PSNR, MSE and SSIM given a range of ratios
% INPUTS
%   data: a 1D input signal 
%   ratio: a 1D array of compression ratio
%   a: approximation coefficient
%   D: object that stores details coefficient at all levels
% OUTPUTS
%   PSNRs: the corresponding PSNR values
%   MSEs: the corresponding PSNR values
%   SSIMs: the corresponding PSNR values
function [PSNRs,MSEs,SSIMs] = ComputePSNRs2D(data,a,h,v,d,ratio)
    PSNRs = zeros(length(ratio),1);
    MSEs = zeros(length(ratio),1);
    SSIMs = zeros(length(ratio),1);
    for i = 1:length(ratio)
       [threshold, ~] = FindThreshold2D(a,h,v,d,ratio(i));
       rec = reconstruction2D(a,h,v,d,threshold);
       MSE = sum(sum((rec - data).^2))./(size(data,1)*size(data,2));
       PSNR = 10*log10(1./MSE);
       SSIM = ssim(rec,data);
       MSEs(i) = MSE;
       PSNRs(i) = PSNR;
       SSIMs(i) = SSIM;
    end
end