image = imread('cheetah.png');
grayimage = rgb2gray(image);
% imshow(grayimage)
% title('threshold = 0')
[a,h,v,d] = haart2(grayimage);

%% Test Different Threshold
% [threshold, data_sort] = FindThreshold2D(a,h,v,d,0.8);
% [rec5,a_rec5,h_rec5,v_rec5,d_rec5] = reconstruction2D(a,h,v,d,threshold);
% ratio5 = Compare(a,h,v,d, a_rec5,h_rec5,v_rec5,d_rec5);
% figure
% imshow(uint8(rec5))
% title('threshold = 80')
% [threshold, data_sort] = FindThreshold2D(a,h,v,d,0.5);
% [rec9,a_rec9,h_rec9,v_rec9,d_rec9] = reconstruction2D(a,h,v,d,threshold);
% ratio9 = Compare(a,h,v,d, a_rec9,h_rec9,v_rec9,d_rec9);
% figure
% imshow(uint8(rec9))
% title('threshold = 50')
% threshold = FindThreshold2D(a,h,v,d,0.9)
% rec70 = reconstruction2D(a,h,v,d,threshold);
% figure
% imshow(rec70)
% title('threshold = 70')

%% Compute PSNR
% ratio = (0.01:0.005:0.9999)';
% [psnr,mse] = ComputePSNRs2D(double(grayimage),a,h,v,d,ratio);
% figure
% plot(ratio,psnr,'x')
% xlabel('compression ratio')
% ylabel('PSNR')
% figure
% plot(ratio,mse,'x')
% xlabel('compression ratio')
% ylabel('MSE')
%% Show compressed images from 0 - 100% ratio
% ratio = (0.01:0.005:0.9999)';
% ratio = [ratio; 0.9960; 0.9970;0.9980;0.9990;0.9995;0.9997;0.9999];
% figure
% pause(5);
% for i = 1:length(ratio)
%     [threshold, data_sort] = FindThreshold2D(a,h,v,d,ratio(i));
%     rec = reconstruction2D(a,h,v,d,threshold);
%     imshow(uint8(rec));
%     titleText = join(['Compressed Ratio =',string(ratio(i)*100),'%']);
%     title(titleText);
%     pause(0);
% end
%% Compute ssim
% ratio = (0.01:0.005:0.9999)';
% ssim = ComputeSSIMs2D(double(grayimage),a,h,v,d,ratio);
% figure
% hold on
% plot(ratio,ssim,'ro')
% plot(ratio,ssim)
% hold off
% xlabel('Compression Ratio(%)')
% ylabel('SSIM')
%% HHLL plot
[a,h,v,d] = haart2(grayimage,1);
[a2,h2,v2,d2] = haart2(a,1);
imshow(uint8(a))

function [rec,a,h,v,d] = reconstruction2D(a,h,v,d,threshold)
    count = 0;
    total = 0; %total nonzero element
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
    count/total
end
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
function [PSNRs,MSEs] = ComputePSNRs2D(data,a,h,v,d,ratio)
    PSNRs = zeros(length(ratio),1);
    MSEs = zeros(length(ratio),1);
    for i = 1:length(ratio)
       [threshold, data_sort] = FindThreshold2D(a,h,v,d,ratio(i));
       rec = reconstruction2D(a,h,v,d,threshold);
       MSE = sum(sum((rec - data).^2))./(size(data,1)*size(data,2));
       PSNR = 10*log10(1./MSE);
       MSEs(i) = MSE;
       PSNRs(i) = PSNR;
    end
end
function SSIMs = ComputeSSIMs2D(data,a,h,v,d,ratio)
    SSIMs = zeros(length(ratio),1);
    for i = 1:length(ratio)
       [threshold, data_sort] = FindThreshold2D(a,h,v,d,ratio(i));
       rec = reconstruction2D(a,h,v,d,threshold);
       SSIM = ssim(rec,data);
       SSIMs(i) = SSIM;
    end
end
function ratio = Compare(a,h,v,d, a_rec,h_rec,v_rec,d_rec)
    n = size(h{1},1)*2;
    coef = zeros(n,n);
    coef_rec = zeros(n,n);
    coef(1) = a;
    coef_rec(1) = a_rec;
    index = 2;
    for i=1:length(d)
        for j = 1:(size(d{i},1)*size(d{i},2))
            coef(index:index+2) = [h{i}(j), v{i}(j),d{i}(j)];
            coef_rec(index:index+2) = [h_rec{i}(j), v_rec{i}(j),d_rec{i}(j)];
            index = index + 3;
        end
    end
    coef_count = sum(coef(:) ~= 0); % number nonzero before compress
    coef_rec_count = sum(coef_rec(:) ~= 0); % number nonzero after compress 
    ratio = (coef_count - coef_rec_count)/coef_count;
end