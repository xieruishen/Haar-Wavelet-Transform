image = imread('cheetah.png');
grayimage = rgb2gray(image);
imshow(grayimage)
title('compression ratio = 0%')
[a,h,v,d] = haart2(grayimage);

%% Plot Different Compression Ratio
[threshold5, data_sort5] = FindThreshold2D(a,h,v,d,0.5);
[rec5,a_rec5,h_rec5,v_rec5,d_rec5] = reconstruction2D(a,h,v,d,threshold5);
ratio5 = Compare(a,h,v,d, a_rec5,h_rec5,v_rec5,d_rec5);
figure
subplot(1,2,1)
imshow(uint8(rec5))
title('compression ratio = 50%')
[threshold8, data_sort8] = FindThreshold2D(a,h,v,d,0.8);
[rec8,a_rec8,h_rec8,v_rec8,d_rec8] = reconstruction2D(a,h,v,d,threshold8);
ratio8 = Compare(a,h,v,d, a_rec8,h_rec8,v_rec8,d_rec8);
subplot(1,2,2)
imshow(uint8(rec8))
title('compression ratio = 80%')

figure
[threshold9, data_sort9] = FindThreshold2D(a,h,v,d,0.9);
[rec9,a_rec9,h_rec9,v_rec9,d_rec9] = reconstruction2D(a,h,v,d,threshold9);
ratio9 = Compare(a,h,v,d, a_rec9,h_rec9,v_rec9,d_rec9);
subplot(1,2,1)
imshow(uint8(rec9))
title('compression ratio = 90%')
[threshold98, data_sort98] = FindThreshold2D(a,h,v,d,0.98);
[rec98,a_rec98,h_rec98,v_rec98,d_rec98] = reconstruction2D(a,h,v,d,threshold98);
ratio98 = Compare(a,h,v,d, a_rec98,h_rec98,v_rec98,d_rec98);
subplot(1,2,2)
imshow(uint8(rec98))
title('compression ratio = 98%')

figure
[threshold99, data_sort99] = FindThreshold2D(a,h,v,d,0.99);
[rec99,a_rec99,h_rec99,v_rec99,d_rec99] = reconstruction2D(a,h,v,d,threshold99);
ratio99 = Compare(a,h,v,d, a_rec99,h_rec99,v_rec99,d_rec99);
subplot(1,2,1)
imshow(uint8(rec99))
title('compression ratio = 99%')
[threshold995, data_sort995] = FindThreshold2D(a,h,v,d,0.995);
[rec995,a_rec995,h_rec995,v_rec995,d_rec995] = reconstruction2D(a,h,v,d,threshold995);
ratio995 = Compare(a,h,v,d, a_rec995,h_rec995,v_rec995,d_rec995);
subplot(1,2,2)
imshow(uint8(rec995))
title('compression ratio = 99.5%')

figure
[threshold999, data_sort999] = FindThreshold2D(a,h,v,d,0.999);
[rec999,a_rec999,h_rec999,v_rec999,d_rec999] = reconstruction2D(a,h,v,d,threshold999);
ratio999 = Compare(a,h,v,d, a_rec999,h_rec999,v_rec999,d_rec999);
subplot(1,2,1)
imshow(uint8(rec999))
title('compression ratio = 99.9%')
[threshold9995, data_sort9995] = FindThreshold2D(a,h,v,d,0.9995);
[rec9995,a_rec9995,h_rec9995,v_rec9995,d_rec9995] = reconstruction2D(a,h,v,d,threshold9995);
ratio9995 = Compare(a,h,v,d, a_rec9995,h_rec9995,v_rec99,d_rec9995);
subplot(1,2,2)
imshow(uint8(rec9995))
title('compression ratio = 99.95%')

figure
[threshold9997, data_sort9997] = FindThreshold2D(a,h,v,d,0.9997);
[rec9997,a_rec9997,h_rec9997,v_rec9997,d_rec9997] = reconstruction2D(a,h,v,d,threshold9997);
ratio9997 = Compare(a,h,v,d, a_rec9997,h_rec9997,v_rec9997,d_rec9997);
subplot(1,2,1)
imshow(uint8(rec9997))
title('compression ratio = 99.97%')
[threshold9999, data_sort9999] = FindThreshold2D(a,h,v,d,0.9999);
[rec9999,a_rec9999,h_rec9999,v_rec9999,d_rec9999] = reconstruction2D(a,h,v,d,threshold9999);
ratio9999 = Compare(a,h,v,d, a_rec9999,h_rec9999,v_rec9999,d_rec9999);
subplot(1,2,2)
imshow(uint8(rec9999))
title('compression ratio = 99.99%')
%% Compute PSNR and MSE
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
    count/total;
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