function [y] = Phase_permute_2D(x)
% find nan values of the input
A=isnan(x);
% change nans to zero
x(A)=0;

% extract FFT of input
x_fft=fftshift(fft2(x));
% extract phase of input matrix using FFT of input.
x_fft_angle=angle(x_fft);
% phase values of positive frequencies to maintain symmetry.
temp=x_fft_angle(1:1+floor(0.5*size(x_fft_angle,1)),:); [s1 s2]=size(temp);

% permute phase
% temp=temp(randperm(numel(temp))); 
temp=2* pi * rand(s1, s2) - pi; 
% temp=reshape(temp,s1,s2);
% assign the permuted phase
x_fft_angle_perm=x_fft_angle;
x_fft_angle_perm(1:s1,:)=temp;
temp=-fliplr(flipud(temp)); temp(1:2*s1-size(x,1),:)=[];
x_fft_angle_perm(1+s1:end,:)=temp;
x_fft_angle_perm(s1,s1)=0;

% apply inverse fft
y_fft=abs(x_fft).*exp(1i*x_fft_angle_perm);
y=ifft2(fftshift(y_fft),'symmetric');
% put nan values back in place
y(A)=nan;
% force y to be symmetric since it is FC (can be commented out if not needed!)
for i = 1 : size( y, 1)
        for j = 1 : size( y, 1)
                if i < j
                        y(j, i) = y(i, j);
                end
        end
end
end

