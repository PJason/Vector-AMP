clear; close all; clc;
rng('default')

%% parameter setting
sR = 0.3;        % Sampling rate
input_snr = 50;  % Noise level 
rho = 1; 
lambda = 0.1;    % trade-off for regularization
damping = 0.6;
Niter = 100;


%% signal
im = imread('lena64.png');
im = im2double(rgb2gray(im));
mask = rand(size(im)) < sR; ind = find(mask==1); % Mask
% Reconstruction from 33% of Fourier measurements
Ma = diag(mask(:));
Ma(all(Ma == 0,2),:) = []; 
Amat = diag(mask(:));
Amat(all(Amat == 0,2),:) = []; 
n = size(im,1);
N = n^2;
r = 2*(n-1)*n; % number of rows in the 2D TV operator
A = @(x) Ma*x(:); % Select 33% of the values in x;
At = @(x) reshape(Ma'*x(:), size(im)); % Adjoint operator + reshape image
y = A(im); 
sigma_noise = 10^(-input_snr/20)*std(im(:)); % Add Gaussian i.i.d. noise
y = y + randn(size(y))*sigma_noise;


% 1D antisometric TV matrix
D = diag(ones(n-1,1),1);
D = D(1:n-1,:);
D(logical(eye(size(D)))) = -1*ones(n-1,1); 


% IFVAMP 
ERRvsRHO = zeros(1,Niter);
% funIFVAMP(A, y, img, D, damping, lambda, Rho, Lipschitz, Niter, N, r);
[resX, resSSIM, resRHO] = funIFVAMP(Amat, y, im, D, damping, lambda, rho, 2, Niter, N, r);
ERRvsRHO(1,:) = resSSIM;
RHOvsITR(1,:) = resRHO;

% plotting
figure(1)
plot(RHOvsITR(1,:),'-', 'linewidth',1,'color',[0.9290 0.6940 0.1250]); hold on
plot(ERRvsRHO(1,:),'-.', 'linewidth',2,'color',[0.9290 0.6940 0.1250]);
xlabel('iterations','Fontname', 'Times New Roman','FontSize',16);
ylabel('SSIM','Fontname', 'Times New Roman','FontSize',16);
legend('\rho^t ','IFVAMP \rho^0 = 1','Fontname', 'Times New Roman','FontSize',16)
hold off 





