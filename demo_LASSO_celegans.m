%% Import hyperspectral image stack as data matrix D
% Width of D is the number of spectral channels
% Height of D is the number of pixels
clear;
clc;
addpath('Data/CElegans');
% Define size of the image
size_x = 200;
size_y = 200;
size_lambda = 100;
y = zeros(size_x,size_y,size_lambda);
D = zeros(size_x*size_y, size_lambda);
load('Data/celegans_spec_init.mat');

for i = 0:(size_lambda - 1)
    if i<10
        temp = ['Stack000',num2str(i),'.tif'];
    else
        temp = ['Stack00',num2str(i),'.tif'];
    end
        img = imread(temp);
        temp_column = reshape(img,size_x*size_y,1);
        y(:,:,i+1) = img;
        D(:,i+1) = temp_column;
end
clear size_x size_y size_lambda temp temp_column i img;
rmpath('Data/CElegans');

%% Run LASSO
lambda = 0.6e1;
copt = C_ML_Lasso( y, celegans_spec_init',lambda);
sopt = celegans_spec_init;

%% Plot the results
% Plot the concentration maps first
%{
% for mcr use
Comp1 = reshape(copt(:,1),200,200);
Comp2 = reshape(copt(:,2),200,200);
Comp3 = reshape(copt(:,3),200,200);
Comp4 = reshape(copt(:,4),200,200);

%}
%{{
% for lasso use
Comp1 = copt(:,:,1);
Comp2 = copt(:,:,2);
Comp3 = copt(:,:,3);
Comp4 = copt(:,:,4);
%}

Comp1 = imadjust(uint8(255.*Comp1/(max(max(Comp1))-min(min(Comp1)))));
Comp2 = imadjust(uint8(255.*Comp2/(max(max(Comp2))-min(min(Comp2)))));
Comp3 = imadjust(uint8(255.*Comp3/(max(max(Comp3))-min(min(Comp3)))));
Comp4 = imadjust(uint8(255.*Comp4/(max(max(Comp4))-min(min(Comp4)))));

rgbComp1 = cat(3, Comp1*0.5, Comp1*0, Comp1*0.5); % magneta for oxidized lipid
rgbComp2 = cat(3, Comp2*1, Comp2*0, Comp2*0); % red for nucleus
rgbComp3 = cat(3, Comp3*0, Comp3*0.8, Comp3*0); % green for protein
rgbComp4 = cat(3, Comp4*0, Comp4*0.3, Comp4*0.3);% cyan for water
C_disp = rgbComp1+rgbComp2+rgbComp3+rgbComp4;

figure;
subplot(1,2,1); 
imshow(C_disp,'InitialMagnification',150);

Spec1 = sopt(1,:);
Spec2 = sopt(2,:);
Spec3 = sopt(3,:);
Spec4 = sopt(4,:);
subplot(1,2,2);
plot(Spec1,'m','Linewidth',1.5);
hold on;
plot(Spec2,'r','Linewidth',1.5);
plot(Spec3,'g','Linewidth',1.5);
plot(Spec4,'c','Linewidth',1.5);
xlabel 'Channel';
ylabel 'Int. (a.u.)';
legend('Oxidized lipid','Fat droplets','lysosome-related organelles','Protein')
pbaspect([2 1 1]);
