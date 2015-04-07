%% 
[retImg, retKernel] = Non_Uniform_L0('Pantheon.jpg',9,2e-3,4e3,0.1);
imwrite(retImg,'Pantheon_9,2e-3,4e3.png');
imwrite(retKernel,'Pantheon_Kernel_9,2e-3,4e3.png');

% [retImg, retKernel] = Non_Uniform_L0('Elephant.jpg',15,6e-3,1e3,0.1);
% imwrite(retImg,'Elephant_15,6e-3,1e31.png');
% imwrite(retKernel,'Elephant_Kernel_15,6e-3,1e31.png');
% 
% [retImg, retKernel] = Non_Uniform_L0('blurred1.jpg',17,6e-3,2e3,0.1);
% imwrite(retImg,'blurred1_17,6e-3,2e3.png');
% imwrite(retKernel,'blurred1_Kernel_17,6e-3,2e3.png');
% 
% [retImg, retKernel] = Non_Uniform_L0('blurred2.jpg',17,6e-3,2e3,0.1);
% imwrite(retImg,'blurred2_17,6e-3,2e3.png');
% imwrite(retKernel,'blurred2_Kernel_17,6e-3,2e3.png');
% 
% [retImg, retKernel] = Non_Uniform_L0('blurred3.jpg',17,6e-3,2e3,0.1);
% imwrite(retImg,'blurred3_17,6e-3,2e3.png');
% imwrite(retKernel,'blurred3_Kernel_17,6e-3,2e3.png');
% 
% 
% [retImg, retKernel] = Non_Uniform_L0('blurred4.jpg',17,6e-3,2e3,0.1);
% imwrite(retImg,'blurred4_17,6e-3,2e3.png');
% imwrite(retKernel,'blurred4_Kernel_17,6e-3,2e3.png');
% 
% [retImg, retKernel] = Non_Uniform_L0('blurred5.jpg',17,6e-3,2e3,0.1);
% imwrite(retImg,'blurred5_17,6e-3,2e3.png');
% imwrite(retKernel,'blurred5_Kernel_17,6e-3,2e3.png');
% 
% 
% [retImg, retKernel] = Non_Uniform_L0('Pertol_Station.jpg',17,6e-3,2e3,0.1);
% imwrite(retImg,'Pertol_Station_17,6e-3,2e3.png');
% imwrite(retKernel,'Pertol_Station_Kernel_17,6e-3,2e3.png');
% 
% [retImg, retKernel] = Non_Uniform_L0('Butcher_Shop.jpg',17,6e-3,4e3,0.1);
% imwrite(retImg,'Butcher_Shop_17,6e-3,4e3.png');
% imwrite(retKernel,'Butcher_Shop_Kernel_17,6e-3,4e3.png');
% 
% [retImg, retKernel] = Non_Uniform_L0('Book.jpg',17,6e-3,2e3,0.1);
% imwrite(retImg,'Book_17,6e-3,2e3.png');
% imwrite(retKernel,'Book_Kernel_17,6e-3,2e3.png');
% 
% [retImg, retKernel] = Non_Uniform_L0('Vintage_Car.jpg',17,4e-3,2e3,0.1);
% imwrite(retImg,'Vintage_Car_17,4e-3,2e3.png');
% imwrite(retKernel,'Vintage_Car_Kernel_17,4e-3,2e3.png');


