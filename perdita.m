vettore_x_lambda = [5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,25,30,35,40,50];
%vettore_y1_prob_blocco_modello = [0.0002275, 0.0006709, 0.001569, 0.003102, 0.0053981, 0.00850898, 0.0124095, 0.017011, 0.0221872, 0.0277946, 0.0336937, 0.0397591, 0.0458846, 0.0519854, 0.0579963, 0.0638694,0.090];
%vettore_y2_prob_blocco_sim = [0, 0, 0, 0, 0, 0, 0.001, 0.002, 0.002, 0.006, 0.012, 0.022, 0.036, 0.053, 0.072, 0.092, 0.2002, 0.304];


%SEED 1, tempo 5000 secondi 5.5 7.5
%vettore_y1_prob_blocco_modello = [0.00027, 0.00087, 0.00217, 0.00456, 0.00842, 0.0140, 0.0215, 0.0310, 0.0423, 0.0553, 0.0697, 0.0852, 0.1016, 0.1186, 0.1360, 0.1535, 0.2395, 0.3166, 0.382, 0.438, 0.526];
vettore_y1_prob_blocco_modello = [0.0001,0.00056,0.0014,0.0030,0.0058,0.0099,0.0156,0.0230, 0.0320,0.0426,0.0547,0.0680,0.0823,0.0973,0.1130,0.1291,0.2103,0.2857,0.3518,0.4086,0.4991];
vettore_y2_prob_blocco_sim = [0, 0, 0, 0, 0.0001, 0.001, 0.001, 0.002, 0.003, 0.007, 0.012, 0.022, 0.035, 0.052, 0.070, 0.091, 0.206, 0.308, 0.395, 0.462, 0.563];
vettore_traffico_perso_sim = [0, 0, 0, 0, 0, 0.001,0.009 ,0.020, 0.034, 0.092, 0.183, 0.354, 0.599, 0.930, 1.332, 1.817, 5.145,9.253 ,13.835,18.469,28.216];

%SEED 2, tempo 1000 secondi 2.8 10.2
%vettore_y3_prob_blocco_modello = [0.0002,0.0007,0.002,0.0044,0.0083,0.014,0.0218,0.0316,0.043,0.0564,0.0708,0.0862,0.1023,0.1188,0.1355,0.1523,0.2331,0.3048,0.3665 , 0.419 , 0.503];
vettore_y3_prob_blocco_modello = [0.00014, 0.00048, 0.0013, 0.0029, 0.0056, 0.0098, 0.0156, 0.0232, 0.0324, 0.0433, 0.0555, 0.0688, 0.0830, 0.0978, 0.1130, 0.1285, 0.2052, 0.2754, 0.3369, 0.3901, 0.4761];
vettore_y4_prob_blocco_sim = [0,0,0,0,0,0,0.0003, 0.0005, 0.001, 0.001, 0.003, 0.006, 0.011, 0.018, 0.028, 0.039, 0.124,0.223,0.311,0.386,0.498];

%SEED 5, tempo 1000 secondi 4.6 8.3
%vettore_y5_prob_blocco_modello = [0.0002, 0.00075, 0.0019, 0.004, 0.0076, 0.0128, 0.02, 0.0291, 0.04, 0.0525, 0.0665, 0.0817, 0.0977, 0.1143, 0.1313, 0.1485,0.2329, 0.3087 ,0.3739 ,0.4293, 0.5169]; 
vettore_y5_prob_blocco_modello = [0.00014, 0.00047, 0.0012, 0.0027,0.0051, 0.0089, 0.0142, 0.0212, 0.0298, 0.0400, 0.0516, 0.0645, 0.0784, 0.0930, 0.1083, 0.1240, 0.2035, 0.2776, 0.3427, 0.3988, 0.4887];
vettore_y6_prob_blocco_sim = [0,0,0,0,0,0,0,0.0001,0.0002,0.001,0.002 ,0.005, 0.010, 0.017, 0.028, 0.040, 0.131, 0.233,0.325, 0.398,0.509];


%SEED tempo 1000 secondi 3.7 9.3
%vettore_y7_prob_blocco_modello =[0.0002,0.0007,0.0018,0.0039,0.0075, 0.0127, 0.0199,0.029,0.04,0.0526, 0.0666,0.0816, 0.0975,0.1140,0.1308,0.1477, 0.2305,0.3045,0.3682,0.4225,0.5089];
vettore_y7_prob_blocco_modello = [0.00013,0.00044, 0.0011, 0.0026, 0.0050, 0.0083, 0.01414, 0.02112, 0.02976, 0.04, 0.0516, 0.0645, 0.0783, 0.0929, 0.10807, 0.12357,0.2017, 0.2741, 0.3377, 0.3926, 0.4810];
vettore_y8_prob_blocco_sim = [0,0,0,0,0,0,0,0,0.0005,0.001,0.002, 0.006, 0.012,0.02,0.032, 0.046, 0.138,0.237,0.326, 0.4,0.510];



% 
plot(vettore_x_lambda,vettore_y1_prob_blocco_modello,'*-','color', 'red');
hold on
plot(vettore_x_lambda,vettore_y2_prob_blocco_sim,'color', 'red');
hold on
% 
plot(vettore_x_lambda,vettore_y3_prob_blocco_modello,'*-','color', 'blue');
hold on
plot(vettore_x_lambda,vettore_y4_prob_blocco_sim,'color', 'blue');
% 
plot(vettore_x_lambda,vettore_y5_prob_blocco_modello,'*-','color', 'green');
hold on
plot(vettore_x_lambda,vettore_y6_prob_blocco_sim,'color', 'green');
% 
% 
plot(vettore_x_lambda,vettore_y7_prob_blocco_modello,'*-','color', 'black');
hold on
plot(vettore_x_lambda,vettore_y8_prob_blocco_sim,'color', 'black');

set(gca,'XTick',[5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,25,30,35]);
set(gca,'YTick',[0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.8 1]);
ylim([0 1])


vettore_traffico_perso_mod = 20*vettore_y1_prob_blocco_modello;
vettore_traffico_perso_sim2 = vettore_x_lambda.*vettore_y2_prob_blocco_sim;

% 
figure;
plot(vettore_x_lambda,vettore_traffico_perso_mod,'*-','color', 'red');
hold on
% plot(vettore_x_lambda,vettore_traffico_perso_sim,'color', 'red');
% hold on
plot(vettore_x_lambda,vettore_traffico_perso_sim2,'+-','color', 'red');
ylim([0 20])
xlim([0 50])
