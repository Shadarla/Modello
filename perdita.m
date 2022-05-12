%seed1
vettore_x_lambda = [5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,25,30];
%vettore_y1_prob_blocco_modello = [0.0002275, 0.0006709, 0.001569, 0.003102, 0.0053981, 0.00850898, 0.0124095, 0.017011, 0.0221872, 0.0277946, 0.0336937, 0.0397591, 0.0458846, 0.0519854, 0.0579963, 0.0638694,0.090];
vettore_y1_prob_blocco_modello = [0.00027, 0.00087, 0.00217, 0.00456, 0.00842, 0.0140, 0.0215, 0.0310, 0.0423, 0.0553, 0.0697, 0.0852, 0.1016, 0.1186, 0.1360, 0.1535, 0.2395, 0.3166];
vettore_y2_prob_blocco_sim = [0, 0, 0, 0, 0, 0, 0.001, 0.002, 0.002, 0.006, 0.012, 0.022, 0.036, 0.053, 0.072, 0.092, 0.2002, 0.304];
vettore_traffico_perso_mod = 20.*vettore_y1_prob_blocco_modello;
vettore_traffico_perso_sim = 20.*vettore_y2_prob_blocco_sim;

plot(vettore_x_lambda,vettore_y1_prob_blocco_modello);
hold on
plot(vettore_x_lambda,vettore_y2_prob_blocco_sim);
set(gca,'XTick',[5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,25,30]);
set(gca,'YTick',[0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.8 1]);
ylim([0 1])
