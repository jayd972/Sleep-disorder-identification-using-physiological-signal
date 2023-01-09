R_n = R_n(:,1:6);
W_n = W_n(:,1:6);
S1_n = S1_n(:,1:6);
S2_n = S2_n(:,1:6);
S3_n = S3_n(:,1:6);
S4_n = S4_n(:,1:6);


R_narco = R_narco(:,1:6);
W_narco = W_narco(:,1:6);
S1_narco = S1_narco(:,1:6);
S2_narco = S2_narco(:,1:6);
S3_narco = S3_narco(:,1:6);
S4_narco = S4_narco(:,1:6);

R_n(:,7) = 0;
W_n(:,7) = 0;
S1_n(:,7) = 0;
S2_n(:,7) = 0;
S3_n(:,7) = 0;
S4_n(:,7) = 0;



R_narco(:,7) = 1;
W_narco(:,7) = 1;
S1_narco(:,7) = 1;
S2_narco(:,7) = 1;
S3_narco(:,7) = 1;
S4_narco(:,7) = 1;

R_n_vs_narco = [R_n;R_narco];
W_n_vs_narco = [W_n;W_narco];
S1_n_vs_narco = [S1_n;S1_narco];
S2_n_vs_narco = [S2_n;S2_narco];
S3_n_vs_narco = [S3_n;S3_narco];
S4_n_vs_narco = [S4_n;S4_narco];
