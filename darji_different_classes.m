W_n = [];
S1_n = [];
S2_n = [];
S3_n = [];
S4_n = [];
R_n = [];


% W_narco = [];
% S1_narco = [];
% S2_narco = [];
% S3_narco = [];
% S4_narco = [];
% R_narco = [];


for i = 1:size(combined_LLE,1)
    if combined_LLE(i,7) == 0
        A = combined_LLE(i,:);
        W_n = cat(1,W_n,A);
%         W_narco = cat(1,W_narco,A);
    elseif combined_LLE(i,7) == 1
        B = combined_LLE(i,:);
        S1_n = cat(1,S1_n,B);
%         S1_narco = cat(1,S1_narco,B);
    elseif combined_LLE(i,7) == 2
        C = combined_LLE(i,:);
        S2_n = cat(1,S2_n,C);
%         S2_narco = cat(1,S2_narco,C);
    elseif combined_LLE(i,7) == 3
        D = combined_LLE(i,:);
        S3_n = cat(1,S3_n,D);
%         S3_narco = cat(1,S3_narco,D);
    elseif combined_LLE(i,7) == 4
        E = combined_LLE(i,:);
        S4_n = cat(1,S4_n,E);
%         S4_narco = cat(1,S4_narco,E);
    elseif combined_LLE(i,7) == 5
        F = combined_LLE(i,:);
         R_n = cat(1,R_n,F);
%         R_narco = cat(1,R_narco,F);
    end
end


 save N_cl W_n S1_n S2_n S3_n S4_n R_n
% save Narco_cl W_narco S1_narco S2_narco S3_narco S4_narco R_narco