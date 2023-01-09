

tic

%Narclepsy files
sdb = {char('sdb3.txt'), 3578880 };
%let's start
all_sdb_combined_Hjorth = [];
for ij=1:size(sdb,1)
        patients = sdb;
        ScoreReader_general
        clearvars -except ij file hyp sdb sdb_results all_sdb_combined_Hjorth nor_features_combined

        edf_file = file;
        edf_file([end-2 end-1 end]) = char('edf');
        disp('Processing...')
        [hdr, record] = edfread(edf_file);

        for i=1:length(hdr.label)
        if contains(hdr.label{1,i},'EMG1EMG2')
        ch=i;
        end
        end

        if sdb{ij,2}~=0
            record = normalize(record(ch,sdb{ij,2}:end));
        else
            record = normalize(record(ch,:));
        end
         
        %segregation
        w = [];s1 = []; s2 = [];s3=[];r=[];
        j = 1;
        for i=1:length(hyp)
        if hyp(i,1)==0
            w(j) = hyp(i,2);
            j = j+1;
        elseif hyp(i,1)==1
            s1(j) = hyp(i,2);
            j=j+1;
        elseif hyp(i,1)==2
            s2(j) = hyp(i,2);
            j = j+1;
        elseif (hyp(i,1)==3) || (hyp(i,1)==4)
            s3(j) = hyp(i,2);
            j = j + 1;
     
        elseif hyp(i,1)==5
            r(j) = hyp(i,2);
            j = j+1;
        else
            continue
        end
        end
        
        fs = hdr.frequency(ch);   %change frequency

        w=w*fs; s1=s1*fs; s2=s2*fs ;s3=s3*fs ;r=r*fs;

        %sleep
        W = [];
        for i = 1:length(w)-3
            if w(i)~=0
                w1 = record(:,((w(i)+1)-(fs*30)):w(i));
                W = cat(2,W,w1);
            end
        end


        S1 = [];
        for i = 1:length(s1)
            if s1(i)~=0
                s11 = record(:,((s1(i)+1)-(fs*30)):s1(i));
                S1 = cat(2,S1,s11);
            end
        end


        S2 = [];
        %if ij==9                            %for n16
        %    for i = 1:length(s2)-1
        %        if s2(i)~=0 
        %            s21 = record(:,((s2(i)+1)-(fs*30)):s2(i));
        %            S2 = cat(2,S2,s21);
        %        end
        %    end
        %else
            for i = 1:length(s2)
                if s2(i)~=0 
                    s21 = record(:,((s2(i)+1)-(fs*30)):s2(i));
                    S2 = cat(2,S2,s21);
                end
            end
        %end
    


        S3 = [];
        for i = 1:length(s3)
            if s3(i)~=0
                s31 = record(:,((s3(i)+1)-(fs*30)):s3(i));
                S3 = cat(2,S3,s31);
            end
        end


 
        R = [];
        for i = 1:length(r)
            if r(i)~=0
                r1 = record(:,((r(i)+1)-(fs*30)):r(i));
                R = cat(2,R,r1);
            end
        end

%Channel_EMG1EMG2_General
    w_ch = W(1,:); 
    s1_ch = S1(1,:); 
    s2_ch = S2(1,:); 
    s3_ch = S3(1,:);
    %s4_ch = S4(1,:); 
    r_ch = R(1,:);

    reshaped_w_ch = reshape(w_ch,30*fs,length(w_ch)/(30*fs));
    reshaped_s1_ch = reshape(s1_ch,30*fs,length(s1_ch)/(30*fs));
    reshaped_s2_ch = reshape(s2_ch,30*fs,length(s2_ch)/(30*fs));
    reshaped_s3_ch = reshape(s3_ch,30*fs,length(s3_ch)/(30*fs));
    %reshaped_s4_ch = reshape(s4_ch,30*fs,length(s4_ch)/(30*fs));
    reshaped_r_ch = reshape(r_ch,30*fs,length(r_ch)/(30*fs));



    
    load('myExample-5');   

    W_EMG1EMG2_Hjorth=[];
    [Lo_D,Hi_D,Lo_R,Hi_R] = biorfilt(h0',f0');  %filter
    for i = 1:length(w_ch)/(30*fs) 
      [C,L] = wavedec(reshaped_w_ch(:,i),7,Lo_D,Hi_D);  
      sb1 = appcoef(C,L,Lo_R,Hi_R,7);
      [sb2,sb3,sb4,sb5,sb6,sb7,sb8] = detcoef(C,L,[1 2 3 4 5 6 7]);
      [a1, m1, c1] = HjorthParameters(sb1);
      [a2, m2, c2] = HjorthParameters(sb2);
      [a3, m3, c3] = HjorthParameters(sb3);
      [a4, m4, c4] = HjorthParameters(sb4);
      [a5, m5, c5] = HjorthParameters(sb5);
      [a6, m6, c6] = HjorthParameters(sb6);
      [a7, m7, c7] = HjorthParameters(sb7);
      [a8, m8, c8] = HjorthParameters(sb8);


     X = [a1 a2 a3 a4 a5 a6 a7 a8 m1 m2 m3 m4 m5 m6 m7 m8 c1 c2 c3 c4 c5 c6 c7 c8];
      W_EMG1EMG2_Hjorth = cat(1,W_EMG1EMG2_Hjorth,X);
    end
    W_EMG1EMG2_Hjorth(:,25)=0;


    S1_EMG1EMG2_Hjorth=[];
    for i = 1:length(s1_ch)/(30*fs) 
      [C,L] = wavedec(reshaped_s1_ch(:,i),7,Lo_D,Hi_D);
      sb1 = appcoef(C,L,Lo_R,Hi_R,7);
      [sb2,sb3,sb4,sb5,sb6,sb7,sb8] = detcoef(C,L,[1 2 3 4 5 6 7]);
      [a1, m1, c1] = HjorthParameters(sb1);
      [a2, m2, c2] = HjorthParameters(sb2);
      [a3, m3, c3] = HjorthParameters(sb3);
      [a4, m4, c4] = HjorthParameters(sb4);
      [a5, m5, c5] = HjorthParameters(sb5);
      [a6, m6, c6] = HjorthParameters(sb6);
      [a7, m7, c7] = HjorthParameters(sb7);
      [a8, m8, c8] = HjorthParameters(sb8);


     X = [a1 a2 a3 a4 a5 a6 a7 a8 m1 m2 m3 m4 m5 m6 m7 m8 c1 c2 c3 c4 c5 c6 c7 c8];
      S1_EMG1EMG2_Hjorth = cat(1,S1_EMG1EMG2_Hjorth,X);
    end
    S1_EMG1EMG2_Hjorth(:,25)=1;


    S2_EMG1EMG2_Hjorth=[];
    for i = 1:length(s2_ch)/(30*fs) 
      [C,L] = wavedec(reshaped_s2_ch(:,i),7,Lo_D,Hi_D);
      sb1 = appcoef(C,L,Lo_R,Hi_R,7);
      [sb2,sb3,sb4,sb5,sb6,sb7,sb8] = detcoef(C,L,[1 2 3 4 5 6 7]);
      [a1, m1, c1] = HjorthParameters(sb1);
      [a2, m2, c2] = HjorthParameters(sb2);
      [a3, m3, c3] = HjorthParameters(sb3);
      [a4, m4, c4] = HjorthParameters(sb4);
      [a5, m5, c5] = HjorthParameters(sb5);
      [a6, m6, c6] = HjorthParameters(sb6);
      [a7, m7, c7] = HjorthParameters(sb7);
      [a8, m8, c8] = HjorthParameters(sb8);


     X = [a1 a2 a3 a4 a5 a6 a7 a8 m1 m2 m3 m4 m5 m6 m7 m8 c1 c2 c3 c4 c5 c6 c7 c8];
      S2_EMG1EMG2_Hjorth = cat(1,S2_EMG1EMG2_Hjorth,X);
    end
    S2_EMG1EMG2_Hjorth(:,25)=2;


    S3_EMG1EMG2_Hjorth=[];
    for i = 1:length(s3_ch)/(30*fs) 
      [C,L] = wavedec(reshaped_s3_ch(:,i),7,Lo_D,Hi_D);
      sb1 = appcoef(C,L,Lo_R,Hi_R,7);
      [sb2,sb3,sb4,sb5,sb6,sb7,sb8] = detcoef(C,L,[1 2 3 4 5 6 7]);
      [a1, m1, c1] = HjorthParameters(sb1);
      [a2, m2, c2] = HjorthParameters(sb2);
      [a3, m3, c3] = HjorthParameters(sb3);
      [a4, m4, c4] = HjorthParameters(sb4);
      [a5, m5, c5] = HjorthParameters(sb5);
      [a6, m6, c6] = HjorthParameters(sb6);
      [a7, m7, c7] = HjorthParameters(sb7);
      [a8, m8, c8] = HjorthParameters(sb8);


     X = [a1 a2 a3 a4 a5 a6 a7 a8 m1 m2 m3 m4 m5 m6 m7 m8 c1 c2 c3 c4 c5 c6 c7 c8];
      S3_EMG1EMG2_Hjorth = cat(1,S3_EMG1EMG2_Hjorth,X);
    end
    S3_EMG1EMG2_Hjorth(:,25)=3;



    R_EMG1EMG2_Hjorth=[];
    for i = 1:length(r_ch)/(30*fs) 
      [C,L] = wavedec(reshaped_r_ch(:,i),7,Lo_D,Hi_D);
      sb1 = appcoef(C,L,Lo_R,Hi_R,7);
      [sb2,sb3,sb4,sb5,sb6,sb7,sb8] = detcoef(C,L,[1 2 3 4 5 6 7]);
      [a1, m1, c1] = HjorthParameters(sb1);
      [a2, m2, c2] = HjorthParameters(sb2);
      [a3, m3, c3] = HjorthParameters(sb3);
      [a4, m4, c4] = HjorthParameters(sb4);
      [a5, m5, c5] = HjorthParameters(sb5);
      [a6, m6, c6] = HjorthParameters(sb6);
      [a7, m7, c7] = HjorthParameters(sb7);
      [a8, m8, c8] = HjorthParameters(sb8);


     X = [a1 a2 a3 a4 a5 a6 a7 a8 m1 m2 m3 m4 m5 m6 m7 m8 c1 c2 c3 c4 c5 c6 c7 c8];
      R_EMG1EMG2_Hjorth = cat(1,R_EMG1EMG2_Hjorth,X);
    end
    R_EMG1EMG2_Hjorth(:,25)=5;

    combined_Hjorth = [W_EMG1EMG2_Hjorth; S1_EMG1EMG2_Hjorth; S2_EMG1EMG2_Hjorth; S3_EMG1EMG2_Hjorth; R_EMG1EMG2_Hjorth ];  


    name = sdb{ij,1};
    name([end-3 end-2 end-1 end])=[];
    sdb_results{ij,1} = name;
    sdb_results{ij,2} = combined_Hjorth;
    all_sdb_combined_Hjorth = cat(1,all_sdb_combined_Hjorth,combined_Hjorth);
    
    sprintf('Completed %s',name)

end
 cd('E:\darji')
 save SDB_Results_HjP sdb sdb_results all_sdb_combined_Hjorth
 

disp('SDBlepsy completed successfully.')
toc    




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('Normal_Results_HjP.mat')
X = all_normal_combined_Hjorth;
W_normal = []; N1_normal = []; N2_normal=[]; N3_normal = []; R_normal=[];
for i=1:length(X)
if X(i,25)==0
W_normal = cat(1,W_normal,X(i,:));
elseif X(i,25)==1
N1_normal = cat(1,N1_normal,X(i,:));
elseif X(i,25)==2
N2_normal = cat(1,N2_normal,X(i,:));
elseif X(i,25)==3
N3_normal = cat(1,N3_normal,X(i,:));
elseif X(i,25)==5
R_normal = cat(1,R_normal,X(i,:));
end
end
W_normal(:,25) = 0; N1_normal(:,25) = 1; N2_normal(:,25)=2;
N3_normal(:,25) = 3; R_normal(:,25)= 5;
%%%%%%%%%%%%%%%%%%  SDBlepsy  %%%%%%%%%%%%%%%%%%%%%%%%%
load('SDB_Results_HjP.mat')
Y = all_sdb_combined_Hjorth;
W_sdb = []; N1_sdb = []; N2_sdb=[]; N3_sdb = []; R_sdb=[];
for i=1:length(Y)
if Y(i,25)==0
W_sdb = cat(1,W_sdb,Y(i,:));
elseif Y(i,25)==1
N1_sdb = cat(1,N1_sdb,Y(i,:));
elseif Y(i,25)==2
N2_sdb = cat(1,N2_sdb,Y(i,:));
elseif Y(i,25)==3
N3_sdb = cat(1,N3_sdb,Y(i,:));
elseif Y(i,25)==5
R_sdb = cat(1,R_sdb,Y(i,:));
end
end
W_sdb(:,25) = 10; N1_sdb(:,25) = 11; N2_sdb(:,25)=12;
N3_sdb(:,25) = 13; R_sdb(:,25)= 15;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
W_normal_vs_sdb = [W_normal; W_sdb];
N1_normal_vs_sdb = [N1_normal; N1_sdb];
N2_normal_vs_sdb = [N2_normal; N2_sdb];
N3_normal_vs_sdb = [N3_normal; N3_sdb];
R_normal_vs_sdb = [R_normal; R_sdb];

temp_LSS_normal=[N1_normal; N2_normal]; temp_LSS_normal(:,25)=0;
temp_LSS_sdb=[N1_sdb; N2_sdb]; temp_LSS_sdb(:,25)=1;
LSS_normal_vs_sdb = [temp_LSS_normal; temp_LSS_sdb];

temp_NREM_normal=[N1_normal; N2_normal; N3_normal]; temp_NREM_normal(:,25)=0;
temp_NREM_sdb=[N1_sdb; N2_sdb; N3_sdb]; temp_NREM_sdb(:,25)=1;
NREM_normal_vs_sdb = [temp_NREM_normal; temp_NREM_sdb];

all_normal_combined_Hjorth(:,25)=0;
all_sdb_combined_Hjorth(:,25)=1;
Allstages_normal_vs_sdb = [all_normal_combined_Hjorth; all_sdb_combined_Hjorth];

cd('E:\darji')
save Normal_vs_SDB_HjP W_normal_vs_sdb N1_normal_vs_sdb N2_normal_vs_sdb N3_normal_vs_sdb R_normal_vs_sdb LSS_normal_vs_sdb NREM_normal_vs_sdb Allstages_normal_vs_sdb



