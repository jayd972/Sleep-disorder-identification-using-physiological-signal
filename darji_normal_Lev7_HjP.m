tic
% 1. A code to stuff all the patient files in a file named normal(a cell).

normal{1,1} = char('n1.txt'); normal{2,1} = char('n2.txt');
normal{3,1} = char('n3.txt'); normal{4,1} = char('n4.txt');
normal{5,1} = char('n5.txt'); normal{6,1} = char('n8.txt');
normal{7,1} = char('n10.txt'); normal{8,1} = char('n11.txt');
normal{9,1} = char('n16.txt');

% 2. Starting points for all normal patients.jjjfjf

normal{1,2} = 107521;   normal{2,2} = 5744641;
normal{3,2} = 1551361;  normal{4,2} = 342001;
normal{5,2} = 614401;   normal{6,2} = 1;
normal{7,2} = 1843201;  normal{8,2} = 15361;
normal{9,2} = 3001;

%let's start
all_normal_combined_Hjorth = [];
for ij=[2 3 5 7 8]
        ScoreReader
        clearvars -except ij file hyp normal normal_results all_normal_combined_Hjorth nor_features_combined

        edf_file = file;
        edf_file([end-2 end-1 end]) = char('edf');
        disp('Processing...')
        [hdr, record] = edfread(edf_file);

        for i=1:length(hdr.label)
        if contains(hdr.label{1,i},'ROCLOC')
        ch=i;
        end
        end

 
        record = normalizeData(record(ch,normal{ij,2}:end));       %starting point
        
        
        %record = downsample(record,4);
        %load('FIR_window_512hz_1_2_30_35.mat');                             %prerequisite
       % filteredsignal = filter(Hd,record);


        %plotting FFT

        % L =length(filteredsignal);
%         n = 2^nextpow2(L);
%         Y = fft(sb1,n);
%         P2 = abs(Y/L);
%         P1 = P2(:, 1:n/2+1);
%         P1(:,2:end-1) = 2*P1(:,2:end-1);
%         Fs = 4;
%         plot(0:(Fs/n):(Fs/2-Fs/n),P1(1,1:n/2))

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

     

    name = normal{ij,1};
    name([end-3 end-2 end-1 end])=[];
    normal_results{ij,1} = name;
    normal_results{ij,2} = combined_Hjorth;
    all_normal_combined_Hjorth = cat(1,all_normal_combined_Hjorth,combined_Hjorth);
    
    sprintf('Completed %s',name)

end
 cd('E:\darji')
 save Normal_Results_Hjorth_512hz normal normal_results all_normal_combined_Hjorth
 

disp('Normal completed successfully.')
toc   





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






tic

%Narclepsy files
narco = {'narco1.txt',15360; 'narco2.txt',0; 'narco3.txt', 15360; 'narco4.txt',61140; 'narco5.txt',15360};

%let's start
all_narco_combined_Hjorth = [];
for ij=1:size(narco,1)
        ScoreReader_narco
        clearvars -except ij file hyp narco narco_results all_narco_combined_Hjorth nor_features_combined

        edf_file = file;
        edf_file([end-2 end-1 end]) = char('edf');
        disp('Processing...')
        [hdr, record] = edfread(edf_file);

        for i=1:length(hdr.label)
        if contains(hdr.label{1,i},'EMG1EMG2')
        ch=i;
        end
        end

        if narco{ij,2}~=0
            record = normalize(record(ch,narco{ij,2}:end));
        else
            record = normalize(record(ch,:));
        end
        
        
        %record = downsample(record,4);
%         load('FIR_window_512hz_1_2_30_35.mat');                             %prerequisite
%         filteredsignal = filter(Hd,record);


        %plotting FFT

        % L =length(filteredsignal);
        % n = 2^nextpow2(L);
        % Y = fft(filteredsignal,n);
        % P2 = abs(Y/L);
        % P1 = P2(:, 1:n/2+1);
        % P1(:,2:end-1) = 2*P1(:,2:end-1);
        % Fs = 128;
        % plot(0:(Fs/n):(Fs/2-Fs/n),P1(1,1:n/2))

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


    name = narco{ij,1};
    name([end-3 end-2 end-1 end])=[];
    narco_results{ij,1} = name;
    narco_results{ij,2} = combined_Hjorth;
    all_narco_combined_Hjorth = cat(1,all_narco_combined_Hjorth,combined_Hjorth);
    
    sprintf('Completed %s',name)

end
 cd('E:\darji')
 save Narco_Results_Hjorth_512hz narco narco_results all_narco_combined_Hjorth
 

disp('Narcolepsy completed successfully.')
toc    




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('Normal_Results_Hjorth_512hz.mat')
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
%%%%%%%%%%%%%%%%%%  Narcolepsy  %%%%%%%%%%%%%%%%%%%%%%%%%
load('Narco_Results_Hjorth_512hz.mat')
Y = all_narco_combined_Hjorth;
W_narco = []; N1_narco = []; N2_narco=[]; N3_narco = []; R_narco=[];
for i=1:length(Y)
if Y(i,25)==0
W_narco = cat(1,W_narco,Y(i,:));
elseif Y(i,25)==1
N1_narco = cat(1,N1_narco,Y(i,:));
elseif Y(i,25)==2
N2_narco = cat(1,N2_narco,Y(i,:));
elseif Y(i,25)==3
N3_narco = cat(1,N3_narco,Y(i,:));
elseif Y(i,25)==5
R_narco = cat(1,R_narco,Y(i,:));
end
end
W_narco(:,25) = 10; N1_narco(:,25) = 11; N2_narco(:,25)=12;
N3_narco(:,25) = 13; R_narco(:,25)= 15;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
W_normal_vs_narco = [W_normal; W_narco];
N1_normal_vs_narco = [N1_normal; N1_narco];
N2_normal_vs_narco = [N2_normal; N2_narco];
N3_normal_vs_narco = [N3_normal; N3_narco];
R_normal_vs_narco = [R_normal; R_narco];
cd('E:\darji')
save Normal_vs_Narco_Hjorth_512hz W_normal_vs_narco N1_normal_vs_narco N2_normal_vs_narco N3_normal_vs_narco R_normal_vs_narco

