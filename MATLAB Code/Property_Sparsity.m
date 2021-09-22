    %% Influence Sparsity Property 

    load(['IEEE2383_Capacity_1_Flow_1_25_InitFail_2_D_SampleNum_5000_Balanced_New_AC.mat'],'D')
    load(['IEEE2383_Capacity_1_Flow_1_25_InitFail_2_A_1_SampleNum_5000_Balanced_New_AC.mat'],'A_1')
    load(['IEEE2383_Capacity_1_Flow_1_25_InitFail_2_A_2_SampleNum_5000_Balanced_New_AC.mat'],'A_2')
    
    for i=1:size(A_2,1)
        A_2(i,i)=0; 
        D(i,i)=0;
    end

    P_1=A_1';
    P_2=A_2';
    M=size(D,1);

    %Calculate M_est,b_est
    M_est=D.*(P_1-P_2);
    sum_tmp_1=zeros(M,1);
    M_1=D.*P_1;
    for i=1:M
        sum_tmp_1(i,1)=sum(M_1(i,:)); 
    end
    
    sum_tmp_2=zeros(M,1);
    M_2=D.*P_2;
    for i=1:M
        sum_tmp_2(i,1)=sum(M_2(i,:)); 
    end
    
    b_est=zeros(M,1);
    for i=1:M
        b_est(i,1)=D(i,:)*P_2(i,:)'; 
    end

    D_new=D>0.01;
    spy(D_new)
