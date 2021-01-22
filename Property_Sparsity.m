function Property_Sparsity

    Test_num=10000;

%     load(['IEEE1354_Capacity_1_D_',num2str(Test_num),'_New.mat'],'D')
%     load(['IEEE1354_Capacity_1_A_1_',num2str(Test_num),'_New.mat'],'A_1')
%     load(['IEEE1354_Capacity_1_A_2_',num2str(Test_num),'_New.mat'],'A_2')

%     load(['IEEE118_Capacity_1_2_InitFail_3_D_SampleNum_',num2str(Test_num),'_Balanced_DiffTrainingSet_2_var.mat'],'D')
%     load(['IEEE118_Capacity_1_2_InitFail_3_A_1_SampleNum_',num2str(Test_num),'_Balanced_DiffTrainingSet_2_var.mat'],'A_1')
%     load(['IEEE118_Capacity_1_2_InitFail_3_A_2_SampleNum_',num2str(Test_num),'_Balanced_DiffTrainingSet_2_var.mat'],'A_2')

%     load(['IEEE118_Capacity_1_5_InitFail_3_D_SampleNum_',num2str(Test_num),'_Balanced_New.mat'],'D')
%     load(['IEEE118_Capacity_1_5_InitFail_3_A_1_SampleNum_',num2str(Test_num),'_Balanced_New.mat'],'A_1')
%     load(['IEEE118_Capacity_1_5_InitFail_3_A_2_SampleNum_',num2str(Test_num),'_Balanced_New.mat'],'A_2')

%     load(['IEEE118_Capacity_1_2_InitFail_3_D_SampleNum_',num2str(Test_num),'_Balanced_DiffTrainingSet_RateFix_0_1.mat'],'D')
%     load(['IEEE118_Capacity_1_2_InitFail_3_A_1_SampleNum_',num2str(Test_num),'_Balanced_DiffTrainingSet_RateFix_0_1.mat'],'A_1')
%     load(['IEEE118_Capacity_1_2_InitFail_3_A_2_SampleNum_',num2str(Test_num),'_Balanced_DiffTrainingSet_RateFix_0_1.mat'],'A_2')

    load(['IEEE2383_Capacity_1_Flow_1_25_InitFail_2_D_SampleNum_5000_Balanced_New_AC.mat'],'D')
    load(['IEEE2383_Capacity_1_Flow_1_25_InitFail_2_A_1_SampleNum_5000_Balanced_New_AC.mat'],'A_1')
    load(['IEEE2383_Capacity_1_Flow_1_25_InitFail_2_A_2_SampleNum_5000_Balanced_New_AC.mat'],'A_2')

%     load(['IEEE2383_Capacity_1_InitFail_2_D_SampleNum_5000_New_2.mat'],'D')
%     load(['IEEE2383_Capacity_1_InitFail_2_A_1_SampleNum_5000_New_2.mat'],'A_1')
%     load(['IEEE2383_Capacity_1_InitFail_2_A_2_SampleNum_5000_New_2.mat'],'A_2')

%     load(['IEEE1354_Capacity_1_InitFail_2_D_SampleNum_',num2str(Test_num),'_Balanced.mat'],'D')
%     load(['IEEE1354_Capacity_1_InitFail_2_A_1_SampleNum_',num2str(Test_num),'_Balanced.mat'],'A_1')
%     load(['IEEE1354_Capacity_1_InitFail_2_A_2_SampleNum_',num2str(Test_num),'_Balanced.mat'],'A_2')
    
    for i=1:size(A_2,1)
        A_2(i,i)=0; 
        D(i,i)=0;
    end

%     load(['IEEE118_Capacity_1_2_Flow_AC_A_1_10000.mat'],'A_1')
%     load(['IEEE118_Capacity_1_2_Flow_AC_A_2_10000.mat'],'A_2')
%     
%     load(['IEEE118_Capacity_1_2_Flow_A_1_10000_Random_0_8_1_0_Prob_0_7.mat'],'A_1')
%     load(['IEEE118_Capacity_1_2_Flow_A_2_10000_Random_0_8_1_0_Prob_0_7.mat'],'A_2') 

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
    sum_tmp_1';
    sum_tmp_2=zeros(M,1);
    M_2=D.*P_2;
    for i=1:M
        sum_tmp_2(i,1)=sum(M_2(i,:)); 
    end
    sum_tmp_2';
    [sum_tmp_1,sum_tmp_2];
    b_est=zeros(M,1);
    for i=1:M
        b_est(i,1)=D(i,:)*P_2(i,:)'; 
    end
    b_est;
    
%     Test_cand=[0.01;0.1];
%     for k=1:2
%         D_new=D>Test_cand(k);
%         subplot(1,2,k)
%         spy(D_new)
%     end

    D_new=D>0.01;
    spy(D_new)
    
end

