        
%     load('IEEE1354_Capacity_1_InitFail_2_New_Train_Tmp.mat','cascade_train')
%     cascade_1=cascade_train;
%     size(cascade_1);
%     load('IEEE2383_Capacity_1_InitFail_2_New_Train_Tmp_2.mat','cascade_train')
%     cascade_2=cascade_train;    
%     load('IEEE2383_Capacity_1_InitFail_2_New_Test_Tmp_2.mat','cascade_test')
%     cascade_3=cascade_test;
%     cascade=[cascade_2,cascade_3];

%     load('IEEE1354_Capacity_1_InitFail_2_Train_Balanced.mat','cascade_train')
%     cascade=cascade_train;

    load('IEEE1354_Capacity_1_Flow_1_5_InitFail_2_Train_Balanced_New_AC.mat','cascade_link_cum')
    cascade_train=cascade_link_cum;
    cascade=cascade_train;
%     load('IEEE2383_Capacity_1_Flow_1_25_InitFail_2_Train_Balanced_New_AC.mat','cascade_train')
%     cascade=cascade_train;

%     load('IEEE118_Capacity_1_8_InitFail_3_Train_Balanced_New.mat','cascade_train')
%     cascade=cascade_train;
%     load('IEEE2383_Capacity_1_InitFail_2_Train_Balanced_AC.mat','cascade_train')
%     cascade=cascade_train;
%     load('IEEE1354_Capacity_1_Flow_1_5_InitFail_2_Train_Balanced_New_AC.mat','cascade_train')
%     cascade=cascade_train;

%     size(cascade)
%     load('IEEE2383_Capacity_1_Flow_1_InitFail_2_Train_Balanced_New_AC.mat','cascade_train')
%     cascade_train=cascade_link_cum;
%     cascade=cascade_train;

%     load('IEEE3012_Capacity_1_Flow_1_5_InitFail_2_Train_Balanced_New.mat','cascade_train')
%     cascade=cascade_train;

%     load('IEEE118_Capacity_1_8_InitFail_2_Train_Balanced_New.mat','cascade_train')
%     cascade=cascade_train;

%     Test_num=10000;
%     load(['IEEE118_Capacity_1_2_InitFail_3_D_SampleNum_',num2str(Test_num),'_Balanced_DiffTrainingSet_2_var.mat'],'D')
%     load(['IEEE118_Capacity_1_2_InitFail_3_A_1_SampleNum_',num2str(Test_num),'_Balanced_DiffTrainingSet_2_var.mat'],'A_1')
%     load(['IEEE118_Capacity_1_2_InitFail_3_A_2_SampleNum_',num2str(Test_num),'_Balanced_DiffTrainingSet_2_var.mat'],'A_2')

%     Test_num=10000;
%     load(['IEEE118_Capacity_1_8_InitFail_3_D_SampleNum_',num2str(Test_num),'_Balanced_New.mat'],'D')
%     load(['IEEE118_Capacity_1_8_InitFail_3_A_1_SampleNum_',num2str(Test_num),'_Balanced_New.mat'],'A_1')
%     load(['IEEE118_Capacity_1_8_InitFail_3_A_2_SampleNum_',num2str(Test_num),'_Balanced_New.mat'],'A_2')

%     Test_num=10000;
%     load(['IEEE118_Capacity_1_8_InitFail_2_D_SampleNum_',num2str(Test_num),'_Balanced_New.mat'],'D')
%     load(['IEEE118_Capacity_1_8_InitFail_2_A_1_SampleNum_',num2str(Test_num),'_Balanced_New.mat'],'A_1')
%     load(['IEEE118_Capacity_1_8_InitFail_2_A_2_SampleNum_',num2str(Test_num),'_Balanced_New.mat'],'A_2')

%     
%     load(['IEEE2383_Capacity_1_InitFail_2_D_SampleNum_5000_New_2.mat'],'D')
%     load(['IEEE2383_Capacity_1_InitFail_2_A_1_SampleNum_5000_New_2.mat'],'A_1')
%     load(['IEEE2383_Capacity_1_InitFail_2_A_2_SampleNum_5000_New_2.mat'],'A_2')
% 
%     load(['IEEE2383_Capacity_1_InitFail_2_D_SampleNum_5000_New_2.mat'],'D')
%     load(['IEEE2383_Capacity_1_InitFail_2_A_1_SampleNum_5000_New_2.mat'],'A_1')
%     load(['IEEE2383_Capacity_1_InitFail_2_A_2_SampleNum_5000_New_2.mat'],'A_2')

%     K=10000;
%     load(['IEEE2383_Capacity_1_Flow_1_5_InitFail_2_D_SampleNum_',num2str(K),'_Balanced_New.mat'],'D')
%     load(['IEEE2383_Capacity_1_Flow_1_5_InitFail_2_A_1_SampleNum_',num2str(K),'_Balanced_New.mat'],'A_1')
%     load(['IEEE2383_Capacity_1_Flow_1_5_InitFail_2_A_2_SampleNum_',num2str(K),'_Balanced_New.mat'],'A_2')

%     K=5000;
%     load(['IEEE2383_Capacity_1_Flow_1_25_InitFail_2_D_SampleNum_',num2str(K),'_Balanced_New_AC.mat'],'D')
%     load(['IEEE2383_Capacity_1_Flow_1_25_InitFail_2_A_1_SampleNum_',num2str(K),'_Balanced_New_AC.mat'],'A_1')
%     load(['IEEE2383_Capacity_1_Flow_1_25_InitFail_2_A_2_SampleNum_',num2str(K),'_Balanced_New_AC.mat'],'A_2')

%     load(['IEEE1354_Capacity_1_InitFail_2_D_SampleNum_5000_Balanced.mat'],'D')
%     load(['IEEE1354_Capacity_1_InitFail_2_A_1_SampleNum_5000_Balanced.mat'],'A_1')
%     load(['IEEE1354_Capacity_1_InitFail_2_A_2_SampleNum_5000_Balanced.mat'],'A_2')

    K=5000;
    load(['IEEE1354_Capacity_1_Flow_1_5_InitFail_2_D_SampleNum_',num2str(K),'_Balanced_New_AC.mat'],'D')
    load(['IEEE1354_Capacity_1_Flow_1_5_InitFail_2_A_1_SampleNum_',num2str(K),'_Balanced_New_AC.mat'],'A_1')
    load(['IEEE1354_Capacity_1_Flow_1_5_InitFail_2_A_2_SampleNum_',num2str(K),'_Balanced_New_AC.mat'],'A_2')

%     K=5000;
%     load(['IEEE1354_Capacity_1_InitFail_2_D_SampleNum_',num2str(K),'_Balanced.mat'],'D')
%     load(['IEEE1354_Capacity_1_InitFail_2_A_1_SampleNum_',num2str(K),'_Balanced.mat'],'A_1')
%     load(['IEEE1354_Capacity_1_InitFail_2_A_2_SampleNum_',num2str(K),'_Balanced.mat'],'A_2')

%     K=5000;
%     load(['IEEE3012_Capacity_1_Flow_1_5_InitFail_2_D_SampleNum_',num2str(K),'_Balanced_New.mat'],'D')
%     load(['IEEE3012_Capacity_1_Flow_1_5_InitFail_2_A_1_SampleNum_',num2str(K),'_Balanced_New.mat'],'A_1')
%     load(['IEEE3012_Capacity_1_Flow_1_5_InitFail_2_A_2_SampleNum_',num2str(K),'_Balanced_New.mat'],'A_2')

    for i=1:size(A_2,1)
        A_2(i,i)=0;
    end

    P_1=A_1';
    P_2=A_2';    
    
    [P_1(:,20),P_2(:,20)];

    K_prev=50000;
    
%     Train_num=40000;
%     load(['IEEE2383_Capacity_1_SampleNum_',num2str(Train_num),'_Initial_state_A_2.mat'],'Initial_state')
%     load(['IEEE2383_Capacity_1_SampleNum_',num2str(Train_num),'_Final_state_A_2.mat'],'Final_state')
%     load(['IEEE2383_Capacity_1_SampleNum_',num2str(Train_num),'_Epsilon_opt_A_2.mat'],'Epsilon_opt')
%     
%     Initial_state_1=Initial_state;
%     Final_state_1=Final_state;
%     Epsilon_opt_1=Epsilon_opt;
%     
%     Train_num=20000;
%     load(['IEEE2383_Capacity_1_SampleNum_',num2str(Train_num),'_Initial_state_A_2_More.mat'],'Initial_state')
%     load(['IEEE2383_Capacity_1_SampleNum_',num2str(Train_num),'_Final_state_A_2_More.mat'],'Final_state')
%     load(['IEEE2383_Capacity_1_SampleNum_',num2str(Train_num),'_Epsilon_opt_A_2_More.mat'],'Epsilon_opt')
%     
%     Initial_state_2=Initial_state;
%     Final_state_2=Final_state;
%     Epsilon_opt_2=Epsilon_opt;
%     
%     Initial_state=[Initial_state_1,Initial_state_2];
%     Final_state=[Final_state_1,Final_state_2];
%     Epsilon_opt=[Epsilon_opt_1,Epsilon_opt_2];

%     Train_num=60000;
%     load(['IEEE1354_Capacity_1_SampleNum_',num2str(Train_num),'_Initial_state_A_Balanced.mat'],'Initial_state')
%     load(['IEEE1354_Capacity_1_SampleNum_',num2str(Train_num),'_Final_state_A_Balanced.mat'],'Final_state')
%     load(['IEEE1354_Capacity_1_SampleNum_',num2str(Train_num),'_Epsilon_opt_A_Balanced.mat'],'Epsilon_opt')
%   
%     Test_num=50000;
%     load(['IEEE118_Capacity_1_2_InitFail_3_SampleNum_',num2str(Test_num),'_Initial_state_Balanced_DiffTrainingSet_2_var.mat'],'Initial_state')
%     load(['IEEE118_Capacity_1_2_InitFail_3_SampleNum_',num2str(Test_num),'_Final_state_Balanced_DiffTrainingSet_2_var.mat'],'Final_state')
%     load(['IEEE118_Capacity_1_2_InitFail_3_SampleNum_',num2str(Test_num),'_Epsilon_opt_Balanced_DiffTrainingSet_2_var.mat'],'Epsilon_opt')
%     load(['IEEE118_Capacity_1_2_InitFail_3_SampleNum_',num2str(Test_num),'_Fail_Index_Set_Balanced_DiffTrainingSet.mat'],'Link_Fail_Index_Set')

%     Test_num=50000;
%     load(['IEEE118_Capacity_1_2_InitFail_3_SampleNum_',num2str(Test_num),'_Initial_state_Balanced_New_AC.mat'],'Initial_state')
%     load(['IEEE118_Capacity_1_2_InitFail_3_SampleNum_',num2str(Test_num),'_Final_state_Balanced_New_AC.mat'],'Final_state')
%     load(['IEEE118_Capacity_1_2_InitFail_3_SampleNum_',num2str(Test_num),'_Epsilon_opt_Balanced_New_AC.mat'],'Epsilon_opt')
%     load(['IEEE118_Capacity_1_2_InitFail_3_SampleNum_',num2str(Test_num),'_Fail_Index_Set_Balanced_DiffTrainingSet.mat'],'Link_Fail_Index_Set')

%     Test_num=50000;
%     load(['IEEE118_Capacity_1_5_InitFail_3_SampleNum_',num2str(Test_num),'_Initial_state_Balanced_New_AC.mat'],'Initial_state')
%     load(['IEEE118_Capacity_1_5_InitFail_3_SampleNum_',num2str(Test_num),'_Final_state_Balanced_New_AC.mat'],'Final_state')
%     load(['IEEE118_Capacity_1_5_InitFail_3_SampleNum_',num2str(Test_num),'_Epsilon_opt_Balanced_New_AC.mat'],'Epsilon_opt')

%     Test_num=50000;
%     load(['IEEE118_Capacity_1_2_InitFail_3_SampleNum_',num2str(Test_num),'_Initial_state_Balanced_SquareRoot.mat'],'Initial_state')
%     load(['IEEE118_Capacity_1_2_InitFail_3_SampleNum_',num2str(Test_num),'_Final_state_Balanced_SquareRoot.mat'],'Final_state')
%     load(['IEEE118_Capacity_1_2_InitFail_3_SampleNum_',num2str(Test_num),'_Epsilon_opt_Balanced_SquareRoot.mat'],'Epsilon_opt')
%     load(['IEEE118_Capacity_1_2_InitFail_3_SampleNum_',num2str(Test_num),'_Fail_Index_Set_Balanced_DiffTrainingSet.mat'],'Link_Fail_Index_Set')

%     Test_num=30000;
%     load(['IEEE2383_Capacity_1_Flow_1_25_InitFail_2_SampleNum_',num2str(Test_num),'_Initial_state_Balanced_AC.mat'],'Initial_state')
%     load(['IEEE2383_Capacity_1_Flow_1_25_InitFail_2_SampleNum_',num2str(Test_num),'_Final_state_Balanced_AC.mat'],'Final_state')
%     load(['IEEE2383_Capacity_1_Flow_1_25_InitFail_2_SampleNum_',num2str(Test_num),'_Epsilon_opt_Balanced_AC.mat'],'Epsilon_opt')
%     load(['IEEE2383_Capacity_1_Flow_1_5_InitFail_2_SampleNum_',num2str(Test_num),'_Beta_Balanced.mat'],'Beta_record')

%     Test_num=50000;
%     load(['IEEE1354_Capacity_1_Flow_1_5_InitFail_2_SampleNum_',num2str(Test_num),'_Initial_state_Balanced_New.mat'],'Initial_state')
%     load(['IEEE1354_Capacity_1_Flow_1_5_InitFail_2_SampleNum_',num2str(Test_num),'_Final_state_Balanced_New.mat'],'Final_state')
%     load(['IEEE1354_Capacity_1_Flow_1_5_InitFail_2_SampleNum_',num2str(Test_num),'_Epsilon_opt_Balanced_New.mat'],'Epsilon_opt')

    Test_num=50000;
    load(['IEEE1354_Capacity_1_Flow_1_5_InitFail_2_SampleNum_',num2str(Test_num),'_Initial_state_Balanced_AC.mat'],'Initial_state')
    load(['IEEE1354_Capacity_1_Flow_1_5_InitFail_2_SampleNum_',num2str(Test_num),'_Final_state_Balanced_AC.mat'],'Final_state')
    load(['IEEE1354_Capacity_1_Flow_1_5_InitFail_2_SampleNum_',num2str(Test_num),'_Epsilon_opt_Balanced_AC.mat'],'Epsilon_opt')
    
%     Test_num=50000;
%     load(['IEEE3012_Capacity_1_Flow_1_5_InitFail_2_SampleNum_',num2str(Test_num),'_Initial_state_Balanced_New.mat'],'Initial_state')
%     load(['IEEE3012_Capacity_1_Flow_1_5_InitFail_2_SampleNum_',num2str(Test_num),'_Final_state_Balanced_New.mat'],'Final_state')
%     load(['IEEE3012_Capacity_1_Flow_1_5_InitFail_2_SampleNum_',num2str(Test_num),'_Epsilon_opt_Balanced_New.mat'],'Epsilon_opt')

%     Initial_index_lf=cell(size(low_freq_link_index,1),1);
%     for q=1:size(low_freq_link_index,1)
%        Initial_index_lf{q}=zeros(3,size(low_freq_sample_index{q},2));
%        for k=1:size(low_freq_sample_index{q},2)
%            Initial_index_lf{q}(:,k)=(find(Initial_state(:,low_freq_sample_index{q}(k))==0))';
%        end
%     end
%     save(['IEEE118_Capacity_1_2_InitFail_3_SampleNum_',num2str(Test_num),'_Initial_index_lf_Balanced.mat'],'Initial_index_lf')
%     load(['IEEE118_Capacity_1_2_InitFail_3_SampleNum_',num2str(Test_num),'_Initial_index_lf_Balanced.mat'],'Initial_index_lf')

%     Train_num=floor(size(cascade,2));
    Train_num=50000;
    Initial_state_prev=Initial_state(:,1:K_prev);
    Initial_state=Initial_state(:,1:Train_num);
    Final_state=Final_state(:,1:Train_num);
    Epsilon_opt=Epsilon_opt(:,1:Train_num);
    
    M=size(Initial_state,1);
    Link_Index=zeros(M,1);
    for i=1:M
        Link_Index(i,1)=i;
    end
    
    fail_frequency=zeros(M,1);
    for i=1:M
        fail_frequency(i,1)=size(find(Initial_state(i,:)==1 & Final_state(i,:)==0),2)/size(find(Initial_state(i,:)==1),2);
    end
    
%     tmp_1=Epsilon_opt(low_freq_link_index(5),:);
%     tmp_2=cascade_train(low_freq_link_index(5),1:Train_num);
%     tmp_3=Final_state(low_freq_link_index(5),:);
%     tmp=[tmp_1;tmp_2;tmp_3]';
%     tmp_sort=sortrows(tmp,1);
    
%     load('IEEE2383_Capa_1_RowColumnIndex.mat','z');
%     load('IEEE2383_Capa_1_Length_Matrix.mat','length_matrix')
    
%     link_index=z;
%     M=size(link_index,1);

%     P_1=rand(M);
%     P_2=rand(M);
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
    
    %Incorporate more information about the training set
%     load('IEEE1888_Capacity_1_InitFail_2_SampleNum_15000_New_Train_Tmp_2.mat','cascade_train');
%     load('IEEE2383_Capacity_1_InitFail_2_New_Train_Tmp.mat','cascade_train');
%     cascade=cascade_train;

%     K=size(cascade_train,2);
    K=Train_num;
    
%     state_series_tensor=cell(K);
%     state_series_tensor=zeros(M,25,K);
    size_state_vector=zeros(K_prev,1);
    failure_size_vector=zeros(K_prev,1);
    weight_vector=zeros(K_prev,1);
%     failure_exist_flag=zeros(M,1);
    for i=1:K_prev
        max_i=max(cascade(:,i));
        tmp_cascade=cascade(:,i);
        tmp_cascade(cascade(:,i)==max_i)=0;
        submax_i=max(tmp_cascade);
        if max_i-submax_i>=2
            size_state=submax_i+1; 
        else
            size_state=submax_i;
        end
%         state_series=ones(M,size_state);
%         for j=1:size_state
%             state_series(find(cascade(:,i)<=j),j)=0;
%         end
%         state_series_tensor(:,1:size_state,i)=state_series;
%         state_series_tensor{i}=state_series;
        size_state_vector(i,1)=size_state;
        failure_size_vector(i,1)=size(find(cascade(:,i)<size_state_vector(i,1) & cascade(:,i)>1),1);
        weight_vector(i,1)=1;
    end
    
    failure_exist_flag=zeros(M,1);
    for i=1:K_prev %We need to fix the set to estimate the link failure frequency
        tmp_failure=zeros(M,1);
        tmp_failure(cascade(:,i)<size_state_vector(i,1) & cascade(:,i)>1,1)=1;
        failure_exist_flag=failure_exist_flag+tmp_failure;
    end

    fail_frequency=zeros(M,1);
    link_index=zeros(M,1);
    link_fail_vector=zeros(M,1);
    for i=1:M
        link_index(i,1)=i;
        link_fail_vector(i,1)=size(find(Initial_state(i,:)==1 & Final_state(i,:)==0),2);
        fail_frequency(i,1)=size(find(Initial_state(i,:)==1 & Final_state(i,:)==0),2)/size(find(Initial_state(i,:)==1),2);
    end

%     [failure_exist_flag,link_fail_vector]
%     [failure_exist_flag/Train_num,fail_frequency]
    
%     size(find(failure_exist_flag==0),1)
%     M
%     size(find(failure_exist_flag==0),1)/M
    
%     Bd_vec=0:0.1:0.9;
%     for i=1:size(Bd_vec,2)
%         eff_index=find(failure_exist_flag>(0.9-Bd_vec(i))*Train_num & failure_exist_flag<(1-Bd_vec(i))*Train_num);
%         size(eff_index)
%     end

%     load('IEEE1888_Capacity_1_InitFail_2_SampleNum_10000_New_Test_Tmp.mat','cascade_test')
%     load('IEEE2383_Capacity_1_InitFail_2_New_Test_Tmp_2.mat','cascade_test')
%     load('IEEE1354_Capacity_1_InitFail_2_Test_Balanced.mat','cascade_test')
%     load('IEEE118_Capacity_1_5_InitFail_3_Test_Balanced_New_AC.mat','cascade_test')
%     load('IEEE2383_Capacity_1_InitFail_2_Test_Balanced_AC.mat','cascade_test')
%     load('IEEE2383_Capacity_1_Flow_1_25_InitFail_2_Test_Balanced_New_AC.mat','cascade_test')

    load('IEEE1354_Capacity_1_Flow_1_5_InitFail_2_Test_Balanced_New_AC.mat','cascade_link_cum')
    cascade_test=cascade_link_cum;
 %    load('IEEE1354_Capacity_1_Flow_1_5_InitFail_2_Test_Balanced_New_AC.mat','cascade_test')
%    load('IEEE2383_Capacity_1_InitFail_2_Test_Balanced_AC.mat','cascade_test')
%     load('IEEE3012_Capacity_1_Flow_1_5_InitFail_2_Test_Balanced_New.mat','cascade_test')
%     load('IEEE118_Capacity_1_2_InitFail_3_Test_Balanced_New_AC.mat','cascade_test')
    cascade=cascade_test;    

    Test_num=500;
    cascade=cascade(:,end-Test_num+1:end);

%     cascade=cascade(:,1:Test_num);
    
    CESR_record=zeros(Test_num,1);
    CESR_record_small_size=zeros(Test_num,1);
    CESR_each_link_cum=zeros(M,1);
    CESR_each_link_cum_small_size=zeros(M,1);
    
    false_alarm_vector=zeros(Test_num,1); 
    false_alarm_vector_small_size=zeros(Test_num,1);
    false_alarm_each_link_cum=zeros(M,1);
    failure_num_vector_cum=zeros(M,1); %Correct false_alarm
    false_alarm_each_link_cum_small=zeros(M,1);
    failure_num_vector_cum_small=zeros(M,1);
    
    normal_alarm_vector=zeros(Test_num,1);
    
    false_alarm_true_each_link_cum=zeros(M,1);
    normal_num_vector_cum=zeros(M,1);
    
    failure_size_record=zeros(Test_num,1);
    relative_failure_size_record=zeros(Test_num,1);
    failure_size_small_record=zeros(Test_num,1);
    
    real_failure_size_record=zeros(Test_num,1);
    predicted_failure_size_record=zeros(Test_num,1);
    
    rate1_record=zeros(Test_num,1);
    rate2_record=zeros(Test_num,1);
    error_record=zeros(Test_num,1);
    error_record_small_size=zeros(Test_num,1);
    error_each_link_cum=zeros(M,1);
    error_each_link_cum_small_size=zeros(M,1);
    
    real_final_state_record=zeros(M,Test_num);
    real_time_record=zeros(M,Test_num);

    %Pure_Direction
    Pure_CESR_record=zeros(Test_num,1);
    Pure_CESR_each_link_cum=zeros(M,1);
    pure_false_alarm_vector=zeros(Test_num,1);
    pure_rate1_record=zeros(Test_num,1);
    pure_rate2_record=zeros(Test_num,1);
    pure_error_record=zeros(Test_num,1);
    pure_error_each_link_cum=zeros(M,1);
    
    Real_Final_Failure_Size=zeros(Test_num,1);
    Est_Final_Failure_Size=zeros(Test_num,1);
    
    %Debugging...
    error_test=zeros(Test_num,1);
    num_index=zeros(Test_num,1);
    count_useless_vector=zeros(Test_num,1);
        
    count_small=0;
    Index_small=[];
    
    fail_size_record=zeros(Test_num,1);
    
    Link_failure_record=zeros(M,1);
    false_alarm_index=zeros(Test_num,1);
    
    time_index_vec_cum=zeros(M,1);

    sum_D_record=zeros(K,1);
    for k=1:K
        sum_D_record(k,1)=sum(sum(D(:,find(Initial_state(:,k)==0)))); 
    end
    
    real_failure_frequency=zeros(M,1);
    predict_failure_frequency=zeros(M,1);
    
    real_failure_count=zeros(M,Test_num);
    predict_failure_count=zeros(M,Test_num);
    count_effective=0;
    
%     Epsilon_opt=round(Epsilon_opt,6);
    Index_Test_Sample=zeros(Test_num,1);
    Test_link_index=find(fail_frequency>=0 & fail_frequency<=1);

    tic
    for k=1:Test_num
        k
        
        real_vector=cascade(:,k);
        tmp_rnd=find(cascade(:,k)==1);
        real_time_record(:,k)=real_vector;

        real_final_state=zeros(M,1);
        max_real=max(real_vector(:,1));
        tmp_real=real_vector(:,1);
        tmp_real(real_vector(:,1)==max_real)=0;
        submax_real=max(tmp_real);
        if max_real-submax_real>=2
            real_final_state(real_vector(:,1)==max_real,1)=1;
        end
        real_final_state_record(:,k)=real_final_state;
        
        count=0;
        state=ones(M,1);
        state(tmp_rnd)=0;
        state_ini=state;
            
        fail_size_record(k,1)=size(find(real_final_state==0),1);
        false_index=find(real_final_state==0 & state_ini==1);
        normal_index=find(real_final_state==1 & state_ini==1);
        test_false_index=intersect(false_index,Test_link_index);

%         tic
                
        tmp_vec=(state)'*(Initial_state);
        index_set=find(tmp_vec>=(max(tmp_vec)));
        if max(tmp_vec)<M-size(find(state_ini==0),1) && max_real>3 %&& size(test_false_index,1)>0
        Index_Test_Sample(k,1)=1;
        Epsilon_Test_1=median(Epsilon_opt(:,index_set),2);
        
%         tmp_2_vec=(state)'*(Initial_state(:,index_set));
%         index_set_2=find(tmp_2_vec==max(tmp_2_vec));

%         Epsilon_Test_1=median(Epsilon_opt(:,index_set(index_set_2)),2);
        
%         Epsilon_Test=zeros(M,1);
%         count_useless=0;
%         for i=1:M
%             Index_Initial_state=Initial_state(i,index_set);
%             if size(find(Index_Initial_state==1),2)>0
%                 Epsilon_Test(i,1)=median(Epsilon_opt(i,index_set(Index_Initial_state==1)));
%             else
%                 count_useless=count_useless+1;
%                 Epsilon_Test(i,1)=0.5;  %0.5 means no information provided to estimate the threshold of the link i.
%             end
%         end
%         
%         Epsilon_Test(eff_index,1)=1;
%         Epsilon_Test(low_freq_link_index,1)=1;

%        Epsilon_Test_2=zeros(M,1);
%        A_diff_0_vec=sum(M_est(tmp_rnd,:),1);
%        A_diff_k=(~Initial_state(:,index_set))'*M_est;
%        A_diff_0=repmat(A_diff_0_vec,[size(index_set,2),1]);
%        A_diff=abs(A_diff_k-A_diff_0);
%        for i=1:M
% % %            if ismember(i,low_freq_link_index)
%               Epsilon_Test_2(i,1)=median(Epsilon_opt(i,index_set(find(A_diff(:,i)==min(A_diff(:,i))))));
% %            end
%        end

%        Epsilon_Test_2=zeros(M,1);
%        A_diff_0_vec=sum(A_1(tmp_rnd,:)-A_2(tmp_rnd,:),1);
%        A_diff_k=(~Initial_state(:,index_set))'*(A_1-A_2);
%        A_diff_0=repmat(A_diff_0_vec,[size(index_set,2),1]);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
%        A_diff=abs(A_diff_k-A_diff_0);
%        for i=1:M
% % %            if ismember(i,low_freq_link_index)
%               Epsilon_Test_2(i,1)=median(Epsilon_opt(i,index_set(find(A_diff(:,i)==min(A_diff(:,i))))));
% %            end
%        end
       
%        Epsilon_Test=(fail_frequency>0.2 & fail_frequency<0.8).*Epsilon_Test_1+(fail_frequency<=0.2 | fail_frequency>=0.8).*Epsilon_Test_2;
        Epsilon_Test = Epsilon_Test_1;

        count=0;
        state=ones(M,1);
        state(tmp_rnd)=0;
        state_ini=state;
        
        state_cum=[];
        state_1_cum=[];
        state_cum=[state_cum,state];
        state_1_cum=[state_1_cum,state];
        state_int=state;
        state_tmp=state;
        Epsilon=zeros(M,1);
        
%         for i=1:M
%             Epsilon(i,1)=1-exp(Beta_record(:,i)'*[1,state_int']')/(1+exp(Beta_record(:,i)'*[1,state_int']'));
%         end

        state_cum_pure=[];
        state_cum_pure=[state_cum_pure,state];
        state_int_pure=state;
        Epsilon_pure=zeros(M,1);
        
%         for i=1:M
%             Epsilon_pure(i,1)=1-exp(Beta_record(:,i)'*[1,state_int_pure']')/(1+exp(Beta_record(:,i)'*[1,state_int_pure']'));
%         end
        
        T=15;
        state_tmp_prev=state;
        
        while count<T
            count=count+1;
            
            state_1 = M_est*state+b_est;
            state_tmp = state_1;
            
            state_tmp = state_tmp>Epsilon_Test;
            state_tmp(state_tmp_prev==0,1)=0;     
            %NOTE: If there is randomness, then we can not break here......
            if sum(abs(state_tmp-state_tmp_prev))==0
                break; 
            end
            
            state_tmp_prev=state_tmp;
            
            state=state_1;
            state_cum=[state_cum,state_tmp];
            state_1_cum=[state_1_cum,state_1];
            
            state_tmp_pure=state_int_pure; 

            state_tmp_pure = state_tmp_pure<0.5;
            state_int_pure=state_tmp_pure;
            state_cum_pure=[state_cum_pure,state_tmp_pure];    
        end       

        predicted_vector=zeros(M,1);
        for i=1:M
            tmp_predicted=find(state_cum(i,:)==0);
            if size(tmp_predicted,2)~=0
                predicted_vector(i,1)=tmp_predicted(1);
            else
                predicted_vector(i,1)=0; 
            end
        end
        max_round=max(predicted_vector)+2;
        for i=1:M
            if predicted_vector(i,1)==0
               predicted_vector(i,1)=max_round; 
            end
        end
                
        predicted_vector_pure=zeros(M,1);
        for i=1:M
            tmp_predicted_pure=find(state_cum_pure(i,:)==0);
            if size(tmp_predicted_pure,2)~=0
                predicted_vector_pure(i,1)=tmp_predicted_pure(1);
            else
                predicted_vector_pure(i,1)=0; 
            end
        end
        max_round_pure=max(predicted_vector_pure)+2;
        for i=1:M
            if predicted_vector_pure(i,1)==0
               predicted_vector_pure(i,1)=max_round_pure; 
            end
        end   
        [state_1_cum,real_vector];
        
        %%
        %Different Criteria
        [Link_Index,real_vector,predicted_vector];
        sum(abs(real_vector-predicted_vector))/M;
        
        %% 1.Correct End State Ratio
        predicted_final_state=zeros(M,1);
        max_predicted=max(predicted_vector(:,1));
        tmp_predicted=predicted_vector(:,1);
        tmp_predicted(predicted_vector(:,1)==max_predicted)=0;
        submax_predicted=max(tmp_predicted);
        if max_predicted-submax_predicted>=2
            predicted_final_state(predicted_vector(:,1)==max_predicted,1)=1;
        end   
        
        %Add a stupid judgement...
%         Final_candidate=Final_state(:,index_set);
%         for i=1:M
%             if Final_fail_rate(i,1)>0.35 && Final_fail_rate(i,1)<0.65
% %                 median(Final_candidate(i,:))
%                 if median(Final_candidate(i,:))~=0.5
%                     predicted_final_state(i,-4)=median(Final_candidate(i,:));
%                 else
%                     predicted_final_state(i,1)=round(rand(1));
%                 end
%             end
%         end
        
        pure_predicted_final_state=zeros(M,1);
        max_predicted_pure=max(predicted_vector_pure(:,1));
        tmp_predicted_pure=predicted_vector_pure(:,1);
        tmp_predicted_pure(predicted_vector_pure(:,1)==max_predicted_pure)=0;
        submax_predicted_pure=max(tmp_predicted_pure);
        if max_predicted_pure-submax_predicted_pure>=2
            pure_predicted_final_state(predicted_vector_pure(:,1)==max_predicted_pure,1)=1;
        end
        
%         pure_predicted_final_state=zeros(M,1);
%         for i=1:M
%             if Epsilon(i,1)<1.5
%                 pure_predicted_final_state(i,1)=1;
%             else
%                 pure_predicted_final_state(i,1)=0;
%             end
%         end
        
        real_failure_count(find(real_final_state==0 & state_ini==1),k)=1;
        predict_failure_count(find(predicted_final_state==0 & state_ini==1),k)=1;
                
        Fail_size=size(find(predicted_final_state==1 & real_final_state==0),1)+size(find(predicted_final_state==0 & real_final_state==1),1);
        [Link_Index,real_final_state,predicted_final_state]';
        Fail_size;
                
        %False_alarm_rate
%         [size(find(predicted_final_state(false_index,1)==1),1),size(false_index,1)];
        false_alarm_vector(k,1)=size(find(predicted_final_state(false_index)==1),1)/size(false_index,1);
        normal_alarm_vector(k,1)=size(find(predicted_final_state(normal_index)==0),1)/size(normal_index,1);
        
        diff=real_final_state-predicted_final_state;
        diff(diff>=0,1)=0;
        diff(diff<0,1)=1;
        false_alarm_each_link_cum=false_alarm_each_link_cum+diff;
        
        %false_alarm_each_link_cum is miss detection
        
        diff=real_final_state-predicted_final_state;
        diff(diff>0,1)=1;
        diff(diff<=0,1)=0;
        false_alarm_true_each_link_cum=false_alarm_true_each_link_cum+diff;
        
        %false_alarm_true_each_link_cum is false alarm
        
        failure_num_vec=zeros(M,1);
        failure_num_vec(real_final_state==0 & state_ini==1)=1;
        failure_num_vector_cum=failure_num_vector_cum+failure_num_vec;
        
        failure_num_vec=zeros(M,1);
        failure_num_vec(real_final_state==1 & state_ini==1)=1;
        normal_num_vector_cum=normal_num_vector_cum+failure_num_vec;       
                
        pure_false_alarm_vector(k,1)=size(find(pure_predicted_final_state(false_index,1)==1),1)/size(false_index,1);
        
        CESR_record(k,1)=sum(abs(real_final_state(Test_link_index)-predicted_final_state(Test_link_index)))/size(Test_link_index,1);
        CESR_each_link_cum=CESR_each_link_cum+abs(real_final_state-predicted_final_state);
        Pure_CESR_record(k,1)=sum(abs(real_final_state-pure_predicted_final_state))/M;
        Pure_CESR_each_link_cum=Pure_CESR_each_link_cum+abs(real_final_state-pure_predicted_final_state);   
        
        time_index_vec=zeros(M,1);
        time_index_vec(false_index)=1;
        
        Time_distance=abs(real_vector-predicted_vector).*time_index_vec;
        Time_distance_pure=abs(real_vector-predicted_vector_pure);

        [real_vector, predicted_vector];
        
        AvgR1=zeros(M,1);
        AvgR1(Time_distance(:,1)<=1,1)=1;
        AvgR2=zeros(M,1);
        AvgR2(Time_distance(:,1)<=2,1)=1;
        
        rate1=sum(AvgR1)/M;
        rate1_record(k,1)=rate1;
        rate2=sum(AvgR2)/M;
        rate2_record(k,1)=rate2;
        
        pure_AvgR1=zeros(M,1);
        pure_AvgR1(Time_distance_pure(:,1)<=1,1)=1;
        pure_AvgR2=zeros(M,1);
        pure_AvgR2(Time_distance_pure(:,1)<=2,1)=1;
        
        rate1_pure=sum(pure_AvgR1)/M;
        pure_rate1_record(k,1)=rate1_pure;
        rate2_pure=sum(pure_AvgR2)/M;
        pure_rate2_record(k,1)=rate2_pure;
        
        %%Error only for failures
        %CESR_record
        error_record(k,1)=sum(Time_distance(test_false_index))/size(test_false_index,1);
        error_each_link_cum=error_each_link_cum+Time_distance.*time_index_vec;
        time_index_vec_cum=time_index_vec_cum+time_index_vec;
        
        failure_size_record(k,1)=abs(size(find(real_final_state(Test_link_index)==0),1)-size(find(predicted_final_state(Test_link_index)==0),1)); %/size(Test_link_index,1);
        real_failure_size_record(k,1)=size(find(real_final_state(Test_link_index)==0),1);
        predicted_failure_size_record(k,1)=size(find(predicted_final_state(Test_link_index)==0),1);
        relative_failure_size_record(k,1)=abs(size(find(real_final_state==0),1)-size(find(predicted_final_state==0),1))/(size(find(real_final_state==0),1)-size(tmp_rnd,1));
        
        pure_error_record(k,1)=sum(Time_distance_pure)/M;
        pure_error_each_link_cum=pure_error_each_link_cum+Time_distance_pure;
                
        if fail_size_record(k,1)<400
%             failure_index=find(real_final_state==0);
            Index_small=[Index_small;k];
            count_small=count_small+1;
            CESR_record_small_size(k,1)=CESR_record(k,1);
            error_record_small_size(k,1)=error_record(k,1);
            false_alarm_vector_small_size(k,1)=false_alarm_vector(k,1);
            false_alarm_index(k,1)=1;
            CESR_each_link_cum_small_size=CESR_each_link_cum_small_size+abs(real_final_state-predicted_final_state);
            error_each_link_cum_small_size=error_each_link_cum_small_size+Time_distance;
           
            diff=real_final_state-predicted_final_state;
            diff(diff>=0,1)=0;
            diff(diff<0,1)=1;
            false_alarm_each_link_cum_small=false_alarm_each_link_cum_small+diff;
            
            failure_num_vec=zeros(M,1);
            failure_num_vec(real_final_state==0 & state_ini==1)=1;
            failure_num_vector_cum_small=failure_num_vector_cum_small+failure_num_vec;      
            
            failure_size_small_record(k,1)=failure_size_record(k,1);
        end
        
        Link_failure_record=Link_failure_record+real_final_state+(1-state_ini);
        count_effective=count_effective+1;
        end
    end
    toc
    
    Test_num=count_effective
    
    %Calculate Failure Frequency
    real_failure_frequency=sum(real_failure_count,2)/Test_num;
    predicted_failure_frequency=sum(predict_failure_count,2)/Test_num;
    failure_frequency_error_record_all=abs(real_failure_frequency-predicted_failure_frequency);
    failure_frequency_error_record=abs(real_failure_frequency(Test_link_index)-predicted_failure_frequency(Test_link_index));
    failure_frequency_avg=mean(failure_frequency_error_record);
    failure_frequency_var=var(failure_frequency_error_record);
    failure_frequency_med=median(failure_frequency_error_record);
    failure_frequency_error_record_sorted=sort(failure_frequency_error_record,'ascend');
    failure_frequency_min=failure_frequency_error_record_sorted(3:-1:1)';
    failure_frequency_max=failure_frequency_error_record_sorted(end:-1:end-2)';
    
    %Calculate Relative Failure Frequency
    relative_failure_frequency_error_record=abs(real_failure_frequency(Test_link_index)-predicted_failure_frequency(Test_link_index))./real_failure_frequency(Test_link_index);
    relative_failure_frequency_avg=mean(relative_failure_frequency_error_record);
    relative_failure_frequency_var=var(relative_failure_frequency_error_record);
    relative_failure_frequency_med=median(relative_failure_frequency_error_record);
    relative_failure_frequency_error_record_sorted=sort(relative_failure_frequency_error_record,'ascend');
    relative_failure_frequency_min=relative_failure_frequency_error_record_sorted(3:-1:1)';
    relative_failure_frequency_max=relative_failure_frequency_error_record_sorted(end:-1:end-2)';   
    
%     sum(false_alarm_each_link_cum)
    
%     z=[error_test,CESR_record,error_record];
%     sortrows(z,-1);
%     Index=zeros(Test_num,1);
%     for i=1:Test_num
%         Index(i,1)=i;
%     end
% 
%     z=[Index,CESR_record,false_alarm_vector,error_record];
%     sortrows(z,-2)
    
    Cov_Real_Est=cov(error_test,CESR_record);
    rho=Cov_Real_Est(1,2)/sqrt(Cov_Real_Est(1,1)*Cov_Real_Est(2,2)); 
    
    [Real_Final_Failure_Size,Est_Final_Failure_Size];
    
    Real_rate=Real_Final_Failure_Size/M;
    Est_rate=Est_Final_Failure_Size/M;
%     Cov_Real_Est=cov(Real_rate,Est_rate);
%     rho=Cov_Real_Est(1,2)/sqrt(Cov_Real_Est(1,1)*Cov_Real_Est(2,2));
    %To add some statistics: expectation & variance
%     save('IEEE39_Capacity_0_7_real_final_state_record_1000.mat','real_final_state_record');
%     save('IEEE39_Capacity_0_7_real_time_record_1000.mat','real_time_record');
    
    rate1_record_avg=sum(rate1_record)/Test_num;
    rate1_record_var=var(rate1_record);
    rate1_record_med=median(rate1_record);
    rate1_record_sorted=sort(rate1_record,'ascend');
    rate1_record_min=rate1_record_sorted(3:-1:1)';
    rate1_record_max=rate1_record_sorted(end:-1:end-2)';
    
    rate2_record_avg=sum(rate2_record)/Test_num;
    rate2_record_var=var(rate2_record);
    rate2_record_med=median(rate2_record);
    rate2_record_sorted=sort(rate2_record,'ascend');
    rate2_record_min=rate2_record_sorted(3:-1:1)';
    rate2_record_max=rate2_record_sorted(end:-1:end-2)';  
    
    error_record_avg=sum(error_record)/Test_num;
    error_record_var=var(error_record);
    error_record_med=median(error_record);
    error_record_sorted=sort(error_record,'ascend');
    error_record_min=error_record_sorted(3:-1:1)';
    error_record_max=error_record_sorted(end:-1:end-2)';        
 
%     error_record_small_size=error_record_small_size(Index_small,1);
%     error_record_small_size_avg=sum(error_record_small_size)/count_small;
%     error_record_small_size_var=var(error_record_small_size);
%     error_record_small_size_med=median(error_record_small_size);
%     error_record_small_size_sorted=sort(error_record_small_size,'ascend');
%     error_record_small_size_min=error_record_small_size_sorted(3:-1:1)';
%     error_record_small_size_max=error_record_small_size_sorted(end:-1:end-2)';   
    
    error_each_link_record=zeros(M,1);
    for i=1:size(time_index_vec_cum,1)
        if time_index_vec_cum(i,1)>0
            error_each_link_record(i,1)=error_each_link_cum(i,1)/time_index_vec_cum(i,1);
        end
    end
    error_each_link_record_avg=sum(error_each_link_record)/M;
    error_each_link_record_var=var(error_each_link_record);
    error_each_link_record_med=median(error_each_link_record);
    error_each_link_record_sorted=sort(error_each_link_record,'ascend');
    error_each_link_record_min=error_each_link_record_sorted(3:-1:1)';
    error_each_link_record_max=error_each_link_record_sorted(end:-1:end-2)';   
    
    error_each_link_small_size_record=error_each_link_cum_small_size/count_small;
    error_each_link_small_size_record_avg=sum(error_each_link_small_size_record)/M;
    error_each_link_small_size_record_var=var(error_each_link_small_size_record);
    error_each_link_small_size_record_med=median(error_each_link_small_size_record);
    error_each_link_small_size_record_sorted=sort(error_each_link_small_size_record,'ascend');
    error_each_link_small_size_record_min=error_each_link_small_size_record_sorted(3:-1:1)';
    error_each_link_small_size_record_max=error_each_link_small_size_record_sorted(end:-1:end-2)';       
    
%     CESR_record
    CESR_record_avg=sum(CESR_record)/Test_num;
    CESR_record_var=var(CESR_record);
    CESR_record_med=median(CESR_record);
    CESR_record_sorted=sort(CESR_record,'ascend');
    CESR_record_min=CESR_record_sorted(3:-1:1)';
    CESR_record_max=CESR_record_sorted(end:-1:end-2)';   
    
%     CESR_record_small_size=CESR_record_small_size(Index_small,1);
%     CESR_record_small_size_avg=sum(CESR_record_small_size)/count_small;
%     CESR_record_small_size_var=var(CESR_record_small_size);
%     CESR_record_small_size_med=median(CESR_record_small_size);
%     CESR_record_small_size_sorted=sort(CESR_record_small_size,'ascend');
%     CESR_record_small_size_min=CESR_record_small_size_sorted(3:-1:1)';
%     CESR_record_small_size_max=CESR_record_small_size_sorted(end:-1:end-2)';
         
    CESR_each_link_record=CESR_each_link_cum/Test_num;
    CESR_each_link_record_avg=sum(CESR_each_link_record)/M;
    CESR_each_link_record_var=var(CESR_each_link_record);
    CESR_each_link_record_med=median(CESR_each_link_record);
    CESR_each_link_record_sorted=sort(CESR_each_link_record,'ascend');
    CESR_each_link_record_min=CESR_each_link_record_sorted(3:-1:1)';
    CESR_each_link_record_max=CESR_each_link_record_sorted(end:-1:end-2)';
    
%     CESR_each_link_small_size_record=CESR_each_link_cum_small_size/count_small;
%     CESR_each_link_small_size_record_avg=sum(CESR_each_link_small_size_record)/M;
%     CESR_each_link_small_size_record_var=var(CESR_each_link_small_size_record);
%     CESR_each_link_small_size_record_med=median(CESR_each_link_small_size_record);
%     CESR_each_link_small_size_record_sorted=sort(CESR_each_link_small_size_record,'ascend');
%     CESR_each_link_small_size_record_min=CESR_each_link_small_size_record_sorted(3:-1:1)';
%     CESR_each_link_small_size_record_max=CESR_each_link_small_size_record_sorted(end:-1:end-2)';   

%     y=[Link_Index,CESR_each_link_record,error_each_link_record];
%     sortrows(y,-3)
%     plot(CESR_each_link_record,error_each_link_record,'*')
    
%     Pure_CESR_record_avg=sum(Pure_CESR_record)/Test_num;
%     Pure_CESR_record_var=var(Pure_CESR_record);
%     Pure_CESR_record_med=median(Pure_CESR_record);
%     Pure_CESR_record_sorted=sort(Pure_CESR_record,'ascend');
%     Pure_CESR_record_min=Pure_CESR_record_sorted(3:-1:1)';
%     Pure_CESR_record_max=Pure_CESR_record_sorted(end:-1:end-2)';
% 
%     Pure_CESR_each_link_record=Pure_CESR_each_link_cum/Test_num;
%     Pure_CESR_each_link_record_avg=sum(Pure_CESR_each_link_record)/M;
%     Pure_CESR_each_link_record_var=var(Pure_CESR_each_link_record);
%     Pure_CESR_each_link_record_med=median(Pure_CESR_each_link_record);
%     Pure_CESR_each_link_record_sorted=sort(Pure_CESR_each_link_record,'ascend');
%     Pure_CESR_each_link_record_min=Pure_CESR_each_link_record_sorted(3:-1:1)';
%     Pure_CESR_each_link_record_max=Pure_CESR_each_link_record_sorted(end:-1:end-2)';   
%     
%     Pure_Binary_error=[Pure_CESR_record_avg,Pure_CESR_record_var,Pure_CESR_record_med,Pure_CESR_record_max,Pure_CESR_record_min]
%     Pure_Binary_each_link_error=[Pure_CESR_each_link_record_avg,Pure_CESR_each_link_record_var,Pure_CESR_each_link_record_med,Pure_CESR_each_link_record_max,Pure_CESR_each_link_record_min]

    false_alarm_avg=sum(false_alarm_vector)/Test_num;
    false_alarm_var=var(false_alarm_vector);
    false_alarm_med=median(false_alarm_vector);
    false_alarm_sorted=sort(false_alarm_vector,'ascend');
    false_alarm_min=false_alarm_sorted(3:-1:1)';
    false_alarm_max=false_alarm_sorted(end:-1:end-2)';
    
    normal_alarm_avg=sum(normal_alarm_vector)/Test_num;
    normal_alarm_var=var(normal_alarm_vector);
    normal_alarm_med=median(normal_alarm_vector);
    normal_alarm_sorted=sort(normal_alarm_vector,'ascend');
    normal_alarm_min=normal_alarm_sorted(3:-1:1)';
    normal_alarm_max=normal_alarm_sorted(end:-1:end-2)';
    
%     false_alarm_vector_small_size=false_alarm_vector_small_size(Index_small,1);
%     false_alarm_small_size_avg=sum(false_alarm_vector_small_size)/count_small;
%     false_alarm_small_size_var=var(false_alarm_vector_small_size);
%     false_alarm_small_size_med=median(false_alarm_vector_small_size);
%     false_alarm_small_size_sorted=sort(false_alarm_vector_small_size,'ascend');
%     false_alarm_small_size_min=false_alarm_small_size_sorted(3:-1:1)';
%     false_alarm_small_size_max=false_alarm_small_size_sorted(end:-1:end-2)';
    
    failure_size_avg=sum(failure_size_record)/Test_num;
    failure_size_var=var(failure_size_record);
    failure_size_med=median(failure_size_record);
    failure_size_sorted=sort(failure_size_record,'ascend');
    failure_size_min=failure_size_sorted(3:-1:1)';
    failure_size_max=failure_size_sorted(end:-1:end-2)';  
    
%     save(['IEEE2383_Capacity_1_Flow_1_25_InitFail_2_Real_Failure_Size_Balanced_New.mat'],'real_failure_size_record')
%     save(['IEEE2383_Capacity_1_Flow_1_25_InitFail_2_Predicted_Failure_Size_Balanced_New.mat'],'predicted_failure_size_record')
    save(['IEEE1354_Capacity_1_Flow_1_5_InitFail_2_Real_Failure_Size_Balanced_New_AC.mat'],'real_failure_size_record')
    save(['IEEE1354_Capacity_1_Flow_1_5_InitFail_2_Predicted_Failure_Size_Balanced_New_AC.mat'],'predicted_failure_size_record')
    
    relative_failure_size_avg=sum(relative_failure_size_record)/Test_num;
    relative_failure_size_var=var(relative_failure_size_record);
    relative_failure_size_med=median(relative_failure_size_record);
    relative_failure_size_sorted=sort(relative_failure_size_record,'ascend');
    relative_failure_size_min=relative_failure_size_sorted(3:-1:1)';
    relative_failure_size_max=relative_failure_size_sorted(end:-1:end-2)';  
    
%     failure_size_small_record=failure_size_small_record(Index_small,1);
%     failure_size_small_record_avg=sum(failure_size_small_record)/count_small;
%     failure_size_small_record_var=var(failure_size_small_record);
%     failure_size_small_record_med=median(failure_size_small_record);
%     failure_size_small_record_sorted=sort(failure_size_small_record,'ascend');
%     failure_size_small_record_min=failure_size_small_record_sorted(3:-1:1)';
%     failure_size_small_record_max=failure_size_small_record_sorted(end:-1:end-2)';
    
%     false_alarm_each_link_record=zeros(M,1);
%     count_false_alarm=0;
%     for i=1:M
%         if failure_num_vector_cum(i,1)>0
%             false_alarm_each_link_record(i,1)=false_alarm_each_link_cum(i,1)/failure_num_vector_cum(i,1);
% %             false_alarm_each_link_cum(i,1);
%             count_false_alarm=count_false_alarm+1;
%         end
%     end

%     false_alarm_each_link_record
%     save('IEEE39_Capacity_0_7_CESR_each_link_20000.mat','CESR_each_link_record')
%     false_alarm_each_link_record_avg=sum(false_alarm_each_link_record)/M;
%     false_alarm_each_link_record_var=var(false_alarm_each_link_record);
%     false_alarm_each_link_record_med=median(false_alarm_each_link_record);
%     false_alarm_each_link_record_sorted=sort(false_alarm_each_link_record,'ascend');
%     false_alarm_each_link_record_min=false_alarm_each_link_record_sorted(3:-1:1)';
%     false_alarm_each_link_record_max=false_alarm_each_link_record_sorted(end:-1:end-2)';   
    
    Binary_error=[CESR_record_avg,CESR_record_var,CESR_record_med,CESR_record_max,CESR_record_min];
%     Binary_each_link_error=[CESR_each_link_record_avg,CESR_each_link_record_var,CESR_each_link_record_med,CESR_each_link_record_max,CESR_each_link_record_min]
    false_alarm_binary_error=[false_alarm_avg,false_alarm_var,false_alarm_med,false_alarm_max,false_alarm_min];
    normal_alarm_binary_error=[normal_alarm_avg,normal_alarm_var,normal_alarm_med,normal_alarm_max,normal_alarm_min]; 
    Time_error=[error_record_avg,error_record_var,error_record_med,error_record_max,error_record_min];
%     Time_each_link_error=[error_each_link_record_avg,error_each_link_record_var,error_each_link_record_med,error_each_link_record_max,error_each_link_record_min]
    Rate1=[rate1_record_avg,rate1_record_var,rate1_record_med,rate1_record_max,rate1_record_min];
    Rate2=[rate2_record_avg,rate2_record_var,rate2_record_med,rate2_record_max,rate2_record_min];
    Failure_size_error=[failure_size_avg,failure_size_var,failure_size_med,failure_size_max,failure_size_min];
    Failure_frequency_error=[failure_frequency_avg,failure_frequency_var,failure_frequency_med,failure_frequency_max,failure_frequency_min];
%     Failure_frequency_error;
    Relative_Failure_size_error=[relative_failure_size_avg,relative_failure_size_var,relative_failure_size_med,relative_failure_size_max,relative_failure_size_min];
    Relative_Failure_frequency_error=[relative_failure_frequency_avg,relative_failure_frequency_var,relative_failure_frequency_med,relative_failure_frequency_max,relative_failure_frequency_min];
    Macro_Result_Record=[Binary_error;false_alarm_binary_error;normal_alarm_binary_error;Time_error;Rate1;Rate2;Failure_size_error;Failure_frequency_error;Relative_Failure_size_error;Relative_Failure_frequency_error]
    
%     LB_vec=[0,0.1,0.2,0.3,0.4];
%     UB_vec=[1,0.9,0.8,0.7,0.6];    
        %%Effective failure testing
%     rate1_record_effective_avg=sum(rate1_record(find(failure_exist_flag>0),1))/size(find(failure_exist_flag>0),1)
%     rate2_record_effective_avg=sum(rate2_record(find(failure_exist_flag>0),1))/size(find(failure_exist_flag>0),1)
%     error_record_effective_avg=sum(error_each_link_record(find(failure_exist_flag>0),1))/size(find(failure_exist_flag>0),1)
%     CESR_record_effective_avg=sum(CESR_each_link_record(find(failure_exist_flag>0),1))/size(find(failure_exist_flag>0),1)
%     false_alarm_effective_avg=sum(false_alarm_each_link_record(find(failure_exist_flag>0 & failure_num_vector_cum>0),1))/count_false_alarm %size(find(failure_exist_flag>0),1)
%     Train_num=5000;

    Initial_normal_vector=zeros(M,1);
    for i=1:M
        Initial_normal_vector(i,1)=size(find(Initial_state_prev(i,:)==1),2);  %Here we need to fix it to be the same.
    end

%     Train_num=5000;
%     for i=1:size(LB_vec,2)
% %         eff_index=find(failure_exist_flag>LB_vec(i)*Initial_normal_vector & failure_exist_flag<UB_vec(i)*Initial_normal_vector);
%         eff_index=find((failure_exist_flag>LB_vec(i)*Initial_normal_vector & failure_exist_flag<(LB_vec(i)+0.1)*Initial_normal_vector)|(failure_exist_flag>(UB_vec(i)-0.1)*Initial_normal_vector & failure_exist_flag<(UB_vec(i))*Initial_normal_vector));
%         error_record_effective_avg=sum(error_each_link_record(eff_index,1))/size(eff_index,1);
%         CESR_each_link_record(eff_index,1)';
%         CESR_record_effective_avg=sum(CESR_each_link_record(eff_index,1))/size(eff_index,1);
%         false_alarm_each_link_cum_sum=sum(false_alarm_each_link_cum(eff_index));
%         failure_num_vector_cum_sum=sum(failure_num_vector_cum(eff_index));
%         false_alarm_effective_avg=false_alarm_each_link_cum_sum/failure_num_vector_cum_sum;
%         [CESR_record_effective_avg,false_alarm_effective_avg,error_record_effective_avg]
%     end
%     
%     Train_num=5000;

    Bd_vec=0:0.1:0.9;
    Result_record=zeros(10,6);
    for i=1:size(Bd_vec,2)
        eff_index=find(failure_exist_flag>(0.9-Bd_vec(i))*Initial_normal_vector & failure_exist_flag<(1-Bd_vec(i))*Initial_normal_vector);
        error_record_effective_avg=sum(error_each_link_record(eff_index,1))/size(eff_index,1);
%         eff_index
        CESR_record_effective_avg=sum(CESR_each_link_record(eff_index,1))/size(eff_index,1);
        false_alarm_each_link_cum_sum=sum(false_alarm_each_link_cum(eff_index));
        failure_num_vector_cum_sum=sum(failure_num_vector_cum(eff_index));
        false_alarm_effective_avg=false_alarm_each_link_cum_sum/failure_num_vector_cum_sum;
        false_alarm_true_each_link_cum_sum=sum(false_alarm_true_each_link_cum(eff_index));
        normal_num_vector_cum_sum=sum(normal_num_vector_cum(eff_index));
        false_alarm_true_effective_avg=false_alarm_true_each_link_cum_sum/normal_num_vector_cum_sum;
        failure_frequency_each_link=mean(failure_frequency_error_record_all(intersect(eff_index,Test_link_index),1));
%         [false_alarm_each_link_cum_sum,failure_num_vector_cum_sum]
        Result_record(i,:)=[Bd_vec(i),CESR_record_effective_avg,false_alarm_effective_avg,false_alarm_true_effective_avg,error_record_effective_avg,failure_frequency_each_link];
%         size(eff_index)
        1;
    end
    Micro_Result_Record=Result_record
    
%     save(['IEEE3012_Capacity_1_Flow_1_5_InitFail_2_SampleNum_',num2str(Train_num),'_Macro_Result_Balanced.mat'],'Macro_Result_Record')
%     save(['IEEE3012_Capacity_1_Flow_1_5_InitFail_2_SampleNum_',num2str(Train_num),'_Micro_Result_Balanced.mat'],'Micro_Result_Record')
    
%     save(['IEEE2383_Capacity_1_Flow_1_25_InitFail_2_SampleNum_',num2str(Train_num),'_Macro_Result_Balanced_AC.mat'],'Macro_Result_Record')
%     save(['IEEE2383_Capacity_1_Flow_1_25_InitFail_2_SampleNum_',num2str(Train_num),'_Micro_Result_Balanced_AC.mat'],'Micro_Result_Record')
    
%     save(['IEEE1354_Capacity_1_Flow_1_5_InitFail_2_SampleNum_',num2str(Train_num),'_Macro_Result_Balanced_AC.mat'],'Macro_Result_Record')
%     save(['IEEE1354_Capacity_1_Flow_1_5_InitFail_2_SampleNum_',num2str(Train_num),'_Micro_Result_Balanced_AC.mat'],'Micro_Result_Record')
    
%     save(['IEEE118_Capacity_1_5_InitFail_3_SampleNum_',num2str(Train_num),'_Macro_Result_Balanced_AC.mat'],'Macro_Result_Record')
%     save(['IEEE118_Capacity_1_5_InitFail_3_SampleNum_',num2str(Train_num),'_Micro_Result_Balanced_AC.mat'],'Micro_Result_Record')

%     save(['IEEE118_Capacity_1_8_InitFail_2_SampleNum_',num2str(Train_num),'_Macro_Result_Balanced_Refined.mat'],'Macro_Result_Record')
%     save(['IEEE118_Capacity_1_8_InitFail_2_SampleNum_',num2str(Train_num),'_Micro_Result_Balanced_Refined.mat'],'Micro_Result_Record')

%     Distribution_Result_Record=[Index_Test_Sample,CESR_record,false_alarm_vector,normal_alarm_vector,error_record,failure_size_record,relative_failure_size_record];
%     save(['IEEE118_Capacity_1_5_InitFail_3_SampleNum_',num2str(Train_num),'_Distribution_Result_Balanced_Refined_AC.mat'],'Distribution_Result_Record')
%     Frequency_Distribution_Record=failure_frequency_error_record;
%     save(['IEEE118_Capacity_1_8_InitFail_3_SampleNum_',num2str(Train_num),'_Frequency_Distribution_Record_Refined.mat'],'Frequency_Distribution_Record')    
%     Relative_Frequency_Distribution_Record=relative_failure_frequency_error_record;
%     save(['IEEE118_Capacity_1_8_InitFail_3_SampleNum_',num2str(Train_num),'_Relative_Frequency_Distribution_Record_Refined.mat'],'Relative_Frequency_Distribution_Record') 

%     Distribution_Result_Record=[Index_Test_Sample,CESR_record,false_alarm_vector,normal_alarm_vector,error_record,failure_size_record,relative_failure_size_record];
%     save(['IEEE118_Capacity_1_8_InitFail_2_SampleNum_',num2str(Train_num),'_Distribution_Result_Balanced_Refined.mat'],'Distribution_Result_Record')
%     Frequency_Distribution_Record=failure_frequency_error_record;
%     save(['IEEE118_Capacity_1_8_InitFail_2_SampleNum_',num2str(Train_num),'_Frequency_Distribution_Record_Refined.mat'],'Frequency_Distribution_Record')    
%     Relative_Frequency_Distribution_Record=relative_failure_frequency_error_record;
%     save(['IEEE118_Capacity_1_8_InitFail_2_SampleNum_',num2str(Train_num),'_Relative_Frequency_Distribution_Record_Refined.mat'],'Relative_Frequency_Distribution_Record')
    
%     Distribution_Result_Record=[Index_Test_Sample,CESR_record,false_alarm_vector,normal_alarm_vector,error_record,failure_size_record,relative_failure_size_record];
%     save(['IEEE1354_Capacity_1_Flow_1_5_InitFail_2_SampleNum_',num2str(Train_num),'_Distribution_Result_Balanced_AC.mat'],'Distribution_Result_Record')
% %     Frequency_Distribution_Record=failure_frequency_error_record;
% %     save(['IEEE1354_Capacity_1_Flow_1_InitFail_3_SampleNum_',num2str(Train_num),'_Frequency_Distribution_Record_Refined_AC.mat'],'Frequency_Distribution_Record')    
% %      Relative_Frequency_Distribution_Record=relative_failure_frequency_error_record;
% %     save(['IEEE1354_Capacity_1_Flow_1_5_InitFail_2_SampleNum_',num2str(Train_num),'_Relative_Frequency_Distribution_Record_Refined.mat'],'Relative_Frequency_Distribution_Record')       

%     Distribution_Result_Record=[Index_Test_Sample,CESR_record,false_alarm_vector,normal_alarm_vector,error_record,failure_size_record,relative_failure_size_record];
%     save(['IEEE2383_Capacity_1_Flow_1_25_InitFail_2_SampleNum_',num2str(Train_num),'_Distribution_Result_Balanced_AC.mat'],'Distribution_Result_Record')
%     Frequency_Distribution_Record=failure_frequency_error_record;
%     save(['IEEE2383_Capacity_1_Flow_1_InitFail_2_SampleNum_',num2str(Train_num),'_Frequency_Distribution_Record_Refined_AC.mat'],'Frequency_Distribution_Record')        
%     Relative_Frequency_Distribution_Record=relative_failure_frequency_error_record;
%     save(['IEEE2383_Capacity_1_Flow_1_InitFail_2_SampleNum_',num2str(Train_num),'_Relative_Frequency_Distribution_Record_Refined.mat'],'Relative_Frequency_Distribution_Record') 
    
%     Distribution_Result_Record=[Index_Test_Sample,CESR_record,false_alarm_vector,normal_alarm_vector,error_record,failure_size_record,relative_failure_size_record];
%     save(['IEEE3012_Capacity_1_Flow_1_25_InitFail_2_SampleNum_',num2str(Train_num),'_Distribution_Result_Balanced_Refined.mat'],'Distribution_Result_Record')
%     Frequency_Distribution_Record=failure_frequency_error_record;
%     save(['IEEE3012_Capacity_1_Flow_1_25_InitFail_2_SampleNum_',num2str(Train_num),'_Frequency_Distribution_Record_Refined.mat'],'Frequency_Distribution_Record') 
%     Relative_Frequency_Distribution_Record=relative_failure_frequency_error_record;
%     save(['IEEE3012_Capacity_1_Flow_1_25_InitFail_2_SampleNum_',num2str(Train_num),'_Relative_Frequency_Distribution_Record_Refined.mat'],'Relative_Frequency_Distribution_Record')    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Small Failure (for AC)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
% %     
%     Binary_small_size_error=[CESR_record_small_size_avg,CESR_record_small_size_var,CESR_record_small_size_med,CESR_record_small_size_max,CESR_record_small_size_min]
%     false_alarm_binary_small_size_error=[false_alarm_small_size_avg,false_alarm_small_size_var,false_alarm_small_size_med,false_alarm_small_size_max,false_alarm_small_size_min]
%     Time_small_size_error=[error_record_small_size_avg,error_record_small_size_var,error_record_small_size_med,error_record_small_size_max,error_record_small_size_min]
%     Failure_size_small_error=[failure_size_small_record_avg,failure_size_small_record_var,failure_size_small_record_med,failure_size_small_record_max,failure_size_small_record_min]
%     
% %     error_record_small_size_effective_avg=sum(error_each_link_small_size_record(find(failure_exist_flag>0),1))/size(find(failure_exist_flag>0),1)
% %     CESR_record_small_size_effective_avg=sum(CESR_each_link_small_size_record(find(failure_exist_flag>0),1))/size(find(failure_exist_flag>0),1)
%     small_eff_index=find(failure_exist_flag_small>LB*Train_num & failure_exist_flag_small<UB*Train_num);
%     error_record_small_effective_avg=sum(error_each_link_small_size_record(small_eff_index,1))/size(small_eff_index,1)
%     CESR_record_small_effective_avg=sum(CESR_each_link_small_size_record(small_eff_index,1))/size(small_eff_index,1)
% %     false_alarm_small_size_effective_avg=sum(false_alarm_each_link_record(find(failure_exist_flag>0 & Link_failure_record>0.1*Test_num & Link_failure_record<0.9*Test_num),1))/size(find(failure_exist_flag>0 & Link_failure_record>0.1*Test_num & Link_failure_record<0.9*Test_num),1) %how to write this???
%     false_alarm_each_link_cum_small_sum=sum(false_alarm_each_link_cum_small(small_eff_index))
%     failure_num_vector_cum_small_sum=sum(failure_num_vector_cum_small(small_eff_index))
%     false_alarm_small_effective_avg=false_alarm_each_link_cum_small_sum/failure_num_vector_cum_small_sum

%     fail_size_record_avg=sum(fail_size_record)/Test_num
%     fail_size_record_max=max(fail_size_record)
%     fail_size_record_min=min(fail_size_record)
%     
%     figure('visible','on')
%     [counts,centers]=hist(fail_size_record,20);
%     fig_tmp=bar(centers,counts/sum(counts));

