    figure

%     load(['IEEE2383_Capacity_1_InitFail_2_SampleNum_',num2str(Test_num),'_Low_Freq_Index.mat'],'low_freq_link_index')
    
%     load('IEEE1354_Capacity_1_InitFail_2_New_Train_Tmp.mat','cascade_train')
%     cascade_1=cascade_train;
%     size(cascade_1);
%     load('IEEE2383_Capacity_1_InitFail_2_New_Train_Tmp_2.mat','cascade_train')
%     cascade_2=cascade_train;    
%     load('IEEE2383_Capacity_1_InitFail_2_New_Test_Tmp_2.mat','cascade_test')
%     cascade_3=cascade_test;
%     cascade=[cascade_2,cascade_3];
    
%     load('IEEE1354_Capacity_1_Flow_1_5_InitFail_2_Train_Balanced_New_AC.mat','cascade_link_cum')
%     cascade_train=cascade_link_cum;
%     cascade=cascade_train;
%     load('IEEE1354_Capacity_1_5_Flow_1_InitFail_2_Train_Balanced_New.mat','cascade_link_cum')
%     cascade=cascade_link_cum;
%     load('IEEE1354_Capacity_1_5_Flow_1_InitFail_2_Train_Balanced_New_AC.mat','cascade_train')
%     cascade=cascade_train;
%     load('IEEE118_Capacity_1_8_InitFail_3_Train_Balanced_New.mat','cascade_train')
%     cascade=cascade_train;
%     load('IEEE2383_Capacity_1_InitFail_2_Train_Balanced_AC.mat','cascade_train')
%     cascade=cascade_train;
%     load('IEEE2383_Capacity_1_InitFail_2_Train_Balanced.mat','cascade_train')
%     cascade=cascade_train; 
%     load('IEEE2383_Capacity_1_5_Flow_1_InitFail_2_Train_Balanced_New_AC.mat','cascade_train')
%     cascade=cascade_train;
%     load('IEEE1354_Capacity_1_InitFail_2_Train_Balanced.mat','cascade_train')
%     cascade=cascade_train;
%     load('IEEE3012_Capacity_1_Flow_1_InitFail_2_Train_Balanced_New.mat','cascade_train')
%     cascade=cascade_train;

%     load('IEEE2383_Capacity_1_InitFail_2_Test_Balanced_AC.mat','cascade_test')
%     cascade=cascade_test;

    %Predicted Failure Size Distribution
    %load(['IEEE2383_Capacity_1_Flow_1_25_InitFail_2_Predicted_Failure_Size_Balanced_New_AC.mat'],'predicted_failure_size_record')
    load(['IEEE1354_Capacity_1_Flow_2_InitFail_2_Real_Failure_Size_Balanced_New.mat'],'real_failure_size_record')
    %Real Failure Size Distribution
    %load(['IEEE2383_Capacity_1_Flow_1_25_InitFail_2_Real_Failure_Size_Balanced_New_AC.mat'],'real_failure_size_record')
    load(['IEEE1354_Capacity_1_Flow_2_InitFail_2_Predicted_Failure_Size_Balanced_New.mat'],'predicted_failure_size_record')

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
  
%     Test_num=50000;
%     load(['IEEE118_Capacity_1_2_InitFail_3_SampleNum_',num2str(Test_num),'_Initial_state_Balanced_DiffTrainingSet_2_var.mat'],'Initial_state')
%     load(['IEEE118_Capacity_1_2_InitFail_3_SampleNum_',num2str(Test_num),'_Final_state_Balanced_DiffTrainingSet_2_var.mat'],'Final_state')
%     load(['IEEE118_Capacity_1_2_InitFail_3_SampleNum_',num2str(Test_num),'_Epsilon_opt_Balanced_DiffTrainingSet_2_var.mat'],'Epsilon_opt')
%     load(['IEEE118_Capacity_1_2_InitFail_3_SampleNum_',num2str(Test_num),'_Fail_Index_Set_Balanced_DiffTrainingSet.mat'],'Link_Fail_Index_Set')

%    Test_num=30000;
%     load(['IEEE118_Capacity_1_8_InitFail_3_SampleNum_',num2str(Test_num),'_Initial_state_Balanced_New.mat'],'Initial_state')
%     load(['IEEE1354_Capacity_1_Flow_1_5_InitFail_2_SampleNum_',num2str(Test_num),'_Initial_state_Balanced_AC.mat'],'Initial_state')
%    load(['IEEE2383_Capacity_1_Flow_1_25_InitFail_2_SampleNum_',num2str(Test_num),'_Initial_state_Balanced_AC.mat'],'Initial_state')
%     load(['IEEE3012_Capacity_1_Flow_1_InitFail_2_SampleNum_',num2str(Test_num),'_Initial_state_Balanced_New.mat'],'Initial_state')

%     load(['IEEE2383_Capacity_1_SampleNum_40000_Initial_state_A_2.mat'],'Initial_state')
%     find(Initial_state(:,1)==0)
%     load(['IEEE1354_Capacity_1_Flow_1_5_InitFail_2_SampleNum_',num2str(Test_num),'_Initial_state_Balanced_AC.mat'],'Initial_state')
%     load(['IEEE1354_Capacity_1_5_Flow_1_InitFail_2_SampleNum_',num2str(Test_num),'_Initial_state_Balanced_New_AC.mat'],'Initial_state')
%     load(['IEEE2383_Capacity_1_5_Flow_1_InitFail_2_SampleNum_',num2str(Test_num),'_Initial_state_Balanced_New_AC.mat'],'Initial_state')

%     load(['IEEE118_Capacity_1_8_InitFail_3_SampleNum_',num2str(Test_num),'_Initial_state_Balanced_DiffTrainingSet.mat'],'Initial_state')
%     load(['IEEE118_Capacity_1_2_InitFail_3_SampleNum_',num2str(Test_num),'_Final_state_Balanced_DiffTrainingSet.mat'],'Final_state')
%     load(['IEEE118_Capacity_1_2_InitFail_3_SampleNum_',num2str(Test_num),'_Epsilon_opt_Balanced_DiffTrainingSet.mat'],'Epsilon_opt')
%     load(['IEEE118_Capacity_1_2_InitFail_3_SampleNum_',num2str(Test_num),'_Fail_Index_Set_Balanced_DiffTrainingSet.mat'],'Link_Fail_Index_Set')

%     Test_num=50000;
%     load(['IEEE2383_Capacity_1_Flow_1_5_InitFail_2_SampleNum_',num2str(Test_num),'_Initial_state_Balanced.mat'],'Initial_state')
%     load(['IEEE2383_Capacity_1_Flow_1_5_InitFail_2_SampleNum_',num2str(Test_num),'_Final_state_Balanced.mat'],'Final_state')
%     load(['IEEE2383_Capacity_1_Flow_1_5_InitFail_2_SampleNum_',num2str(Test_num),'_Epsilon_opt_Balanced.mat'],'Epsilon_opt')
%     load(['IEEE2383_Capacity_1_Flow_1_5_InitFail_2_SampleNum_',num2str(Test_num),'_Beta_Balanced.mat'],'Beta_record')
% 
%     Test_num=50000;
%     load(['IEEE1354_Capacity_1_Flow_2_InitFail_2_SampleNum_',num2str(Test_num),'_Initial_state_Balanced_New.mat'],'Initial_state')
%     load(['IEEE1354_Capacity_1_Flow_2_InitFail_2_SampleNum_',num2str(Test_num),'_Final_state_Balanced_New.mat'],'Final_state')
%     load(['IEEE1354_Capacity_1_Flow_2_InitFail_2_SampleNum_',num2str(Test_num),'_Epsilon_opt_Balanced_New.mat'],'Epsilon_opt')

%     Initial_index_lf=cell(size(low_freq_link_index,1),1);
%     for q=1:size(low_freq_link_index,1)
%        Initial_index_lf{q}=zeros(3,size(low_freq_sample_index{q},2));
%        for k=1:size(low_freq_sample_index{q},2)
%            Initial_index_lf{q}(:,k)=(find(Initial_state(:,low_freq_sample_index{q}(k))==0))';
%        end
%     end
%     save(['IEEE118_Capacity_1_2_InitFail_3_SampleNum_',num2str(Test_num),'_Initial_index_lf_Balanced.mat'],'Initial_index_lf')
%     load(['IEEE118_Capacity_1_2_InitFail_3_SampleNum_',num2str(Test_num),'_Initial_index_lf_Balanced.mat'],'Initial_index_lf')

% %     Train_num=36000;
% %     Initial_state=Initial_state(:,1:Train_num);
%     Train_num=500;
%     Initial_state=Initial_state(:,end-Train_num+1:end);
%     
%     
%     %    cascade=cascade(:,end-Test_num+1:end);
% 
%     M=size(Initial_state,1);
%     Link_Index=zeros(M,1);
%     for i=1:M
%         Link_Index(i,1)=i;
%     end
% 
%     K_prev=Train_num;
% 
%     fail_size_record=zeros(K_prev,1);
%     record_test=zeros(K_prev,3);
%     for k=1:K_prev 
%         real_vector=cascade(:,k);
%         
%         real_final_state=zeros(M,1);
%         max_real=max(real_vector(:,1));
%         tmp_real=real_vector(:,1);
%         tmp_real(real_vector(:,1)==max_real)=0;
%         submax_real=max(tmp_real);
%         if max_real-submax_real>=2
%             real_final_state(real_vector(:,1)==max_real,1)=1;
%         end        
%         
%         fail_size_record(k,1)=size(find(real_final_state==0),1);
%         Initial_failure_index=find(Initial_state(:,k)==0);
%         record_test(k,:)=[fail_size_record(k,1),Initial_failure_index'];
%     end
%     
%     %Compare the small fail size record and the initial failures
%     (sortrows(record_test,1))
    
%     save('IEEE118_Capacity_1_8_InitFail_3_Total_Initial_Outage_Balanced_New.mat','fail_size_record')
%     save('IEEE1354_Capacity_1_Flow_1_5_InitFail_2_Total_Initial_Outage_Balanced_New_AC.mat','fail_size_record')
 %   save('IEEE2383_Capacity_1_Flow_1_25_InitFail_2_Total_Initial_Outage_Balanced_New_AC.mat','fail_size_record')
%     save('IEEE3012_Capacity_1_Flow_1_InitFail_2_Total_Initial_Outage_Balanced_New.mat','fail_size_record')
%     save('IEEE1354_Capacity_1_Flow_1_5_InitFail_2_Total_Initial_Outage_Selective_New_AC.mat','fail_size_record')
%     save('IEEE1354_Capacity_1_5_Flow_1_InitFail_2_Total_Initial_Outage_Balanced_New_AC.mat','fail_size_record')
%     save('IEEE2383_Capacity_1_5_Flow_1_InitFail_2_Total_Initial_Outage_Balanced_New_AC.mat','fail_size_record')

    %%%%
%     load('IEEE118_Capacity_1_8_InitFail_3_Total_Initial_Outage_Balanced_New.mat','fail_size_record')
%     total_fail_size_record=fail_size_record;
%     load('IEEE118_Capacity_1_8_InitFail_3_Max_Initial_Outage_Balanced_New_2.mat','fail_size_record')
%     max_fail_size_record=fail_size_record;
%     load('IEEE118_Capacity_1_8_InitFail_3_Min_Initial_Outage_Balanced_New_2.mat','fail_size_record')
%     min_fail_size_record=fail_size_record;
% 
%     load('IEEE1354_Capacity_1_Flow_1_5_InitFail_2_Total_Initial_Outage_Balanced_New_AC.mat','fail_size_record')
%     total_fail_size_record=fail_size_record;
    load('IEEE1354_Capacity_1_Flow_2_InitFail_2_Max_Initial_Outage_Balanced_New_2.mat','fail_size_record')
    max_fail_size_record=fail_size_record;
%     load('IEEE1354_Capacity_1_Flow_1_5_InitFail_2_Min_Initial_Outage_Balanced_New_2.mat','fail_size_record')
%     min_fail_size_record=fail_size_record;

%     load('IEEE1354_Capacity_1_5_Flow_1_InitFail_2_Total_Initial_Outage_Balanced_New_AC.mat','fail_size_record')
%     total_fail_size_record=fail_size_record;
%     sort(total_fail_size_record','ascend');
%     load('IEEE1354_Capacity_1_5_Flow_1_InitFail_2_Max_Initial_Outage_Balanced_New_AC.mat','fail_size_record')
%     max_fail_size_record=fail_size_record;
%     load('IEEE1354_Capacity_1_5_Flow_1_InitFail_2_Min_Initial_Outage_Balanced_New_AC.mat','fail_size_record')
%     min_fail_size_record=fail_size_record;

%     load('IEEE1354_Capacity_1_Flow_1_5_InitFail_2_Total_Initial_Outage_Selective_New.mat','fail_size_record')
%     total_fail_size_record=fail_size_record;
%     sort(total_fail_size_record','ascend');
%     load('IEEE1354_Capacity_1_Flow_1_5_InitFail_2_Max_Initial_Outage_Selective_New_AC.mat','fail_size_record')
%     max_fail_size_record=fail_size_record;
%     load('IEEE1354_Capacity_1_Flow_1_5_InitFail_2_Min_Initial_Outage_Selective_New.mat','fail_size_record')
%     min_fail_size_record=fail_size_record;

%     load('IEEE2383_Capacity_1_Flow_1_25_InitFail_2_Total_Initial_Outage_Balanced_New_AC.mat','fail_size_record')
%     total_fail_size_record=fail_size_record;
     total_fail_size_record=real_failure_size_record;
%     load('IEEE2383_Capacity_1_Flow_1_25_InitFail_2_Max_Initial_Outage_Balanced_New_AC_2.mat','fail_size_record')
 %   max_fail_size_record=fail_size_record;
%     load('IEEE2383_Capacity_1_Flow_1_InitFail_2_Min_Initial_Outage_Balanced_New_2.mat','fail_size_record')
%     min_fail_size_record=fail_size_record;

%     load('IEEE2383_Capacity_1_5_Flow_1_InitFail_2_Total_Initial_Outage_Balanced_New_AC.mat','fail_size_record')
%     total_fail_size_record=fail_size_record;
%     sort(total_fail_size_record','ascend');
%     load('IEEE2383_Capacity_1_5_Flow_1_InitFail_2_Max_Initial_Outage_Balanced_New_AC.mat','fail_size_record')
%     max_fail_size_record=fail_size_record;
%     load('IEEE2383_Capacity_1_5_Flow_1_InitFail_2_Min_Initial_Outage_Balanced_New_AC.mat','fail_size_record')
%     min_fail_size_record=fail_size_record;

%     load('IEEE3012_Capacity_1_Flow_1_InitFail_2_Total_Initial_Outage_Balanced_New.mat','fail_size_record')
%     total_fail_size_record=fail_size_record;
%     load('IEEE3012_Capacity_1_Flow_1_InitFail_2_Max_Initial_Outage_Balanced_New_2.mat','fail_size_record')
%     max_fail_size_record=fail_size_record;
%     load('IEEE3012_Capacity_1_Flow_1_25_InitFail_2_Min_Initial_Outage_Balanced_New_2.mat','fail_size_record')
%     min_fail_size_record=fail_size_record;

%     xbins=min(total_fail_size_record):max(total_fail_size_record); %max(total_fail_size_record);
    interval=(max(max_fail_size_record)-0*min(max_fail_size_record))/300;
    xbins=(0*min(max_fail_size_record)):interval:max(max_fail_size_record); %max(total_fail_size_record);
%     max(total_fail_size_record)
%     subplot(1,2,1)
%     xbins2=0.025:0.05:0.975;
    [counts,centers] = hist(total_fail_size_record,xbins);
    bar(centers, counts / sum(counts))
    hold on
    [counts,centers] = hist(max_fail_size_record,xbins);
    bar(centers, counts / sum(counts))
    hold on
    [counts,centers] = hist(predicted_failure_size_record,xbins);
    plot(centers,counts / sum(counts),'r*-')
    
    xlabel('Failure Size')
    ylabel('Frequency')
    legend('p_0: from S_{train}','p_{max}: from top-10 links','predicted failure size distribution')
%     axis([0,1,0,0.4]);

    %Ratio that exceeds the median value
    Median_threshold=median(total_fail_size_record);
    Median_threshold
    total_fail_size_record_sort=sort(total_fail_size_record,'ascend');
    Three_over_Four_threshold=quantile(total_fail_size_record_sort,3/4-0.1);
    Three_over_Four_threshold
    max_fail_size_record
    Ratio_surpass_median=size(find(max_fail_size_record>=Median_threshold),1)/size(max_fail_size_record,1);
    Ratio_surpass_3over4=size(find(max_fail_size_record>=Three_over_Four_threshold),1)/size(max_fail_size_record,1);
    [Ratio_surpass_median,Ratio_surpass_3over4]
    
%     text(100,0.005,['Average failure size: ',num2str(mean(total_fail_size_record)),''])
%     axis([0,1,0,0.4]);    
%     subplot(1,2,2)
% %     xbins2=0.025:0.05:0.975;
%     [counts,centers] = hist(max_fail_size_record,xbins);
%     bar(centers, counts / sum(counts))
%     xlabel('Failure Size')
%     ylabel('Frequency')
%     text(100,0.05,['Average failure size: ',num2str(mean(max_fail_size_record)),''])
%     axis([0,1,0,0.4]);
%     subplot(2,3,3)
% %     xbins2=0.025:0.05:0.975;
%     [counts,centers] = hist(min_fail_size_record,xbins);
%     bar(centers, counts / sum(counts))
%     xlabel('Failure Size')
%     ylabel('Frequency')
% %     text(100,0.05,['Average failure size: ',num2str(mean(min_fail_size_record)),''])
% %     axis([0,1,0,0.4]);
%     


