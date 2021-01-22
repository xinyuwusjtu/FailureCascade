function Property_Critical_Component

%     load('IEEE118_Capacity_1_8_InitFail_3_Train_Balanced_New.mat','cascade_train')
%     load('IEEE118_Capacity_1_8_InitFail_3_Test_Balanced_New.mat','cascade_test')
    
%     K=10000;
%     load(['IEEE118_Capacity_1_8_InitFail_3_D_SampleNum_',num2str(K),'_Balanced_New.mat'],'D')
%     load(['IEEE118_Capacity_1_8_InitFail_3_A_1_SampleNum_',num2str(K),'_Balanced_New.mat'],'A_1')
%     load(['IEEE118_Capacity_1_8_InitFail_3_A_2_SampleNum_',num2str(K),'_Balanced_New.mat'],'A_2')

%     K=5000;
%     load(['IEEE1354_Capacity_1_InitFail_2_D_SampleNum_',num2str(K),'_Balanced.mat'],'D')
%     load(['IEEE1354_Capacity_1_InitFail_2_A_1_SampleNum_',num2str(K),'_Balanced.mat'],'A_1')
%     load(['IEEE1354_Capacity_1_InitFail_2_A_2_SampleNum_',num2str(K),'_Balanced.mat'],'A_2')
    
    K=10000;
    load(['IEEE1354_Capacity_1_Flow_1_5_InitFail_2_D_SampleNum_',num2str(K),'_Balanced_New_AC.mat'],'D')
    load(['IEEE1354_Capacity_1_Flow_1_5_InitFail_2_A_1_SampleNum_',num2str(K),'_Balanced_New_AC.mat'],'A_1')
    load(['IEEE1354_Capacity_1_Flow_1_5_InitFail_2_A_2_SampleNum_',num2str(K),'_Balanced_New_AC.mat'],'A_2')

%     K=5000;
%     load(['IEEE2383_Capacity_1_Flow_1_InitFail_2_D_SampleNum_',num2str(K),'_Balanced_New_AC.mat'],'D')
%     load(['IEEE2383_Capacity_1_Flow_1_InitFail_2_A_1_SampleNum_',num2str(K),'_Balanced_New_AC.mat'],'A_1')
%     load(['IEEE2383_Capacity_1_Flow_1_InitFail_2_A_2_SampleNum_',num2str(K),'_Balanced_New_AC.mat'],'A_2')

%     K=5000;
%     load(['IEEE2383_Capacity_1_InitFail_2_D_SampleNum_',num2str(K),'_New.mat'],'D')
%     load(['IEEE2383_Capacity_1_InitFail_2_A_1_SampleNum_',num2str(K),'_New.mat'],'A_1')
%     load(['IEEE2383_Capacity_1_InitFail_2_A_2_SampleNum_',num2str(K),'_New.mat'],'A_2')
    
    M=size(D,1);
  
    for i=1:size(D,1)
        A_2(i,i)=0; 
        D(i,i)=0;
    end
    
    P_1=A_1';
    P_2=A_2';
    
    M_est=D.*(P_1-P_2);
    b_est=zeros(M,1);
    for i=1:M
        b_est(i,1)=D(i,:)*P_2(i,:)'; 
    end    
    
    Link_index=zeros(M,1);
    for i=1:M
        Link_index(i,1)=i;
    end
    
    M_tmp=M_est;
    b_tmp=repmat(b_est,[1,M]);
    
    M_vec=sum((M_tmp),1);
%     M_vec=max(M_tmp,[],1);
%     M_vec(M_vec<-0.1)=[];

    M_profile=[Link_index,M_vec'];
    M_vec_sort=sortrows(M_profile,2);
%     M_vec_sort(:,1)
       
%     M_profile
%     M_vec_sort(1:100,:)
%     sort(M_tmp(:,437),'ascend')
    
%     save('IEEE118_Capacity_1_8_InitFail_3_M_vec_sort_Balanced_New_2.mat','M_vec_sort')
    save('IEEE1354_Capacity_1_Flow_1_5_InitFail_2_M_vec_sort_Balanced_New_AC_2.mat','M_vec_sort')    
%     save('IEEE2383_Capacity_1_Flow_1_InitFail_2_M_vec_sort_Balanced_New_AC_2.mat','M_vec_sort')    

end

