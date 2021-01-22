function LearningMonteCarlo_Matpower_IPOPT_Large_Variant                              

%    load('italyNet.mat')
%     mpc = case118; %case1354pegase
%     mpc = case3012wp; %case1888rte;case3012wp;case1354pegase;case2383wp
    
%     load('IEEE118_Capacity_1_2_Flow_InitFail_1_New_Train.mat','cascade_train');
%     cascade_1=cascade_train;
%     load('IEEE118_Capacity_1_2_Flow_InitFail_2_New_Train.mat','cascade_train');
%     cascade_2=cascade_train;
%     load('IEEE118_Capacity_1_2_Flow_InitFail_3_New_Train.mat','cascade_train');
%     cascade_3=cascade_train;
%     load('IEEE118_Capacity_1_2_Flow_InitFail_4_New_Train.mat','cascade_train');
%     cascade_4=cascade_train;
%     load('IEEE118_Capacity_1_2_Flow_InitFail_5_New_Train.mat','cascade_train');
%     cascade_5=cascade_train;
%     cascade=[cascade_1,cascade_2,cascade_3,cascade_4,cascade_5];
% 
%     shuffle=randperm(size(cascade,2));
%     cascade=cascade(:,shuffle);
    
%     load('IEEE118_Capacity_1_2_Flow_AC_Train_Set.mat','cascade')
%     load('IEEE1354_Capacity_1_InitFail_2_Train_Balanced.mat','cascade_train')
%     cascade=cascade_train;
%     load('IEEE118_Capacity_1_5_InitFail_3_Train_Balanced.mat','cascade_train')
%     cascade=cascade_train;
    load('IEEE2383_Capacity_1_Flow_1_InitFail_2_Train_Balanced_New.mat','cascade_train')
    cascade=cascade_train;
%     load('IEEE1354_Capacity_1_Flow_2_InitFail_2_Train_Balanced_New_AC.mat','cascade_train')
%     cascade=cascade_train;
%     load('IEEE1354_Capacity_1_5_Flow_1_InitFail_2_Train_Balanced_New_AC.mat','cascade_train')
%     cascade=cascade_train;
%     load('IEEE2383_Capacity_1_5_Flow_1_InitFail_2_Train_Balanced_New_AC.mat','cascade_train')
%     cascade=cascade_train;
%     load('IEEE1354_Capacity_1_Flow_1_5_InitFail_2_Train_Balanced_New_Retry.mat','cascade_link_cum')
%     cascade_train=cascade_link_cum;
%     cascade=cascade_train;
% %     cascade_train=cascade_link_cum;

%     load('IEEE1354_Capacity_1_Flow_1_5_InitFail_2_Train_Balanced_New_AC.mat','cascade_link_cum')
%     cascade_train=cascade_link_cum;
%         cascade=cascade_train;

%     load('IEEE3012_Capacity_1_Flow_1_InitFail_2_Train_Balanced_New.mat','cascade_train')
%     cascade=cascade_train;

%     load('IEEE3012_Capacity_1_InitFail_2_New_Train_Tmp_1.mat','cascade_train');
%     cascade_1=cascade_train;
%     load('IEEE3012_Capacity_1_InitFail_2_New_Train_Tmp_2.mat','cascade_train');
%     cascade_2=cascade_train;
%     load('IEEE3012_Capacity_1_InitFail_2_New_Train_Tmp_3.mat','cascade_train');
%     cascade_3=cascade_train;
%     load('IEEE3012_Capacity_1_InitFail_2_New_Train_Tmp_4.mat','cascade_train');
%     cascade_4=cascade_train;
%     cascade=[cascade_1,cascade_2,cascade_3,cascade_4];
%     shuffle=randperm(size(cascade,2));
%     cascade=cascade(:,shuffle);
     
%     cascade_train=cascade;

%     save('IEEE3012_Capacity_1_InitFail_2_New_Train_Tmp.mat','cascade_train')    
%     load('IEEE3012_Capacity_1_InitFail_2_New_Train_Tmp.mat','cascade_train')

%     load('IEEE2383_Capacity_1_InitFail_2_Train_Balanced_AC.mat','cascade_train')
%     cascade=cascade_train;
    
    M=size(cascade,1);
    K=5000;
    
    state_series_tensor=zeros(M,25,K); 
    size_state_vector=zeros(K,1);
    weight_vector=zeros(M,1);
    failure_exist_flag=zeros(M,1);
    for i=1:K
        max_i=max(cascade(:,i));
        tmp_cascade=cascade(:,i);
        tmp_cascade(find(cascade(:,i)==max_i))=0;
        submax_i=max(tmp_cascade);
        if max_i-submax_i>=2
            size_state=submax_i+1; 
        else
            size_state=submax_i;
        end
        state_series=ones(M,size_state);
        for j=1:size_state
            state_series(find(cascade(:,i)<=j),j)=0;
        end
        state_series_tensor(:,1:size_state,i)=state_series;
        size_state_vector(i,1)=size_state;
    end
    
    state_series_sparse=cell(size(state_series_tensor,3));
    for k=1:K
        state_series_matrix=reshape(state_series_tensor(:,1:size_state_vector(k,1),k),[size(state_series_tensor,1),size_state_vector(k,1)]);
        state_series_sparse{k}=(state_series_matrix);
    end
    
    for i=1:K
        tmp_failure=zeros(M,1);
        tmp_failure(find(cascade(:,i)<size_state_vector(i,1) & cascade(:,i)>1),1)=1;
        failure_exist_flag=failure_exist_flag+tmp_failure;
    end
    
    size(find(failure_exist_flag==0),1)
    M
    size(find(failure_exist_flag==0),1)/M
    
%     Test_num=50000;
%     load(['IEEE118_Capacity_1_2_InitFail_3_SampleNum_',num2str(Test_num),'_Initial_state_Balanced_1.mat'],'Initial_state')
%     load(['IEEE118_Capacity_1_2_InitFail_3_SampleNum_',num2str(Test_num),'_Final_state_Balanced_1.mat'],'Final_state')
%     
%     fail_frequency=zeros(M,1);
%     link_index=zeros(M,1);
%     for i=1:M
%         link_index(i,1)=i;
%         fail_frequency(i,1)=size(find(Initial_state(i,:)==1 & Final_state(i,:)==0),2)/size(find(Initial_state(i,:)==1),2);
%     end    
%     link_fail_profile=[link_index,fail_frequency];
%     link_fail_low_freq=link_fail_profile(link_fail_profile(:,2)<0.1 & link_fail_profile(:,2)>0,:);
%     low_freq_link_index=link_fail_low_freq(:,1);        

    D_ini=rand(M);
    A_1_ini=zeros(M);
    A_2_ini=zeros(M);
    count_1=zeros(M,1);
    count_2=zeros(M,1);
    
    for i=1:M
        D_ini(i,i)=0;
        D_ini(i,:)=D_ini(i,:)/sum(D_ini(i,:)); 
    end

    %MonteCarlo for A1 and A2
    
    %Dealing with A1 & A2 (Eliminate Monte Carlo)
     parfor l=1:M
         l
%          tic
         count_l=0;
        for k=1:M
            k;
            if failure_exist_flag(k,1)>0
                count_l=count_l+1;
                Count_lk_1=0;
                Count_lk_2=0;
                for i=1:K
                    for j=1:size_state_vector(i,1)-1
                        if state_series_tensor(l,j,i)==1
                            Count_lk_1=Count_lk_1+1;
                            if state_series_tensor(k,j+1,i)==1 
                                A_1_ini(k,l)=A_1_ini(k,l)+1;
                            else
                               break; 
                            end
                        end
                        if state_series_tensor(l,j,i)==0
                            Count_lk_2=Count_lk_2+1;
                            if state_series_tensor(k,j+1,i)==1 
                                A_2_ini(k,l)=A_2_ini(k,l)+1;
                            else
                               break; 
                            end
                        end
                    end
                end
                if Count_lk_1~=0
                    A_1_ini(k,l)=A_1_ini(k,l)/Count_lk_1;
                else
                    A_1_ini(k,l)=0;
                end
                if Count_lk_2~=0
                    A_2_ini(k,l)=A_2_ini(k,l)/Count_lk_2;
                else
                    A_2_ini(k,l)=0;
                end
            else
                A_1_ini(k,l)=1;
                A_2_ini(k,l)=1;
            end
        end
%         toc
     end    
    
%     failure_exist_flag = failure_exist_flag & failure_exist_flag;
%     failure_exist_flag_rep=repmat(failure_exist_flag,[1,M]);
%     A_1_ini(find(failure_exist_flag_rep==0))=1;
%     A_2_ini(find(failure_exist_flag_rep==0))=1;
%     tmp_flag=find(failure_exist_flag_rep==1);  %failure_exist_flag和k有关系

% %      for l=1:M
% %          l
% %         for k=1:M
% %             k;
% %             if failure_exist_flag(k,1)>0
% % tic
%                 Count_lk_1=zeros(M,1);
%                 Count_lk_2=zeros(M,1);
%                 for i=1:K
%                     i
%                     tic
% %                     tmp_size_1=size(find(state_series_tensor(l,1:(size_state_vector(i,1)-1),i)==1),2);
% %                     tmp_size_2=size_state_vector(i,1)-1-tmp_size_1;
% %                     tmp_size_A_1=size(find(state_series_tensor(l,1:(size_state_vector(i,1)-1),i)==1 & state_series_tensor(k,2:(size_state_vector(i,1)),i)==1),2);
% %                     tmp_size_A_2=size(find(state_series_tensor(l,1:(size_state_vector(i,1)-1),i)==0 & state_series_tensor(k,2:(size_state_vector(i,1)),i)==1),2);
% %                     
% %                     Count_lk_1=Count_lk_1+tmp_size_1;
% %                     A_1_ini(k,l)=A_1_ini(k,l)+tmp_size_A_1;
% %                     Count_lk_2=Count_lk_2+tmp_size_2; %size(find(state_series_tensor(l,1:(size_state_vector(i,1)-1),i)==0),2);
% %                     A_2_ini(k,l)=A_2_ini(k,l)+tmp_size_A_2;%size(find(state_series_tensor(l,1:(size_state_vector(i,1)-1),i)==0 & state_series_tensor(k,2:(size_state_vector(i,1)),i)==1),2);
% %                     tic
%                     for j=1:size_state_vector(i,1)-1
% %                         if state_series_tensor(l,j,i)==1
%                             Count_lk_1=Count_lk_1+state_series_tensor(:,j,i);%Count_lk_1+repmat(state_series_tensor(:,j,i),[1,M]);
%                             Tmp_1=state_series_tensor(:,j+1,i)*state_series_tensor(:,j,i)';
%                             A_1_ini(tmp_flag)=A_1_ini(tmp_flag)+Tmp_1(tmp_flag); %Count_lk_1*state_series_tensor(:,j+1,i)';
% %                             if state_series_tensor(k,j+1,i)==1 
% %                                 A_1_ini(k,l)=A_1_ini(k,l)+1;
% %                             end
% %                         end
% %                         if state_series_tensor(l,j,i)==0
%                             Count_lk_2=Count_lk_2+1-state_series_tensor(:,j,i);%Count_lk_2+1;
%                             Tmp_2=state_series_tensor(:,j+1,i)*(1-state_series_tensor(:,j,i))';
%                             A_2_ini(tmp_flag)=A_2_ini(tmp_flag)+Tmp_2(tmp_flag); %Count_lk_2*state_series_tensor(:,j+1,i)';
% %                             if state_series_tensor(k,j+1,i)==1 
% %                                 A_2_ini(k,l)=A_2_ini(k,l)+1;
% %                             end
% %                         end
%                     end
%                     toc
% %                     [A_1_ini(6,1)/Count_lk_1(1),A_2_ini(6,1)/Count_lk_2(1)]
%                     1;
%                 end
%                 Count_lk_1_rep=repmat(Count_lk_1,[1,M]);
%                 Count_lk_2_rep=repmat(Count_lk_2,[1,M]);
%                 Count_lk_1_rep=Count_lk_1_rep';
%                 Count_lk_2_rep=Count_lk_2_rep';
%                 tmp_1=find(Count_lk_1_rep~=0 & failure_exist_flag_rep==1);
%                 A_1_ini(tmp_1)=A_1_ini(tmp_1)./Count_lk_1_rep(tmp_1);
%                 tmp_2=find(Count_lk_2_rep~=0 & failure_exist_flag_rep==1);
%                 A_2_ini(tmp_2)=A_2_ini(tmp_2)./Count_lk_2_rep(tmp_2);    
                
% %                 toc
% %                 A_1_ini=A_1_ini./repmat(Count_lk_1,[1,M]);
% %                 A_2_ini=A_2_ini./repmat(Count_lk_2,[1,M]);
% %                 if Count_lk_1~=0
% %                     A_1_ini(k,l)=A_1_ini(k,l)/Count_lk_1;
% %                 else
% %                     A_1_ini(k,l)=0;
% %                 end
% %                 if Count_lk_2~=0
% %                     A_2_ini(k,l)=A_2_ini(k,l)/Count_lk_2;
% %                 else
% %                     A_2_ini(k,l)=0;
% %                 end
% %             else
% %                 A_1_ini(k,l)=1;
% %                 A_2_ini(k,l)=1;
% %             end
% %         end
% %      end
%  
%      parfor l=1:M
%          l
%         for k=1:M
%             Ratio_1=zeros(K,1);
%             Ratio_2=zeros(K,2);
%             for i=1:K
%                 Count_lk_1=0;
%                 Count_lk_2=0;
%                 A_1_ini_tmp=0;
%                 A_2_ini_tmp=0;
%                 for j=1:size_state_vector(i,1)-1
%                     if state_series_tensor(l,j,i)==1
%                         Count_lk_1=Count_lk_1+1;
%                         if state_series_tensor(k,j+1,i)==1 
%                             A_1_ini_tmp=A_1_ini_tmp+1;
%                         end
%                     end
%                     if state_series_tensor(l,j,i)==0
%                         Count_lk_2=Count_lk_2+1;
%                         if state_series_tensor(k,j+1,i)==1 
%                             A_2_ini_tmp=A_2_ini_tmp+1;
%                         end
%                     end
%                 end
%                 if Count_lk_1~=0
%                     Ratio_1(i,1)=A_1_ini_tmp/Count_lk_1;
%                 else
%                     Ratio_1(i,1)=0;
%                 end
%                 if Count_lk_2~=0
%                     Ratio_2(i,1)=A_2_ini_tmp/Count_lk_2;
%                 else
%                     Ratio_2(i,1)=0;
%                 end
%             end
%             A_1_ini(k,l)=sum(sum(Ratio_1))/K;
%             A_2_ini(k,l)=sum(sum(Ratio_2))/K;
%         end
%      end
%      
     %Currently the P_1/P_2 is reverse compared with theoretical analysis!!!
    
     P_1=A_1_ini;
     P_2=A_2_ini;
     
     A_1=P_1';
     A_2=P_2';
    
    count=0;
    A_k=[ones(1,M);zeros(1,M)];
    b_k=[1;0];       
%  
%     save(['IEEE39_Capacity_0_7_AC_A_1_',num2str(K),'_New.mat'],'A_1')
%     save(['IEEE39_Capacity_0_7_AC_A_2_',num2str(K),'_New.mat'],'A_2')
%     
%     load(['IEEE118_Capacity_1_2_Flow_AC_A_1_',num2str(K),'_New.mat'],'A_1')
%     load(['IEEE118_Capacity_1_2_Flow_AC_A_2_',num2str(K),'_New.mat'],'A_2')

%     save(['IEEE1354_Capacity_1_InitFail_2_A_1_SampleNum_',num2str(K),'_Balanced.mat'],'A_1')
%     save(['IEEE1354_Capacity_1_InitFail_2_A_2_SampleNum_',num2str(K),'_Balanced.mat'],'A_2')
%      
%     load(['IEEE1354_Capacity_1_InitFail_2_A_1_SampleNum_',num2str(K),'_Balanced.mat'],'A_1')
%     load(['IEEE1354_Capacity_1_InitFail_2_A_2_SampleNum_',num2str(K),'_Balanced.mat'],'A_2')

%     save(['IEEE1354_Capacity_1_Flow_1_5_InitFail_2_A_1_SampleNum_',num2str(K),'_Balanced_New_Retry.mat'],'A_1')
%     save(['IEEE1354_Capacity_1_Flow_1_5_InitFail_2_A_2_SampleNum_',num2str(K),'_Balanced_New_Retry.mat'],'A_2')
% %      
%     load(['IEEE1354_Capacity_1_Flow_1_5_InitFail_2_A_1_SampleNum_',num2str(K),'_Balanced_New_Retry.mat'],'A_1')
%     load(['IEEE1354_Capacity_1_Flow_1_5_InitFail_2_A_2_SampleNum_',num2str(K),'_Balanced_New_Retry.mat'],'A_2')

%     save(['IEEE2383_Capacity_1_5_Flow_1_InitFail_2_A_1_SampleNum_',num2str(K),'_Balanced_New_AC.mat'],'A_1')
%     save(['IEEE2383_Capacity_1_5_Flow_1_InitFail_2_A_2_SampleNum_',num2str(K),'_Balanced_New_AC.mat'],'A_2')
% %      
%     load(['IEEE2383_Capacity_1_5_Flow_1_InitFail_2_A_1_SampleNum_',num2str(K),'_Balanced_New_AC.mat'],'A_1')
%     load(['IEEE2383_Capacity_1_5_Flow_1_InitFail_2_A_2_SampleNum_',num2str(K),'_Balanced_New_AC.mat'],'A_2')

%     save(['IEEE118_Capacity_1_5_InitFail_3_A_1_SampleNum_',num2str(K),'_Balanced.mat'],'A_1')
%     save(['IEEE118_Capacity_1_5_InitFail_3_A_2_SampleNum_',num2str(K),'_Balanced.mat'],'A_2')
%      
%     load(['IEEE118_Capacity_1_5_InitFail_3_A_1_SampleNum_',num2str(K),'_Balanced.mat'],'A_1')
%     load(['IEEE118_Capacity_1_5_InitFail_3_A_2_SampleNum_',num2str(K),'_Balanced.mat'],'A_2')

    save(['IEEE2383_Capacity_1_Flow_1_InitFail_2_A_1_SampleNum_',num2str(K),'_Balanced_New.mat'],'A_1')
    save(['IEEE2383_Capacity_1_Flow_1_InitFail_2_A_2_SampleNum_',num2str(K),'_Balanced_New.mat'],'A_2')
     
    load(['IEEE2383_Capacity_1_Flow_1_InitFail_2_A_1_SampleNum_',num2str(K),'_Balanced_New.mat'],'A_1')
    load(['IEEE2383_Capacity_1_Flow_1_InitFail_2_A_2_SampleNum_',num2str(K),'_Balanced_New.mat'],'A_2')

%     save(['IEEE3012_Capacity_1_Flow_1_InitFail_2_A_1_SampleNum_',num2str(K),'_Balanced_New.mat'],'A_1')
%     save(['IEEE3012_Capacity_1_Flow_1_InitFail_2_A_2_SampleNum_',num2str(K),'_Balanced_New.mat'],'A_2')
% %      
%     load(['IEEE3012_Capacity_1_Flow_1_InitFail_2_A_1_SampleNum_',num2str(K),'_Balanced_New.mat'],'A_1')
%     load(['IEEE3012_Capacity_1_Flow_1_InitFail_2_A_2_SampleNum_',num2str(K),'_Balanced_New.mat'],'A_2')
    
    for i=1:size(A_2,1)
        A_2(i,i)=0; 
    end

%     load(['IEEE1888_Capacity_1_InitFail_2_A_1_SampleNum_5000_New.mat'],'A_1')
%     load(['IEEE1888_Capacity_1_InitFail_2_A_2_SampleNum_5000_New.mat'],'A_2')

    P_1=A_1;
    P_2=A_2;

    f=zeros(M,1);
    
    %%Adjustment of D: Eliminate the always safe ones...
    Always_safe_index=zeros(M,1);
    Always_safe_index(find(failure_exist_flag==0),1)=1;
    size(find(Always_safe_index==1),1)/M;
    
    size_concatenate=0;
    for q=1:K
        size_concatenate=size_concatenate+(size_state_vector(q)-1);
    end
    
        for k=1:M %k imp
            if Always_safe_index(k,1)==0
                tic
                k

                A_tmp=A_k;
                A_tmp(1,k)=0;
                A_tmp(2,k)=1;

                A_1_vec=P_1(:,k);
                A_2_vec=P_2(:,k);

                state_t_plus_1=zeros(1,size_concatenate);
                count_t=0;
                for q=1:K
                    state_t_plus_1((count_t+1):(count_t+(size_state_vector(q)-1)))=state_series_sparse{q}(k,2:size_state_vector(q)); 
                    count_t=count_t+(size_state_vector(q)-1);
                end
                
                state_t=zeros(M,size_concatenate);
                count_t=0;
                for q=1:K
                    A_1_vec_rep=repmat(A_1_vec,[1,(size_state_vector(q)-1)]);
                    A_2_vec_rep=repmat(A_2_vec,[1,(size_state_vector(q)-1)]);
                    state_t(:,(count_t+1):(count_t+(size_state_vector(q)-1)))=A_1_vec_rep.*state_series_sparse{q}(:,1:(size_state_vector(q)-1))+A_2_vec_rep.*(1-state_series_sparse{q}(:,1:(size_state_vector(q)-1)));
                    count_t=count_t+(size_state_vector(q)-1);
                end
                
                state_t_prod=state_t*state_t';
                state_t_plus_1_prod=state_t*state_t_plus_1';

%                 Hess=2*state_t_prod;
%                 Hess=sparse(tril(Hess));
%                 save('Hess_tmp.mat','Hess');

%                 Hess = zeros(size(A_1_vec,1)); 
%                 diff=A_1_vec(:,1)-A_2_vec(:,1);
%                 for l=1:K
%                     for t=1:size_state_vector(l,1)-1
%                         tmp=state_series_tensor(:,t,l).*(diff)+A_2_vec(:,1);
%                         Hess=Hess+2*weight_vector(l,1)*(tmp*tmp');
%                     end
%                 end
%     %             K*sum(size_state_vector(:,1))
%                 Hess=tril(sparse(1/(K)*(sigma*(Hess)+lambda(1)*zeros(M)+lambda(2)*zeros(M))));
%                 save('Hess_tmp.mat','Hess');

                lambda=0.4;
                % x_0
                x0=(D_ini(k,:))';

                %options
                options.lb = zeros(M,1);
                options.ub = ones(M,1);
                options.cl = [b_k];
                options.cu = [b_k];
    %             options.zl = ones(M,1);
    %             options.zu = ones(M,1);
    %             options.lambda= ones(2,1);

                % Set the IPOPT options.
                options.ipopt.jac_c_constant        = 'yes';
                options.ipopt.hessian_constant      = 'yes';
                options.ipopt.jac_d_constant        = 'yes';
                options.ipopt.hessian_approximation = 'limited-memory'; %'limited-memory';
                options.ipopt.mu_strategy           = 'adaptive';
                options.ipopt.tol                   = 1e-7;
                options.ipopt.print_level           = 0;   

                % funcs
                funcs.objective=@(Coeff)fun_full_D_weighted_function_vector(Coeff,state_t_plus_1,state_t,K);
                funcs.constraints=@(Coeff)[A_tmp*Coeff]; %;Coeff'*Coeff
                funcs.gradient=@(Coeff)fun_full_D_weighted_gradient_vector(Coeff,state_t_plus_1_prod,state_t_prod,K);
                funcs.jacobian=@(Coeff)jacobian(Coeff,A_tmp);
                funcs.jacobianstructure=@()sparse([A_tmp]); %;ones(1,M)
%                 funcs.hessian=@fun_full_D_weighted_hessian; %(Coeff,sigma,lambda)%,A_1_vec,A_2_vec,K,size_state_vector,state_series_tensor,weight_vector);
%                 funcs.hessianstructure=@()sparse(tril(ones(M)));
                
                %options=optimoptions('fmincon','Algorithm','interior-point','GradObj','on','HessianFcn',@(Coeff,lambda)HessianD(Coeff,lambda,Hess));
                %options=optimoptions('fmincon','MaxFunEvals',2000,'Algorithm','trust-region-reflective','GradObj','on','SpecifyObjectiveGradient',true,'HessianFcn','objective');
    %             [x,fval,exitflag,info] = opti_ipopt(nlprob,x0)
                nlprob.funcs=funcs;
                nlprob.options=options;
                [x_tmp,fval,exitflag,info]=opti_ipopt(nlprob,x0); 
    %             Opt=opti(nlprob,x0);
    %             [x_tmp,fval,exitflag,info] = solve(Opt);
                f(k,1)=fval;
                fval
                D_ini(k,:)=x_tmp';

                toc
                1;
            else
               D_ini(k,:)=zeros(1,M);
               D_ini(k,k)=1;
            end
        end
        
        D=D_ini;
%         save(['IEEE118_Capacity_1_5_InitFail_3_D_SampleNum_',num2str(K),'_Balanced.mat'],'D')
        save(['IEEE2383_Capacity_1_Flow_1_InitFail_2_D_SampleNum_',num2str(K),'_Balanced_New.mat'],'D')
%         save(['IEEE1354_Capacity_1_Flow_1_5_InitFail_2_D_SampleNum_',num2str(K),'_Balanced_New_Retry.mat'],'D')
%         save(['IEEE2383_Capacity_1_5_Flow_1_InitFail_2_D_SampleNum_',num2str(K),'_Balanced_New_AC.mat'],'D')
%         save(['IEEE3012_Capacity_1_Flow_1_InitFail_2_D_SampleNum_',num2str(K),'_Balanced_New.mat'],'D')

end

function J = jacobian (Coeff,A_k) 
    J = sparse([A_k]); %;2*Coeff'
end
% 
% function Hess = fun_full_D_weighted_hessian(Coeff,A_1_vec,A_2_vec,K,size_state_vector,state_series_tensor,weight_vector,sigma,lambda)
%     Hess = zeros(size(A_1_vec,1)); 
%     diff=A_1_vec(:,1)-A_2_vec(:,1);
%     for l=1:K
%         for t=1:size_state_vector(l,1)-1
%             tmp=state_series_tensor(:,t,l).*(diff)+A_2_vec(:,1);
%             Hess=Hess+2*weight_vector(l,1)*(tmp*tmp');
%         end
%     end
%     Hess=sparse(sigma*tril(Hess)+lambda(1)*zeros(M)+lambda(2)*zeros(M));
% end









