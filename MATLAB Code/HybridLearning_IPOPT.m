    %% Learning Framework: Monte Carlo and Quadratic Optimization
    load('IEEE2383_Capacity_1_Flow_1_InitFail_2_Train_Balanced_New.mat','cascade_train')
    cascade=cascade_train;
    
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

    D_ini=rand(M);
    A_1_ini=zeros(M);
    A_2_ini=zeros(M);
    count_1=zeros(M,1);
    count_2=zeros(M,1);
    
    for i=1:M
        D_ini(i,i)=0;
        D_ini(i,:)=D_ini(i,:)/sum(D_ini(i,:)); 
    end

     parfor l=1:M
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
     end    
    
     P_1=A_1_ini;
     P_2=A_2_ini;
     
     A_1=P_1';
     A_2=P_2';
    
    count=0;
    A_k=[ones(1,M);zeros(1,M)];
    b_k=[1;0];       

    save(['IEEE2383_Capacity_1_Flow_1_InitFail_2_A_1_SampleNum_',num2str(K),'_Balanced_New.mat'],'A_1')
    save(['IEEE2383_Capacity_1_Flow_1_InitFail_2_A_2_SampleNum_',num2str(K),'_Balanced_New.mat'],'A_2')
     
    load(['IEEE2383_Capacity_1_Flow_1_InitFail_2_A_1_SampleNum_',num2str(K),'_Balanced_New.mat'],'A_1')
    load(['IEEE2383_Capacity_1_Flow_1_InitFail_2_A_2_SampleNum_',num2str(K),'_Balanced_New.mat'],'A_2')
  
    for i=1:size(A_2,1)
        A_2(i,i)=0; 
    end
    
    P_1=A_1;
    P_2=A_2;

    f=zeros(M,1);
    
    Always_safe_index=zeros(M,1);
    Always_safe_index(find(failure_exist_flag==0),1)=1;

    size_concatenate=0;
    for q=1:K
        size_concatenate=size_concatenate+(size_state_vector(q)-1);
    end
    
        for k=1:M
            if Always_safe_index(k,1)==0

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

                lambda=0.4;
                x0=(D_ini(k,:))';

                %options
                options.lb = zeros(M,1);
                options.ub = ones(M,1);
                options.cl = [b_k];
                options.cu = [b_k];

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
                funcs.constraints=@(Coeff)[A_tmp*Coeff];
                funcs.gradient=@(Coeff)fun_full_D_weighted_gradient_vector(Coeff,state_t_plus_1_prod,state_t_prod,K);
                funcs.jacobian=@(Coeff)jacobian(Coeff,A_tmp);
                funcs.jacobianstructure=@()sparse([A_tmp]);

                nlprob.funcs=funcs;
                nlprob.options=options;
                [x_tmp,fval,exitflag,info]=opti_ipopt(nlprob,x0); 
                f(k,1)=fval;
                D_ini(k,:)=x_tmp';
            else
               D_ini(k,:)=zeros(1,M);
               D_ini(k,k)=1;
            end
        end
        
        D=D_ini;
    save(['IEEE2383_Capacity_1_Flow_1_InitFail_2_D_SampleNum_',num2str(K),'_Balanced_New.mat'],'D')
        
function J = jacobian (Coeff,A_k) 
    J = sparse([A_k]); 
end