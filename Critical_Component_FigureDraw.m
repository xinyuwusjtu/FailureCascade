    %% Critical Component Figure Draw
    figure

    load(['IEEE1354_Capacity_1_Flow_2_InitFail_2_Real_Failure_Size_Balanced_New.mat'],'real_failure_size_record')
    load(['IEEE1354_Capacity_1_Flow_2_InitFail_2_Predicted_Failure_Size_Balanced_New.mat'],'predicted_failure_size_record')

    load('IEEE1354_Capacity_1_Flow_2_InitFail_2_Max_Initial_Outage_Balanced_New_2.mat','fail_size_record')
    max_fail_size_record=fail_size_record;

    total_fail_size_record=real_failure_size_record;

    interval=(max(max_fail_size_record)-0*min(max_fail_size_record))/300;
    xbins=(0*min(max_fail_size_record)):interval:max(max_fail_size_record); %max(total_fail_size_record);
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

%     %Ratio that exceeds the median value
%     Median_threshold=median(total_fail_size_record);
%     Median_threshold
%     total_fail_size_record_sort=sort(total_fail_size_record,'ascend');
%     Three_over_Four_threshold=quantile(total_fail_size_record_sort,3/4-0.1);
%     Three_over_Four_threshold
%     max_fail_size_record
%     Ratio_surpass_median=size(find(max_fail_size_record>=Median_threshold),1)/size(max_fail_size_record,1);
%     Ratio_surpass_3over4=size(find(max_fail_size_record>=Three_over_Four_threshold),1)/size(max_fail_size_record,1);
%     [Ratio_surpass_median,Ratio_surpass_3over4]


