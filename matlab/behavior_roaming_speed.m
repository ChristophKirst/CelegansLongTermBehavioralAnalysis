clear all;
load ;

history=30;

for j=1:
    display(j);
    a=individual_speed_AV{j};
    [s1 s2]=size(a);
    Speed_AV=zeros(s1,3);
    for i=1+history/2:s1-history/2
        
        
        Speed_AV(i,1)=mean(a(i-history/2:i+history/2,1));
        Speed_AV(i,2)=mean(a(i-history/2:i+history/2,2));
        Speed_AV(i,3)=a(i,3);
    end;
    individual_Speed_AV2{j}=Speed_AV;
    clear Speed_AV;
end;




%%%%%%%%%%%%%%%%%%%




clear all;
load  
bins_num_speed=50;
bins_num_AV=50;


 speed_min=0;
 speed_max=30;
 
 
min_hist=0.05;
step_hist=0.1;
max_hist=1;

min_hist_speed=0;
step_hist_speed=0.5;
max_hist_speed=5;



AV_min=0;
AV_max=181;
num_of_windows=;
num_individuals=;
num_iterations_random=;
history=15;



diagonal_factor=[5,2.5,2.3,2,1.5];

for k=1:num_individuals
    display(k);
    Speed_AV=individual_Speed_AV{k};
    f=find(Speed_AV(:,1)>speed_min & Speed_AV(:,1)<speed_max);
    temp=Speed_AV(:,1);
    temp2=Speed_AV(:,2);
    temp3=Speed_AV(:,3);
    Speed_AV2(:,1)=temp(f);
    Speed_AV2(:,2)=temp2(f);
    Speed_AV2(:,5)=temp3(f);

    [s1_2 s2_2]=size(Speed_AV2);
    bin_size_speed=(speed_max-speed_min)/bins_num_speed;
    bin_size_AV=(AV_max-AV_min)/bins_num_AV;
    array=zeros(bins_num_AV+1,bins_num_speed+1);
    for i=1+history:s1_2-history
        Speed_AV2(i,3)=round((Speed_AV2(i,1)-speed_min)/bin_size_speed)+1;
        Speed_AV2(i,4)=round((Speed_AV2(i,2)-AV_min)/bin_size_AV)+1;
        array(Speed_AV2(i,4),Speed_AV2(i,3))=array(Speed_AV2(i,4),Speed_AV2(i,3))+1;
        if Speed_AV2(i,4)<diagonal_factor(1,Speed_AV2(i,5))*Speed_AV2(i,3)
            Speed_AV2(i,6)=1;
        else
            Speed_AV2(i,6)=0;
        end;
    end;

    array2=array./sum(sum(array));
    
    
    individuals_array{k}=array2;

    for j=1:1
        array_stage=zeros(bins_num_AV+1,bins_num_speed+1);
        f=find(Speed_AV2(:,5)==j);
        Speed_AV_stage=Speed_AV2(f,3:4);
        [s1_3 s2_3]=size(Speed_AV_stage);
        for i=1+history:s1_3-history
            array_stage(Speed_AV_stage(i,2),Speed_AV_stage(i,1))=array_stage(Speed_AV_stage(i,2),Speed_AV_stage(i,1))+1;
        end;
        array_stage=array_stage./sum(sum(array_stage));
        
        
    end;
    
    individuals_array_stage{k}=array_stage;
    individual_Speed_AV2{k}=Speed_AV2;
    clear Speed_AV2;
    
   
end;
temp_array=zeros(bins_num_AV+1,bins_num_speed+1);
temp_array_stage=zeros(bins_num_AV+1,bins_num_speed+1);
for k=1:num_individuals
    display(k);
    temp_array=temp_array+individuals_array{k};
    temp_array_stage=temp_array_stage+individuals_array_stage{k};
    
end;




for k=1:num_individuals
    ind=1;
    display(k);
    a=individual_Speed_AV2{k};
    for j=1:5
        f=find(a(:,5)==j);
        [s1 s2]=size(f);
        temp=a(f,3:6);
        temp_s=a(f,1);
        win_size=floor(s1/num_of_windows);
        for i=1:win_size:num_of_windows*win_size;
            time_speed(1:win_size,ind)=temp(i:i+win_size-1,1);
            time_speed2=temp_s(i:i+win_size-1,1);
            time_AV(1:win_size,ind)=temp(i:i+win_size-1,2);
            temp_RD=temp(i:i+win_size-1,4);
            f2=find(temp_RD==0);
            f3=find(temp_RD==1);
            [s1_D s2_D]=size(f2);
            [s1_R s2_R]=size(f3);
            fraction_roaming=s1_R/(s1_D+s1_R);
            time_speed3=time_speed2(f3);
            time_RD(1:win_size,ind)=temp(i:i+win_size-1,4);
            time_fraction_roaming(1,ind)=fraction_roaming;
            time_speed4(1,ind)=mean(time_speed3);
            time_speed5(1,ind)=std(time_speed3)/mean(time_speed3);
            ind=ind+1;
            clear temp_RD;
            clear time_speed2;
        end;
    end;
    individual_time_speed{k}=time_speed;
    individual_time_AV{k}=time_AV;
    individual_time_RD{k}=time_RD;
    individual_time_roaming(k,:)=time_fraction_roaming;
    individual_time_average_speed(k,:)=time_speed4;
    individual_time_CV(k,:)=time_speed5;
     
    
    [s1 s2]=size(time_speed);
    for j=1:s2
        array=zeros(bins_num_AV+1,bins_num_speed+1);
        temp_speed=time_speed(:,j);
        temp_AV=time_AV(:,j);
        temp_speed2=temp_speed(find(temp_speed>0));
        temp_AV2=temp_AV(find(temp_AV>0));
        [s1_2 s2_2]=size(temp_speed2);
        for i=1:s1_2
            array(temp_AV2(i,1),temp_speed2(i,1))=array(temp_AV2(i,1),temp_speed2(i,1))+1;
        end;

        array2=array./sum(sum(array));
        states{j}=array2;
    end;
    individual_states{k}=states;
    clear states;
    clear time_speed;
    clear time_AV;
    clear time_RD;
    clear time_fraction_roaming;
end;

for i=1:s2
    average_array_states{i}=zeros(bins_num_AV+1,bins_num_speed+1);
end;
for i=1:num_individuals
    c=individual_states{i};
    for j=1:s2
        average_array_states{j}=average_array_states{j}+c{j};
    end;
end;
for i=1:s2
    average_array_states_norm{i}=average_array_states{i}./sum(sum(average_array_states{i}));
end;



 
 for i=1:num_individuals
    
    c=individual_states{i};
    for j=1:s2
        sum1=0;
        sum2=0;
        d=c{j};
        e=average_array_states_norm{j};
        for m=1:bins_num_AV+1
            for n=1:bins_num_speed+1
                if (d(m,n)>0 && e(m,n)>0)
                mi=(d(m,n)+e(m,n))/2;
                sum1=sum1+(d(m,n)*log2(d(m,n)/mi));
                sum2=sum2+(e(m,n)*log2(e(m,n)/mi));
                end;
            end;
        end;
        distance=0.5*(sum1+sum2);
        individuals_distance_arrays(i,j)=distance;
    end;
    
end;
 
 
 
 


cor=corr(individual_time_roaming,'type','spearman');

for i=1:5
    max_roaming_stage(((i-1)*num_of_windows)+1:i*num_of_windows)=max(max(individual_time_roaming(:,((i-1)*num_of_windows)+1:i*num_of_windows)));
    mean_roaming_stage(((i-1)*num_of_windows)+1:i*num_of_windows)=mean(mean(individual_time_roaming(:,((i-1)*num_of_windows)+1:i*num_of_windows)));
    
    
end;

for i=1:5
    max_speed_stage(((i-1)*num_of_windows)+1:i*num_of_windows)=max(max(individual_time_average_speed(:,((i-1)*num_of_windows)+1:i*num_of_windows)));
    mean_speed_stage(((i-1)*num_of_windows)+1:i*num_of_windows)=mean(mean(individual_time_average_speed(:,((i-1)*num_of_windows)+1:i*num_of_windows)));
    
end;

individual_time_roaming=individual_time_average_speed;
individual_time_roaming(isnan(individual_time_roaming))=0;

for i=1:num_of_windows*5
    [B I]=sort(individual_time_roaming(:,i));
    for j=1:num_individuals
    individual_time_roaming_rank(I(j),i)=j;
    end;
    
end;

for i=1:num_of_windows*5
    for j=1:num_individuals
    individual_time_roaming_rank_norm(j,i)=(individual_time_roaming_rank(j,i)-1)*(1/(num_individuals-1));
    end;
    
    
end;

average_individual_time_roaming=mean(individual_time_roaming);


[s1 s2]=size(individual_time_roaming);

for i=1:num_individuals
    for j=1:s2
        individual_time_roaming_norm(i,j)=individual_time_roaming_rank_norm(i,j);
    end;
end;
individual_time_roaming_original=individual_time_roaming;

individual_time_roaming=individual_time_roaming_norm;
average_individual_time_roaming=mean(individual_time_roaming);


for i=1:num_individuals
    individual_time_roaming_hist(i,:)=hist(individual_time_roaming(i,:),min_hist:step_hist:max_hist);
    individual_time_average_speed_hist(i,:)=hist(individual_time_average_speed(i,:),min_hist_speed:step_hist_speed:max_hist_speed);
    individual_time_roaming_hist2(i,:)=individual_time_roaming_hist(i,:)./sum(individual_time_roaming_hist(i,:));
    individual_time_average_speed_hist2(i,:)=individual_time_average_speed_hist(i,:)./sum(individual_time_average_speed_hist(i,:));
    
    individual_time_average_speed_norm(i,:)=individual_time_average_speed(i,:)./max_speed_stage;
    
    individual_time_average_speed_norm_hist(i,:)=hist(individual_time_average_speed_norm(i,:),min_hist:step_hist:max_hist);
    individual_time_average_speed_norm_hist2(i,:)=individual_time_average_speed_norm_hist(i,:)./sum(individual_time_average_speed_norm_hist(i,:));
end;

[s1_3 s2_3]=size(individual_time_roaming_hist2);

for i=1:num_individuals
    for j=1:s2_3
        individual_time_roaming_hist2_integrated(i,j)=sum(individual_time_roaming_hist2(i,1:j));
        individual_time_average_speed_hist2_integrated(i,j)=sum(individual_time_average_speed_hist2(i,1:j));
    end;
end;

average_individual_time_roaming_hist2_integrated=mean(individual_time_roaming_hist2_integrated);
A=reshape(individual_time_roaming,[1,num_of_windows*5*num_individuals]);

for i=1:num_individuals
    
    distance=sum(individual_time_roaming_hist2_integrated(i,:)-average_individual_time_roaming_hist2_integrated);
    [h p]=kstest2(individual_time_roaming(i,:),A);
    individuals_ks_distance(i,1)=p;
    individuals_ks_distance(i,2)=distance;
    [p h]=ranksum(individual_time_roaming(i,:),A);
    individuals_ks_distance(i,6)=p;
    
    
end;

individuals_ks_distance(:,3)=mafdr(individuals_ks_distance(:,1),'BHFDR', true);
individuals_ks_distance(:,7)=mafdr(individuals_ks_distance(:,6),'BHFDR', true);


% 
for i=1:num_iterations_random
    for j=1:s2
        rand_individual_time_roaming(i,j)=individual_time_roaming(randperm(s1,1),j);
        rand_individual_time_average_speed(i,j)=individual_time_average_speed(randperm(s1,1),j);
        
    end;
    rand_individual_time_roaming_hist(i,:)=hist(rand_individual_time_roaming(i,:),min_hist:step_hist:max_hist);
    rand_individual_time_average_speed_hist(i,:)=hist(rand_individual_time_average_speed(i,:),min_hist_speed:step_hist_speed:max_hist_speed);
    rand_individual_time_roaming_hist2(i,:)=rand_individual_time_roaming_hist(i,:)./sum(rand_individual_time_roaming_hist(i,:));
    rand_individual_time_average_speed_hist2(i,:)=rand_individual_time_average_speed_hist(i,:)./sum(rand_individual_time_average_speed_hist(i,:));
end;
% 
for i=1:num_iterations_random
    for j=1:s2_3
        rand_individual_time_roaming_hist2_integrated(i,j)=sum(rand_individual_time_roaming_hist2(i,1:j));
        rand_individual_time_average_speed_hist2_integrated(i,j)=sum(rand_individual_time_average_speed_hist2(i,1:j));
 
    end;
end;

for i=1:num_iterations_random
    
    distance=sum(rand_individual_time_roaming_hist2_integrated(i,:)-average_individual_time_roaming_hist2_integrated);
    
    rand_individuals_distance(i,1)=distance;
    
    
end;

for i=1:num_individuals
    if(individuals_ks_distance(i,2))>0
        f=find(rand_individuals_distance>individuals_ks_distance(i,2));
        [s1_6 s2_6]=size(f);
        individuals_ks_distance(i,4)=s1_6/num_iterations_random;
    else
        f=find(rand_individuals_distance<individuals_ks_distance(i,2));
        [s1_6 s2_6]=size(f);
        individuals_ks_distance(i,4)=s1_6/num_iterations_random;
    end;

end;

individuals_ks_distance(:,5)=mafdr(individuals_ks_distance(:,4),'BHFDR', true);

for i=1:num_individuals
    h=hist(individual_time_roaming(i,:),0:0.125:1);h=h./sum(h);
    individual_time_roaming_hist2_view(i,:)=h;
end;

[s1_5 s2_5]=size(find(individuals_ks_distance(:,5)<0.05 & abs(individuals_ks_distance(:,2))>1))

cor_time_roaming=corr(individual_time_roaming,'type','spearman');


s=(sum(sum(cor_time_roaming))-num_of_windows)/2;

[f_positive s]=size(find(cor_time_roaming>0.3));
[f_negative s]=size(find(cor_time_roaming<-0.3));
index_corr_time_roaming=(f_positive-f_negative)/(f_positive+f_negative);


hist_distances=hist(individuals_ks_distance(:,2),-3.5:0.5:3.5);
hist_distances_norm=hist_distances./sum(hist_distances);
hist_distances_rand=hist(rand_individuals_distance,-3.5:0.5:3.5);
hist_distances_rand_norm=hist_distances_rand./sum(hist_distances_rand);

std_MAD_IQR_distances(1,1)=std(individuals_ks_distance(:,2));
std_MAD_IQR_distances(2,1)=std(rand_individuals_distance);
std_MAD_IQR_distances(1,2)=mad(individuals_ks_distance(:,2));
std_MAD_IQR_distances(2,2)=mad(rand_individuals_distance);
std_MAD_IQR_distances(1,3)=iqr(individuals_ks_distance(:,2));
std_MAD_IQR_distances(2,3)=iqr(rand_individuals_distance);

average_cor=mean(cor(find(cor<1)));
cronbach=(num_individuals*average_cor)/(1+((num_individuals-1)*average_cor));


for i=1:s1
    [f10 f11]=size(find(individual_time_roaming(i,:)>=0.5));
    bias_bin(i,1)=(f11/s2)/((1-f11/s2));
    for j=1:s2
        if individual_time_roaming(i,j)>=0.5
            individual_time_roaming_bin(i,j)=1;
        else
            individual_time_roaming_bin(i,j)=-1;
        end;
    end;
end;



sum_bin=sum(individual_time_roaming_bin')';

bias_weights=[0.9 0.7 0.5 0.3 0.1 0.1 0.3 0.5 0.7 0.9];
for i=1:s1
    
    bias_bin(i,1)=skewness(individual_time_roaming(i,:));
end;

for i=1:num_iterations_random
    
    rand_bias_bin(i,1)=skewness(rand_individual_time_roaming(i,:));
end;

for i=1:s1
    if(bias_bin(i,1))>0
        f14=find(rand_bias_bin>bias_bin(i,1));
        [f15 f16]=size(f14);
        bias_bin(i,2)=f15/num_iterations_random;
    else
        f14=find(rand_bias_bin<bias_bin(i,1));
        [f15 f16]=size(f14);
        bias_bin(i,2)=f15/num_iterations_random;
    end;

end;





for i=1:num_individuals
    q=find(individual_time_roaming_norm(i,1:num_of_windows*5)>0.5);
    [s1 s2]=size(q);
    count_above_median(i,1)=s2;
end;

for i=1:num_iterations_random
    q=find(rand_individual_time_roaming(i,1:num_of_windows*5)>0.5);
    [s1 s2]=size(q);
    rand_count_above_median(i,1)=s2;
end;


for i=1:num_individuals
    
    middle_hist(i,1)=sum(individual_time_roaming_hist2(i,6:10))/sum(individual_time_roaming_hist2(i,1:5));
    
end;

for i=1:num_iterations_random
    
    rand_middle_hist(i,1)=sum(rand_individual_time_roaming_hist2(i,6:10))/sum(rand_individual_time_roaming_hist2(i,1:5));
    
end;

for i=1:num_individuals
    if(middle_hist(i,1)>1)
        q=find(middle_hist(i,1)<=rand_middle_hist(1:num_iterations_random,1));
        [s1 s2]=size(q);
        temp4(i,1)=s1/num_iterations_random;
    end;
    if(middle_hist(i,1)<=1)
        q=find(middle_hist(i,1)>=rand_middle_hist(1:num_iterations_random,1));
        [s1 s2]=size(q);
        temp4(i,1)=s1/num_iterations_random;
    end;
    
end;
%temp4(:,2)=mafdr(temp4(:,1),'BHFDR', true);

temp4(:,2)=mafdr(temp4(:,1),'BHFDR', true);


