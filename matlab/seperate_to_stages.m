clear all;


corrFiles = dir('*corrd*.mat');
sizeFiles = dir('*Size*.mat');
numcorr = length(corrFiles);
num_errors=0;

ind=1;
coord=zeros(1000000,6);
coord2=zeros(1000000,6);
ang_round=100;
time_window=30;
bin=1000;


file_start=5;
file_end=numcorr;

        
for i=file_start:file_end
    load(corrFiles(i).name);
    [x,y]=size(x_y_coor);
    display(i);
    
  
    for j=1:x
        coord(ind,1)=i;
        coord(ind,2)=j;
        coord(ind,3)=x_y_coor(j,1);
        coord(ind,4)=x_y_coor(j,2);
        ind=ind+1;
    end;
end;


    for k=(1+time_window):(ind-time_window-1)
        display(k);
        if coord(k+time_window,3)==1 || coord(k,3)==1 || coord(k-time_window,3)==1
            display('error');
            num_errors=num_errors+1;
             coord(k,5)=0;
             coord(k,6)=0;
             
        else
             coord(k,5)=sqrt((coord(k+time_window,3)-coord(k,3))^2+(coord(k+time_window,4)-coord(k,4))^2);
             
             DirVector1=[coord(k,3),coord(k,4)]-[coord(k-time_window,3),coord(k-time_window,4)];
             DirVector2=[coord(k+time_window,3),coord(k+time_window,4)]-[coord(k,3),coord(k,4)];
         
            ang = -180/pi*atan2( DirVector2(1)*DirVector1(2)-DirVector2(2)*DirVector1(1) , DirVector2(1)*DirVector1(1)+DirVector2(2)*DirVector1(2));
            
            coord(k,6)=ang;
            coord(k,7)=round(coord(k,3)*1)/1;
            coord(k,8)=round(coord(k,4)*1)/1;
            coord(k,9)=abs(ang);
            
        end;
    
    end;

    
    
    figure;
    
    tic;smooth_speed=smooth(coord(1:ind-time_window-1,5),300);toc;
    plot(1/10800:1/10800:(ind-time_window-1)/10800,smooth_speed);
    ylim([0 200]);
    
  
    set(gca,'XMinorTick','on')
    grid on
    grid minor
    
    
    min_matrix_speed=1;
    max_matrix_speed=200;
    min_matrix_ang=-181;
    max_matrix_ang=181;
    numofcells=25;
    
    matrix_time_window=10800;
    
    ind5=1
   for a=time_window+1:matrix_time_window:(ind-matrix_time_window)
        matrix=zeros(numofcells,numofcells);
        ind6=1;
        for b=a:a+matrix_time_window-1
            if coord(b,5)<max_matrix_speed && coord(b,5)>min_matrix_speed
        cell1=floor((coord(b,5)-min_matrix_speed)/((max_matrix_speed-min_matrix_speed)/numofcells))+1;
        coord(b,10)=cell1;
        cell2=floor((coord(b,6)-min_matrix_ang)/((max_matrix_ang-min_matrix_ang)/numofcells))+1;
        coord(b,11)=cell2;
        matrix(cell1,cell2)=matrix(cell1,cell2)+1;
        ave_speed(ind6)=coord(b,5);
        ind6=ind6+1;
            end;    
        end;
        array_matrix{ind5}=matrix./sum(sum(matrix));
        counts_matrix{ind5}=sum(sum(matrix));
        if ind6==1
            array_ave_speed(ind5)=50;
        else
        array_ave_speed(ind5)=mean(ave_speed);
        end;
        ind5=ind5+1;
        clear ave_speed;
    end;
 
 flag=0;   
 for i=1:ind5-1
    if counts_matrix{i}>1000 && flag==0
        start_traj=i;
        flag=1;
    end;
    
 end;
 flag2=0;
 for i=start_traj+1:ind5-1
     if counts_matrix{i}<1000 && counts_matrix{i-1}<1000 && counts_matrix{i-2}<1000 && flag2==0 && counts_matrix{i-3}<1000 && counts_matrix{i-4}<1000
        end_traj=i-5;
        flag2=1;
     end;
 end;
 if flag2==0
     end_traj=i;
 end;
 

 

 
 ind=1;
 W_size=zeros(k,3);
 for i=file_start:file_end
    load(sizeFiles(i).name);
    [x,y]=size(big);
    display(i);
    
    
    for j=1:y
        
        W_size(ind,1)=i;
        W_size(ind,2)=j;
        if(big(j)>1)
        W_size(ind,3)=big(j);
        else
        W_size(ind,3)=0;
        end;
        
        ind=ind+1;
    end;
end;



 
for k=start_traj:end_traj
    display(k);
    for l=start_traj:end_traj
        
        array1=array_matrix{k};

        [size2 s]=size(array1);
       
        array2=array_matrix{l};

        [size2 s2]=size(array2);

        
        
                mat1_2=array1;
                mat2_2=array2;
                mat1=mat1_2;
                mat2=mat2_2;

                 mat1(find(mat1_2==0))=10^-12;
                mat2(find(mat2_2==0))=10^-12;
                for i=1:s
                    for j=1:s
                        mi=(mat1(i,j)+mat2(i,j))/2;    
                        mat3(i,j)=mat1(i,j)*log2(mat1(i,j)/mi);
                        mat4(i,j)=mat2(i,j)*log2(mat2(i,j)/mi);
                        sum_mat3=sum(mat3(:));
                        sum_mat4=sum(mat4(:));
                        D=0.5*(sum_mat3+sum_mat4);


                    end;
                end;
                if counts_matrix{k}>1000 && counts_matrix{l}>1000
                    similarity_mat(k,l)=D;
                else
                    similarity_mat(k,l)=0;
                end;
                
                
    end;
end;    

for i=start_traj:end_traj-1
    similarity_mat2=similarity_mat(i,i+1:end_traj);
    stat_corr(i)=mean(similarity_mat2(find(similarity_mat(i,i+1:end_traj)>0)));
end;

% 
 L1_s=start_traj;
     L1_e=input('L1 ends: ');
     L2_s=input('L2 starts: ');
     L2_e=input('L2 ends: ');
     L3_s=input('L3 starts: ');
     L3_e=input('L3 ends: ');
     L4_s=input('L4 starts: ');
     L4_e=input('L4 ends: ');
     A_s=input('A starts: ');
     A_e=A_s+16;
     
     filename2='/Users/shaystern/Desktop/results/';
     filename3=[filename2 corrFiles(i).name];
     filename3=filename3(1:end-4);
      filename4=[filename3 '_L1'];
      coord_L1=coord(L1_s*10800:L1_e*10800,:);
      save(filename4,'coord_L1');
%     
      filename4=[filename3 '_L2'];
      coord_L2=coord(L2_s*10800:L2_e*10800,:);
      save(filename4,'coord_L2');
%      
      filename4=[filename3 '_L3'];
      coord_L3=coord(L3_s*10800:L3_e*10800,:);
      save(filename4,'coord_L3');
%      
      filename4=[filename3 '_L4'];
      coord_L4=coord(L4_s*10800:L4_e*10800,:);
      save(filename4,'coord_L4');
%      
      filename4=[filename3 '_A'];
     coord_A=coord(A_s*10800:A_e*10800,:);
      save(filename4,'coord_A');
     
     
     
     coord_L1_2=coord_L1(:,5);
     coord_L2_2=coord_L2(:,5);
     coord_L3_2=coord_L3(:,5);
     coord_L4_2=coord_L4(:,5);
     coord_A_2=coord_A(:,5);
     
     ave_speed_stages=[mean(coord_L1_2(find(coord_L1_2>1))) mean(coord_L2_2(find(coord_L2_2>1))) mean(coord_L3_2(find(coord_L3_2>1))) mean(coord_L4_2(find(coord_L4_2>1))) mean(coord_A_2(find(coord_A_2>1)))]; 
     max_speed_stages=[max(coord_L1_2(find(coord_L1_2<50))) max(coord_L2_2(find(coord_L2_2>1))) max(coord_L3_2(find(coord_L3_2>1))) max(coord_L4_2(find(coord_L4_2>1))) max(coord_A_2(find(coord_A_2<350)))]; 
     all_speed_stages=[ave_speed_stages max_speed_stages];
     
     ang_L1_2=coord_L1(:,9);
     ang_L2_2=coord_L2(:,9);
     ang_L3_2=coord_L3(:,9);
     ang_L4_2=coord_L4(:,9);
     ang_A_2=coord_A(:,9);
     
     ave_ang_stages=[mean(ang_L1_2(find(ang_L1_2>1))) mean(ang_L2_2(find(ang_L2_2>1))) mean(ang_L3_2(find(ang_L3_2>1))) mean(ang_L4_2(find(ang_L4_2>1))) mean(ang_A_2(find(ang_A_2>1)))]; 
     
     
     
     figure;
    
     subplot(5,1,1);
     [s1,s2]=size(coord_L1);
     plot(1:s1,coord_L1(1:s1,5));
     title('L1');
     ylim([1 50]);
     
     subplot(5,1,2);
     [s1,s2]=size(coord_L2);
     plot(1:s1,coord_L2(1:s1,5));
     title('L2');
     ylim([1 100]);
     
     subplot(5,1,3);
     [s1,s2]=size(coord_L3);
     plot(1:s1,coord_L3(1:s1,5));
     title('L3');
     ylim([1 150]);
     
     subplot(5,1,4);
     [s1,s2]=size(coord_L4);
     plot(1:s1,coord_L4(1:s1,5));
     title('L4');
     ylim([1 200]);
     
     subplot(5,1,5);
     [s1,s2]=size(coord_A);
     plot(1:s1,coord_A(1:s1,5));
     title('A');
     ylim([1 300]);
     
     
    min_matrix_speed=2;
    max_matrix_speed=225;
    min_matrix_ang=-181;
    max_matrix_ang=181;
    numofcells_x=75;
    numofcells_y=25;
    matrix_time_window=(L1_e-L1_s)*10800;
    
    ind5=1
   for a=1:matrix_time_window:matrix_time_window
        matrix_L1=zeros(numofcells_x,numofcells_y);
        ind6=1;
        for b=a:a+matrix_time_window-1
            if coord_L1(b,5)<max_matrix_speed && coord_L1(b,5)>min_matrix_speed
        cell1=floor((coord_L1(b,5)-min_matrix_speed)/((max_matrix_speed-min_matrix_speed)/numofcells_x))+1;
        cell2=floor((coord_L1(b,9)-min_matrix_ang)/((max_matrix_ang-min_matrix_ang)/numofcells_y))+1;
        matrix_L1(cell1,cell2)=matrix_L1(cell1,cell2)+1;

            end;    
        end;
        array_matrix_stages{1}=matrix_L1./sum(sum(matrix_L1));

    end;
    
    min_matrix_speed=2;
    max_matrix_speed=225;
    min_matrix_ang=-181;
    max_matrix_ang=181;
    numofcells_x=75;
    numofcells_y=25;
    matrix_time_window=(L2_e-L2_s)*10800;
    
    ind5=1
   for a=1:matrix_time_window:matrix_time_window
        matrix_L2=zeros(numofcells_x,numofcells_y);
        ind6=1;
        for b=a:a+matrix_time_window-1
            if coord_L2(b,5)<max_matrix_speed && coord_L2(b,5)>min_matrix_speed
        cell1=floor((coord_L2(b,5)-min_matrix_speed)/((max_matrix_speed-min_matrix_speed)/numofcells_x))+1;
        cell2=floor((coord_L2(b,9)-min_matrix_ang)/((max_matrix_ang-min_matrix_ang)/numofcells_y))+1;
        matrix_L2(cell1,cell2)=matrix_L2(cell1,cell2)+1;

            end;    
        end;
        array_matrix_stages{2}=matrix_L2./sum(sum(matrix_L2));

    end;
    
    min_matrix_speed=2;
    max_matrix_speed=225;
    min_matrix_ang=-181;
    max_matrix_ang=181;
    numofcells_x=75;
    numofcells_y=25;
    matrix_time_window=(L3_e-L3_s)*10800;
    
    ind5=1
   for a=1:matrix_time_window:matrix_time_window
        matrix_L3=zeros(numofcells_x,numofcells_y);
        ind6=1;
        for b=a:a+matrix_time_window-1
            if coord_L3(b,5)<max_matrix_speed && coord_L3(b,5)>min_matrix_speed
        cell1=floor((coord_L3(b,5)-min_matrix_speed)/((max_matrix_speed-min_matrix_speed)/numofcells_x))+1;
        cell2=floor((coord_L3(b,9)-min_matrix_ang)/((max_matrix_ang-min_matrix_ang)/numofcells_y))+1;
        matrix_L3(cell1,cell2)=matrix_L3(cell1,cell2)+1;

            end;    
        end;
        array_matrix_stages{3}=matrix_L3./sum(sum(matrix_L3));

    end;
    
    min_matrix_speed=2;
    max_matrix_speed=225;
    min_matrix_ang=-181;
    max_matrix_ang=181;
    numofcells_x=75;
    numofcells_y=25;
    matrix_time_window=(L4_e-L4_s)*10800;
    
    ind5=1
   for a=1:matrix_time_window:matrix_time_window
        matrix_L4=zeros(numofcells_x,numofcells_y);
        ind6=1;
        for b=a:a+matrix_time_window-1
            if coord_L4(b,5)<max_matrix_speed && coord_L4(b,5)>min_matrix_speed
        cell1=floor((coord_L4(b,5)-min_matrix_speed)/((max_matrix_speed-min_matrix_speed)/numofcells_x))+1;
        cell2=floor((coord_L4(b,9)-min_matrix_ang)/((max_matrix_ang-min_matrix_ang)/numofcells_y))+1;
        matrix_L4(cell1,cell2)=matrix_L4(cell1,cell2)+1;

            end;    
        end;
        array_matrix_stages{4}=matrix_L4./sum(sum(matrix_L4));

    end;
    
    min_matrix_speed=2;
    max_matrix_speed=225;
    min_matrix_ang=-181;
    max_matrix_ang=181;
    numofcells_x=75;
    numofcells_y=25;
    matrix_time_window=(A_e-A_s)*10800;
    
    ind5=1;
   for a=1:matrix_time_window:matrix_time_window
        matrix_A=zeros(numofcells_x,numofcells_y);
        ind6=1;
        for b=a:a+matrix_time_window-1
            if coord_A(b,5)<max_matrix_speed && coord_A(b,5)>min_matrix_speed
        cell1=floor((coord_A(b,5)-min_matrix_speed)/((max_matrix_speed-min_matrix_speed)/numofcells_x))+1;
        cell2=floor((coord_A(b,9)-min_matrix_ang)/((max_matrix_ang-min_matrix_ang)/numofcells_y))+1;
        matrix_A(cell1,cell2)=matrix_A(cell1,cell2)+1;

            end;    
        end;
        array_matrix_stages{5}=matrix_A./sum(sum(matrix_A));

    end;
    
    
    filename4=[filename3 '_stages_matrix_array'];
     
     save(filename4,'array_matrix_stages');
