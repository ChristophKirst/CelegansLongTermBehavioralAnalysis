clear all;


corrFiles = dir('*_L1.mat');
sizeFiles = dir('*Size*.mat');
numcorr = length(corrFiles);
num_errors=0;

ind=1;
coord=zeros(1000000,6);
coord2=zeros(1000000,6);
ang_round=100;
time_window=30;
time_window_2=3;
time_bin=10800;

thresh=4;
upper_limit=500;


for n=1:numcorr
ind=1;
coord=zeros(1000000,6);
coord2=zeros(1000000,6);
file_start=n;
file_end=n;

        
for i=file_start:file_end
    load(corrFiles(i).name);
    [x,y]=size(coord_L1);
    display(i);
    
  
    for j=1:x
        coord(ind,1)=i;
        coord(ind,2)=j;
        coord(ind,3)=coord_L1(j,3);
        coord(ind,4)=coord_L1(j,4);
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
    
    for k=(1+time_window_2):(ind-time_window_2-1)
        display(k);
        if coord(k+time_window_2,3)==1 || coord(k,3)==1 || coord(k-time_window_2,3)==1
            display('error');
            num_errors=num_errors+1;
             coord(k,5)=0;
             coord(k,6)=0;
             
        else
             coord(k,12)=sqrt((coord(k+time_window_2,3)-coord(k,3))^2+(coord(k+time_window_2,4)-coord(k,4))^2);
             
             DirVector1=[coord(k,3),coord(k,4)]-[coord(k-time_window_2,3),coord(k-time_window_2,4)];
             DirVector2=[coord(k+time_window_2,3),coord(k+time_window_2,4)]-[coord(k,3),coord(k,4)];
         
            ang = -180/pi*atan2( DirVector2(1)*DirVector1(2)-DirVector2(2)*DirVector1(1) , DirVector2(1)*DirVector1(1)+DirVector2(2)*DirVector1(2));
            
            
            coord(k,13)=ang;
            
            
        end;
    
    end;
    
    
    coord(1:ind-time_window-1,10)=smooth(coord(1:ind-time_window-1,5),360);
    for i=(1+time_window):(ind-time_window-1)
        
        if(coord(i,10)>thresh && coord(i,5)< upper_limit)
                coord(i,11)=1;
        else
                coord(i,11)=0;
        end;
    end;
    

     indx=find(coord(:,11)==1);
     total_indx=find(coord(:,3)>0);

  
    
    size_indx=size(indx);
    size_total_indx=size(total_indx);
    bursts_statistics(n,1)=size_indx(1,1)/size_total_indx(1,1);
    sum_speed=0;
    for i=1:ind-361
        sum_speed=sum_speed+coord(i,5)*coord(i,11);
    end;        
     
    bursts_statistics(n,2)=sum_speed/size_indx(1,1);
    clear sum_speed;
    clear size_indx;
    
    ind7=1
    
    for i=1:time_bin:size_total_indx(1,1)-2
        average_activity(n,ind7)=mean(coord(i:i+time_bin-1,11));
        speed_during_activity_10s=times(coord(i:i+time_bin-1,11),coord(i:i+time_bin-1,5));
        speed_during_activity_1s=times(coord(i:i+time_bin-1,11),coord(i:i+time_bin-1,12));
        if sum(speed_during_activity_10s)>0
            average_speed_10s(n,ind7)=mean(speed_during_activity_10s(find(speed_during_activity_10s>0)));
            average_speed_1s(n,ind7)=mean(speed_during_activity_1s(find(speed_during_activity_1s>0)));
        else
            average_speed_10s(n,ind7)=NaN;
            average_speed_1s(n,ind7)=NaN;
        end;
        ind7=ind7+1;
    end;
    average_activity(n,ind7)=NaN;
    average_speed_10s(n,ind7)=0;
    average_speed_1s(n,ind7)=0;
    
     clear speed_during_activity_10s;
     clear speed_during_activity_1s
     filename2='/Users/shaystern/Desktop/results/';
     filename3=[filename2 corrFiles(n).name];
     filename3=filename3(1:end-4);
     filename4=[filename3 '_bursts'];
     save(filename4,'coord');
end;





clear all;


corrFiles = dir('*_L2.mat');
sizeFiles = dir('*Size*.mat');
numcorr = length(corrFiles);
num_errors=0;

ind=1;
coord=zeros(1000000,6);
coord2=zeros(1000000,6);
ang_round=100;
time_window=30;
time_window_2=3;
time_bin=10800;

thresh=8;
upper_limit=500;


for n=1:numcorr
ind=1;
coord=zeros(1000000,6);
coord2=zeros(1000000,6);
file_start=n;
file_end=n;

        
for i=file_start:file_end
    load(corrFiles(i).name);
    [x,y]=size(coord_L2);
    display(i);
    
  
    for j=1:x
        coord(ind,1)=i;
        coord(ind,2)=j;
        coord(ind,3)=coord_L2(j,3);
        coord(ind,4)=coord_L2(j,4);
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
    
    for k=(1+time_window_2):(ind-time_window_2-1)
        display(k);
        if coord(k+time_window_2,3)==1 || coord(k,3)==1 || coord(k-time_window_2,3)==1
            display('error');
            num_errors=num_errors+1;
             coord(k,5)=0;
             coord(k,6)=0;
             
        else
             coord(k,12)=sqrt((coord(k+time_window_2,3)-coord(k,3))^2+(coord(k+time_window_2,4)-coord(k,4))^2);
             
             DirVector1=[coord(k,3),coord(k,4)]-[coord(k-time_window_2,3),coord(k-time_window_2,4)];
             DirVector2=[coord(k+time_window_2,3),coord(k+time_window_2,4)]-[coord(k,3),coord(k,4)];
         
            ang = -180/pi*atan2( DirVector2(1)*DirVector1(2)-DirVector2(2)*DirVector1(1) , DirVector2(1)*DirVector1(1)+DirVector2(2)*DirVector1(2));
            
            
            coord(k,13)=ang;
            
            
        end;
    
    end;
    
    
    coord(1:ind-time_window-1,10)=smooth(coord(1:ind-time_window-1,5),360);
    for i=(1+time_window):(ind-time_window-1)
        
        if(coord(i,10)>thresh && coord(i,5)< upper_limit)
                coord(i,11)=1;
        else
                coord(i,11)=0;
        end;
    end;
    

    

     indx=find(coord(:,11)==1);
     total_indx=find(coord(:,3)>0);

  
    
    size_indx=size(indx);
    size_total_indx=size(total_indx);
    bursts_statistics(n,1)=size_indx(1,1)/size_total_indx(1,1);
    sum_speed=0;
    for i=1:ind-361
        sum_speed=sum_speed+coord(i,5)*coord(i,11);
    end;        
     
    bursts_statistics(n,2)=sum_speed/size_indx(1,1);
    clear sum_speed;
    clear size_indx;
    
    ind7=1
    
    for i=1:time_bin:size_total_indx(1,1)-2
        average_activity(n,ind7)=mean(coord(i:i+time_bin-1,11));
        speed_during_activity_10s=times(coord(i:i+time_bin-1,11),coord(i:i+time_bin-1,5));
        speed_during_activity_1s=times(coord(i:i+time_bin-1,11),coord(i:i+time_bin-1,12));
        if sum(speed_during_activity_10s)>0
            average_speed_10s(n,ind7)=mean(speed_during_activity_10s(find(speed_during_activity_10s>0)));
            average_speed_1s(n,ind7)=mean(speed_during_activity_1s(find(speed_during_activity_1s>0)));
        else
            average_speed_10s(n,ind7)=NaN;
            average_speed_1s(n,ind7)=NaN;
        end;
        ind7=ind7+1;
    end;
    average_activity(n,ind7)=NaN;
    average_speed_10s(n,ind7)=0;
    average_speed_1s(n,ind7)=0;
    
     clear speed_during_activity_10s;
     clear speed_during_activity_1s
     filename2='/Users/shaystern/Desktop/results/';
     filename3=[filename2 corrFiles(n).name];
     filename3=filename3(1:end-4);
     filename4=[filename3 '_bursts'];
     save(filename4,'coord');
end;




clear all;


corrFiles = dir('*_L3.mat');
sizeFiles = dir('*Size*.mat');
numcorr = length(corrFiles);
num_errors=0;

ind=1;
coord=zeros(1000000,6);
coord2=zeros(1000000,6);
ang_round=100;
time_window=30;
time_window_2=3;
time_bin=10800;

thresh=12;
upper_limit=500;


for n=1:numcorr
ind=1;
coord=zeros(1000000,6);
coord2=zeros(1000000,6);
file_start=n;
file_end=n;
        
for i=file_start:file_end
    load(corrFiles(i).name);
    [x,y]=size(coord_L3);
    display(i);
    
  
    for j=1:x
        coord(ind,1)=i;
        coord(ind,2)=j;
        coord(ind,3)=coord_L3(j,3);
        coord(ind,4)=coord_L3(j,4);
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
    
    for k=(1+time_window_2):(ind-time_window_2-1)
        display(k);
        if coord(k+time_window_2,3)==1 || coord(k,3)==1 || coord(k-time_window_2,3)==1
            display('error');
            num_errors=num_errors+1;
             coord(k,5)=0;
             coord(k,6)=0;
             
        else
             coord(k,12)=sqrt((coord(k+time_window_2,3)-coord(k,3))^2+(coord(k+time_window_2,4)-coord(k,4))^2);
             
             DirVector1=[coord(k,3),coord(k,4)]-[coord(k-time_window_2,3),coord(k-time_window_2,4)];
             DirVector2=[coord(k+time_window_2,3),coord(k+time_window_2,4)]-[coord(k,3),coord(k,4)];
         
            ang = -180/pi*atan2( DirVector2(1)*DirVector1(2)-DirVector2(2)*DirVector1(1) , DirVector2(1)*DirVector1(1)+DirVector2(2)*DirVector1(2));
            
            coord(k,13)=ang;
            
            
        end;
    
    end;
    
    
    coord(1:ind-time_window-1,10)=smooth(coord(1:ind-time_window-1,5),360);
    for i=(1+time_window):(ind-time_window-1)
        
        if(coord(i,10)>thresh && coord(i,5)< upper_limit)
                coord(i,11)=1;
        else
                coord(i,11)=0;
        end;
    end;
    


     indx=find(coord(:,11)==1);
     total_indx=find(coord(:,3)>0);

  
    
    size_indx=size(indx);
    size_total_indx=size(total_indx);
    bursts_statistics(n,1)=size_indx(1,1)/size_total_indx(1,1);
    sum_speed=0;
    for i=1:ind-361
        sum_speed=sum_speed+coord(i,5)*coord(i,11);
    end;        
     
    bursts_statistics(n,2)=sum_speed/size_indx(1,1);
    clear sum_speed;
    clear size_indx;
    
    ind7=1
    
    for i=1:time_bin:size_total_indx(1,1)-2
        average_activity(n,ind7)=mean(coord(i:i+time_bin-1,11));
        speed_during_activity_10s=times(coord(i:i+time_bin-1,11),coord(i:i+time_bin-1,5));
        speed_during_activity_1s=times(coord(i:i+time_bin-1,11),coord(i:i+time_bin-1,12));
        if sum(speed_during_activity_10s)>0
            average_speed_10s(n,ind7)=mean(speed_during_activity_10s(find(speed_during_activity_10s>0)));
            average_speed_1s(n,ind7)=mean(speed_during_activity_1s(find(speed_during_activity_1s>0)));
        else
            average_speed_10s(n,ind7)=NaN;
            average_speed_1s(n,ind7)=NaN;
        end;
        ind7=ind7+1;
    end;
    average_activity(n,ind7)=NaN;
    average_speed_10s(n,ind7)=0;
    average_speed_1s(n,ind7)=0;
    
     clear speed_during_activity_10s;
     clear speed_during_activity_1s
     filename2='/Users/shaystern/Desktop/results/';
     filename3=[filename2 corrFiles(n).name];
     filename3=filename3(1:end-4);
     filename4=[filename3 '_bursts'];
     save(filename4,'coord');
end;




clear all;


corrFiles = dir('*_L4.mat');
sizeFiles = dir('*Size*.mat');
numcorr = length(corrFiles);
num_errors=0;

ind=1;
coord=zeros(1000000,6);
coord2=zeros(1000000,6);
ang_round=100;
time_window=30;
time_window_2=3;
time_bin=10800;

thresh=16;
upper_limit=500;


for n=1:numcorr
ind=1;
coord=zeros(1000000,6);
coord2=zeros(1000000,6);
file_start=n;
file_end=n;

        
for i=file_start:file_end
    load(corrFiles(i).name);
    [x,y]=size(coord_L4);
    display(i);
    
  
    for j=1:x
        coord(ind,1)=i;
        coord(ind,2)=j;
        coord(ind,3)=coord_L4(j,3);
        coord(ind,4)=coord_L4(j,4);
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
    
    for k=(1+time_window_2):(ind-time_window_2-1)
        display(k);
        if coord(k+time_window_2,3)==1 || coord(k,3)==1 || coord(k-time_window_2,3)==1
            display('error');
            num_errors=num_errors+1;
             coord(k,5)=0;
             coord(k,6)=0;
             
        else
             coord(k,12)=sqrt((coord(k+time_window_2,3)-coord(k,3))^2+(coord(k+time_window_2,4)-coord(k,4))^2);
             
             DirVector1=[coord(k,3),coord(k,4)]-[coord(k-time_window_2,3),coord(k-time_window_2,4)];
             DirVector2=[coord(k+time_window_2,3),coord(k+time_window_2,4)]-[coord(k,3),coord(k,4)];
         
            ang = -180/pi*atan2( DirVector2(1)*DirVector1(2)-DirVector2(2)*DirVector1(1) , DirVector2(1)*DirVector1(1)+DirVector2(2)*DirVector1(2));
            
            
            coord(k,13)=ang;
            
            
        end;
    
    end;
    
    
    coord(1:ind-time_window-1,10)=smooth(coord(1:ind-time_window-1,5),360);
    for i=(1+time_window):(ind-time_window-1)
        
        if(coord(i,10)>thresh && coord(i,5)< upper_limit)
                coord(i,11)=1;
        else
                coord(i,11)=0;
        end;
    end;
    

 
     indx=find(coord(:,11)==1);
     total_indx=find(coord(:,3)>0);

  
    
    size_indx=size(indx);
    size_total_indx=size(total_indx);
    bursts_statistics(n,1)=size_indx(1,1)/size_total_indx(1,1);
    sum_speed=0;
    for i=1:ind-361
        sum_speed=sum_speed+coord(i,5)*coord(i,11);
    end;        
     
    bursts_statistics(n,2)=sum_speed/size_indx(1,1);
    clear sum_speed;
    clear size_indx;
    
    ind7=1
    
    for i=1:time_bin:size_total_indx(1,1)-2
        average_activity(n,ind7)=mean(coord(i:i+time_bin-1,11));
        speed_during_activity_10s=times(coord(i:i+time_bin-1,11),coord(i:i+time_bin-1,5));
        speed_during_activity_1s=times(coord(i:i+time_bin-1,11),coord(i:i+time_bin-1,12));
        if sum(speed_during_activity_10s)>0
            average_speed_10s(n,ind7)=mean(speed_during_activity_10s(find(speed_during_activity_10s>0)));
            average_speed_1s(n,ind7)=mean(speed_during_activity_1s(find(speed_during_activity_1s>0)));
        else
            average_speed_10s(n,ind7)=NaN;
            average_speed_1s(n,ind7)=NaN;
        end;
        ind7=ind7+1;
    end;
    average_activity(n,ind7)=NaN;
    average_speed_10s(n,ind7)=0;
    average_speed_1s(n,ind7)=0;
    
     clear speed_during_activity_10s;
     clear speed_during_activity_1s
     filename2='/Users/shaystern/Desktop/results/';
     filename3=[filename2 corrFiles(n).name];
     filename3=filename3(1:end-4);
     filename4=[filename3 '_bursts'];
     save(filename4,'coord');
end;



clear all;


corrFiles = dir('*_A.mat');
sizeFiles = dir('*Size*.mat');
numcorr = length(corrFiles);
num_errors=0;

ind=1;
coord=zeros(1000000,6);
coord2=zeros(1000000,6);
ang_round=100;
time_window=30;
time_window_2=3;
time_bin=10800;

thresh=20;
upper_limit=500;


for n=1:numcorr
ind=1;
coord=zeros(1000000,6);
coord2=zeros(1000000,6);
file_start=n;
file_end=n;

        
for i=file_start:file_end
    load(corrFiles(i).name);
    [x,y]=size(coord_A);
    display(i);
    
  
    for j=1:x
        coord(ind,1)=i;
        coord(ind,2)=j;
        coord(ind,3)=coord_A(j,3);
        coord(ind,4)=coord_A(j,4);
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
    
    for k=(1+time_window_2):(ind-time_window_2-1)
        display(k);
        if coord(k+time_window_2,3)==1 || coord(k,3)==1 || coord(k-time_window_2,3)==1
            display('error');
            num_errors=num_errors+1;
             coord(k,5)=0;
             coord(k,6)=0;
             
        else
             coord(k,12)=sqrt((coord(k+time_window_2,3)-coord(k,3))^2+(coord(k+time_window_2,4)-coord(k,4))^2);
             
             DirVector1=[coord(k,3),coord(k,4)]-[coord(k-time_window_2,3),coord(k-time_window_2,4)];
             DirVector2=[coord(k+time_window_2,3),coord(k+time_window_2,4)]-[coord(k,3),coord(k,4)];
         
            ang = -180/pi*atan2( DirVector2(1)*DirVector1(2)-DirVector2(2)*DirVector1(1) , DirVector2(1)*DirVector1(1)+DirVector2(2)*DirVector1(2));
            
            
            coord(k,13)=ang;
            
            
        end;
    
    end;
    
    
    coord(1:ind-time_window-1,10)=smooth(coord(1:ind-time_window-1,5),360);
    for i=(1+time_window):(ind-time_window-1)
        
        if(coord(i,10)>thresh && coord(i,5)< upper_limit)
                coord(i,11)=1;
        else
                coord(i,11)=0;
        end;
    end;
    
    
    figure;
    subplot(4,1,1);
    
    plot(1/10800:1/10800:(ind-time_window-1)/10800,coord(1:ind-time_window-1,5));
    ylim([0 100]);
    
    subplot(4,1,2);
   
    plot(1/10800:1/10800:(ind-time_window-1)/10800,coord(1:ind-time_window-1,10));
    ylim([0 100]);
    
    subplot(4,1,3);
    
    plot(1/10800:1/10800:(ind-time_window-1)/10800,coord(1:ind-time_window-1,11));
    ylim([0 2]);
    
    subplot(4,1,4);
    
    plot(1/10800:1/10800:(ind-time_window-1)/10800,coord(1:ind-time_window-1,6));
    ylim([-180 180]);
    
%     
     indx=find(coord(:,11)==1);
     total_indx=find(coord(:,3)>0);

  
    
    size_indx=size(indx);
    size_total_indx=size(total_indx);
    bursts_statistics(n,1)=size_indx(1,1)/size_total_indx(1,1);
    sum_speed=0;
    for i=1:ind-361
        sum_speed=sum_speed+coord(i,5)*coord(i,11);
    end;        
     
    bursts_statistics(n,2)=sum_speed/size_indx(1,1);
    clear sum_speed;
    clear size_indx;
    
    ind7=1
    
    for i=1:time_bin:size_total_indx(1,1)-2
        average_activity(n,ind7)=mean(coord(i:i+time_bin-1,11));
        speed_during_activity_10s=times(coord(i:i+time_bin-1,11),coord(i:i+time_bin-1,5));
        speed_during_activity_1s=times(coord(i:i+time_bin-1,11),coord(i:i+time_bin-1,12));
        if sum(speed_during_activity_10s)>0
            average_speed_10s(n,ind7)=mean(speed_during_activity_10s(find(speed_during_activity_10s>0)));
            average_speed_1s(n,ind7)=mean(speed_during_activity_1s(find(speed_during_activity_1s>0)));
        else
            average_speed_10s(n,ind7)=NaN;
            average_speed_1s(n,ind7)=NaN;
        end;
        ind7=ind7+1;
    end;
    average_activity(n,ind7)=NaN;
    average_speed_10s(n,ind7)=0;
    average_speed_1s(n,ind7)=0;
    
     clear speed_during_activity_10s;
     clear speed_during_activity_1s;
     filename2='/Users/shaystern/Desktop/results/';
     filename3=[filename2 corrFiles(n).name];
     filename3=filename3(1:end-4);
     filename4=[filename3 '_bursts'];
     save(filename4,'coord');
end;






