clear all;

corrFiles_L1 = dir('*_L1*');
corrFiles_L2 = dir('*_L2*');
corrFiles_L3 = dir('*_L3*');
corrFiles_L4 = dir('*_L4*');
corrFiles_A = dir('*_A*');

numcorr = length(corrFiles_L1);

for n=1:numcorr
    display(n);
    load(corrFiles_L1(n).name);
    total_indx=find(coord(:,3)>0);
    size_total_indx_L1=size(total_indx);
    coord_L1=coord(4:size_total_indx_L1(1,1)-3,12);
    coord_L1(:,2)=abs(coord(4:size_total_indx_L1(1,1)-3,13));
    coord_L1(:,3)=1;
    
    load(corrFiles_L2(n).name);
    total_indx=find(coord(:,3)>0);
    size_total_indx_L2=size(total_indx);
    coord_L2=coord(4:size_total_indx_L2(1,1)-3,12);
    coord_L2(:,2)=abs(coord(4:size_total_indx_L2(1,1)-3,13));
    coord_L2(:,3)=2;
    
    load(corrFiles_L3(n).name);
    total_indx=find(coord(:,3)>0);
    size_total_indx_L3=size(total_indx);
    coord_L3=coord(4:size_total_indx_L3(1,1)-3,12);
    coord_L3(:,2)=abs(coord(4:size_total_indx_L3(1,1)-3,13));
    coord_L3(:,3)=3;
    
    load(corrFiles_L4(n).name);
    total_indx=find(coord(:,3)>0);
    size_total_indx_L4=size(total_indx);
    coord_L4=coord(4:size_total_indx_L4(1,1)-3,12);
    coord_L4(:,2)=abs(coord(4:size_total_indx_L4(1,1)-3,13));
    coord_L4(:,3)=4;
    
    load(corrFiles_A(n).name);
    total_indx=find(coord(:,3)>0);
    size_total_indx_A=size(total_indx);
    coord_A=coord(4:size_total_indx_A(1,1)-3,12);
    coord_A(:,2)=abs(coord(4:size_total_indx_A(1,1)-3,13));
    coord_A(:,3)=5;
    
    coord_all=[coord_L1;coord_L2;coord_L3;coord_L4;coord_A];
    individual_speed_AV{n}=coord_all;
    
end;