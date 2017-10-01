clear all;

mkdir('E:\Results\CAM207A1');

aviFiles = dir('*CAM207*.avi');
file_start=1;
file_end=272;
numavi = length(aviFiles);
shuttleVideo2 = VideoReader(aviFiles(file_start).name);
[corrd rect] = imcrop(read(shuttleVideo2,1));


Mask=uint8(roipoly(imcrop(read(shuttleVideo2,1),[rect(1),rect(2),rect(3),rect(4)])));
result = input('How many deleted objects?');

for a=1:result
Mask_2=roipoly(imcrop(read(shuttleVideo2,1),[rect(1),rect(2),rect(3),rect(4)]));
Mask_2=uint8(times(Mask_2,-1)+1);

Mask=uint8(times(Mask_2,Mask));
end;


Mask2(:,:,1)=Mask;
Mask2(:,:,2)=Mask;
Mask2(:,:,3)=Mask;
frame_window=30;
frame_window2=1;






for d=file_start:file_end+8
    shuttleVideo = VideoReader(aviFiles(d).name);
    Mean_zeros(1:4000,1:4000,1:3)=uint8(0);
    Mean_mat2=imcrop(Mean_zeros,[rect(1),rect(2),rect(3),rect(4)]);
    numf=shuttleVideo.NumberOfFrames;
    numf2=round(numf/8-1);

    for p=1:8
    
    Mean_mat{p}=imcrop(read(shuttleVideo,p*numf2),[rect(1),rect(2),rect(3),rect(4)]);
    Mean_mat2=Mean_mat2+((Mean_mat{p}-Mean_mat2)./p);
    Mean_mat2=times(Mask2,Mean_mat2);
    end;

    
    Norm_array{d}=Mean_mat2;
    display(d);
    
    
end;

tic;
for n = file_start:file_end
if n<100
    thresh=0.37;
    WormMin=30;
    WormMax=1000;
else
    thresh=0.35;
    WormMin=100;
    WormMax=1500;
end;
    
ind=1;
shuttleVideo = VideoReader(aviFiles(n).name);
display(aviFiles(n).name);


numfiles=shuttleVideo.NumberOfFrames;
m=mod(numfiles,frame_window);
numfiles2=(numfiles-m)/frame_window;
mydata = cell(1, numfiles);


Norm=Norm_array{n+8};





   
for i = 1:numfiles2
    
    flag=0;
    display(aviFiles(n).name);
    
    img = imcrop(read(shuttleVideo,i*frame_window-frame_window+1),[rect(1),rect(2),rect(3),rect(4)]);
    display(i);
    img2=times(img,Mask2)+100-Norm;
    
    J=imcrop(img2,[1,1,4096,2160]);
    
    I3=imcrop(img,[1,1,4096,2160]);
  
  
    Q=imfill(((im2bw(J, thresh))-1).*-1,'holes');
    CC = bwconncomp(Q);
    S = regionprops(CC,'Centroid');
    numPixels = cellfun(@numel,CC.PixelIdxList);
    [biggest,idx] = max(numPixels);
    if CC.NumObjects==0 || biggest<WormMin || biggest>WormMax 
        display('error');
        S2(1)=350;
        S2(2)=350;
        Norm2=imcrop(Norm,[S2(1)-350,S2(2)-350,700,700]);
        flag=1;
            
    else
    S = regionprops(CC,'Centroid');
    numPixels = cellfun(@numel,CC.PixelIdxList);
    [biggest,idx] = max(numPixels);
    S2=S(idx).Centroid;
    
    if S2(1)<350
        left_upper_corner(1)=rect(1);
    else
        left_upper_corner(1)=rect(1)+S2(1)-350;
    end;
    
    if S2(2)<350
        left_upper_corner(2)=rect(2);
    else
        left_upper_corner(2)=rect(2)+S2(2)-350;
    end;
  
    Norm2=imcrop(Norm,[S2(1)-350,S2(2)-350,700,700]);
  end;
    st=i*frame_window-frame_window+1;
    en=i*frame_window;
    
    for q=st:frame_window2:en
        
    img = imcrop(times(imcrop(read(shuttleVideo,q),[rect(1),rect(2),rect(3),rect(4)]),Mask2),[S2(1)-350,S2(2)-350,700,700]);
    
    img2=img+100-Norm2;
    
    J=imcrop(img2,[1,1,4096,2160]);
    
    I3=imcrop(img,[1,1,4096,2160]);
  
  
    Q=imfill(((im2bw(J, thresh))-1).*-1,'holes');
    CC = bwconncomp(Q);
    S = regionprops(CC,'Centroid');
    numPixels = cellfun(@numel,CC.PixelIdxList);
    [biggest,idx] = max(numPixels);
    if CC.NumObjects==0 || biggest<WormMin || biggest>WormMax || flag==1
   
        display('error');
        x_y_coor(ind,1:2)=[1,1];
        big(ind)=1;
        ind=ind+1;

    else
    S = regionprops(CC,'Centroid');
    numPixels = cellfun(@numel,CC.PixelIdxList);
    [biggest,idx] = max(numPixels);
    S2_2=S(idx).Centroid;
    mydata{ind}=imcrop(img2,[S2_2(1)-75,S2_2(2)-75,150,150]); 
    
    
    
  big(ind)=biggest;
  
    
    
    x_y_coor(ind,1:2)= [left_upper_corner(1)+S2_2(1),left_upper_corner(2)+S2_2(2)];
    ind=ind+1;

 
  display(q);

   
    end;
    end;
    
    
end;
   
   
    st=(i*frame_window)+1;
    en=numfiles;
    
    
    for q=st:frame_window2:en
    img = imcrop(times(imcrop(read(shuttleVideo,q),[rect(1),rect(2),rect(3),rect(4)]),Mask2),[S2(1)-350,S2(2)-350,700,700]);
    
    Norm2=imcrop(Norm,[S2(1)-350,S2(2)-350,700,700]);
    img2=img+100-Norm2;
        
   
    J=imcrop(img2,[1,1,4096,2160]);
    
    I3=imcrop(img,[1,1,4096,2160]);
  
  
    Q=imfill(((im2bw(J, thresh))-1).*-1,'holes');
    CC = bwconncomp(Q);
    S = regionprops(CC,'Centroid');
    numPixels = cellfun(@numel,CC.PixelIdxList);
    [biggest,idx] = max(numPixels);
    if CC.NumObjects==0 || biggest<WormMin || biggest>WormMax || flag==1
    
        display('error');
        x_y_coor(ind,1:2)=[1,1];
        big(ind)=1;
        ind=ind+1;
    else
    S = regionprops(CC,'Centroid');
    numPixels = cellfun(@numel,CC.PixelIdxList);
    [biggest,idx] = max(numPixels);
    S2_2=S(idx).Centroid;
    
    
    mydata{ind}=imcrop(img2,[S2_2(1)-75,S2_2(2)-75,150,150]);  
    
    
    
    big(ind)=biggest;
  
    
    x_y_coor(ind,1:2)= [left_upper_corner(1)+S2_2(1),left_upper_corner(2)+S2_2(2)];
    %x_y_coor(ind,1:2)= [rect(1)+S2(1)-350+S2_2(1),rect(2)+S2(2)-350+S2_2(2)];
    ind=ind+1;
   
 
  display(q);
    end;
    end;


filename2='E:\Results\CAM207A1\corrdCAM207A1';
filename3=[filename2 aviFiles(n).name];
filename3=filename3(1:end-4);
save(filename3,'x_y_coor');

Trace=read(shuttleVideo,1);
for i=1:ind-1
    Trace(round(x_y_coor(i,2)),round(x_y_coor(i,1)),1)=255;
end;

%filename2=sprintf('/Users/shaystern/Desktop/short/newillumination/Trajectory%d',n);
filename2='E:\Results\CAM207A1\TrajectoryCAM207A1';
filename3=[filename2 aviFiles(n).name];
filename3=filename3(1:end-4);
save(filename3,'Trace');

%filename=sprintf('/Users/shaystern/Desktop/short/newillumination/short%d.avi',n);
filename2='E:\Results\CAM207A1\shortCAM207A1';
filename3=[filename2 aviFiles(n).name];
filename3=filename3(1:end-4);
save(filename3,'mydata');

filename2='E:\Results\CAM207A1\WormSizeCAM207A1';
filename3=[filename2 aviFiles(n).name];
filename3=filename3(1:end-4);
save(filename3,'big');


clear big;
clear mydata;
clear x_y_coor;

end;

%close(outputVideo);

%close(outputVideo_bw);
toc;


















