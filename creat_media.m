clear;close all;clc;
nx=500;   
fid=fopen('media.dat','w');
fmt_media='%15.7f %15.7f';
fmt_media=[fmt_media '\n'];
for i=1:nx
    p(i)=2700+500*(1-2*rand(1));
    c(i)=5000+1000*(1-2*rand(1));
     fprintf(fid,fmt_media ,p(i),c(i));
end
fclose(fid)