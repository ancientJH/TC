
for h=1:7  

     fprintf('Mesh : %g\n',h);
     j_count=mod(h-1,10);

     meshType = strcat('Edgecrack6-h',char(48+j_count),'.inp');
     fid=fopen(meshType, 'r');
     [p,t]=readinp(fid);
     filename=strcat('Edgecrack6pt-h',char(48+j_count),'.mat');
     save(filename,'p','t');
     
end