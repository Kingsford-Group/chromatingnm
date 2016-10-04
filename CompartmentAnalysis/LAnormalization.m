% Lieberman-Aiden et al normalization method

res = 5000;
resname = '5kb';
sumvec = zeros(50000,1);
countvec = zeros(length(sumvec),1);

for chrnum = 23:-1:1
    if chrnum == 23
        chrstr = 'X';
    else
        chrstr = num2str(chrnum);
    end
    %     filename = strcat('data/Hi-C/chr',num2str(chrnum),'_10kb.RAWobserved');
    filename = strcat('data/',resname,'_resolution_intrachromosomal/chr',chrstr,'/MAPQGE30/chr',chrstr,'_',resname,'.RAWobserved');
    
    currdata = load(filename);
    
    currdata(:,1:2) = currdata(:,1:2)/res + 1;
    
    maxrow = max(currdata(:,1));
    maxcol = max(currdata(:,2));
    if maxrow ~= maxcol
        matsize = max(maxrow, maxcol);
        currdata = [currdata; matsize matsize 0];
    end
    
    currdata = spconvert(currdata);
    currdata = currdata + currdata' - diag(diag(currdata));
    
    currdata = full(currdata);
    currdata(isnan(currdata)) = 0;
    
    %     figure; imagesc(currdata)
    
    newfilename = strcat('data/Hi-C/RaoChr', chrstr,'_res',resname,'RAW.txt');
    %     if exist(newfilename, 'file') ~= 2
    dlmwrite(newfilename,currdata,'\t')
    %     end
    
    for i=1:length(currdata)
        for j=i:length(currdata)
            sumvec(j-i+1) = sumvec(j-i+1) + currdata(i,j);
            countvec(j-i+1) = countvec(j-i+1) + 1;
        end
    end
    
    fprintf(strcat('chr ',chrstr,'_',resname,' finished\n'));
    
    clear currdata
end
distavgs = sumvec./countvec;
dlmwrite(strcat('data/Hi-C/LAnormvals_',resname,'.txt'),distavgs,'\t');

fprintf('\n distance averages calcualated\n')

for chrnum = 23:-1:1
    if chrnum == 23
        chrstr = 'X';
    else
        chrstr = num2str(chrnum);
    end
    rawfilename = strcat('data/Hi-C/RaoChr', chrstr,'_res',resname,'RAW.txt');
    currdata = load(rawfilename);
    normdata = currdata;
    
    for i=1:length(currdata)
        for j=1:length(currdata)
            normdata(i,j) = normdata(i,j)/distavgs(abs(j-i)+1);
        end
    end
    
    normdata = normdata + normdata' - diag(diag(normdata));
    
    newfilename = strcat('data/Hi-C/RaoChr', chrstr,resname,'LAnorm.txt');
    %     if exist(newfilename, 'file') ~= 2
    dlmwrite(newfilename,normdata,'\t')
    %     end
    
    fprintf(strcat('chr ',chrstr,'_',resname,' normalized\n'));
    
    clear normdata currdata
end
