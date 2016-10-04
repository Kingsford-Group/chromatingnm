% script to find and save compartment boundaries for all chromosomes

resname = '5kb';

for chrnum = 23:-1:1
    if chrnum == 23
        chrstr = 'X';
    else
        chrstr = num2str(chrnum);
    end
    
    sprintf('Finding compartments for chr%s...\n',chrstr)
    
    filename = strcat('data/Hi-C/RaoChr', chrstr,resname,'LAnorm.txt');
    
    if exist(filename, 'file') == 0
        LAnormalization;
    end
    
    [compartments, ~,~] = HiCcompartments(filename);
    
    newfilename = strcat('data/CompartmentBoundaries/RaoChr', chrstr,'_',resname,'_compartments.txt');
    dlmwrite(newfilename, compartments, '\t');
    
end