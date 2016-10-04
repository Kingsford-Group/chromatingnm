function [ chrom, starts, ends, name, scores, strand, signal, post, qvalue ] = readBroadPk( filename, split2mat, outputpath )
%READBROADPK parse broad peak files and save data to a mat file
% filename   - the path to the broad peak file
% split2mat  - save peaks of different chromosomes to
%              different mat files
% outputpath - the path where the files are saved

data = importdata(filename);
raw = data.textdata;

chrom_raw  = raw(:, 1);
starts_raw = raw(:, 2);
ends_raw   = raw(:, 3);
name_raw   = raw(:, 4);
score_raw  = raw(:, 5);
strand_raw = raw(:, 6);

reads_raw = data.data(:, 1);
post_raw = data.data(:, 2);
qvalue_raw = data.data(:, 3);

% convert cell array to list
starts_raw = cellfun(@str2num, starts_raw);
ends_raw   = cellfun(@str2num, ends_raw);
score_raw = cellfun(@str2num, score_raw);

% split the BED file into parts for each chromosome if required
if nargin > 1 && split2mat 
    [pathstr,fname, ext] = fileparts(filename); %#ok<ASGLU>
    if nargin == 3
        pathstr = outputpath;
    end
    unichrom = unique(chrom_raw);  
    n_rows = length(chrom_raw);
    n_chrom = length(unichrom);
    
    for n = 1:n_chrom
        index = [];
        curr_chrom = char(unichrom(n));
        outfile = fullfile(pathstr, [fname, '.', curr_chrom, '.mat']);
        for i = 1:n_rows
            if strcmp(char(chrom_raw(i)), curr_chrom)   % if 'chrom' matches current chrome name
                index = [index i]; %#ok<AGROW>
            end
        end
        chrom = chrom_raw(index); %#ok<*NASGU>
        starts = starts_raw(index);
        ends = ends_raw(index);
        name = name_raw(index);
        scores = score_raw(index);
        strand = strand_raw(index);
        
        signal = reads_raw(index);
        post = post_raw(index);
        qvalue = qvalue_raw(index);
        save(outfile, 'chrom', 'starts', 'ends', 'name', 'scores', 'strand', 'signal', 'post', 'qvalue')
    end
end

chrom = chrom_raw;
starts = starts_raw;
ends = ends_raw;
name = name_raw;
scores = score_raw;
strand = strand_raw;

signal = reads_raw;
post = post_raw;
qvalue = qvalue_raw;
end

