%% PROGRAM FILE FOR COMPARING BOS READS TO GENES (FROM FASTA) TO MAKE ID
%
%   File          : BOCS_Simulation.m
%
%   Description   : Simulation for the BOCS algorithm comparing BOS reads
%                   to known gene databases for biomarker identification.
%
%   Created by    : Lee E. Korshoj (lee.korshoj@colorado.edu)
%

close all
clear
clc
warning('off', 'MATLAB:mir_warning_maybe_uninitialized_temporary');

%% USER SPECIFICATIONS ===================================================
% ========================================================================

% Choose database --------------------------------------------------------
% Specify the (1) database type and (2) name of the file. If deviating from
% the 3 built-in database types, coding changes must be made. The file must
% be in the location 'Data/{database_name}.fasta', and the file must be
% in the .fasta format.
%
% Select one of the following for 'g_database'...
% 'resistance'
% 'genetic'
% 'cancer'

g_database='resistance';
database_name='megares_database_v1.01';

% Set location for output file -------------------------------------------
% Specify the folder location for the output .txt file to be written.

file_output_loc='';

% Length of k-mer --------------------------------------------------------
% Specify (1) the k-mer splitting method and (2) the k-mer length.
%
% Select one of the following for 'kmer_split_method'...
% 'constant' (for k-mers all of constant length)
% 'variable' (for k-mers of varying length, centered around the avg 
%             specified by 'kmer_length', normal distribution 
%             with stdev = 2)
% 
% Note the first and last blocks will deviate from this set length or the
% length distribution in order to account for the entire sequence. Also,
% the deviation is capped at +/- length of 4 from the set 'kmer_length'.

kmer_split_method='constant';
kmer_length=10;

% Coverage per nucleotide ------------------------------------------------
% Specify the coverage at which each nucleotide in the sequence is seen in 
% the blocks. Breaks are made in different locations for each additional 
% +1X coverage. The value must be an integer.

gene_coverage=1;

% Number of genes and select genes ---------------------------------------
% Specify the (1) number of genes from which the blocks are comprised and
% (2) the number(s) within the database of the specific genes (if any) to
% use. The genes will be split into blocks and randomized together in a 
% batch with the blocks from all genes. The value for 'num_genes' must 
% be an integer, and the number of entries in 'sel_genes' must match. If
% 'sel_genes' is left blank, random genes will be selected.

num_genes=1;
sel_genes=[];

% Errors -----------------------------------------------------------------
% Specify (1) whether random errors should be inserted and (2) the rate at
% which they are seen. Note that the specified error rate corresponds to 
% the number of random point errors, which is actually only half of the 
% error rate observed in content-based sequencing. So, the actual error 
% rate in the block optical method is double the entered value (enter this 
% value accordingly).
%
% Choose one of the following options for 'err_mode'...
% 'on'  (with errors)
% 'off' (no errors - err_rate variable is not used)

err_mode='off';
err_rate=0.01;

% Penalty score ----------------------------------------------------------
% Specify the score given to genes when no mathces are found for a specific
% block. It is suggested that a value of 0.1 is used for starting and all
% normal analyses.

penalty_score=0.1;

% Thresholding -----------------------------------------------------------
% Specify (1) the multiplier to be multiplied to each of the standard 
% thresholding trends and (2) which of the probability factors to use for 
% thresholding. 
%
% The 'thresh_multiplier' can be thought of as a sensitivity, where the
% values are as follows...
% >1 - This corresponds to a LESS SENSITIVE state (i.e., more genes remain
% in consideration after each block is analyzed --> better for analyses
% with multiple genes and errors).
% <1 This corresponds to a MORE SENSITIVE state (i.e., fewer genes remain
% in consideration after each block is analyzed --> better for analyses
% with a single gene and no errors).
%
% The 'thresh_prob_facts_X' variables are an on/off (1/0) toggle for 
% whether a specific probability factor is used for thresholding or not. 
% There are seven options - the content score (X=CS) and each of the
% six factors comprising the content score (X=F1-F6).

thresh_multiplier=1;
thresh_prob_facts_CS=1;
thresh_prob_facts_F1=1;
thresh_prob_facts_F2=1;
thresh_prob_facts_F3=1;
thresh_prob_facts_F4=1;
thresh_prob_facts_F5=1;
thresh_prob_facts_F6=1;

% Entropy screening ------------------------------------------------------
% Specify (1) the entropy screening mode and (2) the threshold for what is
% considered 'high entropy'. It is suggested to use 10000 as the marker for
% high entropy since there is a natural break around this value.
%
% Select one of the following options for 'entropy_screening_mode'...
% 'rand' (for random screening)
% 'ideal' (for screening idealized from lowest to highest)
% 'none' (for no entropy screening)

entropy_screening_mode='rand';
perms_thresh=10000;

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Analysis/Troubleshooting/Output options --------------------------------
% Specify (1) what kind of analysis is being done, (2) if factor analysis
% is needed (i.e., figures displayed), and (3) the level at which to track
% gene class.
% 
% Select one of the following options for 'analysis_type'...
% 'standard' (for normal operation and output)
% 'benchmarking' (for extra output including all of the factor values for 
%                 all of the genes in the database - for establishing new
%                 thresholding trends)
%
% Select one of the following options for 'disp_fact_figs...
% 'yes' (display)
% 'no' (do not display)
%
% Set the level of output. This 'tracking_level' variable is the number of
% unique subCLASSes of genes with the top content score after each block is
% analyzed. This number should be increased as more genes are combined, and
% helps in analyzing the level of identification of the selected genes
% (i.e., positive and false positive identifications).

analysis_type='standard';
disp_fact_figs='yes';
tracking_level=10;

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% ========================================================================
% ========================================================================
%% READING IN GENE DATABASE FROM FASTA FILE

full_gene_database=fastaread(['Data/',database_name,'.fasta']);
  
%% SELECTION OF GENE

% Gene selection ---------------------------------------------------------
if isempty(sel_genes)
    gene_selection=randperm(length(full_gene_database),num_genes);
else
    if length(sel_genes)~=num_genes
        error('The value entered for number of genes in the analysis (num_genes) is not consistent with the number of selected genes entered (sel_genes).')
    else
        gene_selection=sel_genes;
    end
end

% Collect name and sequence ----------------------------------------------
genes=struct('Header',[],'Sequence',[]);
for i=1:num_genes
    genes(i).Header=full_gene_database(gene_selection(i)).Header;
    genes(i).Sequence=upper(full_gene_database(gene_selection(i)).Sequence);
end

% Find total seq length of all seqs --------------------------------------
tot_seqs_len=0;
for i=1:num_genes
    tot_seqs_len=tot_seqs_len+length(genes(i).Sequence);
end

% Collect subCLASS/CLASS -------------------------------------------------
select_gene_subCLASS=cell(num_genes,1);
select_gene_CLASS=cell(num_genes,1);
for i=1:num_genes
    if strcmpi(g_database,'resistance')
        [select_gene_subCLASS{i}, select_gene_CLASS{i}]=ResExtract(genes(i).Header);
    elseif strcmpi(g_database,'genetic')
        select_gene_subCLASS{i}=GeneticExtract(genes(i).Header);
    elseif strcmpi(g_database,'cancer')
        select_gene_subCLASS{i}=CancerExtract(genes(i).Header);
    end
end

%% CHOPPING INTO BLOCKS

% Random start locations for each gene -----------------------------------
rand_start=zeros(num_genes,gene_coverage);
for i=1:num_genes
    rand_start(i,:)=randperm(kmer_length/2,gene_coverage)+kmer_length/2;
end

% Chop into blocks based on start locations ------------------------------
% Pre-allocate structure
blocks=struct('Header',[],'Blocks',[]);
for i=1:num_genes
    blocks(i).Header=genes(i).Header;
end
% Loop through all genes
for i=1:num_genes
    blocks_hold_out=cell(1,gene_coverage);
    for j=1:gene_coverage
        seq_split=genes(i).Sequence;
        blocks_hold_in=cell(1,1);
        split_its=1;
        while length(seq_split)>kmer_length
            if split_its==1
                blocks_hold_in{split_its,1}=seq_split(1:rand_start(i,j));
                seq_split(1:rand_start(i,j))='';
                split_its=split_its+1;
            else
                if strcmpi(kmer_split_method,'variable')
                    next_split=round(random('normal',kmer_length,2));
                    while next_split<=0 || next_split<(kmer_length-4) || next_split>(kmer_length+4) || next_split>length(seq_split)
                        next_split=round(random('Normal',kmer_length,2));
                    end
                else
                    next_split=kmer_length;
                end
                blocks_hold_in{split_its,1}=seq_split(1:next_split);
                seq_split(1:next_split)='';
                split_its=split_its+1;
            end
        end
        if length(seq_split)>=1
            blocks_hold_in{split_its,1}=seq_split;
        end
        blocks_hold_out{j}=blocks_hold_in;
    end
    blocks(i).Blocks=blocks_hold_out;
end

%% INTRODUCING RANDOM ERRORS AT SPECIFIED RATE

if strcmpi(err_mode,'on')

    % Generate random error points in genes ------------------------------
    % Copy genes structure
    genes_err=genes;
    % Pre-allocate error tracking
    g_len=zeros(num_genes,1);
    n_err=zeros(num_genes,1);
    err_points=cell(num_genes,1);
    % Insert random errors (loop through all genes)
    for i=1:num_genes
        g_len(i)=length(genes_err(i).Sequence);
        n_err(i)=ceil(g_len(i)*err_rate);
        err_points{i}=randperm(g_len(i),n_err(i));
        for j=1:n_err(i)
            pnt_err=randperm(3,1);
            if genes_err(i).Sequence(err_points{i}(j))=='A'
                err_choice='GCT';
                genes_err(i).Sequence(err_points{i}(j))=err_choice(pnt_err);
            elseif genes_err(i).Sequence(err_points{i}(j))=='G'
                err_choice='ACT';
                genes_err(i).Sequence(err_points{i}(j))=err_choice(pnt_err);
            elseif genes_err(i).Sequence(err_points{i}(j))=='C'
                err_choice='AGT';
                genes_err(i).Sequence(err_points{i}(j))=err_choice(pnt_err);
            else
                err_choice='AGC';
                genes_err(i).Sequence(err_points{i}(j))=err_choice(pnt_err);
            end
        end
    end
    
    % Use random error points in blocks ----------------------------------
    % Pre-allocate structure
    blocks_err=struct('Header',[],'Blocks',[]);
    for i=1:num_genes
        blocks_err(i).Header=genes_err(i).Header;
    end
    % Loop through all genes
    for i=1:num_genes
        blocks_hold_out=cell(1,gene_coverage);
        for j=1:gene_coverage
            seq_split=genes_err(i).Sequence;
            blocks_hold_in=cell(1,1);
            split_its=1;
            while length(seq_split)>kmer_length
                if split_its==1
                    blocks_hold_in{split_its,1}=seq_split(1:rand_start(i,j));
                    seq_split(1:rand_start(i,j))='';
                    split_its=split_its+1;
                else
                    if strcmpi(kmer_split_method,'variable')
                        next_split=round(random('normal',kmer_length,2));
                        while next_split<=0 || next_split<(kmer_length-4) || next_split>(kmer_length+4) || next_split>length(seq_split)
                            next_split=round(random('Normal',kmer_length,2));
                        end
                    else
                        next_split=kmer_length;
                    end
                    blocks_hold_in{split_its,1}=seq_split(1:next_split);
                    seq_split(1:next_split)='';
                    split_its=split_its+1;
                end
            end
            if length(seq_split)>=1
                blocks_hold_in{split_its,1}=seq_split;
            end
            blocks_hold_out{j}=blocks_hold_in;
        end
        blocks_err(i).Blocks=blocks_hold_out;
    end
    
    % Replace blocks with blocks_err -------------------------------------
    blocks=blocks_err;

end

%% GENERATING BLOCK CONTENT AND RANDOMIZING BLOCKS FOR TESTING

% Generate block content -------------------------------------------------
% Note that some genes use special letters to denote multiple possible
% letters at that nucleotide position. These special letters are accounted
% for in the mapping part of the content scoring algorithm. For the initial
% content calculations here, only the standard ATGC letters are analyzed.
% Pre-allocate structure
blocks_content=struct('Header',[],'Blocks',[]);
for i=1:num_genes
    blocks_content(i).Header=genes(i).Header;
end
% Loop through all genes
for i=1:num_genes
    for j=1:gene_coverage
        for k=1:length(blocks(i).Blocks{j})
            tot_A=zeros(1,kmer_length);
            tot_G=zeros(1,kmer_length);
            tot_C=zeros(1,kmer_length);
            tot_T=zeros(1,kmer_length);
            for m=1:length(blocks(i).Blocks{j}{k})
                if blocks(i).Blocks{j}{k}(m)=='A'
                    tot_A(m)=1;
                elseif blocks(i).Blocks{j}{k}(m)=='G'
                    tot_G(m)=1;
                elseif blocks(i).Blocks{j}{k}(m)=='C'
                    tot_C(m)=1;
                elseif blocks(i).Blocks{j}{k}(m)=='T'
                    tot_T(m)=1;
                end
            end
            cont_A=sum(tot_A)/length(blocks(i).Blocks{j}{k});
            cont_G=sum(tot_G)/length(blocks(i).Blocks{j}{k});
            cont_C=sum(tot_C)/length(blocks(i).Blocks{j}{k});
            cont_T=sum(tot_T)/length(blocks(i).Blocks{j}{k});
            blocks_content(i).Blocks{j}{k}=[cont_A cont_G cont_C cont_T];
        end
    end
end

% Combining blocks from all genes and randomizing ------------------------
% Pre-allocate structures
blocks_content_test=struct('Header',[],'Error',[],'Blocks',[],'BlocksCont',[],'BlocksGT',[],'BlocksComb',[],'BlocksContComb',[],'BlocksCombGT',[],'Perms',[],'PermsScreen',[],'BlocksCombScreen',[],'BlocksContCombScreen',[],'BlocksCombScreenGT',[]);
% Set fields initially
blocks_content_test(1).Header=cell(1,num_genes);
blocks_content_test(1).Error=cell(1,num_genes);
cov_track=1;
% Loop through genes to get combined Header, Error, and blocks (randomized for each gene)
for i=1:num_genes
    blocks_content_test.Header{i}=blocks_content(i).Header;
    if strcmpi(err_mode,'on')
        blocks_content_test.Error{i}=err_rate;
    else
        blocks_content_test.Error{i}=0;
    end
    for j=1:gene_coverage
        b_len=length(blocks_content(i).Blocks{j});
        b_rand=randperm(b_len,b_len);
        blocks_content_test.Blocks{cov_track}=blocks(i).Blocks{j}(b_rand);
        blocks_content_test.BlocksCont{cov_track}=blocks_content(i).Blocks{j}(b_rand);
        blocks_content_test.BlocksGT{cov_track}=ones(1,b_len)+(i-1);
        cov_track=cov_track+1;
    end
end
% Adding all randomized blocks to one cell matrix (BlocksComb and BlocksContComb) for testing
bcthold=[blocks_content_test.BlocksCont{:}];
bcGTthold=[blocks_content_test.BlocksGT{:}];
blocks_per_cov=zeros(1,length(blocks_content_test.Blocks));
for y=1:length(blocks_content_test.Blocks)
    blocks_per_cov(y)=length(blocks_content_test.Blocks{y});
end
blocks_tot=sum(blocks_per_cov);
bthold=cell(1,blocks_tot);
track=1;
for w=1:length(blocks_content_test.Blocks)
    bthold(track:track+blocks_per_cov(w)-1)=blocks_content_test.Blocks{w};
    track=track+blocks_per_cov(w);
end
b_len_rand=length(bcthold);
final_rand=randperm(b_len_rand,b_len_rand);
blocks_content_test.BlocksComb=bthold(final_rand);
blocks_content_test.BlocksContComb=bcthold(final_rand);
blocks_content_test.BlocksCombGT=bcGTthold(final_rand);

% Entropy screening adjustments ------------------------------------------
% Pre-allocates and/or copies substructures in prep for entropy screening
blocks_content_test.Perms=cell(length(blocks_content_test.BlocksComb),1);
if strcmpi(entropy_screening_mode,'rand')==1 || strcmpi(entropy_screening_mode,'ideal')==1
    blocks_content_test.PermsScreen=cell(length(blocks_content_test.BlocksComb),1);
    blocks_content_test.BlocksCombScreen=blocks_content_test.BlocksComb;
    blocks_content_test.BlocksContCombScreen=blocks_content_test.BlocksContComb;
    blocks_content_test.BlocksCombScreenGT=blocks_content_test.BlocksCombGT;
end

%% CREATING MAPPING STRUCTURE FOR CONTENT SCORING ALGORITHM

% Find how many max iterations (how many total blocks)
tot_blocks=length(blocks_content_test.BlocksComb);
tot_genes_database=length(full_gene_database);
% Pre-allocate structure
gene_database_map=struct('Header',[],'Status',[],'Sequence',[],'ConsensusMap',[],'MatchLoc',[],'BlockLen',[],'ERM',[],'ProbRaw',[],'ProbCumRaw',[],'Prob',[],'ProbCum',[],'ProbCumNorm',[],'D',[],'DCum',[],'PD',[],'PDCum',[],'PDCumSc',[],'PDSlope',[],'Fact1',[],'Fact2',[],'Fact3',[],'Fact4',[],'Fact5',[],'Fact6',[],'CSUN',[],'CS',[],'CSSlope',[]);
% Loop to creat entry for each gene in database
for i=1:tot_genes_database
    gene_database_map(i).Header=full_gene_database(i).Header;
    gene_database_map(i).Status='Possible';
    gene_database_map(i).Sequence=upper(full_gene_database(i).Sequence);
    gene_database_map(i).ConsensusMap=zeros(1,length(full_gene_database(i).Sequence));
    gene_database_map(i).MatchLoc=cell(tot_blocks,1);
    gene_database_map(i).BlockLen=zeros(tot_blocks,1);
    gene_database_map(i).ERM=zeros(tot_blocks,1);
    gene_database_map(i).ProbRaw=zeros(tot_blocks,1);
    gene_database_map(i).ProbCumRaw=zeros(tot_blocks,1);
    gene_database_map(i).Prob=zeros(tot_blocks,1);
    gene_database_map(i).ProbCum=zeros(tot_blocks,1);
    gene_database_map(i).ProbCumNorm=zeros(tot_blocks,1);
    gene_database_map(i).D=zeros(tot_blocks,1);
    gene_database_map(i).DCum=zeros(tot_blocks,1);
    gene_database_map(i).PD=zeros(tot_blocks,1);
    gene_database_map(i).PDCum=zeros(tot_blocks,1);
    gene_database_map(i).PDCumSc=zeros(tot_blocks,1);
    gene_database_map(i).PDSlope=zeros(tot_blocks,1);
    gene_database_map(i).Fact1=zeros(tot_blocks,1);
    gene_database_map(i).Fact2=zeros(tot_blocks,1);
    gene_database_map(i).Fact3=zeros(tot_blocks,1);
    gene_database_map(i).Fact4=zeros(tot_blocks,1);
    gene_database_map(i).Fact5=zeros(tot_blocks,1);
    gene_database_map(i).Fact6=zeros(tot_blocks,1);
    gene_database_map(i).CSUN=zeros(tot_blocks,1);
    gene_database_map(i).CS=zeros(tot_blocks,1);
    gene_database_map(i).CSSlope=zeros(tot_blocks,1);
end

%% ENTROPY CALCULATIONS AND SCREENING

% Calculate entropy (or number of permutations) --------------------------
for i=1:tot_blocks
    blen=length(blocks_content_test.BlocksComb{i});
    num_As=blocks_content_test.BlocksContComb{i}(1)*blen;
    num_Gs=blocks_content_test.BlocksContComb{i}(2)*blen;
    num_Cs=blocks_content_test.BlocksContComb{i}(3)*blen;
    num_Ts=blocks_content_test.BlocksContComb{i}(4)*blen;
    perms=factorial(blen)/(factorial(num_As)*factorial(num_Gs)*factorial(num_Cs)*factorial(num_Ts));
    blocks_content_test.Perms{i}=perms;
end

% Entropy screening - random case ----------------------------------------
if strcmpi(entropy_screening_mode,'rand')==1
    % Loop through all blocks
    mark_end='go';
    blocks_content_test.PermsScreen=blocks_content_test.Perms;
    for i=1:tot_blocks
        % Go to shift if number of permutations is greater than the threshold
        perms=blocks_content_test.PermsScreen{i};
        if perms > perms_thresh
            % Shift all blocks until the end or a block with permutations less than the threshold is reached
            perms_cycle=perms;
            count_perms_cycle=0;
            while (perms_cycle > perms_thresh) && (count_perms_cycle < (tot_blocks-(i-1))) && strcmp(mark_end,'go')
                block_hold=blocks_content_test.BlocksCombScreen{i};
                cont_hold=blocks_content_test.BlocksContCombScreen{i};
                perms_hold=blocks_content_test.PermsScreen{i};
                GT_hold=blocks_content_test.BlocksCombScreenGT(i);
                for j=i:tot_blocks
                    if j==tot_blocks
                        blocks_content_test.BlocksCombScreen{j}=block_hold;
                        blocks_content_test.BlocksContCombScreen{j}=cont_hold;
                        blocks_content_test.BlocksCombScreenGT(j)=GT_hold;
                        blocks_content_test.PermsScreen{j}=perms_hold;
                    else
                        blocks_content_test.BlocksCombScreen{j}=blocks_content_test.BlocksCombScreen{j+1};
                        blocks_content_test.BlocksContCombScreen{j}=blocks_content_test.BlocksContCombScreen{j+1};
                        blocks_content_test.BlocksCombScreenGT(j)=blocks_content_test.BlocksCombScreenGT(j+1);
                        blocks_content_test.PermsScreen{j}=blocks_content_test.PermsScreen{j+1};
                    end
                end
                perms_cycle=blocks_content_test.PermsScreen{i};
                count_perms_cycle=count_perms_cycle+1;
                if count_perms_cycle == (tot_blocks-(i-1))
                    mark_end='done';
                end
            end
        end
    end
end

% Entropy screening - ideal case, ordered low to high --------------------
if strcmpi(entropy_screening_mode,'ideal')==1
    [perms_sort, perms_sort_ind]=sort(cell2mat(blocks_content_test.Perms));
    for i=1:tot_blocks
        blocks_content_test.BlocksCombScreen{i}=blocks_content_test.BlocksComb{perms_sort_ind(i)};
        blocks_content_test.BlocksContCombScreen{i}=blocks_content_test.BlocksContComb{perms_sort_ind(i)};
        blocks_content_test.BlocksCombScreenGT(i)=blocks_content_test.BlocksCombGT(perms_sort_ind(i));
        blocks_content_test.PermsScreen{i}=perms_sort(i);
    end
end

% Selecting test blocks to use -------------------------------------------
blocks_final_test=struct('BlocksComb',[],'BlocksContComb',[],'BlocksCombGT',[],'Perms',[]);
if strcmpi(entropy_screening_mode,'rand')==1 || strcmpi(entropy_screening_mode,'ideal')==1
    blocks_final_test.BlocksComb=blocks_content_test.BlocksCombScreen;
    blocks_final_test.BlocksContComb=blocks_content_test.BlocksContCombScreen;
    blocks_final_test.BlocksCombGT=blocks_content_test.BlocksCombScreenGT;
    blocks_final_test.Perms=blocks_content_test.PermsScreen;
else
    blocks_final_test.BlocksComb=blocks_content_test.BlocksComb;
    blocks_final_test.BlocksContComb=blocks_content_test.BlocksContComb;
    blocks_final_test.BlocksCombGT=blocks_content_test.BlocksCombGT;
    blocks_final_test.Perms=blocks_content_test.Perms;
end

%% MAPPING BLOCKS TO GENE DATABASE - IMPLEMENTING CONTENT SCORING ALGORITHM

% Pre-allocate tracking variables used throughout algorithm --------------
matches=zeros(tot_blocks, tot_genes_database);
ConsensusMap_results=zeros(tot_blocks,tot_genes_database);
analyzed_genes=zeros(tot_blocks,1);
remaining_genes=zeros(tot_blocks,1);
prob_nc_eblock_max_raw=zeros(tot_blocks,1);
prob_nc_eblock_avg=zeros(tot_blocks,1);
prob_nc_eblock_std=zeros(tot_blocks,1);
pdslope_eblock_avg=zeros(tot_blocks,1);
csslope_eblock_avg=zeros(tot_blocks,1);
specificity_tracker=zeros(tot_blocks,1);
cov_tracker=zeros(tot_blocks,num_genes);
cov_tracker_overall=zeros(tot_blocks,1);
tot_blen_cov_tracker=zeros(tot_blocks,num_genes);
tot_blen_cov_tracker_overall=zeros(tot_blocks,1);
ID_tracker=zeros(tot_blocks,num_genes);
Log_TOP_eblock_subCLASS_probsnorm=zeros(tot_blocks,tracking_level);
Log_TOP_eblock_CLASS_probsnorm=zeros(tot_blocks,tracking_level);
Log_TOP_eblock_subCLASS=cell(tot_blocks,tracking_level);
Log_TOP_eblock_CLASS=cell(tot_blocks,tracking_level);

% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
% Content scoring algorithm ----------------------------------------------

% Start timer over the algorithm ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
start_cputime=cputime;

% Loop through all blocks ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for i=1:tot_blocks
    
    % Visual output tracker ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    disp(['Block number: ',num2str(i)])
    disp(['Total blocks: ',num2str(tot_blocks)])
    
    % Determine block and length ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    b_frag=blocks_final_test.BlocksComb{i};
    blen=length(blocks_final_test.BlocksComb{i});
    
    % Reset possible_genes count ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    possible_genes=0;
    track_poss_genes=zeros(1,tot_genes_database);
    
    % Loop through all genes~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    for j=1:tot_genes_database
        
        % Check if status is still 'possible' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if strcmp(gene_database_map(j).Status,'Possible')==1
            
            % Log length of block ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            gene_database_map(j).BlockLen(i)=blen;
            
            % Find length of gene ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            glen=length(gene_database_map(j).Sequence);
            
            % Reset special characters analysis indicator ~~~~~~~~~~~~~~~~
            spec_char_analysis_ind=0;
            
            % Loop through all possible k-mers ~~~~~~~~~~~~~~~~~~~~~~~~~~~
            for w=1:(glen-blen+1)
                
                % Identify k-mer, set match condition ~~~~~~~~~~~~~~~~~~~~
                g_frag=gene_database_map(j).Sequence(w:(w+blen-1));
                match_yn='no';
                
                % Special characters analysis ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                if any(b_frag=='R' | b_frag=='Y' | b_frag=='K' | b_frag=='M' | b_frag=='S' | b_frag=='W' | b_frag=='B' | b_frag=='D' | b_frag=='H' | b_frag=='V' | b_frag=='N')...
                        || any(g_frag=='R' | g_frag=='Y' | g_frag=='K' | g_frag=='M' | g_frag=='S' | g_frag=='W' | g_frag=='B' | g_frag=='D' | g_frag=='H' | g_frag=='V' | g_frag=='N')
                    % Turn on special character analysis indicator
                    spec_char_analysis_ind=1;
                    % Pre-allocate counting character matrices
                    tot_A=zeros(2,blen);
                    tot_G=zeros(2,blen);
                    tot_C=zeros(2,blen);
                    tot_T=zeros(2,blen);
                    % --------------------
                    tot_R=zeros(2,blen);
                    tot_Y=zeros(2,blen);
                    tot_K=zeros(2,blen);
                    tot_M=zeros(2,blen);
                    tot_S=zeros(2,blen);
                    tot_W=zeros(2,blen);
                    tot_B=zeros(2,blen);
                    tot_D=zeros(2,blen);
                    tot_H=zeros(2,blen);
                    tot_V=zeros(2,blen);
                    tot_N=zeros(2,blen);
                    % Find special characters in block and gene
                    for h=1:2
                        if h==1
                            t_frag=g_frag;
                        else
                            t_frag=b_frag;
                        end
                        for m=1:length(t_frag)
                            if t_frag(m)=='A'
                                tot_A(h,m)=1;
                            elseif t_frag(m)=='G'
                                tot_G(h,m)=1;
                            elseif t_frag(m)=='C'
                                tot_C(h,m)=1;
                            elseif t_frag(m)=='T'
                                tot_T(h,m)=1;
                            % -----------------------
                            elseif t_frag(m)=='R'
                                tot_R(h,m)=1;
                            elseif t_frag(m)=='Y'
                                tot_Y(h,m)=1;
                            elseif t_frag(m)=='K'
                                tot_K(h,m)=1;
                            elseif t_frag(m)=='M'
                                tot_M(h,m)=1;
                            elseif t_frag(m)=='S'
                                tot_S(h,m)=1;
                            elseif t_frag(m)=='W'
                                tot_W(h,m)=1;
                            elseif t_frag(m)=='B'
                                tot_B(h,m)=1;
                            elseif t_frag(m)=='D'
                                tot_D(h,m)=1;
                            elseif t_frag(m)=='H'
                                tot_H(h,m)=1;
                            elseif t_frag(m)=='V'
                                tot_V(h,m)=1;
                            elseif t_frag(m)=='N'
                                tot_N(h,m)=1;
                            end
                        end
                    end
                    % Pre-allocate cell matrices for totals, type, and meaning
                    tot_all_let=cell(1,11);
                    typ_all_let=cell(1,11);
                    lib_all_let=cell(1,11);
                    % Fill cell matrices - totals
                    tot_all_let{1}=tot_R;
                    tot_all_let{2}=tot_Y;
                    tot_all_let{3}=tot_K;
                    tot_all_let{4}=tot_M;
                    tot_all_let{5}=tot_S;
                    tot_all_let{6}=tot_W;
                    tot_all_let{7}=tot_B;
                    tot_all_let{8}=tot_D;
                    tot_all_let{9}=tot_H;
                    tot_all_let{10}=tot_V;
                    tot_all_let{11}=tot_N;
                    % Fill cell matrices - type
                    typ_all_let{1}='R';
                    typ_all_let{2}='Y';
                    typ_all_let{3}='K';
                    typ_all_let{4}='M';
                    typ_all_let{5}='S';
                    typ_all_let{6}='W';
                    typ_all_let{7}='B';
                    typ_all_let{8}='D';
                    typ_all_let{9}='H';
                    typ_all_let{10}='V';
                    typ_all_let{11}='N';
                    % Fill cell matrices - meaning
                    lib_all_let{1}='AG';
                    lib_all_let{2}='CT';
                    lib_all_let{3}='GT';
                    lib_all_let{4}='AC';
                    lib_all_let{5}='GC';
                    lib_all_let{6}='AT';
                    lib_all_let{7}='CGT';
                    lib_all_let{8}='AGT';
                    lib_all_let{9}='ACT';
                    lib_all_let{10}='ACG';
                    lib_all_let{11}='ACGT';
                    % Find all permutations in block and gene
                    for h=1:2
                        if h==1
                            t_frag=g_frag;
                        else
                            t_frag=b_frag;
                        end
                        perms_seq=cell(1,1);
                        perms_track=1;
                        for y=1:11
                            if sum(tot_all_let{y}(h,:))>=1
                                k=sum(tot_all_let{y}(h,:));
                                values=lib_all_let{y};
                                n = numel(values);
                                combs = bsxfun(@minus, nchoosek(1:n+k-1,k), 0:k-1);
                                combs = reshape(values(combs),[],k);
                                [num_combs, ~]=size(combs);
                                spec_char_loc=find(t_frag==typ_all_let{y});
                                if isempty(perms_seq{1})
                                    t_frag_log={t_frag};
                                else
                                    t_frag_log=perms_seq;
                                end
                                for p=1:length(t_frag_log)
                                    t_frag_log_sel=t_frag_log{p};
                                    for m=1:num_combs
                                        t_frag_c=t_frag_log_sel;
                                        for q=1:length(spec_char_loc)
                                            t_frag_c(spec_char_loc(q))=combs(m,q);
                                        end
                                        perms_seq{perms_track}=t_frag_c;
                                        perms_track=perms_track+1;
                                    end
                                end
                            end
                        end
                        if isempty(perms_seq{1})
                            perms_seq={t_frag};
                        end
                        if h==1
                            g_frag_perms=perms_seq;
                        else
                            b_frag_perms=perms_seq;
                        end
                    end
                    % Remove any combinations with special letters still remaining
                    for m=1:2
                        if m==1
                            frag_perms=g_frag_perms;
                        else
                            frag_perms=b_frag_perms;
                        end
                        for h=1:length(frag_perms)
                            rev_frag=frag_perms{h};
                            if any(rev_frag=='R' | rev_frag=='Y' | rev_frag=='K' | rev_frag=='M' | rev_frag=='S' | rev_frag=='W' | rev_frag=='B' | rev_frag=='D' | rev_frag=='H' | rev_frag=='V' | rev_frag=='N')
                                frag_perms(h)=[];
                            end
                        end
                        if m==1
                            g_frag_perms=frag_perms;
                        else
                            b_frag_perms=frag_perms;
                        end
                    end
                    % Convert fragments to content
                    for k=1:2
                        if k==1
                            g_frag_perms_cont=cell(size(g_frag_perms));
                            frag_perms=g_frag_perms;
                        else
                            b_frag_perms_cont=cell(size(b_frag_perms));
                            frag_perms=b_frag_perms;
                        end
                        for h=1:length(frag_perms)
                            frag_t=frag_perms{h};
                            tot_A=zeros(1,blen);
                            tot_G=zeros(1,blen);
                            tot_C=zeros(1,blen);
                            tot_T=zeros(1,blen);
                            for m=1:length(frag_t)
                                if frag_t(m)=='A'
                                    tot_A(m)=1;
                                elseif frag_t(m)=='G'
                                    tot_G(m)=1;
                                elseif frag_t(m)=='C'
                                    tot_C(m)=1;
                                elseif frag_t(m)=='T'
                                    tot_T(m)=1;
                                end
                            end
                            cont_A=sum(tot_A)/length(frag_t);
                            cont_G=sum(tot_G)/length(frag_t);
                            cont_C=sum(tot_C)/length(frag_t);
                            cont_T=sum(tot_T)/length(frag_t);
                            test_con=[cont_A cont_G cont_C cont_T];
                            if k==1
                                g_frag_perms_cont{h}=test_con;
                            else
                                b_frag_perms_cont{h}=test_con;
                            end
                        end
                    end
                    % Determine if a match is made
                    for h=1:length(b_frag_perms_cont)
                        for m=1:length(g_frag_perms_cont)
                            if all(b_frag_perms_cont{h}==g_frag_perms_cont{m})
                                match_yn='yes';
                            end
                        end
                    end   
                % Standard AGCT analysis ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                else
                    % Find bases present in gene
                    tot_A=zeros(1,blen);
                    tot_G=zeros(1,blen);
                    tot_C=zeros(1,blen);
                    tot_T=zeros(1,blen);
                    for m=1:length(g_frag)
                        if g_frag(m)=='A'
                            tot_A(m)=1;
                        elseif g_frag(m)=='G'
                            tot_G(m)=1;
                        elseif g_frag(m)=='C'
                            tot_C(m)=1;
                        elseif g_frag(m)=='T'
                            tot_T(m)=1;
                        end
                    end
                    % Calculate content
                    cont_A=sum(tot_A)/length(g_frag);
                    cont_G=sum(tot_G)/length(g_frag);
                    cont_C=sum(tot_C)/length(g_frag);
                    cont_T=sum(tot_T)/length(g_frag);
                    test_con=[cont_A cont_G cont_C cont_T];
                    % Make note of matches
                    if all(test_con==blocks_final_test.BlocksContComb{i})
                        match_yn='yes';
                    end
                end
                
                % Make note of match ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                if strcmp(match_yn,'yes')
                    gene_database_map(j).ConsensusMap(w:(w+blen-1))=1;
                    matches(i,j)=matches(i,j)+1;
                    gene_database_map(j).MatchLoc{i}(matches(i,j))=w;
                end
                
            end
            
            % Eliminations for case of no errors AND single gene ~~~~~~~~~
            if strcmpi(err_mode,'off') && num_genes==1
                % Check if no matches made in gene
                if matches(i,j)==0
                    gene_database_map(j).Status='Out';
                % Checks if only 1 match is made
                elseif matches(i,j)==1 && i>1 && gene_coverage==1
                    % Loop through all previous blocks
                    corr_y_n=0;
                    for k=1:i-1
                        % If one of previous blocks only had one match
                        if matches(k,j)==1
                            i_len=gene_database_map(j).BlockLen(i);
                            i_st=gene_database_map(j).MatchLoc{i};
                            k_len=gene_database_map(j).BlockLen(k);
                            k_st=gene_database_map(j).MatchLoc{k};
                            if i_st - k_st <= 0
                                if (i_st + i_len - 1) >= k_st
                                    gene_database_map(j).Status='Out';
                                    matches(i,j)=0;
                                    corr_y_n=1;
                                end
                            elseif i_st - k_st > 0
                                if (k_st + k_len - 1) >= i_st
                                    gene_database_map(j).Status='Out';
                                    matches(i,j)=0;
                                    corr_y_n=1;
                                end
                            end
                        % If one of previous blocks had more than one match
                        elseif matches(k,j)>1
                            corr_track=0;
                            i_len=gene_database_map(j).BlockLen(i);
                            i_st=gene_database_map(j).MatchLoc{i};
                            k_len=gene_database_map(j).BlockLen(k);
                            k_st=gene_database_map(j).MatchLoc{k};
                            k_st_corr=k_st;
                            k_st_corr_l=ones(1,length(k_st));
                            for r=1:matches(k,j)
                                if i_st - k_st(r) <= 0
                                    if (i_st + i_len - 1) >= k_st(r)
                                        gene_database_map(j).ConsensusMap(k_st(r):(k_st(r)+k_len-1))=0;
                                        corr_track=corr_track+1;
                                        k_st_corr_l(r)=0;
                                        corr_y_n=1;
                                    end  
                                elseif i_st - k_st(r) > 0
                                    if (k_st(r) + k_len - 1) >= i_st
                                        gene_database_map(j).ConsensusMap(k_st(r):(k_st(r)+k_len-1))=0;
                                        corr_track=corr_track+1;
                                        k_st_corr_l(r)=0;
                                        corr_y_n=1;
                                    end
                                end
                            end
                            % Correct number of match locations
                            k_st_corr_l=logical(k_st_corr_l);
                            k_st_corr=k_st_corr(k_st_corr_l);
                            gene_database_map(j).MatchLoc{k}=zeros(1,length(k_st_corr));
                            for u=1:length(k_st_corr)
                                gene_database_map(j).MatchLoc{k}(u)=k_st_corr(u);
                            end
                            % Correct number of matches and eliminate gene if zero
                            matches(k,j)=matches(k,j)-corr_track;
                            if matches(k,j)==0
                                gene_database_map(j).Status='Out';
                            end
                        end
                    end
                    % Re-fill in ConsensusMap
                    for k=1:i
                        k_st=gene_database_map(j).MatchLoc{k};
                        k_len=gene_database_map(j).BlockLen(k);
                        for r=1:matches(k,j)
                            gene_database_map(j).ConsensusMap(k_st(r):(k_st(r)+k_len-1))=1;
                        end
                    end
                end
            end
            
            % Use ConsensusMap to make eliminations ~~~~~~~~~~~~~~~~~~~~~~
            if strcmpi(err_mode,'off') && num_genes==1 && gene_coverage==1
                nucleotides_seen=sum(gene_database_map(j).BlockLen(1:i));
                nucleotides_mapped=sum(gene_database_map(j).ConsensusMap);
                if nucleotides_seen > nucleotides_mapped
                    gene_database_map(j).Status='Out';
                end
            end
            
            % Raw probabilities ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            % Calculate permutations
            if spec_char_analysis_ind==1
                block_cont_variations=b_frag_perms_cont;
            else
                block_cont_variations={blocks_final_test.BlocksContComb{i}};
            end
            perms_cont_variations=zeros(1,length(block_cont_variations));
            for h=1:length(block_cont_variations)
                num_As=block_cont_variations{h}(1)*blen;
                num_Gs=block_cont_variations{h}(2)*blen;
                num_Cs=block_cont_variations{h}(3)*blen;
                num_Ts=block_cont_variations{h}(4)*blen;
                perms_cont_variations(h)=factorial(blen)/(factorial(num_As)*factorial(num_Gs)*factorial(num_Cs)*factorial(num_Ts));
            end
            perms=sum(perms_cont_variations);
            % Calculate expected matches
            rand_matches=(glen-blen+1)/(4^blen);
            % Calculate expected random matches from permutations and expected matches
            ERM=perms*rand_matches;
            gene_database_map(j).ERM(i)=ERM;
            % Calculate raw probability for this block - number of observed matches / expected random matches (or the penalty score if no matches)
            if matches(i,j) > 0
                gene_database_map(j).ProbRaw(i)=matches(i,j)/gene_database_map(j).ERM(i);
            else
                gene_database_map(j).ProbRaw(i)=penalty_score;
            end
            % Calculate raw cumulative probability - sum all raw probs
            gene_database_map(j).ProbCumRaw(i)=sum(gene_database_map(j).ProbRaw(1:i));
            
            % Tracking genes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            possible_genes=possible_genes+1;
            track_poss_genes(j)=j;
            
        end
        
    end
    
    % Add up genes that were looked at this round ~~~~~~~~~~~~~~~~~~~~~~~~
    if strcmpi(err_mode,'off') && num_genes==1
        analyzed_genes(i)=sum(matches(i,:)>0);
    else
        analyzed_genes(i)=possible_genes;
    end
    
    % Check for case of no genes remaining and exit program ~~~~~~~~~~~~~~
    if analyzed_genes(i)==0
        disp('--------------------------------------------------------------------------------------------')
        disp('--------------------------------------------------------------------------------------------')
        disp('Analysis terminated: All genes have been eliminated. Adjust simulation settings and re-run.')
        disp('--------------------------------------------------------------------------------------------')
        disp('--------------------------------------------------------------------------------------------')
        break
    end
    
    % Compile genes remaining ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if strcmpi(err_mode,'off') && num_genes==1
        rem_vals_eblock=find(matches(i,:));
    else
        rem_vals_eblock=track_poss_genes(track_poss_genes>0);
    end
    
    % ConsensusMap analysis ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    for j=1:analyzed_genes(i)
        ConsensusMap_results(i,rem_vals_eblock(j))=sum(gene_database_map(rem_vals_eblock(j)).ConsensusMap)/length(gene_database_map(rem_vals_eblock(j)).ConsensusMap);
    end
    
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % CALCULATE BASE PROBS AND SCALING/NORMALIZING ~~~~~~~~~~~~~~~~~~~~~~~
    % Collect ProbRaw values to get max for normalization
    prob_nc_eblock_all=zeros(1,analyzed_genes(i));
    for q=1:analyzed_genes(i)
        prob_nc_eblock_all(1,q)=gene_database_map(rem_vals_eblock(q)).ProbRaw(i);
    end
    prob_nc_eblock_max_raw(i)=max(prob_nc_eblock_all);
    % Scale ProbRaw to Probs and calculate ProbCum
    for q=1:analyzed_genes(i)
        gene_database_map(rem_vals_eblock(q)).Prob(i)=gene_database_map(rem_vals_eblock(q)).ProbRaw(i)/prob_nc_eblock_max_raw(i);
        gene_database_map(rem_vals_eblock(q)).ProbCum(i)=sum(gene_database_map(rem_vals_eblock(q)).Prob(1:i));
    end
    % Collect all Probs and ProbCum to get avg, std, and total
    prob_nc_eblock_all_stats=zeros(1,analyzed_genes(i));
    probcum_c_eblock_all_norm=zeros(1,analyzed_genes(i));
    for q=1:analyzed_genes(i)
        prob_nc_eblock_all_stats(1,q)=gene_database_map(rem_vals_eblock(q)).Prob(i);
        probcum_c_eblock_all_norm(1,q)=gene_database_map(rem_vals_eblock(q)).ProbCum(i);
    end
    prob_nc_eblock_avg(i)=mean(prob_nc_eblock_all_stats);
    prob_nc_eblock_std(i)=std(prob_nc_eblock_all_stats);
    tot_probs_c=sum(probcum_c_eblock_all_norm);
    % Calculate ProbCumNorm, PD, and PDCum from cumulative and normalized trends
    pdcum_hold_fact=zeros(1,analyzed_genes(i));
    for q=1:analyzed_genes(i)
        gene_database_map(rem_vals_eblock(q)).ProbCumNorm(i)=gene_database_map(rem_vals_eblock(q)).ProbCum(i)/tot_probs_c;
        gene_database_map(rem_vals_eblock(q)).D(i)=gene_database_map(rem_vals_eblock(q)).Prob(i)-prob_nc_eblock_avg(i);
        gene_database_map(rem_vals_eblock(q)).DCum(i)=sum(gene_database_map(rem_vals_eblock(q)).D(1:i));
        gene_database_map(rem_vals_eblock(q)).PD(i)=(gene_database_map(rem_vals_eblock(q)).Prob(i)-prob_nc_eblock_avg(i))/prob_nc_eblock_avg(i);
        gene_database_map(rem_vals_eblock(q)).PDCum(i)=sum(gene_database_map(rem_vals_eblock(q)).PD(1:i));
        pdcum_hold_fact(q)=gene_database_map(rem_vals_eblock(q)).PDCum(i);
    end
    pdcum_hold_fact_max=max(pdcum_hold_fact);
    % Calculate slopes for PDCum and PDCumSc
    pdslope_hold=zeros(1,analyzed_genes(i));
    for q=1:analyzed_genes(i)
        gene_database_map(rem_vals_eblock(q)).PDCumSc(i)=sum(gene_database_map(rem_vals_eblock(q)).PD(1:i))/pdcum_hold_fact_max;
        if i>9
            p9=polyfit((i-9:i)',gene_database_map(rem_vals_eblock(q)).PDCumSc(i-9:i),1);
            pdslope_hold(q)=p9(1);
            gene_database_map(rem_vals_eblock(q)).PDSlope(i)=p9(1);
        end
    end
    pdslope_eblock_avg(i)=mean(pdslope_hold);
    
    % FACTORS CALCULATION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % 1 - (Cumulative percent diff from avg) * (Cumulative normalized probs)
    % 2 - Cumulative number of blocks having at least one match
    % 3 - Cumulative probs multiplied together (taken as log2)
    % 4 - Exponential of gene coverage
    % 5 - Cumulative of the percent diff from avg slope
    % 6 - Cumulative diff from avg
    %---------------------------------------------------------------------
    fact1_c_eblock_all=zeros(1,analyzed_genes(i));
    %-----
    fact2_c_eblock_all=zeros(1,analyzed_genes(i));
    %-----
    fact3_c_eblock_all=zeros(1,analyzed_genes(i));
    %-----
    fact4_c_eblock_all=zeros(1,analyzed_genes(i));
    %-----
    fact5_c_eblock_all=zeros(1,analyzed_genes(i));
    %-----
    fact6_c_eblock_all=zeros(1,analyzed_genes(i));
    %-----
    for q=1:analyzed_genes(i)
        %-----
        fact1_c_eblock_all(q)=gene_database_map(rem_vals_eblock(q)).PDCum(i)*gene_database_map(rem_vals_eblock(q)).ProbCumNorm(i);
        %-----
        fact2_c_eblock_all(q)=sum(gene_database_map(rem_vals_eblock(q)).ProbRaw(1:i)>penalty_score);
        %-----
        fact3_c_eblock_all(q)=sum(log2(gene_database_map(rem_vals_eblock(q)).Prob(1:i)));
        %-----
        fact4_c_eblock_all(q)=ConsensusMap_results(i,rem_vals_eblock(q));
        %-----
        fact5_c_eblock_all(q)=sum(gene_database_map(rem_vals_eblock(q)).PDSlope(1:i));
        %-----
        fact6_c_eblock_all(q)=gene_database_map(rem_vals_eblock(q)).DCum(i);
        %-----
    end
    %---------------------------------------------------------------------
    if analyzed_genes(i)>1
        fact1_c_eblock_all_min=min(fact1_c_eblock_all);
        fact1_c_eblock_all=fact1_c_eblock_all+((-1)*fact1_c_eblock_all_min);
        %-----
        fact5_c_eblock_all_min=min(fact5_c_eblock_all);
        fact5_c_eblock_all=fact5_c_eblock_all+((-1)*fact5_c_eblock_all_min);
        %-----
        fact6_c_eblock_all_min=min(fact6_c_eblock_all);
        fact6_c_eblock_all=fact6_c_eblock_all+((-1)*fact6_c_eblock_all_min);
    end
    %---------------------------------------------------------------------
    fact1_c_eblock_max=max(fact1_c_eblock_all);
    %-----
    fact2_c_eblock_max=max(fact2_c_eblock_all);
    %-----
    fact3_hold_log_max=max(abs(fact3_c_eblock_all));
    fact3_hold_log_max_sub=fact3_hold_log_max-abs(fact3_c_eblock_all);
    fact3_c_eblock_max=max(fact3_hold_log_max_sub);
    %-----
    fact4_hold=exp(500*fact4_c_eblock_all)/exp(500);
    fact4_c_eblock_max=max(fact4_hold);
    %-----
    fact5_c_eblock_max=max(fact5_c_eblock_all);
    %-----
    fact6_c_eblock_max=max(fact6_c_eblock_all);
    %-----
    for q=1:analyzed_genes(i)
        %-----
        if fact1_c_eblock_max==0
            gene_database_map(rem_vals_eblock(q)).Fact1(i)=1;
        else
            gene_database_map(rem_vals_eblock(q)).Fact1(i)=fact1_c_eblock_all(q)/fact1_c_eblock_max;
        end
        %-----
        gene_database_map(rem_vals_eblock(q)).Fact2(i)=fact2_c_eblock_all(q)/fact2_c_eblock_max;
        %-----
        if analyzed_genes(i)==1 || fact3_c_eblock_max==0
            gene_database_map(rem_vals_eblock(q)).Fact3(i)=1;
        else
            gene_database_map(rem_vals_eblock(q)).Fact3(i)=fact3_hold_log_max_sub(q)/fact3_c_eblock_max;
        end
        %-----
        gene_database_map(rem_vals_eblock(q)).Fact4(i)=fact4_hold(q)/fact4_c_eblock_max;
        %-----
        if i>9
            if fact5_c_eblock_max==0
                gene_database_map(rem_vals_eblock(q)).Fact5(i)=1;
            else
                gene_database_map(rem_vals_eblock(q)).Fact5(i)=fact5_c_eblock_all(q)/fact5_c_eblock_max;
            end
        end
        %-----
        if fact6_c_eblock_max==0
            gene_database_map(rem_vals_eblock(q)).Fact6(i)=1;
        else
            gene_database_map(rem_vals_eblock(q)).Fact6(i)=fact6_c_eblock_all(q)/fact6_c_eblock_max;
        end
        %-----
    end
    %---------------------------------------------------------------------
    
    % CONTENT SCORE CALCULATION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Adding all 6 factors together and normalizing
    CSUN_hold=zeros(1,analyzed_genes(i));
    for q=1:analyzed_genes(i)
        gene_database_map(rem_vals_eblock(q)).CSUN(i)=gene_database_map(rem_vals_eblock(q)).Fact1(i) + gene_database_map(rem_vals_eblock(q)).Fact2(i) + gene_database_map(rem_vals_eblock(q)).Fact3(i) + gene_database_map(rem_vals_eblock(q)).Fact4(i) + gene_database_map(rem_vals_eblock(q)).Fact5(i) + gene_database_map(rem_vals_eblock(q)).Fact6(i);
        CSUN_hold(q)=gene_database_map(rem_vals_eblock(q)).CSUN(i);
    end
    CSUN_eblock_sum=abs(sum(CSUN_hold));
    for q=1:analyzed_genes(i)
        gene_database_map(rem_vals_eblock(q)).CS(i)=CSUN_hold(q)/CSUN_eblock_sum;
    end
    % Slopes for CS
    csslope_hold=zeros(1,analyzed_genes(i));
    for q=1:analyzed_genes(i)
        if i>9
            p9=polyfit((i-9:i)',gene_database_map(rem_vals_eblock(q)).CS(i-9:i),1);
            csslope_hold(q)=p9(1);
            gene_database_map(rem_vals_eblock(q)).CSSlope(i)=p9(1);
        end
    end
    csslope_eblock_avg(i)=mean(csslope_hold);
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    % Compile data on genes remaining ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ALL_eblock_genes=cell(analyzed_genes(i),1);
    ALL_eblock_genes_probsnorm=zeros(analyzed_genes(i),1);
    for q=1:analyzed_genes(i)
        ALL_eblock_genes{q}=full_gene_database(rem_vals_eblock(q)).Header;
        ALL_eblock_genes_probsnorm(q)=gene_database_map(rem_vals_eblock(q)).CS(i);
    end
    
    % Check if largest prob is rand_select gene(s) ~~~~~~~~~~~~~~~~~~~~~~~
    [ALL_eblock_genes_probsnorm_sort, ALL_eblock_genes_probsnorm_sort_inds]=sort(ALL_eblock_genes_probsnorm,'descend');
    for y=1:num_genes
        if any(ALL_eblock_genes_probsnorm_sort(1:num_genes)==gene_database_map(gene_selection(y)).CS(i))
            ID_tracker(i,y)=1;
        end
    end
    
    % Order genes remaining ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % See if 10 are available to rank
    if length(ALL_eblock_genes_probsnorm_sort)>=10
        sort_num=10;
    else
        sort_num=length(ALL_eblock_genes_probsnorm_sort);
    end
    % Collect indices of genes with top probs (multiples if ties)
    top_gprobs_track=1;
    top_gprobs_ind=zeros();
    for y=1:sort_num
        ind=find(ALL_eblock_genes_probsnorm==ALL_eblock_genes_probsnorm_sort(y));
        for q=1:length(ind)
            if ~ismember(ind(q),top_gprobs_ind)
                top_gprobs_ind(top_gprobs_track)=ind(q);
                top_gprobs_track=top_gprobs_track+1;
            end
        end
    end
    % Compile top genes from sorting and indices analysis
    top_gprobs_tot=top_gprobs_track-1;
    ALL_eblock_genes_top=cell(top_gprobs_tot,1);
    ALL_eblock_probsnorm_top=zeros(top_gprobs_tot,1);
    for y=1:top_gprobs_tot
        ALL_eblock_genes_top{y}=ALL_eblock_genes{top_gprobs_ind(y)};
        ALL_eblock_probsnorm_top(y)=ALL_eblock_genes_probsnorm(top_gprobs_ind(y));
    end
    
    % Check for different subCLASSes in top genes and correct ~~~~~~~~~~~~
    % Collect subCLASS for all top genes
    check_subCLASS=cell(length(ALL_eblock_genes_top),1);
    for y=1:length(ALL_eblock_genes_top)
        if strcmpi(g_database,'resistance')
            [ext, ~]=ResExtract(ALL_eblock_genes_top{y});
        elseif strcmpi(g_database,'genetic')
            ext=GeneticExtract(ALL_eblock_genes_top{y});
        elseif strcmpi(g_database,'cancer')
            ext=CancerExtract(ALL_eblock_genes_top{y});
        end
        check_subCLASS{y}=upper(ext);
    end
    % Correct top genes list to have set number of different subCLASSes
    top_num_genes=length(ALL_eblock_genes_top);
    tot_diff_subCLASS=length(unique(check_subCLASS));
    while tot_diff_subCLASS<tracking_level && top_num_genes<analyzed_genes(i)
        ALL_eblock_genes_top{top_num_genes+1}=ALL_eblock_genes{ALL_eblock_genes_probsnorm_sort_inds(top_num_genes+1)};
        ALL_eblock_probsnorm_top(top_num_genes+1)=ALL_eblock_genes_probsnorm(ALL_eblock_genes_probsnorm_sort_inds(top_num_genes+1));
        if strcmpi(g_database,'resistance')
            [ext, ~]=ResExtract(ALL_eblock_genes_top{top_num_genes+1});
        elseif strcmpi(g_database,'genetic')
            ext=GeneticExtract(ALL_eblock_genes_top{top_num_genes+1});
        elseif strcmpi(g_database,'cancer')
            ext=CancerExtract(ALL_eblock_genes_top{top_num_genes+1});
        end
        check_subCLASS{top_num_genes+1}=upper(ext);
        top_num_genes=length(ALL_eblock_genes_top);
        tot_diff_subCLASS=length(unique(check_subCLASS));
    end
    
    % Get correct loop settings and prefill cells ~~~~~~~~~~~~~~~~~~~~~~~~
    % Set loop size (tracking_level or tot_diff_subCLASS - whichever is less)
    if tracking_level<tot_diff_subCLASS
        class_loop_sz=tracking_level;
    else
        class_loop_sz=tot_diff_subCLASS;
    end
    % Fill slots in Log_TOP_eblock_subCLASS/CLASS that will be empty
    fill_spaces=tracking_level-class_loop_sz;
    if fill_spaces>0
        for y=1:fill_spaces
            Log_TOP_eblock_subCLASS{i,tracking_level-fill_spaces+y}='-';
            Log_TOP_eblock_CLASS{i,tracking_level-fill_spaces+y}='-';
        end
    end
        
    % Extract class of top gene(s) and find probs ~~~~~~~~~~~~~~~~~~~~~~~~
    % Loop through tracking level setpoint
    y_track=1;
    for y=1:class_loop_sz
        % Find subCLASS/CLASS with next highest probability
        uniq_subCLASS=false;
        while uniq_subCLASS==false
            if y_track>length(ALL_eblock_genes_top) 
                header=ALL_eblock_genes_top{length(ALL_eblock_genes_top)};
                y_track=length(ALL_eblock_genes_top);
            else
                header=ALL_eblock_genes_top{y_track};
            end
            if strcmpi(g_database,'resistance')
                [Log_TOP_eblock_subCLASS{i,y}, Log_TOP_eblock_CLASS{i,y}]=ResExtract(header);
                Log_TOP_eblock_subCLASS_probsnorm(i,y)=ALL_eblock_probsnorm_top(y_track);
                Log_TOP_eblock_CLASS_probsnorm(i,y)=ALL_eblock_probsnorm_top(y_track);
            elseif strcmpi(g_database,'genetic')
                Log_TOP_eblock_subCLASS{i,y}=GeneticExtract(header);
                Log_TOP_eblock_subCLASS_probsnorm(i,y)=ALL_eblock_probsnorm_top(y_track);
            elseif strcmpi(g_database,'cancer')
                Log_TOP_eblock_subCLASS{i,y}=CancerExtract(header);
                Log_TOP_eblock_subCLASS_probsnorm(i,y)=ALL_eblock_probsnorm_top(y_track);
            end
            if y>1 && y_track<length(ALL_eblock_genes_top) && ismember(upper(Log_TOP_eblock_subCLASS(i,y)),upper(Log_TOP_eblock_subCLASS(i,1:y-1)))
                y_track=y_track+1;
            else
                uniq_subCLASS=true;
                y_track=y_track+1;
            end
        end
    end
    
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Use thresholding to make eliminations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Collect prob values for each gene
    thresh_cs=zeros(length(full_gene_database),1);
    thresh_F1=zeros(length(full_gene_database),1);
    thresh_F2=zeros(length(full_gene_database),1);
    thresh_F3=zeros(length(full_gene_database),1);
    thresh_F4=zeros(length(full_gene_database),1);
    thresh_F5=zeros(length(full_gene_database),1);
    thresh_F6=zeros(length(full_gene_database),1);
    for y=1:length(full_gene_database)
        thresh_cs(y)=gene_database_map(y).CS(i);
        thresh_F1(y)=gene_database_map(y).Fact1(i);
        thresh_F2(y)=gene_database_map(y).Fact2(i);
        thresh_F3(y)=gene_database_map(y).Fact3(i);
        thresh_F4(y)=gene_database_map(y).Fact4(i);
        thresh_F5(y)=gene_database_map(y).Fact5(i);
        thresh_F6(y)=gene_database_map(y).Fact6(i);
    end
    % Rank prob values
    [thresh_cs_sort, thresh_cs_ind]=sort(thresh_cs,'descend');
    [thresh_F1_sort, thresh_F1_ind]=sort(thresh_F1,'descend');
    [thresh_F2_sort, thresh_F2_ind]=sort(thresh_F2,'descend');
    [thresh_F3_sort, thresh_F3_ind]=sort(thresh_F3,'descend');
    [thresh_F4_sort, thresh_F4_ind]=sort(thresh_F4,'descend');
    [thresh_F5_sort, thresh_F5_ind]=sort(thresh_F5,'descend');
    [thresh_F6_sort, thresh_F6_ind]=sort(thresh_F6,'descend');
    % Use thresholding trends to identify at which genes (ranked) to remove
    if i<=25
        cutoff_rank=ceil((1.03-0.02575*i)*length(full_gene_database)*thresh_multiplier);
        if cutoff_rank>length(full_gene_database)
            cutoff_rank=length(full_gene_database);
        elseif cutoff_rank<5
            cutoff_rank=5;
        end
    else
        cutoff_rank=ceil((52.22*i^(-1.525))*length(full_gene_database)*thresh_multiplier);
        if cutoff_rank>length(full_gene_database)
            cutoff_rank=length(full_gene_database);
        elseif cutoff_rank<5
            cutoff_rank=5;
        end
    end
    % Look at ties in probability
    cutoff_tie_cs=sum(thresh_cs_sort(cutoff_rank+1:length(full_gene_database))==thresh_cs_sort(cutoff_rank));
    cutoff_tie_F1=sum(thresh_F1_sort(cutoff_rank+1:length(full_gene_database))==thresh_F1_sort(cutoff_rank));
    cutoff_tie_F2=sum(thresh_F2_sort(cutoff_rank+1:length(full_gene_database))==thresh_F2_sort(cutoff_rank));
    cutoff_tie_F3=sum(thresh_F3_sort(cutoff_rank+1:length(full_gene_database))==thresh_F3_sort(cutoff_rank));
    cutoff_tie_F4=sum(thresh_F4_sort(cutoff_rank+1:length(full_gene_database))==thresh_F4_sort(cutoff_rank));
    cutoff_tie_F5=sum(thresh_F5_sort(cutoff_rank+1:length(full_gene_database))==thresh_F5_sort(cutoff_rank));
    cutoff_tie_F6=sum(thresh_F6_sort(cutoff_rank+1:length(full_gene_database))==thresh_F6_sort(cutoff_rank));
    % Remove genes
    % CS ----------
    if thresh_prob_facts_CS==1
        if cutoff_rank+cutoff_tie_cs<length(full_gene_database)
            thresh_elim=thresh_cs_ind(cutoff_rank+cutoff_tie_cs+1:length(full_gene_database));
            for y=1:length(thresh_elim)
                gene_database_map(thresh_elim(y)).Status='Out';
            end
        end
    end
    % F1 ----------
    if thresh_prob_facts_F1==1
        if cutoff_rank+cutoff_tie_F1<length(full_gene_database)
            thresh_elim=thresh_F1_ind(cutoff_rank+cutoff_tie_F1+1:length(full_gene_database));
            for y=1:length(thresh_elim)
                gene_database_map(thresh_elim(y)).Status='Out';
            end
        end
    end
    % F2 ----------
    if thresh_prob_facts_F2==1
        if cutoff_rank+cutoff_tie_F2<length(full_gene_database)
            thresh_elim=thresh_F2_ind(cutoff_rank+cutoff_tie_F2+1:length(full_gene_database));
            for y=1:length(thresh_elim)
                gene_database_map(thresh_elim(y)).Status='Out';
            end
        end
    end
    % F3 ----------
    if thresh_prob_facts_F3==1
        if cutoff_rank+cutoff_tie_F3<length(full_gene_database)
            thresh_elim=thresh_F3_ind(cutoff_rank+cutoff_tie_F3+1:length(full_gene_database));
            for y=1:length(thresh_elim)
                gene_database_map(thresh_elim(y)).Status='Out';
            end
        end
    end
    % F4 ----------
    if thresh_prob_facts_F4==1
        if cutoff_rank+cutoff_tie_F4<length(full_gene_database)
            thresh_elim=thresh_F4_ind(cutoff_rank+cutoff_tie_F4+1:length(full_gene_database));
            for y=1:length(thresh_elim)
                gene_database_map(thresh_elim(y)).Status='Out';
            end
        end
    end
    % F5 ----------
    if thresh_prob_facts_F5==1
        if cutoff_rank+cutoff_tie_F5<length(full_gene_database)
            thresh_elim=thresh_F5_ind(cutoff_rank+cutoff_tie_F5+1:length(full_gene_database));
            for y=1:length(thresh_elim)
                gene_database_map(thresh_elim(y)).Status='Out';
            end
        end
    end
    % F6 ----------
    if thresh_prob_facts_F6==1
        if cutoff_rank+cutoff_tie_F6<length(full_gene_database)
            thresh_elim=thresh_F6_ind(cutoff_rank+cutoff_tie_F6+1:length(full_gene_database));
            for y=1:length(thresh_elim)
                gene_database_map(thresh_elim(y)).Status='Out';
            end
        end
    end
    % Determine genes remaining
    rems=0;
    for y=1:length(full_gene_database)
        if strcmpi(gene_database_map(y).Status,'Out')
            rems=rems+1;
        end
    end
    remaining_genes(i)=length(full_gene_database)-rems;
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    % Track coverage of gene(s) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Overall coverage (for all genes comprising blocks)
    if i==1
        tot_blen_cov_tracker_overall(i)=length(blocks_final_test.BlocksComb{i});
    else
        tot_blen_cov_tracker_overall(i)=tot_blen_cov_tracker_overall(i-1)+length(blocks_final_test.BlocksComb{i});
    end
    cov_tracker_overall(i)=tot_blen_cov_tracker_overall(i)/tot_seqs_len;
    % Coverage for individual genes
    GT_select=blocks_final_test.BlocksCombGT(i);
    for j=1:num_genes
        if GT_select==j
            if i==1
                tot_blen_cov_tracker(i,j)=length(blocks_final_test.BlocksComb{i});
            else
                tot_blen_cov_tracker(i,j)=tot_blen_cov_tracker(i-1,j)+length(blocks_final_test.BlocksComb{i});
            end
            cov_tracker(i,j)=tot_blen_cov_tracker(i,j)/length(gene_database_map(gene_selection(j)).Sequence);
        else
            if i==1
                tot_blen_cov_tracker(i,j)=0;
                cov_tracker(i,j)=0;
            else
                tot_blen_cov_tracker(i,j)=tot_blen_cov_tracker(i-1,j);
                cov_tracker(i,j)=cov_tracker(i-1,j);
            end
        end
    end
    
    % Track specificity ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    specificity_tracker(i)=(length(full_gene_database)-remaining_genes(i))/length(full_gene_database);
    
    % Visual output tracker ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    status_all=cell(1,num_genes);
    status_string=char();
    cum_prob_all=zeros(1,num_genes);
    strack=1;
    for j=1:num_genes
        status_all{j}=gene_database_map(gene_selection(j)).Status;
        status_string(strack:strack+length(status_all{j})-1)=status_all{j};
        if j~=num_genes
            status_string(strack+length(status_all{j}))=' ';
        end
        strack=strack+length(status_all{j})+1;
        cum_prob_all(j)=gene_database_map(gene_selection(j)).CS(i);
    end
    disp(['gene_selection status: ',status_string])
    disp(['gene_selection CS: ',num2str(cum_prob_all,4)])
    disp(['Genes remaining: ',num2str(remaining_genes(i))])
    fprintf('\n')
    
end

% End timer (in minutes) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
end_cputime=(cputime-start_cputime)/60;
disp(['Runtime (min): ',num2str(end_cputime)])
fprintf('\n')

% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------

%% COMPILE RESULTS FROM CONTENT SCORING ALGORITHM

% Compile cumulative probabilities - all ---------------------------------
CS_all=zeros(tot_blocks,length(full_gene_database));
Fact1_all=zeros(tot_blocks,length(full_gene_database));
Fact2_all=zeros(tot_blocks,length(full_gene_database));
Fact3_all=zeros(tot_blocks,length(full_gene_database));
Fact4_all=zeros(tot_blocks,length(full_gene_database));
Fact5_all=zeros(tot_blocks,length(full_gene_database));
Fact6_all=zeros(tot_blocks,length(full_gene_database));
PDSlope_all=zeros(tot_blocks,length(full_gene_database));
CSSlope_all=zeros(tot_blocks,length(full_gene_database));
for w=1:length(full_gene_database)
    CS_all(:,w)=gene_database_map(w).CS(:);
    Fact1_all(:,w)=gene_database_map(w).Fact1(:);
    Fact2_all(:,w)=gene_database_map(w).Fact2(:);
    Fact3_all(:,w)=gene_database_map(w).Fact3(:);
    Fact4_all(:,w)=gene_database_map(w).Fact4(:);
    Fact5_all(:,w)=gene_database_map(w).Fact5(:);
    Fact6_all(:,w)=gene_database_map(w).Fact6(:);
    PDSlope_all(:,w)=gene_database_map(w).PDSlope(:);
    CSSlope_all(:,w)=gene_database_map(w).CSSlope(:);
end

% Compile cumulative probabilities - select genes ------------------------
CS_select=zeros(tot_blocks,num_genes);
for w=1:num_genes
    CS_select(:,w)=gene_database_map(gene_selection(w)).CS(:);
end

%% OUTPUT TO TXT FILE FOR ANALYSIS

% Set filename -----------------------------------------------------------
select_string=char();
rstrack=1;
for i=1:num_genes
    rs_len=length(num2str(gene_selection(i)));
    select_string(rstrack:rstrack+rs_len-1)=num2str(gene_selection(i));
    if i~=num_genes
        select_string(rstrack+rs_len)='_';
    end
    rstrack=rstrack+rs_len+1;
end
fid = fopen([file_output_loc,'/',g_database,' ',select_string,' ','k-mers_',kmer_split_method,'_',num2str(kmer_length),' ',num2str(gene_coverage),' ',num2str(num_genes),' ','error_',err_mode,'_',num2str(err_rate),' ',num2str(penalty_score),' ',num2str(thresh_multiplier),' ','entropy_',entropy_screening_mode,'_',num2str(perms_thresh),' ',analysis_type,' ',num2str(tracking_level),' ','BOS Simulation.txt'],'wt');

% Runtime ----------------------------------------------------------------
fprintf(fid,'%s\n','RUNTIME--------------------------------------------------------------------------------');
fprintf(fid,'%f %s\n',end_cputime,'min');

% Inputs -----------------------------------------------------------------
fprintf(fid,'\n%s\n','INPUTS--------------------------------------------------------------------------------');
fprintf(fid,'%s\n',g_database);
fprintf(fid,'%s\n',database_name);
fprintf(fid,'%s\n',file_output_loc);
fprintf(fid,'%s\n',kmer_split_method);
fprintf(fid,'%.0f\n',kmer_length);
fprintf(fid,'%.0f\n',gene_coverage);
fprintf(fid,'%.0f\n',num_genes);
if isempty(sel_genes)
    fprintf(fid,'%s\n','random selection');
else
    fprintf(fid,'%.0f\n',sel_genes);
end
fprintf(fid,'%s\n',err_mode);
fprintf(fid,'%f\n',err_rate);
fprintf(fid,'%f\n',penalty_score);
fprintf(fid,'%f\n',thresh_multiplier);
fprintf(fid,'%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\n',[thresh_prob_facts_CS thresh_prob_facts_F1 thresh_prob_facts_F2 thresh_prob_facts_F3 thresh_prob_facts_F4 thresh_prob_facts_F5 thresh_prob_facts_F6]);
fprintf(fid,'%s\n',entropy_screening_mode);
fprintf(fid,'%.0f\n',perms_thresh);
fprintf(fid,'%s\n',analysis_type);
fprintf(fid,'%.0f\n\n',tracking_level);

% Select Gene(s) ---------------------------------------------------------
fprintf(fid,'%s\n','SELECT_GENES--------------------------------------------------------------------------------');
for i=1:num_genes
    fprintf(fid,'%.0f\n',gene_selection(i));
    fprintf(fid,'%s\n',select_gene_subCLASS{i});
    if strcmpi(g_database,'resistance')
        fprintf(fid,'%s\n',select_gene_CLASS{i});
    end
    fprintf(fid,'%s\n',genes(i).Header);
    fprintf(fid,'%s\n',genes(i).Sequence);
    if strcmpi(err_mode,'on')
        fprintf(fid,'%s\n',genes_err(i).Sequence);
    end
    fprintf(fid,'\n');
end

% Select Block(s) --------------------------------------------------------
% Ordered
for i=1:num_genes
    for j=1:gene_coverage
        fprintf(fid,'%s\n',['BLOCKS_ORDERED_',num2str(i),'_',num2str(j),'--------------------------------------------------------------------------------']);
        write=[num2cell(1:length(blocks(i).Blocks{j}))' blocks(i).Blocks{j} blocks_content(i).Blocks{j}'];
        for k=1:length(write)
            fprintf(fid,'%-5.0f\t%-10s\t%-6f\t%-6f\t%-6f\t%-6f\n',write{k,:});
        end
        fprintf(fid,'\n');
    end
end
% Random
fprintf(fid,'%s\n','BLOCKS_RAND--------------------------------------------------------------------------------');
write=[num2cell((1:tot_blocks))' num2cell(blocks_content_test.BlocksCombGT)' blocks_content_test.BlocksComb' blocks_content_test.Perms blocks_content_test.BlocksContComb'];
for i=1:length(write)
    fprintf(fid,'%-5.0f\t%-4.0f\t%-10s\t%-6.0f\t%-6f\t%-6f\t%-6f\t%-6f\n',write{i,:});
end
% Final order
fprintf(fid,'\n%s\n','BLOCKS_RAND_FINAL--------------------------------------------------------------------------------');
write=[num2cell((1:tot_blocks))' num2cell(blocks_final_test.BlocksCombGT)' blocks_final_test.BlocksComb' blocks_final_test.Perms blocks_final_test.BlocksContComb'];
for i=1:length(write)
    fprintf(fid,'%-5.0f\t%-4.0f\t%-10s\t%-6.0f\t%-6f\t%-6f\t%-6f\t%-6f\n',write{i,:});
end

% Results ----------------------------------------------------------------
% Coverage
fprintf(fid,'\n%s\n','COVERAGE--------------------------------------------------------------------------------');
write=[(1:tot_blocks)' cov_tracker cov_tracker_overall];
fmtspec_cov=repmat('%-8f\t',1,num_genes);
for i=1:length(write)
    fprintf(fid,['%-5.0f\t',fmtspec_cov,'%-8f\n'],write(i,:));
end
% Remaining genes and specificity
fprintf(fid,'\n%s\n','REM_GENES_SPECIFICITY--------------------------------------------------------------------------------');
write=[(1:tot_blocks)' remaining_genes specificity_tracker];
for i=1:length(write)
    fprintf(fid,'%-5.0f\t%-8.0f\t%-8f\n',write(i,:));
end
% subCLASS/CLASS analysis
fprintf(fid,'\n%s\n','CLASS--------------------------------------------------------------------------------');
if strcmpi(g_database,'resistance')
    write=[num2cell((1:tot_blocks))' num2cell(ID_tracker) Log_TOP_eblock_subCLASS num2cell(Log_TOP_eblock_subCLASS_probsnorm) Log_TOP_eblock_CLASS num2cell(Log_TOP_eblock_CLASS_probsnorm)];
else
    write=[num2cell((1:tot_blocks))' num2cell(ID_tracker) Log_TOP_eblock_subCLASS num2cell(Log_TOP_eblock_subCLASS_probsnorm)];
end
fmtspec_ID=repmat('%-5.0f\t',1,num_genes);
fmtspec_subCLASS=repmat('%-8s\t',1,tracking_level);
fmtspec_subCLASS_pr=repmat('%-8f\t',1,tracking_level);
fmtspec_CLASS=repmat('%-65s\t',1,tracking_level);
fmtspec_CLASS_pr=repmat('%-8f\t',1,tracking_level);
if strcmpi(g_database,'resistance')
    fmtspec_CLASS_pr=fmtspec_CLASS_pr(1:length(fmtspec_CLASS_pr)-2);
else
    fmtspec_subCLASS_pr=fmtspec_subCLASS_pr(1:length(fmtspec_subCLASS_pr)-2);
end
[r_write, ~]=size(write);
if strcmpi(g_database,'resistance')
    for i=1:r_write
        fprintf(fid,['%-5.0f\t',fmtspec_ID,fmtspec_subCLASS,fmtspec_subCLASS_pr,fmtspec_CLASS,fmtspec_CLASS_pr,'\n'],write{i,:});
    end
else
    for i=1:r_write
        fprintf(fid,['%-5.0f\t',fmtspec_ID,fmtspec_subCLASS,fmtspec_subCLASS_pr,'\n'],write{i,:});
    end
end
% Factors (only output for 'benchmarking' mode
if strcmpi(analysis_type,'benchmarking')
    fmtspec_probs_all=repmat('%-10f\t',1,tot_blocks);
    fmtspec_probs_all=fmtspec_probs_all(1:length(fmtspec_probs_all)-2);
    % 1 -----------------
    fprintf(fid,'\n%s\n','FACT_1--------------------------------------------------------------------------------');
    Fact1_all_t=Fact1_all';
    for i=1:tot_genes_database
        fprintf(fid,['%-.0f\t',fmtspec_probs_all,'\n'],[i Fact1_all_t(i,:)]);
    end
    % 2 -----------------
    fprintf(fid,'\n%s\n','FACT_2--------------------------------------------------------------------------------');
    Fact2_all_t=Fact2_all';
    for i=1:tot_genes_database
        fprintf(fid,['%-.0f\t',fmtspec_probs_all,'\n'],[i Fact2_all_t(i,:)]);
    end
    % 3 -----------------
    fprintf(fid,'\n%s\n','FACT_3--------------------------------------------------------------------------------');
    Fact3_all_t=Fact3_all';
    for i=1:tot_genes_database
        fprintf(fid,['%-.0f\t',fmtspec_probs_all,'\n'],[i Fact3_all_t(i,:)]);
    end
    % 4 -----------------
    fprintf(fid,'\n%s\n','FACT_4--------------------------------------------------------------------------------');
    Fact4_all_t=Fact4_all';
    for i=1:tot_genes_database
        fprintf(fid,['%-.0f\t',fmtspec_probs_all,'\n'],[i Fact4_all_t(i,:)]);
    end
    % 5 -----------------
    fprintf(fid,'\n%s\n','FACT_5--------------------------------------------------------------------------------');
    Fact5_all_t=Fact5_all';
    for i=1:tot_genes_database
        fprintf(fid,['%-.0f\t',fmtspec_probs_all,'\n'],[i Fact5_all_t(i,:)]);
    end
    % 6 -----------------
    fprintf(fid,'\n%s\n','FACT_6--------------------------------------------------------------------------------');
    Fact6_all_t=Fact6_all';
    for i=1:tot_genes_database
        fprintf(fid,['%-.0f\t',fmtspec_probs_all,'\n'],[i Fact6_all_t(i,:)]);
    end
    % PDSlope avg -------
    fprintf(fid,'\n%s\n','PDSlope_Avg--------------------------------------------------------------------------------');
    fprintf(fid,[fmtspec_probs_all,'\n'],pdslope_eblock_avg');
    % PDSlope -----------
    fprintf(fid,'\n%s\n','PDSlope--------------------------------------------------------------------------------');
    PDSlope_all_t=PDSlope_all';
    for i=1:tot_genes_database
        fprintf(fid,['%-.0f\t',fmtspec_probs_all,'\n'],[i PDSlope_all_t(i,:)]);
    end
    % CSSlope avg -------
    fprintf(fid,'\n%s\n','CSSlope_Avg--------------------------------------------------------------------------------');
    fprintf(fid,[fmtspec_probs_all,'\n'],csslope_eblock_avg');
    % CSSlope -----------
    fprintf(fid,'\n%s\n','CSSlope--------------------------------------------------------------------------------');
    CSSlope_all_t=CSSlope_all';
    for i=1:tot_genes_database
        fprintf(fid,['%-.0f\t',fmtspec_probs_all,'\n'],[i CSSlope_all_t(i,:)]);
    end
end
% All final normalized probs
fprintf(fid,'\n%s\n','CS_ALL--------------------------------------------------------------------------------');
fmtspec_probs_all=repmat('%-10f\t',1,tot_blocks);
fmtspec_probs_all=fmtspec_probs_all(1:length(fmtspec_probs_all)-2);
CS_all_t=CS_all';
for i=1:tot_genes_database
    fprintf(fid,['%-.0f\t',fmtspec_probs_all,'\n'],[i CS_all_t(i,:)]);
end
% Probs for selected gene(s)
fprintf(fid,'\n%s\n','CS_SELECT--------------------------------------------------------------------------------');
write=[(1:tot_blocks)' CS_select];
fmtspec_probs_sel=repmat('%-10f\t',1,num_genes);
fmtspec_probs_sel=fmtspec_probs_sel(1:length(fmtspec_probs_sel)-2);
for i=1:length(write)
    fprintf(fid,['%-5.0f\t',fmtspec_probs_sel,'\n'],write(i,:));
end

% Close file -------------------------------------------------------------
fclose(fid);

%% PLOTS FOR FACTOR ANALYSIS

% Figures for Factors ----------------------------------------------------
% 1 - (Cumulative precent diff from avg) * (Cumulative normalized probs)
% 2 - Cumulative number of blocks having at least one match
% 3 - Cumulative probs multiplied together (taken as log2)
% 4 - Exponential of gene coverage
% 5 - Cumulative of the percent diff from avg slope
% 6 - Cumulative diff from avg
% ------------------------------------------------------------------------

if strcmpi(disp_fact_figs,'yes')
    
    % 1 ------------------------------------------------------------------
    figure(1)
    plot((1:tot_blocks)',Fact1_all,'k',(1:tot_blocks)',Fact1_all(:,gene_selection),'r')
    xlabel('Blocks')
    ylabel('(Cumulative precent diff from avg) * (Cumulative normalized probs)')
    title('Factor 1')
    hold on

    % 2 ------------------------------------------------------------------
    figure(2)
    plot((1:tot_blocks)',Fact2_all,'k',(1:tot_blocks)',Fact2_all(:,gene_selection),'r')
    xlabel('Blocks')
    ylabel('Cumulative number of blocks having at least one match')
    title('Factor 2')
    hold on

    % 3 ------------------------------------------------------------------
    figure(3)
    plot((1:tot_blocks)',Fact3_all,'k',(1:tot_blocks)',Fact3_all(:,gene_selection),'r')
    xlabel('Blocks')
    ylabel('Cumulative probs multiplied together (taken as log2)')
    title('Factor 3')
    hold on

    % 4 ------------------------------------------------------------------
    figure(4)
    plot((1:tot_blocks)',Fact4_all,'k',(1:tot_blocks)',Fact4_all(:,gene_selection),'r')
    xlabel('Blocks')
    ylabel('Exponential of gene coverage')
    title('Factor 4')
    hold on

    % 5 ------------------------------------------------------------------
    figure(5)
    plot((1:tot_blocks)',Fact5_all,'k',(1:tot_blocks)',Fact5_all(:,gene_selection),'r')
    xlabel('Blocks')
    ylabel('Cumulative of the percent diff from avg slope')
    title('Factor 5')
    hold on

    % 6 ------------------------------------------------------------------
    figure(6)
    plot((1:tot_blocks)',Fact6_all,'k',(1:tot_blocks)',Fact6_all(:,gene_selection),'r')
    xlabel('Blocks')
    ylabel('Cumulative diff from avg')
    title('Factor 6')
    hold on

    % Content Score (CS) ------------------------------------------------
    figure(7)
    plot((1:tot_blocks)',CS_all,'k',(1:tot_blocks)',CS_all(:,gene_selection),'r')
    xlabel('Blocks')
    ylabel('Content score')
    title('Content Score (CS)')
    hold on
    
end

%% FUNCTIONS USED THROUGHOUT

% Resistance database ----------------------------------------------------
function [subCLASS_res, CLASS_res] = ResExtract(header)
    
    inds=find(header=='|',3,'last');
    subCLASS_res=header(inds(3)+1:length(header));
    CLASS_res=header(inds(2)+1:inds(3)-1);
    if strcmpi('RequiresSNPConfirmation',subCLASS_res)
        subCLASS_res=header(inds(2)+1:inds(3)-1);
        CLASS_res=header(inds(1)+1:inds(2)-1);
    end
    
end

% Genetic database -------------------------------------------------------
function subCLASS_genetic = GeneticExtract(header)
    
    ind_s=find(header=='(',1);
    ind_e=find(header==')',1);
    subCLASS_genetic=header(ind_s+1:ind_e-1);
    
end

% Cancer database --------------------------------------------------------
function subCLASS_cancer = CancerExtract(header)
    
    ind_und=find(header=='_',1);
    ind_sp=find(header==' ',1);
    if ind_und < ind_sp
        subCLASS_cancer=header(1:ind_und-1);
    else
        subCLASS_cancer=header(1:ind_sp-1);
    end
    
end
