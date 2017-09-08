#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;
use FindBin;
use Cwd;

use lib ("$FindBin::Bin/PerlLib");

use Getopt::Long qw(:config no_ignore_case bundling pass_through);


my $trinotate_conf_file = "$FindBin::Bin/conf.txt";

my $usage = <<__EOUSAGE__;

#################################################################################
#
#  Optional:
#
#  --autopilot         automatically run the pipeline end-to-end
#
#################################################################################


__EOUSAGE__

    ;


my $help_flag = 0;
my $AUTO_MODE = 0;


&GetOptions( 'help|h' => \$help_flag,
             'autopilot' => \$AUTO_MODE,
             
             
    );

if ($help_flag) {
    die $usage;
}

my $OS_type = `uname`;

my $workdir = cwd();

## first check for tools needed.

my @tools = qw (Trinity
    bowtie2
    bowtie
    samtools
    blastx
);


{
    my $missing_tool_flag = 0;
    foreach my $tool (@tools) {
        my $path = `which $tool`;
        unless ($path =~ /\w/) {
            print STDERR "Error, cannot find path to: $tool\n";
            $missing_tool_flag = 1;
        }
    }
    
    if ($missing_tool_flag) {
        die "\n\nTools must be in PATH setting before proceeding.\n\n";
    }

}

my $trinity_dir = dirname(`which Trinity`);


my $checkpoints_dir = $workdir . "/__TrinDemo_checkpoints_dir";
unless (-d $checkpoints_dir) {
    mkdir $checkpoints_dir or die "Error, cannot mkdir $checkpoints_dir";
}



##############
# run Trinity.


my $run_Trinity_cmd = "Trinity --seqType fq "
    . " --samples_file samples.txt "
    . " --CPU 2 --max_memory 2G "; 


&process_cmd($run_Trinity_cmd, "$checkpoints_dir/trinity.ok");

# Examine top of Trinity.fasta file
&process_cmd("head trinity_out_dir/Trinity.fasta", "$checkpoints_dir/head_trinity.ok");

# count the number of transcripts assembled.
&process_cmd("grep '>' trinity_out_dir/Trinity.fasta | wc -l ", "$checkpoints_dir/count_trans.ok");


# Get Trinity stats:
&process_cmd("$trinity_dir/util/TrinityStats.pl trinity_out_dir/Trinity.fasta", "$checkpoints_dir/trin_stats.ok");


## representation of reads by the assembly
&process_cmd("bowtie2-build trinity_out_dir/Trinity.fasta trinity_out_dir/Trinity.fasta", "$checkpoints_dir/bowtie2_build_read_assess.ok");

&process_cmd("bowtie2 --local --no-unal -x trinity_out_dir/Trinity.fasta -q -1 STND-A_1.fq.gz -2 STND-A_2.fq.gz | samtools view -Sb - | samtools sort -o - - > bowtie2.bam", "$checkpoints_dir/bowtie2_align_reads_assess.ok");


## Examine read alignments in IGV

&process_cmd("samtools index bowtie2.bam", "index_bowtie2_bam.ok");

my $igv_cmd = "igv.sh -g trinity_out_dir/Trinity.fasta bowtie2.bam";
if ($AUTO_MODE) {
    $igv_cmd .= " & ";
}
&process_cmd($igv_cmd, "$checkpoints_dir/igv_trinity_reads.ok");




###########################################
## assess number of full-length transcripts

&process_cmd("makeblastdb -in mini_sprot.pep -dbtype prot", "$checkpoints_dir/makeblastdb.ok");

&process_cmd("blastx -query trinity_out_dir/Trinity.fasta -db mini_sprot.pep -out blastx.outfmt6 -evalue 1e-20 -num_threads 2 -max_target_seqs 1 -outfmt 6", "$checkpoints_dir/blastx_for_full_length.ok");

&process_cmd("$trinity_dir/util/analyze_blastPlus_topHit_coverage.pl blastx.outfmt6 trinity_out_dir/Trinity.fasta mini_sprot.pep", "$checkpoints_dir/tophat_blast_cov_stats.ok");


####################################
## Abundance estimation using Salmon
####################################


        
my $align_estimate_command = "$trinity_dir/util/align_and_estimate_abundance.pl --seqType fq "
    . " --samples_file samples.txt "
    . " --transcripts trinity_out_dir/Trinity.fasta "
    . " --est_method salmon "
    . " --trinity_mode --prep_reference ";

&process_cmd($align_estimate_command, "$checkpoints_dir/salmon.ok");


&process_cmd("ls -1 */quant.sf | tee  quant_files.txt", "$checkpoints_dir/quant_files_list.ok");

## generate matrix of counts and perform TMM normalization
&process_cmd("$trinity_dir/util/abundance_estimates_to_matrix.pl --est_method salmon --quant_files quant_files.txt --name_sample_by_basedir" , "$checkpoints_dir/trans_matrices.ok");


## Look at the counts matrix
&process_cmd("head -n20 salmon.counts.matrix", "$checkpoints_dir/head.counts.matrix.ok");

## Look at the expression matrix:
&process_cmd("head -n20 salmon.TMM.EXPR.matrix", "$checkpoints_dir/head.expr.matrix.ok");


## Examine the E90N50 statistic
&process_cmd("$trinity_dir//util/misc/contig_ExN50_statistic.pl  salmon.TMM.EXPR.matrix trinity_out_dir/Trinity.fasta > ExN50.stats", "$checkpoints_dir/ExNstats.ok");

&process_cmd("cat ExN50.stats", "$checkpoints_dir/cat_ExNstats.ok");

## Plot the values
&process_cmd("$trinity_dir/util/misc/plot_ExN50_statistic.Rscript ExN50.stats", "$checkpoints_dir/plot_ExN50.ok");

&show("ExN50.stats.plot.pdf");


exit(0);

####
sub process_cmd {
    my ($cmd, $checkpoint) = @_;

    unless ($checkpoint) {
        die "Error, need checkpoint file defined";
    }
    
    
    if (-e $checkpoint) { return; }

    
    unless ($AUTO_MODE) {
        
        my $response = "";
        while ($response !~ /^[YN]/i) {
            print STDERR "\n\n"
                . "###############################################\n"
                . "CMD: $cmd\n"
                . "###############################################\n\n"
                . "Execute (Y/N)? ";

            $response = <STDIN>;
        }

        if ($response =~ /^N/i) {
            print STDERR "\t *** Exiting on demand. ****\n\n"
                . "Goodbye. \n\n\n";
            exit(0);
        }
    }
    
    print STDERR "\tNow running:\n\t\t$cmd\n\n\n";
    
    my $ret = system($cmd);
    if ($ret) {
        die "Error, cmd: $cmd died with ret $ret";
    }

    system("touch $checkpoint");
    
    return;
}


sub show {
    my ($image) = @_;

    my $cmd;

    if ($OS_type =~ /linux/i) {
        ## use xpdf
        $cmd = "xpdf $image";
    }
    else {
        ## rely on ImageMagick:
        $cmd = "open $image";
    }
    
    if ($AUTO_MODE) {
        $cmd .= " & ";
    }
    
    &process_cmd($cmd, "$checkpoints_dir/view." . basename($image) . ".ok");

    return;
}



####
sub get_fq_files_listings {
    my (%samples) = @_;
    
    my @left_fq_files;
    my @right_fq_files;

    foreach my $condition_fq_lists_href (values %samples) {
        
        foreach my $replicates_fq_lists_aref (values %$condition_fq_lists_href) {
            
            my ($left_fq, $right_fq) = @$replicates_fq_lists_aref;
            push (@left_fq_files, "$left_fq");
            push (@right_fq_files, "$right_fq");
        }
    }

    return(\@left_fq_files, \@right_fq_files);
}


####
sub change_dir {
    my ($dest_dir, $checkpoint) = @_;

    
    eval {
        &process_cmd("cd $dest_dir", $checkpoint);
    };
    if ($@) {
        print STDERR "\n ** Note: if you see an error message above about not being able to cd, just ignore it... it's a weird behavior of this demo script. Rest assured we've \'cd $dest_dir\' just fine.   :)\n\n";
        system("touch $checkpoint");
    }

    # now do it in the script. :)
    chdir("$dest_dir") or die "Error, could not cd to $dest_dir"; 
 
    return;
}


####
sub substitute_tokens {
    my ($cmd, $globals_href) = @_;

    my %token_templates;
    while ($cmd =~ /(\{__\S+__\})/g) {
        my $token_template = $1;
        
        $token_templates{$token_template}++;
    }

    if (%token_templates) {
        foreach my $token_template (keys %token_templates) {
            $token_template =~ /\{__(\S+)__\}/ or die "Error, not able to parse token template: $token_template";
            my $token_name = $1;

            my $replacement_val = $globals_href->{$token_name};
            unless (defined $replacement_val) {
                die "Error, unable to identify global value for token name: $token_name of cmd: $cmd";
            }
            $cmd =~ s/$token_template/$replacement_val/g;
        }
    }
    
    return($cmd);
}
