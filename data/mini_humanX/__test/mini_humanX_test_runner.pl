#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib("$FindBin::Bin/../../__lib");
use Pipeliner;

main: {

    my $pipeliner = new Pipeliner(-verbose => 2);

    ## run hisat
    $pipeliner->add_commands(new Command("hisat2_extract_splice_sites.py minigenome.gtf > minigenome.gtf.ss", "hisat_ss.ok"));

    $pipeliner->add_commands(new Command("hisat2_extract_exons.py  minigenome.gtf > minigenome.gtf.exons", "hisat_exons.ok"));

    $pipeliner->add_commands(new Command("hisat2-build --exon minigenome.gtf.exons --ss minigenome.gtf.ss minigenome.fa minigenome.fa",
                                         "hisat_build.ok"));

    $pipeliner->add_commands(new Command("hisat2 --dta -x minigenome.fa --max-intronlen 5000 -1 reads_1.fq.gz -2 reads_2.fq.gz | samtools view -Sb | samtools sort -o alignments.hisat2.bam", "hisat_align.ok"));

    ## run stringtie to reconstruct transcripts
    $pipeliner->add_commands(new Command("stringtie -o stringtie.gtf alignments.hisat2.bam", "stringtie.ok"));

    
    $pipeliner->add_commands(new Command("Trinity --genome_guided_bam alignments.hisat2.bam --max_memory 1G --genome_guided_max_intron 5000 --CPU 2", "trinityGG.ok"));
    

    $pipeliner->run();

    exit(0);
}


    
    
