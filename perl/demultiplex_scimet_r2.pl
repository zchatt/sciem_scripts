#!/usr/bin/perl

$die = "
This script is for demultiplexing of HiSeq Runs for sciMET. Upstream processing is as follows:

# module load bcl2fastq/1.8.4
# cd /home/groups/oroaklab/fastq/170407_Illumina_sciMET_HiSeq/${run_one}/RunDir
# configureBclToFastq.pl --input-dir /home/groups/oroaklab/seq/madbum/170407_Illumina_sciMET_HiSeq/${run_one}/Data/Intensities/BaseCalls --output-dir /home/groups/oroaklab/fastq/170407_Illumina_sciMET_HiSeq/${run_one}/fastq_index_test --use-bases-mask Y*,Y*,Y*,Y* --with-failed-reads
# cd /home/groups/oroaklab/fastq/170407_Illumina_sciMET_HiSeq/${run_one}/fastq
# make -j 20
# #For Read 1
# for n in {1..8}; do cat ./Sample_lane${n}/lane${n}_NoIndex_L00${n}_R2_*.fastq.gz > Undetermined.Lane${n}.R2.fq.gz; done &
# #For Index 1
# for n in {1..8};  do cat ./Sample_lane${n}/lane${n}_NoIndex_L00${n}_R2_*.fastq.gz  > Undetermined.Lane${n}.I1.fq.gz; done &
# #For Index 2
# for n in {1..8};  do cat ./Sample_lane${n}/lane${n}_NoIndex_L00${n}_R3_*.fastq.gz > Undetermined.Lane${n}.I2.fq.gz; done &
# #For Index 3
# for n in {1..8};  do cat ./Sample_lane${n}/lane${n}_NoIndex_L00${n}_R4_*.fastq.gz > Undetermined.Lane${n}.I3.fq.gz; done &


ARGV0 = directory with Undetermined fastq files generated from HiSeq Runs with bcl2fastq version 1.8.4 with following commands:
      OR colon separated fastq files R2:I1:I2:I3
         with \"R2\" etc being a comma separated list of fastq files
ARGV1 = index file
ARGV2 = output prefix

Will split with a hamming distance of 2

NOTE: for methylCPT, with 3 indexes, format:
   I1 = index read 1, only index in read, PCR i7 index
   I2 = index read 2, methylNext i5 index
   I3 = index read 2, concatenated after I2, PCR i5 index

";

if (!defined $ARGV[2]) {die $die}

open IN, "$ARGV[1]";
while ($l = <IN>) {
        chomp $l;
        ($id,$pos,$seq) = split(/\t/, $l);
        $POS_SEQ_seq{$pos}{$seq} = $seq;
        $POS_length{$pos} = length($seq);
        } close IN;

# make all-1-away hash
foreach $pos (keys %POS_SEQ_seq) {
        foreach $seq (keys %{$POS_SEQ_seq{$pos}}) {
                @TRUE = split(//, $seq);
                for ($i = 0; $i < @TRUE; $i++) {
                        if ($TRUE[$i] =~ /A/i) {
                                @NEW = @TRUE; $NEW[$i] = "C"; $new = join("", @NEW); $POS_SEQ_seq{$pos}{$new} = $seq;
                                @NEW = @TRUE; $NEW[$i] = "G"; $new = join("", @NEW); $POS_SEQ_seq{$pos}{$new} = $seq;
                                @NEW = @TRUE; $NEW[$i] = "T"; $new = join("", @NEW); $POS_SEQ_seq{$pos}{$new} = $seq;
                                @NEW = @TRUE; $NEW[$i] = "N"; $new = join("", @NEW); $POS_SEQ_seq{$pos}{$new} = $seq;
                        } elsif ($TRUE[$i] =~ /C/i) {
                                @NEW = @TRUE; $NEW[$i] = "A"; $new = join("", @NEW); $POS_SEQ_seq{$pos}{$new} = $seq;
                                @NEW = @TRUE; $NEW[$i] = "G"; $new = join("", @NEW); $POS_SEQ_seq{$pos}{$new} = $seq;
                                @NEW = @TRUE; $NEW[$i] = "T"; $new = join("", @NEW); $POS_SEQ_seq{$pos}{$new} = $seq;
                                @NEW = @TRUE; $NEW[$i] = "N"; $new = join("", @NEW); $POS_SEQ_seq{$pos}{$new} = $seq;
                        } elsif ($TRUE[$i] =~ /G/i) {
                                @NEW = @TRUE; $NEW[$i] = "C"; $new = join("", @NEW); $POS_SEQ_seq{$pos}{$new} = $seq;
                                @NEW = @TRUE; $NEW[$i] = "A"; $new = join("", @NEW); $POS_SEQ_seq{$pos}{$new} = $seq;
                                @NEW = @TRUE; $NEW[$i] = "T"; $new = join("", @NEW); $POS_SEQ_seq{$pos}{$new} = $seq;
                                @NEW = @TRUE; $NEW[$i] = "N"; $new = join("", @NEW); $POS_SEQ_seq{$pos}{$new} = $seq;
                        } elsif ($TRUE[$i] =~ /T/i) {
                                @NEW = @TRUE; $NEW[$i] = "C"; $new = join("", @NEW); $POS_SEQ_seq{$pos}{$new} = $seq;
                                @NEW = @TRUE; $NEW[$i] = "G"; $new = join("", @NEW); $POS_SEQ_seq{$pos}{$new} = $seq;
                                @NEW = @TRUE; $NEW[$i] = "A"; $new = join("", @NEW); $POS_SEQ_seq{$pos}{$new} = $seq;
                                @NEW = @TRUE; $NEW[$i] = "N"; $new = join("", @NEW); $POS_SEQ_seq{$pos}{$new} = $seq;
                        }
                }
        }
}

# make all-2-away hash
foreach $pos (keys %POS_SEQ_seq) {
        foreach $id_seq (keys %{$POS_SEQ_seq{$pos}}) {
                $seq = $POS_SEQ_seq{$pos}{$id_seq};
                @TRUE = split(//, $seq);
                for ($i = 0; $i < @TRUE; $i++) {
                        if ($TRUE[$i] =~ /A/i) {
                                @NEW = @TRUE; $NEW[$i] = "C"; $new = join("", @NEW); if (!defined $POS_SEQ_seq{$pos}{$new}) {$POS_SEQ_seq{$pos}{$new} = $seq};
                                @NEW = @TRUE; $NEW[$i] = "G"; $new = join("", @NEW); if (!defined $POS_SEQ_seq{$pos}{$new}) {$POS_SEQ_seq{$pos}{$new} = $seq};
                                @NEW = @TRUE; $NEW[$i] = "T"; $new = join("", @NEW); if (!defined $POS_SEQ_seq{$pos}{$new}) {$POS_SEQ_seq{$pos}{$new} = $seq};
                                @NEW = @TRUE; $NEW[$i] = "N"; $new = join("", @NEW); if (!defined $POS_SEQ_seq{$pos}{$new}) {$POS_SEQ_seq{$pos}{$new} = $seq};
                        } elsif ($TRUE[$i] =~ /C/i) {
                                @NEW = @TRUE; $NEW[$i] = "A"; $new = join("", @NEW); if (!defined $POS_SEQ_seq{$pos}{$new}) {$POS_SEQ_seq{$pos}{$new} = $seq};
                                @NEW = @TRUE; $NEW[$i] = "G"; $new = join("", @NEW); if (!defined $POS_SEQ_seq{$pos}{$new}) {$POS_SEQ_seq{$pos}{$new} = $seq};
                                @NEW = @TRUE; $NEW[$i] = "T"; $new = join("", @NEW); if (!defined $POS_SEQ_seq{$pos}{$new}) {$POS_SEQ_seq{$pos}{$new} = $seq};
                                @NEW = @TRUE; $NEW[$i] = "N"; $new = join("", @NEW); if (!defined $POS_SEQ_seq{$pos}{$new}) {$POS_SEQ_seq{$pos}{$new} = $seq};
                        } elsif ($TRUE[$i] =~ /G/i) {
                                @NEW = @TRUE; $NEW[$i] = "C"; $new = join("", @NEW); if (!defined $POS_SEQ_seq{$pos}{$new}) {$POS_SEQ_seq{$pos}{$new} = $seq};
                                @NEW = @TRUE; $NEW[$i] = "A"; $new = join("", @NEW); if (!defined $POS_SEQ_seq{$pos}{$new}) {$POS_SEQ_seq{$pos}{$new} = $seq};
                                @NEW = @TRUE; $NEW[$i] = "T"; $new = join("", @NEW); if (!defined $POS_SEQ_seq{$pos}{$new}) {$POS_SEQ_seq{$pos}{$new} = $seq};
                                @NEW = @TRUE; $NEW[$i] = "N"; $new = join("", @NEW); if (!defined $POS_SEQ_seq{$pos}{$new}) {$POS_SEQ_seq{$pos}{$new} = $seq};
                        } elsif ($TRUE[$i] =~ /T/i) {
                                @NEW = @TRUE; $NEW[$i] = "C"; $new = join("", @NEW); if (!defined $POS_SEQ_seq{$pos}{$new}) {$POS_SEQ_seq{$pos}{$new} = $seq};
                                @NEW = @TRUE; $NEW[$i] = "G"; $new = join("", @NEW); if (!defined $POS_SEQ_seq{$pos}{$new}) {$POS_SEQ_seq{$pos}{$new} = $seq};
                                @NEW = @TRUE; $NEW[$i] = "A"; $new = join("", @NEW); if (!defined $POS_SEQ_seq{$pos}{$new}) {$POS_SEQ_seq{$pos}{$new} = $seq};
                                @NEW = @TRUE; $NEW[$i] = "N"; $new = join("", @NEW); if (!defined $POS_SEQ_seq{$pos}{$new}) {$POS_SEQ_seq{$pos}{$new} = $seq};
                        }
                }
        }
}


if ($ARGV[0] =~ /:/) {
    ($r2list,$i1list, $i2list, $i3list) = split(/:/, $ARGV[0]);
    @R2 = split(/,/, $r2list);
    @I1 = split(/,/, $i1list);
    @I2 = split(/,/, $i2list);
    @I3 = split(/,/, $i3list);

} else {
        $ARGV[0] =~ s/\/$//;
        @R2 = ("$ARGV[0]/Undetermined_S0_L00${num}_R2_001.fastq.gz");
        @I1 = ("$ARGV[0]/Undetermined_S0_L00${num}_I1_001.fastq.gz");
        @I2 = ("$ARGV[0]/Undetermined_S0_L00${num}_I2_001.fastq.gz");
        @I3 = ("$ARGV[0]/Undetermined_S0_L00${num}_I3_001.fastq.gz");}

open R2OUT, "| gzip > $ARGV[2].L00${num}.1.fq.gz";
open R2FAIL, "| gzip > $ARGV[2].L00${num}.fail.1.fq.gz";

$totalCT = 0; $failCT = 0;

for ($i = 0; $i < @R2; $i++) {
        open R2, "zcat $R2[$i] |";
        open I1, "zcat $I1[$i] |";
        open I2, "zcat $I2[$i] |";
        open I3, "zcat $I3[$i] |";

        while ($i1tag = <I1>){
                $i2tag = <I2>;
                $i3tag = <I3>;

                chomp $i1tag; $i1seq = <I1>; chomp $i1seq;

                chomp $i2tag; $i2seq = <I2>; chomp $i2seq;

                chomp $i3tag; $i3seq = <I3>; chomp $i3seq;

                $null = <R2>;$r2seq = <R2>; chomp $r2seq;
                $null = <R2>;$r2qual = <R2>; chomp $r2qual;

                $null = <I1>; $null = <I2>;$null = <I3>;
                $null = <I1>; $null = <I2>;$null = <I3>;
                #need to offset line reads to keep them all a the same pace with $null


                $ix1 = $i1seq;
                $ix2 = $i2seq;
                $ix3 = $i3seq;

                #print STDERR "$ix1\t$ix2\t$ix3\n";

                if (defined $POS_SEQ_seq{'1'}{$ix1} &&
                        defined $POS_SEQ_seq{'2'}{$ix2} &&
                        defined $POS_SEQ_seq{'3'}{$ix3}) {

                        $barc = $POS_SEQ_seq{'1'}{$ix1}.$POS_SEQ_seq{'2'}{$ix2}.$POS_SEQ_seq{'3'}{$ix3};

                        #print R2OUT "\@$ix1:$ix2:$ix3\n";

                        $totalCT++;

                        print R2OUT "\@$barc:$totalCT#0/1\n$r2seq\n\+\n$r2qual\n";

                } else {

                        $barc = $ix1.$ix2.$ix3;

                        #print R2FAIL "\@$ix1:$ix2:$ix3\n";

                        print R2FAIL "\@$barc:F_$failCT#0/1\n$r2seq\n\+\n$r2qual\n";

                        $failCT++;
                }

        }}

        close R2; close I1; close I2; close I3;
