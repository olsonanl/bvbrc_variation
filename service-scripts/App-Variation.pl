#
# The Variation Analysis application.
#

use strict;
use Carp;
use Cwd qw(abs_path getcwd);
use Data::Dumper;
use File::Temp;
use File::Copy;
use File::Basename;
use File::Slurp;
use JSON;
use P3DataAPI;
use File::Which;
use IPC::Run 'run';
use List::Util qw(min max);

#
# Arrange for the proper samtools/bcftools.
# Assume freebayes in default bin.
#

my $path_prefix;

BEGIN {
    my $rt = $ENV{KB_RUNTIME};
    $path_prefix = "$rt/samtools-1.17/bin:$rt/bcftools-1.17/bin";
    $ENV{PATH} = "$path_prefix:$ENV{PATH}";
}

#
# Find right snpeff name to use
#
my @snpEffOptions = qw(snpEff.sh snpEff);
my $snpEffPath;
for my $snpEffTest (@snpEffOptions)
{
    my $vers;
    print "Try '$snpEffTest'\n";
    my $ok = eval { run([$snpEffTest, "-version"], ">", \$vers); };
    if ($ok)
    {
	my($path) = File::Which::which($snpEffTest);
	warn "Using $path version $vers\n";
	$snpEffPath = $path;;
	last;
    }
}
$snpEffPath or die "Cannot find snpEff in $ENV{PATH}; options are @snpEffOptions\n";

use Bio::KBase::AppService::AppConfig;
use Bio::KBase::AppService::AppScript;
use Bio::KBase::AppService::ReadSet;

my $data_url = Bio::KBase::AppService::AppConfig->data_api_url;
my $script = Bio::KBase::AppService::AppScript->new(\&process_variation_data, \&preflight);

my $rc = $script->run(\@ARGV);
exit $rc;

our $global_ws;
our $global_token;

sub preflight
{
    my($app, $app_def, $raw_params, $params) = @_;



    my $token = $app->token();
    my $ws = $app->workspace();

    my $readset;
    eval {
	$readset = Bio::KBase::AppService::ReadSet->create_from_asssembly_params($params);
    };
    if ($@)
    {
	die "Error parsing assembly parameters: $@";
    }

    my($ok, $errs, $comp_size, $uncomp_size) = $readset->validate($ws);

    if (!$ok)
    {
	die "Reads as defined in parameters failed to validate. Errors:\n\t" . join("\n\t", @$errs) . "\n";
    }
    print STDERR "comp=$comp_size uncomp=$uncomp_size\n";

    my $est_comp = $comp_size + 0.75 * $uncomp_size;
    $est_comp /= 1e6;

    my $est_storage = int(1.3e6 * $est_comp / 0.75);

    my $time = 86400 * 2;

    my $pf = {
	cpu => 4,
	memory => "48000M",
	runtime => $time,
	storage => $est_storage,
	is_control_task => 0,
    };
    return $pf;
}



sub process_variation_data {
    my ($app, $app_def, $raw_params, $params) = @_;

    system("samtools", "--version") == 0 or die "samtools --version failed";
    system("bcftools", "--version") == 0 or die "bcftools --version failed";
    system("freebayes", "--version") == 0 or die "freebayes --version failed";
    system("snippy", "--version") == 0 or die "snippy --version failed";

    #
    # Redirect tmp to large NFS if more than 4 input files.
    # (HACK)
    #

    my $file_count = count_params_files($params);
    print STDERR "File count: $file_count\n";
    my $bigtmp = "/vol/patric3/tmp";
    if ($file_count > 4 && -d $bigtmp)
    {
	print STDERR "Changing tmp from $ENV{TEMPDIR} to $bigtmp\n";
	$ENV{TEMPDIR} = $ENV{TMPDIR} = $bigtmp;
    }

    print "Proc variation data ", Dumper($app_def, $raw_params, $params);
    my $time1 = `date`;

    $global_token = $app->token();
    $global_ws = $app->workspace;

    my $output_folder = $app->result_folder();

    my $cleanup = 1;
    if ($params->{debug})
    {
	print STDERR "Running in debug mode, keeping tmpdir in place after run\n";
	$cleanup = 0;
    }

    my $run_dir = getcwd();

    my $tmpdir = File::Temp->newdir( CLEANUP => $cleanup );
    
    # my $tmpdir = "/tmp/oIGe_LLBbt";
    # my $tmpdir = "/disks/tmp/var_bam";
    # my $tmpdir = "/disks/tmp/var_bam1";
    # my $tmpdir = "/disks/tmp/var_debug";

    chmod(0755, "$tmpdir");
    print STDERR "tmpdir=$tmpdir\n";

    print STDERR '$params = '. Dumper($params);
    $params = localize_params($tmpdir, $params);
    $params = parse_bam_files($tmpdir, $params);
    print STDERR '$params = '. Dumper($params);
    # exit;

    my $ref_id = $params->{reference_genome_id} or die "Reference genome is required for variation analysis\n";

    my $has_gbk = prepare_ref_data($ref_id, $tmpdir);

    my $mapper = $params->{mapper} || 'bwa_mem';
    my $caller = $params->{caller} || 'freebayes';

    my $map = "var-map";

    my $threads = $ENV{P3_ALLOCATED_CPU} // 2;

    my @basecmd = ($map, "--path-prefix", $path_prefix);
    push @basecmd, ("-a", $mapper);
    push @basecmd, ("--vc", $caller);
    push @basecmd, ("--threads", $threads);
    push @basecmd, "$tmpdir/$ref_id/$ref_id.fna";

    my $lib_txt = "$tmpdir/libs.txt";
    open(LIBS, ">$lib_txt") or die "Could not open $lib_txt";
    print LIBS "Library\tReads\n";
    my ($pe_no, $se_no);
    my @libs;
    for (@{$params->{paired_end_libs}}) {
        my $lib = "PE". ++$pe_no;
        my $outdir = "$tmpdir/$lib";
        push @libs, $lib;
        my @cmd = @basecmd;
        push @cmd, ("-o", $outdir);
        push @cmd, ($_->{read1}, $_->{read2});
        print LIBS $lib."\t".join(",", basename($_->{read1}), basename($_->{read2}))."\n";
        run_cmd(\@cmd, 1);
    }
    for (@{$params->{single_end_libs}}) {
        my $lib = "SE". ++$se_no;
        my $outdir = "$tmpdir/$lib";
        push @libs, $lib;
        my @cmd = @basecmd;
        push @cmd, ("-o", $outdir);
        push @cmd, $_->{read};
        print LIBS $lib."\t".basename($_->{read})."\n";
        run_cmd(\@cmd, 1);
    }
   # SRR is a special case; we will need to pull metadata.
    # my $total_comp_size = 0;
    for my $srr (@{$params->{srr_ids}})
    	{
	    my $tmp = File::Temp->new();
	    close($tmp);
	    # my $file_name = "$srr" . "_metadata.txt";
	    # my $tmp = "$tmpdir/$file_name";
	    print STDERR "Downloading $srr\n";
	    my $rc = system("p3-sra", "--metaonly", "--metadata-file", "$tmp", "--id", $srr);
	    if ($rc != 0)
	    	{
	    	die "p3-sra failed: $rc";
	    	}
	    else
		    {
			my $mtxt = read_file("$tmp");
			my $meta = eval { JSON::decode_json($mtxt); };
			$meta or die "Error loading or evaluating json metadata: $mtxt";
			print Dumper(MD => $meta);
			my($me) = grep { $_->{accession} eq $srr } @$meta;
			my $rc = system("p3-sra", "--out", "$tmpdir", "--id", $srr);
			# print STDERR lc $me->{library_layout} . "\n";
			if ((lc $me->{library_layout}) eq "paired") {
				my $lib = "PE". ++$pe_no;
        		my $outdir = "$tmpdir/$lib";
        		push @libs, $lib;
        		my @cmd = @basecmd;
        		push @cmd, ("-o", $outdir);
        		push @cmd, ($tmpdir . "/" . $srr . "_1.fastq", $tmpdir . "/" . $srr . "_2.fastq");
        		print LIBS $lib."\t".join(",", basename($srr . "_1.fastq"), basename($srr . "_2.fastq"))."\n";
        		run_cmd(\@cmd, 1);
			} else {
				my $lib = "SE". ++$se_no;
        		my $outdir = "$tmpdir/$lib";
        		push @libs, $lib;
				my @cmd = @basecmd;
        		push @cmd, ("-o", $outdir);
        		push @cmd, $tmpdir . "/" . $srr . ".fastq";
        		print LIBS $lib."\t".basename($srr . ".fastq")."\n";
        		run_cmd(\@cmd, 1);
			}
			# $total_comp_size += $me->{size};
		    }
	    next;
	    }

    close(LIBS);

	# All done getting alignments and variants.
    for (@libs) {
    	print STDERR "GBK existence: $has_gbk \n";
        run_snpeff($tmpdir, $ref_id, $_) if $has_gbk;
        run_var_annotate($tmpdir, $ref_id, $_) if $has_gbk;
        link_snpeff_annotate($tmpdir, $ref_id, $_) if $has_gbk;

        # system("ln -s $tmpdir/$_/aln.bam $tmpdir/$_.aln.bam") if -s "$tmpdir/$_/aln.bam";
        # system("ln -s $tmpdir/$_/aln.bam.bai $tmpdir/$_.aln.bam.bai") if -s "$tmpdir/$_/aln.bam.bai";

	    link_if_present("$tmpdir/$_/aln.bam", "$tmpdir/$_.aln.bam");
        link_if_present("$tmpdir/$_/aln.bam.bai", "$tmpdir/$_.aln.bam.bai");

        # system("cp $tmpdir/$_/var.vcf $tmpdir/$_.var.vcf") if -s "$tmpdir/$_/var.vcf";
        # system("cp $tmpdir/$_/var.vcf.gz $tmpdir/$_.var.vcf.gz") if -s "$tmpdir/$_/var.vcf.gz";
        # system("cp $tmpdir/$_/var.vcf.gz.tbi $tmpdir/$_.var.vcf.gz.tbi") if -s "$tmpdir/$_/var.vcf.gz.tbi";
        # system("cp $tmpdir/$_/var.annotated.tsv $tmpdir/$_.var.annotated.tsv") if -s "$tmpdir/$_/var.annotated.tsv";
        # system("cp $tmpdir/$_/var.annotated.raw.tsv $tmpdir/$_.var.annotated.tsv") if ! -s "$tmpdir/$_/var.annotated.tsv" && -s "$tmpdir/$_/var.annotated.raw.tsv";
        # system("cp $tmpdir/$_/var.snpEff.vcf $tmpdir/$_.var.snpEff.vcf") if -s "$tmpdir/$_/var.snpEff.vcf";

        cp_if_present("$tmpdir/$_/cov.bigwig", "$tmpdir/$_.cov.bigwig");
        cp_if_present("$tmpdir/$_/var.vcf", "$tmpdir/$_.var.vcf");
        cp_if_present("$tmpdir/$_/var.vcf.gz", "$tmpdir/$_.var.vcf.gz");
        cp_if_present("$tmpdir/$_/var.vcf.gz.tbi", "$tmpdir/$_.var.vcf.gz.tbi");
        cp_if_present("$tmpdir/$_/var.annotated.tsv", "$tmpdir/$_.var.annotated.tsv");
        cp_if_present("$tmpdir/$_/var.annotated.raw.tsv", "$tmpdir/$_.var.annotated.tsv");
        cp_if_present("$tmpdir/$_/var.snpEff.vcf", "$tmpdir/$_.var.snpEff.vcf");
        cp_if_present("$tmpdir/$_/Nonsyn_var_circular_view_input.txt", "$tmpdir/$_.Nonsyn_var_circular_view_input.txt");
        cp_if_present("$tmpdir/$_/Deletion_var_circular_view_input.txt", "$tmpdir/$_.Deletion_var_circular_view_input.txt");
        cp_if_present("$tmpdir/$_/Insertion_var_circular_view_input.txt", "$tmpdir/$_.Insertion_var_circular_view_input.txt");
        cp_if_present("$tmpdir/$_/High_var_circular_view_input.txt", "$tmpdir/$_.High_var_circular_view_input.txt");
        cp_if_present("$tmpdir/$_/Moderate_var_circular_view_input.txt", "$tmpdir/$_.Moderate_var_circular_view_input.txt");
        cp_if_present("$tmpdir/$_/Low_var_circular_view_input.txt", "$tmpdir/$_.Low_var_circular_view_input.txt");
        cp_if_present("$tmpdir/$_/Modifier_var_circular_view_input.txt", "$tmpdir/$_.Modifier_var_circular_view_input.txt");

        run(["sed", "s/^>/>$_./g"],
          "<", "$tmpdir/$_/consensus",
          ">", "$tmpdir/$_.consensus.fa") if -s "$tmpdir/$_/consensus";

        # system("cat $tmpdir/$_/consensus | sed 's/^>/>$_./g' > $tmpdir/$_.consensus.fa") if -s "$tmpdir/$_/consensus";
    }

    run_var_combine($tmpdir, \@libs);
    summarize($tmpdir, \@libs, $mapper, $caller);

    # create a combined var txt file that can be used as an input for the circular viewer
    my $sequences = parse_fasta("$tmpdir/$ref_id/$ref_id.fna");
    my %seq_length;
    for my $name (keys %$sequences) {
        $seq_length{$name} = $sequences->{$name};
        print STDERR "contig " . $name . " " . $seq_length{$name} . "\n";	    
    } 

    my $all_var = "$tmpdir/all.var.tsv";
    my $all_var_table = `cat $all_var`;
    
    print STDERR 'all_var_table = $all_var_table';  

    my $circular_file_Nonsyn = "$tmpdir/all.Nonsyn_var_circular_view_input.txt";
    my $circular_file_Deletion = "$tmpdir/all.Deletion_var_circular_view_input.txt";
    my $circular_file_Insertion = "$tmpdir/all.Insertion_var_circular_view_input.txt";
    my $circular_file_High = "$tmpdir/all.High_var_circular_view_input.txt";
    my $circular_file_Moderate = "$tmpdir/all.Moderate_var_circular_view_input.txt";
    my $circular_file_Low = "$tmpdir/all.Low_var_circular_view_input.txt";
    my $circular_file_Modifier = "$tmpdir/all.Modifier_var_circular_view_input.txt";
    open(CIR1, ">$circular_file_Nonsyn") or die "Could not open $circular_file_Nonsyn";
    open(CIR2, ">$circular_file_Deletion") or die "Could not open $circular_file_Deletion";
    open(CIR3, ">$circular_file_Insertion") or die "Could not open $circular_file_Insertion";
    open(CIR4, ">$circular_file_High") or die "Could not open $circular_file_High";
    open(CIR5, ">$circular_file_Moderate") or die "Could not open $circular_file_Moderate";
    open(CIR6, ">$circular_file_Low") or die "Could not open $circular_file_Low";
    open(CIR7, ">$circular_file_Modifier") or die "Could not open $circular_file_Modifier";
    my $len = 2000;
    my @lines = split /\n/, $all_var_table;
    for my $i (0..$#lines) {
        next if $i == 0;
        my $line = $lines[$i];
        # print STDERR "$line\n" ;
        my @col = split(/\t/, $line);
        my ($contig, $pos, $end) = ($col[1], $col[2], $col[2]+$len);
        # print STDERR "tsv file  " . $contig . " " . $pos . "\n";
        my $max_end = $seq_length{$contig};
        # print STDERR "max_end  " . $max_end . "\n";
        $end = min($end, $max_end);
        $pos = max($end - $len, 1);
        if ($col[8] =~ /Nonsyn/) {
            print CIR1 "$contig\t$pos\t$end\n";
        } elsif ($col[8] =~ /Deletion/) {
            print CIR2 "$contig\t$pos\t$end\n";
        } elsif ($col[8] =~ /Insertion/) {
            print CIR3 "$contig\t$pos\t$end\n";
        }
        if ($col[21] =~ /HIGH/) {
            print CIR4 "$contig\t$pos\t$end\n";
        } elsif ($col[21] =~ /MODERATE/) {
            print CIR5 "$contig\t$pos\t$end\n";
        } elsif ($col[21] =~ /LOW/) {
            print CIR6 "$contig\t$pos\t$end\n";
        } elsif ($col[21] =~ /MODIFIER/) {
            print CIR7 "$contig\t$pos\t$end\n";
        }
    }
    close(CIR1);    
    close(CIR2);    
    close(CIR3);    
    close(CIR4);    
    close(CIR5);    
    close(CIR6);    
    close(CIR7);     
    
    my @outputs;
    push @outputs, map { [ $_, 'txt' ] } glob("$tmpdir/*.txt");
    push @outputs, map { [ $_, 'tsv' ] } glob("$tmpdir/*.tsv");
    push @outputs, map { [ $_, 'vcf' ] } glob("$tmpdir/*.vcf");
    push @outputs, map { [ $_, 'html'] } glob("$tmpdir/*.html");
    push @outputs, map { [ $_, 'bam' ] } glob("$tmpdir/*.bam");
    push @outputs, map { [ $_, 'bigwig' ] } glob("$tmpdir/*.bigwig");
    push @outputs, map { [ $_, 'contigs' ] } glob("$tmpdir/*.consensus.fa");
    push @outputs, map { [ $_, 'unspecified' ] } glob("$tmpdir/*.tbi $tmpdir/*.vcf.gz $tmpdir/*.bam.bai");

    print STDERR '\@outputs = '. Dumper(\@outputs);
    # return @outputs;

    for (@outputs) {
        my ($ofile, $type) = @$_;
        if (-f "$ofile") {
            my $filename = basename($ofile);
            print STDERR "Output folder = $output_folder\n";
            print STDERR "Saving $ofile => $output_folder/$filename ...\n";
            if ($filename =~ /circular_view_input.txt/) {
                my $ws_path = "$output_folder/Text_Files_Circular_Viewer";
                if (!$app->workspace->exists($ws_path)) {
                    print "Create folder $ws_path\n";
                    $app->workspace->create({objects => [[$ws_path, 'folder']]});
                }
        
                $app->workspace->save_file_to_file("$ofile", {}, "$ws_path/$filename", $type, 1,
                       (-s "$ofile" > 10_000 ? 1 : 0), # use shock for larger files
                       # (-s "$ofile" > 20_000_000 ? 1 : 0), # use shock for larger files
                       $global_token);            
            } else {
                $app->workspace->save_file_to_file("$ofile", {}, "$output_folder/$filename", $type, 1,
                       (-s "$ofile" > 10_000 ? 1 : 0), # use shock for larger files
                       # (-s "$ofile" > 20_000_000 ? 1 : 0), # use shock for larger files
                       $global_token);            
            }
        } else {
            warn "Missing desired output file $ofile\n";
        }
    }

    my $time2 = `date`;
    write_output("Start: $time1"."End:   $time2", "$tmpdir/DONE");

    chdir($run_dir);
}

sub run_var_annotate {
    my ($tmpdir, $ref_id, $lib) = @_;
    my $annotate = "var-annotate";

    my $fna = "$tmpdir/$ref_id/$ref_id.fna";
    my $gff = "$tmpdir/$ref_id/$ref_id.gff";
    my $dir = "$tmpdir/$lib";
    my @cmd = split(' ', "$annotate --header $fna $gff $dir/var.snpEff.raw.vcf");
    my ($out) = run_cmd(\@cmd, 0);
    write_output($out, "$dir/var.annotated.raw.tsv");

    # create a txt file that can be used as an input for the circular viewer
    my $sequences = parse_fasta($fna);
    my %seq_length;
    for my $name (keys %$sequences) {
        $seq_length{$name} = $sequences->{$name};
        print STDERR "contig " . $name . " " . $seq_length{$name} . "\n";
    } 

    my $circular_file_Nonsyn = "$dir/Nonsyn_var_circular_view_input.txt";
    my $circular_file_Deletion = "$dir/Deletion_var_circular_view_input.txt";
    my $circular_file_Insertion = "$dir/Insertion_var_circular_view_input.txt";
    my $circular_file_High = "$dir/High_var_circular_view_input.txt";
    my $circular_file_Moderate = "$dir/Moderate_var_circular_view_input.txt";
    my $circular_file_Low = "$dir/Low_var_circular_view_input.txt";
    my $circular_file_Modifier = "$dir/Modifier_var_circular_view_input.txt";
    open(CIR1, ">$circular_file_Nonsyn") or die "Could not open $circular_file_Nonsyn";
    open(CIR2, ">$circular_file_Deletion") or die "Could not open $circular_file_Deletion";
    open(CIR3, ">$circular_file_Insertion") or die "Could not open $circular_file_Insertion";
    open(CIR4, ">$circular_file_High") or die "Could not open $circular_file_High";
    open(CIR5, ">$circular_file_Moderate") or die "Could not open $circular_file_Moderate";
    open(CIR6, ">$circular_file_Low") or die "Could not open $circular_file_Low";
    open(CIR7, ">$circular_file_Modifier") or die "Could not open $circular_file_Modifier";
    my $len = 2000;
    my @lines = split /\n/, $out;
    for my $i (0..$#lines) {
        next if $i == 0;
        my $line = $lines[$i];
        # print STDERR "$line\n" ;
        my @col = split(/\t/, $line);
        my ($contig, $pos, $end) = ($col[1], $col[2], $col[2]+$len);
        # print STDERR "tsv file  " . $contig . " " . $pos . "\n";
        my $max_end = $seq_length{$contig};
        # print STDERR "max_end  " . $max_end . "\n";
        $end = min($end, $max_end);
        $pos = max($end - $len, 1);
        if ($col[8] =~ /Nonsyn/) {
            print CIR1 "$contig\t$pos\t$end\n";
        } elsif ($col[8] =~ /Deletion/) {
            print CIR2 "$contig\t$pos\t$end\n";
        } elsif ($col[8] =~ /Insertion/) {
            print CIR3 "$contig\t$pos\t$end\n";
        }
        if ($col[21] =~ /HIGH/) {
            print CIR4 "$contig\t$pos\t$end\n";
        } elsif ($col[21] =~ /MODERATE/) {
            print CIR5 "$contig\t$pos\t$end\n";
        } elsif ($col[21] =~ /LOW/) {
            print CIR6 "$contig\t$pos\t$end\n";
        } elsif ($col[21] =~ /MODIFIER/) {
            print CIR7 "$contig\t$pos\t$end\n";
        }
    }
    close(CIR1);    
    close(CIR2);    
    close(CIR3);    
    close(CIR4);    
    close(CIR5);    
    close(CIR6);    
    close(CIR7);      
}

sub run_snpeff {
    my ($tmpdir, $ref_id, $lib) = @_;
    my $genome_name = get_genome_name($ref_id);
    my @accessions = get_genome_accessions($ref_id);
    my $dir = "$tmpdir/$lib";
    my $config = "$dir/snpEff.config";
    open(F, ">$config") or die "Could not open $config";
    print F "data.dir = $tmpdir\n\n";
    print F "$ref_id.genome : $genome_name\n";
    print F "  $ref_id.chromosomes : ". join(",", @accessions)."\n";
    close(F);
    my $here = getcwd();
    chdir($dir);
    my @cmd = ($snpEffPath, "build", "-c", $config, "-genbank", "-v", $ref_id);
    run_cmd(\@cmd, 1);
    @cmd = ($snpEffPath, "eff", "-no-downstream", "-no-upstream", "-no-utr", "-o", "vcf", "-c", $config,  $ref_id, "var.vcf");
    my ($out) = run_cmd(\@cmd, 0);
    write_output($out, "var.snpEff.raw.vcf");
    chdir($here);
}

sub link_snpeff_annotate {
    my ($tmpdir, $ref_id, $lib) = @_;

    my $url = "$data_url/genome_feature/?and(eq(genome_id,$ref_id),eq(annotation,PATRIC),or(eq(feature_type,CDS),eq(feature_type,tRNA),eq(feature_type,rRNA)))&sort(+accession,+start,+end)&select(patric_id,alt_locus_tag)&http_accept=application/json&limit(25000)";
    my $data = curl_json($url);
    my %vbi2peg = map { $_->{alt_locus_tag} => $_->{patric_id} } @$data;

    my $dir = "$tmpdir/$lib";
    my $eff_raw = "$dir/var.snpEff.raw.vcf";
    my $ann_raw = "$dir/var.annotated.raw.tsv";
    return unless -s $eff_raw && -s $ann_raw;
    my %var2eff;
    my $eff = "$dir/var.snpEff.vcf";
    my $ann = "$dir/var.annotated.tsv";
    my @lines = `cat $eff_raw`;
    open(EFF, ">$eff") or die "Could not open $eff";
    for my $line (@lines) {
        if ($line =~ /EFF=\S+\|(VBI\w+)/) {
            my $vbi = $1;
            my $peg = $vbi2peg{$vbi}; $peg =~ s/fig\|//;
            $line =~ s/\|$vbi/\|$peg/;
        }
        print EFF $line;
        my ($ctg, $pos) = split(/\t/, $line);
        if ($line =~ /EFF=(\w+)?\((\S+?)\|/) {
            $var2eff{"$ctg,$pos"} = [$1, $2]; # EFF=missense_variant(MODERATE|MISSENSE|tCt/tTt...)
        }
    }
    close(EFF);

    cp_if_present($ann_raw, $ann);
}

sub run_var_combine {
    my ($tmpdir, $libs) = @_;
    my $combine = "var-combine";
    my @files = map { -s "$tmpdir/$_/var.annotated.tsv" ? "$tmpdir/$_/var.annotated.tsv" :
                      -s "$tmpdir/$_/var.annotated.raw.tsv" ? "$tmpdir/$_/var.annotated.raw.tsv" :
                      undef } @$libs;
    my $bool = 1;
    for (@files) {
    	if (defined $_) {
    		$bool = 0;
    		last;
    	}
    }
    if (@files == 0 or $bool) {
        print STDERR "No var.annotated.[raw.]tsv files to combine.  The GBK file could be empty.\n";
        return;
    }
    my $cmd = join(" ", @files);
    sysrun("cat $cmd | $combine --header > $tmpdir/all.var.tsv");
    sysrun("cat $cmd | $combine --html > $tmpdir/all.var.html");
}

sub summarize {
    my ($tmpdir, $libs, $mapper, $caller) = @_;
    my $summary;
    my $total = 0;
    my $combined = "$tmpdir/all.var.tsv";
    if (-s $combined) {
        $total = `wc -l $combined`; chomp($total);
        $total--;
    }
    $summary .= "A total of $total variants have been found.";
    my $n_libs = scalar @$libs;
    if ($n_libs > 1) {
        my $shared = `grep "^$n_libs" $combined |wc -l`; chomp($shared);
        $summary .= " $shared of these variants are identified in all read libraries.";
    }
    $summary .= "\n\nOutput summary:\n";
    $summary .= "all.var.tsv             - combined sheet of annotated variants from all libraries\n";
    $summary .= "all.var.html            - sortable HTML table of all annotated variants\n";
    $summary .= "libs.txt                - mapping from read library to file(s)\n";
    $summary .= "<LIB>.aln.bam           - $mapper aligned reads for a read library\n";
    $summary .= "<LIB>.var.vcf           - $caller generated high quality variants for a read library\n";
    $summary .= "<LIB>.var.snpEff.vcf    - snpEff augmented VCF with variant effect prediction\n";
    $summary .= "<LIB>.var.annotated.tsv - PATRIC annotated variants in a spreadsheet\n";
    $summary .= "<LIB> denotes the name of a read library (e.g., PE1 for paired end 1, SE3 for single end 3).\n";
    $summary .= "\nResults by read library:\n";

    my $lib_txt = "$tmpdir/libs.txt";
    my @lines = `tail -n $n_libs $lib_txt`;
    my $i;
    for (@$libs) {
        my ($lib, $reads) = split(/\t/, $lines[$i]); chomp($reads);
        $i++;
        $summary .= "\n[$i] $lib ($reads)\n\n";
        my $sumfile = "$tmpdir/$lib/summary.txt";
        if (-s $sumfile) {
            $summary .= "    $_" for `cat $sumfile`;
        }
    }
    $summary .= "\n";
    write_output($summary, "$tmpdir/summary.txt");
}

sub get_genome_name {
    my ($gid) = @_;
    ($gid) = $gid =~ /(\d+\.\d+)/;
    my $url = "$data_url/genome/?eq(genome_id,$gid)&select(genome_name)&http_accept=application/json";
    my $data = curl_json($url);
    return $data->[0]->{genome_name};
}

sub get_genome_accessions {
    my ($gid) = @_;
    ($gid) = $gid =~ /(\d+\.\d+)/;
    my $url = "$data_url/genome_sequence/?eq(genome_id,$gid)&select(accession)&sort(+accession)&http_accept=application/json";
    my $data = curl_json($url);
    my @accs = map { $_->{accession} } @$data;
    wantarray ? @accs : \@accs;
}

sub prepare_ref_data {
    my ($gid, $basedir) = @_;
    $gid or die "Missing reference genome id\n";

    my $api = P3DataAPI->new();
    my @res = $api->query("genome", ["eq", "genome_id",$gid], ["select", "cds"]);

    if (!@res)
    {
	die "Could not query data api for genome $gid\n";
    }
    my $cds_count = $res[0]->{cds};
    if ($cds_count > 0)
    {
	print STDERR "Genome $gid has $cds_count patric CDSs\n";
    }

    my $dir = "$basedir/$gid";
    sysrun("mkdir -p $dir");

    my $api_url = "$data_url/genome_sequence/?eq(genome_id,$gid)&http_accept=application/sralign+dna+fasta&limit(25000)";
    my $ftp_url = "ftp://ftp.patricbrc.org/genomes/$gid/$gid.fna";

    my $url = $api_url;
    # $url = $ftp_url;
    my $out = curl_text($url);

    if (!$out)
    {
	die "Error retrieving fasta reference data for $gid\n";
    }

    # $out = break_fasta_lines($out."\n");
    $out =~ s/\n+/\n/g;
    write_output($out, "$dir/$gid.fna");

    if ($cds_count == 0)
    {
	print STDERR "Reference genome $gid has no PATRIC CDSs; skipping annotation requiring that data\n";
	return 0;
    }

    #Generate genbank file
    {
	  my $api = P3DataAPI->new();
	  my $gto = $api->gto_of($gid);
	  $gto or die "Could not retreive GTO for $gid\n";
	  $gto->destroy_to_file("$dir/$gid.gto");
	  -f "$dir/$gid.gto" or die "Could not create $dir/$gid.gto from gto\n";
    }
    sysrun("rast_export_genome", "-i", "$dir/$gid.gto", "-o", "$dir/genes.gbk", "genbank");

    my $api_url = "$data_url/genome_feature/?and(eq(genome_id,$gid),eq(annotation,PATRIC),or(eq(feature_type,CDS),eq(feature_type,tRNA),eq(feature_type,rRNA)))&sort(+accession,+start,+end)&http_accept=application/gff&limit(25000)";
    my $ftp_url = "ftp://ftp.patricbrc.org/genomes/$gid/$gid.PATRIC.gff";

    my $url = $api_url;
    # my $url = $ftp_url;
    my $out = curl_text($url);

    if (!$out)
    {
	die "Error retrieving GFF reference data for $gid\n";
    }
    $out =~ s/^accn\|//gm;

    write_output($out, "$dir/$gid.gff");

    my $has_gbk = 0;
    $has_gbk = 1 if -s "$dir/genes.gbk";

    return $has_gbk;
}

sub curl_text {
    my ($url) = @_;
    my @cmd = ("curl", curl_options(), $url);
    print STDERR join(" ", @cmd)."\n";
    my $out;
    my $ok = run(\@cmd, '>',  \$out);
    $ok or die "Error running curl @cmd: $?\n";
    return $out;
}

sub curl_json {
    my ($url) = @_;
    my $out = curl_text($url);
    my $hash = JSON::decode_json($out);
    return $hash;
}

sub curl_options {
    my @opts;
    my $token = get_token()->token;
    push(@opts, "-H", "Authorization: $token");
    push(@opts, "-H", "Content-Type: multipart/form-data");
    return @opts;
}

sub run_cmd {
    my ($cmd, $verbose) = @_;
    my ($out, $err);
    print STDERR "cmd = ", join(" ", @$cmd) . "\n\n" if $verbose;

    run($cmd, '>', \$out,
	# '2>', \$err,
	)
        or die "Error running cmd=@$cmd, stdout:\n$out\nstderr:\n$err\n";
    print STDERR "STDOUT:\n$out\n" if $verbose;
    print STDERR "STDERR:\n$err\n" if $verbose;
    return ($out, $err);
}

sub parse_bam_files {
    my ($tmpdir, $params) = @_;
    for my $lib (@{$params->{single_end_libs}}) {
        if ($lib->{read} =~ /(.*)\.bam$/) {
            my $input = $lib->{read};
            my $pe1 = "$1_R1.fq";
            my $pe2 = "$1_R2.fq";
            sysrun("samtools sort -n $input | samtools fastq -1 $pe1 -2 $pe2 -");
            if (-s $pe1 && -s $pe2) {
                push @{$params->{paired_end_libs}}, { read1 => $pe1, read2 => $pe2 };
            }
            $lib->{read} = undef;
        }
    }
    $params->{single_end_libs} = [ grep { defined($_->{read}) } @{$params->{single_end_libs}} ];
    return $params;
}

sub localize_params {
    my ($tmpdir, $params) = @_;
    for (@{$params->{paired_end_libs}}) {
        $_->{read1} = get_ws_file($tmpdir, $_->{read1}) if $_->{read1};
        $_->{read2} = get_ws_file($tmpdir, $_->{read2}) if $_->{read2};
    }
    for (@{$params->{single_end_libs}}) {
        $_->{read} = get_ws_file($tmpdir, $_->{read}) if $_->{read};
    }
    return $params;
}


sub count_params_files {
    my ($params) = @_;
    my $count = 0;
    if (ref($params->{paired_end_libs}))
    {
	$count += 2 * @{$params->{paired_end_libs}};
    }
    if (ref($params->{single_end_libs}))
    {
	$count += @{$params->{single_end_libs}};
    }
    if (ref($params->{srr_ids}))
    {
    $count += @{$params->{srr_ids}};
    }
    return $count;
}

sub get_ws {
    return $global_ws;
}

sub get_token {
    return $global_token;
}

sub get_ws_file {
    my ($tmpdir, $id) = @_;
    # return $id; # DEBUG
    my $ws = get_ws();
    my $token = get_token();

    my $base = basename($id);
    $base =~ s/[\s()]/_/g;
    my $file = "$tmpdir/$base";
    # return $file; # DEBUG

    my $fh;
    open($fh, ">", $file) or die "Cannot open $file for writing: $!";

    print STDERR "GET WS => $tmpdir $base $id\n";

    eval {
	$ws->copy_files_to_handles(1, $token, [[$id, $fh]]);
    };
    if ($@)
    {
	die "ERROR getting file $id\n$@\n";
    }
    close($fh);

    return $file;
}

sub write_output {
    my ($string, $ofile) = @_;
    open(F, ">$ofile") or die "Could not open $ofile";
    print F $string;
    close(F);
}

sub verify_cmd {
    my ($cmd) = @_;
    print STDERR "verifying executable $cmd ...\n";
    system("which $cmd >/dev/null") == 0 or die "Command not found: $cmd\n";
}

sub sysrun { system(@_) == 0 or confess("Command FAILED: ". join(" ", @_)); }

sub link_if_present
{
    my($from, $to) = @_;

    if (-s $from)
    {
	symlink($from, $to);
    }
}

sub cp_if_present
{
    my($from, $to) = @_;

    if (-s $from)
    {
	copy($from, $to);
    }
}


sub parse_fasta {
    my ($file) = @_;
    my %sequences;
    my $current_header;

    open my $fh, '<', $file or die "Could not open $file: $!";

    while (my $line = <$fh>) {
        chomp $line;
        $line =~ s/\r//g;
        $line =~ s/\t//g;
        if ($line =~ /^>(\S+)/) {
            $current_header = $1;
            $sequences{$current_header} = 0;
        }
        elsif ($current_header) {
            $sequences{$current_header} += length $line;
        }
    }
    close $fh;
    return \%sequences;
}
