#!/usr/bin/env perl
# 
# Creates HHblits databases from HMMER-files, HMM-files or A3M-files

################################################################################################################################
# Update the following variables

use lib $ENV{HHLIB}."/scripts"
use HHPaths;

################################################################################################################################

$|= 1; # Activate autoflushing on STDOUT

# Default values:
our $v=2;              # verbose mode
my $a3mext = "a3m";       # default A3M-file extension
my $hhmext = "hhm";       # default HHM-file extension
my $x = 0.3;
my $c = 4;

my $append = 0;

my $help="
Creates HHblits databases from HMMER-files, HHM-files or A3M-files.

The recommended way to use this script is to start with a directory
of A3M-files (-a3mdir <DIR>) and let this script generates an A3M-
database (-oa3m <FILE>) and an HHM-database (-ohhm <FILE>).
If you already have HHM-models for your A3M-files, you can use them 
as additional input (-hhmdir <DIR>).
If you don't need the A3M-database, you can also start this script
with an directory of HHM-files (-hhmdir <DIR>) and as output only
the HHM-database (-ohhm <FILE>).

Usage: perl create_db.pl -i <dir> [options]

Options:
  -a3mdir <dir>  Input directory (directories) with A3M-files
  -hhmdir <dir>  Input directory (directories) with HHM- or HMMER-files 
                 (WARNING! Using HMMER databases could result in a decreased sensitivity!)

  -oa3m  <FILE>  Output basename for the A3M database (output will be BASENAME_a3m_db)
                 (if not given, no A3M database will be build)
  -ohhm  <FILE>  Output basename for the HHM database (output will be BASENAME_hhm_db)
                 (if not given, no HHM database will be build)

  -a3mext        Extension of A3M-files (default: $a3mext)
  -hhmext        Extension of HHM- or HMMER-files (default: $hhmext)

  -append        If the output file exists, append new files (default: overwrite)

  -v [0-5]       verbose mode (default: $v)

Examples:

   perl create_db.pl -a3mdir /databases/scop_a3ms -oa3m /databases/scop -ohhm /databases/scop

   perl create_db.pl -a3mdir /databases/scop_a3ms -hhmdir /databases/scop_hhms -oa3m /databases/scop -ohhm /databases/scop

   perl create_db.pl -hhmdir /databases/scop_hhms -ohhm /databases/scop
\n";

# Variable declarations
my $line;
my $command;
my $a3mdir  = "";
my $hhmdir  = "";
my $a3mfile = "";
my $hhmfile = "";
my @a3mfiles;
my @hhmfiles;

my $file;
my @files;
my $dir;
my @dirs;

###############################################################################################
# Processing command line input
###############################################################################################

if (@ARGV<1) {die ($help);}

for (my $i=0; $i<@ARGV; $i++) {

    if ($ARGV[$i] eq "-a3mdir") {
	if (++$i<@ARGV) {
	    $a3mdir=$ARGV[$i];
	} else {
	    die ("$help\n\nERROR! Missing directory after -a3mdir option!\n");
	}
    } elsif ($ARGV[$i] eq "-hhmdir") {
	if (++$i<@ARGV) {
	    $hhmdir=$ARGV[$i];
	} else {
	    die ("$help\n\nERROR! Missing directory after -hhmdir option!\n");
	}
    } elsif ($ARGV[$i] eq "-oa3m") {
	if (++$i<@ARGV) {
	    $a3mfile=$ARGV[$i] . "_a3m_db";
	} else {
	    die ("$help\n\nERROR! Missing filename after -oa3m option!\n");
	}
    } elsif ($ARGV[$i] eq "-ohhm") {
	if (++$i<@ARGV) {
	    $hhmfile=$ARGV[$i] . "_hhm_db";
	} else {
	    die ("$help\n\nERROR! Missing filename after -ohhm option!\n");
	}
    } elsif ($ARGV[$i] eq "-a3mext") {
	if (++$i<@ARGV) {
	    $a3mext=$ARGV[$i];
	} else {
	    die ("$help\n\nERROR! Missing extension after -a3mext option!\n");
	}
    } elsif ($ARGV[$i] eq "-hhmext") {
	if (++$i<@ARGV) {
	    $hhmext=$ARGV[$i];
	} else {
	    die ("$help\n\nERROR! Missing extension after -hhmext option!\n");
	}
    } elsif ($ARGV[$i] eq "-v") {
	if (++$i<@ARGV) {
	    $v=$ARGV[$i];
	} else {
	    $v = 2;
	}
    } elsif ($ARGV[$i] eq "-append") {
	$append=1;
    } else {
	print "WARNING! Unknown option $ARGV[$i]!\n";
    }
}

# Check input
if ($a3mdir eq "" && $hhmdir eq "") {
    print($help); print "ERROR! At least one input directory must be given!\n"; exit(1);
}

if ($a3mfile eq "" && $hhmfile eq "") {
    print($help); print "ERROR! At least one database (-ao3m or -ohhm) must be created!\n"; exit(1);
}

if ($a3mfile ne "") {
    if ($a3mdir eq "") {
	print($help); print "ERROR! Input directory with A3M-files needed for A3M database!\n"; exit(1);
    }
}

# Create tmp directory (plus path, if necessary)
my $tmpdir="/tmp/$ENV{USER}/$$";  # directory where all temporary files are written: /tmp/UID/PID
my $suffix=$tmpdir;
while ($suffix=~s/^\/[^\/]+//) {
    $tmpdir=~/(.*)$suffix/;
    if (!-d $1) {mkdir($1,0777);}
} 

##############################################################################################
# Main part
##############################################################################################

if ($a3mfile ne "") {
    print "Creating A3M database $a3mfile ...\n";

    open (OUT, ">$tmpdir/a3m.filelist");

    @dirs = glob($a3mdir);
    foreach $dir (@dirs) {
	my @a3mfiles = glob("$dir/*.$a3mext");
	foreach my $file (@a3mfiles) {
	    print OUT "$file\n";
	}
    }

    close OUT;

    if ($append) {
	$command = "ffindex_build -as -f $tmpdir/a3m.filelist $a3mfile $a3mfile.index";
    } else {
	$command = "ffindex_build -s -f $tmpdir/a3m.filelist $a3mfile $a3mfile.index";
    }
    if (&System($command) != 0) {
	print "WARNING! Error with command $command!\n";
    }
} 
if ($hhmfile ne "") {

    open (OUT, ">$tmpdir/hhm.filelist");

    if ($hhmdir eq "") {
	# Build HHMs from A3Ms
	print "Generate HHMs from A3Ms ...\n";
	$hhmdir = "$tmpdir/hhms";
	$hhmext = "hhm";
	&System("mkdir $hhmdir");
	
	@dirs = glob($a3mdir);
	foreach $dir (@dirs) {
	    @files = glob("$dir/*.$a3mext");
	    foreach $file (@files) {
		$file =~ /^\S+\/(\S+?)\.$a3mext$/;
		$command = "hhmake -i $file -o $hhmdir/$1.hhm";
		if (&System($command) != 0) {
		    print "WARNING! Error with command $command!\n";
		}
		print OUT "$hhmdir/$1.hhm\n";
	    }
	}
    } else {
	@dirs = glob($hhmdir);
	foreach $dir (@dirs) {
	    my @hhmfiles = glob("$dir/*.$hhmext");
	    foreach my $file (@hhmfiles) {
		print OUT "$file\n";
	    }
	}
    }
    
    close OUT;

    print "Creating HHM database $hhmfile ...\n";

    if ($append) {
	$command = "ffindex_build -as -f $tmpdir/hhm.filelist $hhmfile $hhmfile.index";
    } else {
	$command = "ffindex_build -s -f $tmpdir/hhm.filelist $hhmfile $hhmfile.index";
    }
    if (&System($command) != 0) {
	print "WARNING! Error with command $command!\n";
    }

}

if ($v < 4) {
    $command = "rm -rf $tmpdir";
    &System($command);
}

exit;


################################################################################################
### System command with return value parsed from output
################################################################################################
sub System()
{
    if ($v>2) {printf("\$ %s\n",$_[0]);} 
    return system($_[0])/256;
}
