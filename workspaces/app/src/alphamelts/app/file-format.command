#!/usr/bin/perl -w

use strict;

my ($windows, $infile, $outfile, $savefile, @argv2);

$windows = 0;
if(exists($ENV{'OS'})) {
  if($ENV{'OS'} =~ /Windows/i) {
    $windows = 1;
  }
}

$infile = '';
@argv2 = ();
until (@argv2) {

    if (@ARGV) {

 	@argv2 = @ARGV;

    }
    else {

	print ("Usage: file-format.command list_of_files\n");
	print ("Please enter list of files (separated by spaces) now.\n");

	$_ = <STDIN>;
	chomp;

	if (/\</) {

	    $infile = $_;
	    $infile =~ s/^(.*)(\<)(\s*)(.+)$/$4/;
	    $infile =~ s/\s+$//;
	    $infile =~ s/^(\"|\')//; # Trim any leftover quotes from start and end of file name
	    $infile =~ s/(\"|\')$//;
	    # Fix any apostrophes in the name (if had single quotes or no quotes originally)
	    $infile =~ s/\'\\\'\'/\'/g || $outfile =~ s/\\\'/\'/g;
	    $infile =~ s/\\//g unless $windows; # Remove escaping from any other characters (e.g. '(' and ')')

	    unless (open(INFILE, "<$infile")) {
		warn "Error: could not open file \"$infile\" : $!\n";
		next;
	    }
	    $_ = <INFILE>;

	}

	s/\s+$//; # Trim any trailing white space
	s/\'/\"/g; # Windows can only deal with name enclosed in double quotes; *nix can take single or double
	# For a list of files inside quotes it follows that the white space is inside quotes (e.g. '" "' in '"file 1" "file 2"')
	# Otherwise split by spaces and adjust later if any are escaped (e.g. ' ' versus '\ ' in 'file\ 1 file\ 2')
	@argv2 = (/[^\\]\"/) ? split /\"\s+\"/ : split /\s+/; # ... but don't match escaped apostrophes (currently '\"')

    }

    $savefile = '';
    print "\n";
    FOREACH: foreach $infile (@argv2) {

    	$infile =~ s/^\"//; # Trim any leftover quotes from start and end of list
	$infile =~ s/\"$//;
	# Fix any apostrophes in the name (if had single quotes or no quotes or double quotes originally)
 	$infile =~ s/\"\\\"\"/\'/g || $infile =~ s/\\\"/\'/g || $infile =~ s/\"/\'/g;

        # A trailing backslash indicates an escaped space (except on Windows where it means a folder, not a file)
	$infile = $savefile.' '.$infile if ($savefile);
        $savefile = ($infile =~ s/\\$//) ? $infile : '';
        next if ($savefile);

	$infile =~ s/\\//g unless $windows; # Remove escaping from any other characters (e.g. '(' and ')')

	$outfile = $infile;
	$infile .= '~';

	# Don't die as error message will not be seen if user double-clicked
	unless (rename($outfile, $infile)) {
	    warn "Could not back up file : $!\n";
	    warn "Check that \"$outfile\" exists and \"$infile\" is not write protected or open elsewhere!\n";
	    next FOREACH;
	}
	unless (open(IN, "<$infile")) {
	    warn "Could not open file \"$infile\" : $!\n";
	    next FOREACH;
	}
	unless (open(OUT, ">$outfile")) {
	    warn "Could not open file \"$outfile\" : $!\n";
	    next FOREACH;
	}

	while (<IN>) {

	    ($windows) ? (/\r\n/g || s/\r/\r\n/g || s/\n/\r\n/g)
		: (s/\r\n/\n/g || s/\r/\n/g);
	    print OUT;

	}

	close(IN);
	close(OUT);

	print "Processed $outfile\n";

    }

    print "\n";
    foreach my $infile (@argv2) {
      (unlink $infile) || (warn "Could not delete backup file \"$infile\" : $!\n");
    }

}
continue {
 
    unless (@ARGV) {
	$infile = '';
	print ("Format more files (y or n)?\n");
	$_ = <STDIN>;
	chomp;
	# Mostly stays open if user accidently drags and drops before answering...	
	push @argv2, 'y';
	@argv2 = () unless (/^$/ || /^n/i); # shuts if no answer or begins with n/N.
    }

}

