#!/usr/bin/perl -w

#use strict;

my ($savefile, @argv2, $path, $indelim, $outdelim, $padding, $header, $file, $table, $tablename, @index);
my (@labels, $nlabels, @columns, $ncolumns, %values, @output, $lineno, @lines, $windows, $infile, $outfile);

$| = 1; # print buffer immediately
$indelim = "\\s"; # default delimiter is space
$outdelim = "\t"; # default is tab delimited
$padding = ''; 
$header = 'excel'; # default format is two header lines

$file = $table = $tablename = '';
@index = ();

$lineno = 2;
@lines = ();

$windows = 0;
if(exists($ENV{'OS'})) {
  if($ENV{'OS'} =~ /Windows/i) {
    $windows = 1;
  }
}

$infile = '';
$outfile = '';
until ($outfile && !$infile) {

    if (@ARGV) {

 	@argv2 = @ARGV;
	$_ = $argv2[$#argv2];
	if (/\>./) {
	    $outfile = '>'; # Redirect to Matlab
	    pop @argv2;
	}
	else {
	    $outfile = @ARGV;
	}

    }
    else {

	if ($infile) {
	    $_ = <INFILE>;
	}
	else {
	    print ("Usage: column_pick.pl column_list_files > table_file\n");
	    print ("Please enter column_list_files, separated by spaces, and / or\n");
	    print ("type '>' and the table_file name to finish and write output.\n");
	    print ("First processed table must have an entry for each recorded step.\n");
	    $_ = <STDIN>;
	}
	
	chomp;

	unless (/.+/) { # Reached empty line or end of file
	    $infile = '';
	    next;
	}

	if (/\>/) {
	    $outfile = $_;
	    $outfile =~ s/^(.*)(\>)(\s*)(.+)$/$4/;
	    $outfile =~ s/\s+$//;
	}
	s/^(.*)(\>)(.*)$/$1/; # Trim output file name and any trailing white space

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
		$infile = '';
		next;
	    }

	}

	s/\s+$//;
	s/\'/\"/g; # Windows can only deal with name enclosed in double quotes; *nix can take single or double
	# For a list of files inside quotes it follows that the white space is inside quotes (e.g. '" "' in '"file 1" "file 2"')
	# Otherwise split by spaces and adjust later if any are escaped (e.g. ' ' versus '\ ' in 'file\ 1 file\ 2')
	@argv2 = (/[^\\]\"/) ? split /\"\s+\"/ : split /\s+/; # ... but don't match escaped apostrophes (currently '\"')

    }

  FOREACH: foreach my $colfile (@argv2) {

      my $thisfile = '';

      $colfile =~ s/^\"//; # Trim any leftover quotes from start and end of list
      $colfile =~ s/\"$//;
      # Fix any apostrophes in the name (if had single quotes or no quotes or double quotes originally)
      $colfile =~ s/\"\\\"\"/\'/g || $colfile =~ s/\\\"/\'/g || $colfile =~ s/\"/\'/g;

      # A trailing backslash indicates an escaped space (except on Windows where it means a folder, not a file)
      $colfile = $savefile.' '.$colfile if ($savefile);
      $savefile = ($colfile =~ s/\\$//) ? $colfile : '';
      next if ($savefile);

      $colfile =~ s/\\//g unless $windows; # Remove escaping from any other characters (e.g. '(' and ')')

      unless (open(COLUMNLIST, "<$colfile") || ($outfile eq '>')) {
	  warn "Error: could not open file \"$colfile\" : $!\n";
	  last FOREACH;
      }
      ($outfile eq '>') || warn "*** \"$colfile\" ***\n\n";

      $path = '';

      MAIN: while(<COLUMNLIST>) {

	  # Process the column_list_file triplets
	  s/!.*//;
	  chomp;
	  if(/^Name/) {
	      if (!$thisfile) {
		  $path = $colfile;
		  $colfile = $_;
		  $colfile =~ s/Name(:*)(\s+)//;
		  $colfile =~ s/^[\'\"]//; # Trim any quotes from start and end
		  $colfile =~ s/[\'\"]$//;
		  ($outfile eq '>') || warn ("Column file name set to: \"$colfile\"\n");
		  $thisfile = $colfile;
	      }
	  }
	  elsif(/^Header/i) {

	      if (/(excel|matlab)/i) {
		  $header = lc ((split)[1]);
	      }
	      else {
		  s/^Header:\s+//i;
		  $header = $_;
	      }
	      $header =~ s/[\'\"]//g; # strip quotes
	      ($outfile eq '>') || warn ("Header type set to: '$header'\n\n");
	      
	  }
	  elsif(/^Format/i) {

	      $outdelim = (split)[1];
	      $outdelim =~ s/[\'\"]//g; # strip quotes
	      if(lc $outdelim eq 'txt') {
		  ($outfile eq '>') || warn ("Output format set to: 'txt' (tab delimited)\n\n");
		  $outdelim = "\t";
	      }
	      elsif(lc $outdelim eq 'prn') {
		  ($outfile eq '>') || warn ("Output format set to: 'prn' (space delimited)\n\n");
		  $outdelim = ' ';
	      }
	      elsif(lc $outdelim eq 'csv') {
		  ($outfile eq '>') || warn ("Output format set to: 'csv' (comma delimited)\n\n");
		  $outdelim = ',' ;
	      }

	  }
	  elsif(/^Delimiter/i) {

	      $indelim = (split)[1];
	      $indelim =~ s/[\'\"]//g; # strip quotes
	      if(lc $indelim eq 'space') {
		  $indelim = ' ';
	      }
	      elsif(lc $indelim eq 'comma') {
		  $indelim = ',' ;
	      }
	      ($outfile eq '>') || warn ("Delimiter set to: '$indelim'\n\n");

	  }
	  elsif(/^Padding/i) {

	      $padding = (split)[1];
	      $padding =~ s/[\'\"]//g; # strip quotes
	      ($outfile eq '>') || warn ("Padding set to: '$padding'\n\n");

	  }
	  elsif(/^File/i) {
	      
	      s/^File(s*)(\:*)(\s+)//; # The alphaMELTS (or other) output file
	      s/^[\'\"]//; # Trim any quotes from start and end
	      s/[\'\"]$//;
	      $file = $_;

	      if ($path) {
		  $_ = $path;
		  s/$colfile/$file/;
		  while (/\.\./) {
		      ($windows) ? (s/\\[^\.]+\\\.\.\\/\\/ || s/[^\.]+\\\.\.\\//) : 
			  (s/\/[^\.]+\/\.\.\//\// || s/[^\.]+\/\.\.\///);
		  }
		  $file = $_;
	      }
	      ($outfile eq '>') || warn ("Input file set to: \"$file\"\n");

	  }
	  elsif(/^Table/i) {

	      s/[a-z\:\s]+[\'\"]//i; # Strip up to first quote, if present
	      $table = (/[\'\"]/) ? $_ : (split)[1]; # First word(s) of the table name in that file
	      $table =~ s/[\'\"]//g; # strip the other quote
	      ($outfile eq '>') || warn ("Table name set to: '$table'\n");

	  }
	  elsif(/^Column/i) {

	      s/^Column([s:]*)\s*//ig; # remove the 'Column:' part

	      # in future strip quoted spaces in here
	      # strip double quotes	      
	      s/\'//g; # strip single quotes for now
	      @columns = split /$indelim+/; # list of column names

	      # Process text file to pick up matching data
	      unless (open(FILE, "<$file") || ($outfile eq '>')){ 
		  warn ("Error: could not open file \"$file\" at line $. : $!\n");
		  last FOREACH;
	      }
	      ($outfile eq '>') || warn "... now reading \"$file\" to look for columns\n"; 

	    COLUMNS: while(<FILE>) {

		chomp;
		if ($table && !$tablename) {

		    # Skip forward to the table whose title begins with the given word(s),
		    # with or without trailing underscore zero for first appearance of phases.
		    next unless (/^$table(\_0|1)/i || ("$_\_0" =~ /^$table/i) || ("$_\1" =~ /^$table/i) || /^$table/i);

		    my $i = @columns;
		    if ($header eq 'excel') {
			$tablename = (/^$table(\_0|1)/i) ? $table : $_; # for output print '_0' only if given
			
			# encase full table name in quotes so it will be imported into a single spreadsheet cell
			$lines[0] .= "\"$tablename\"";
		    }
		    else {
			$_ = $table;
			$tablename = (split)[0]; # use just first word of table name in printed columnn names
		    }
		    $lines[0] .= "$padding$outdelim" x $i; # pad so that next table name will be correctly aligned

		}
		elsif (/^\s*[0-9-]/) {

		    # For alphaMELTS output, no column names begin with a number or minus sign,
		    # whereas a line of data will always have at least one number or '---'
		    my $line = $_;
		    my $i = 0;
		    # Can eventually have a delimiter for split as well...
		    map {$values{$labels[$i++]} = $_} split /$indelim+/; # one value for each label (even if not matched in columns)
		    $values{$labels[$i-1]} = '' if ($line =~ /---/); # put in empty value instead of '---'

		    my $match = '';
		    until ($match) { # evaluate loop at least once

			push @index, $values{'index'} if (@index <= $lineno); # first time through
			$match = ($values{'index'} eq $index[$lineno]) if exists $values{'index'};
			$match = ($values{'Index'} eq $index[$lineno]) if exists $values{'Index'};

			if($match) {
			    for (my $i = 0; $i < $ncolumns; $i++) { # all columns, including the missing ones
				push @output, $values{$columns[$i]};
			    }			   
			    for (my $i = 0; $i < @columns - $ncolumns; $i++) {
				my $res = eval "$labels[$nlabels+$i];";				
				warn ("Invalid expressions: $!") if ($@);
				
				$values{$columns[$ncolumns+$i]} = $res;
				push @output, $res;
			    }
			    $lines[$lineno++] .= (join $outdelim, @output).$outdelim;

			    #map {warn "$_ $values{$_}\n";} keys %values;
			    
			    foreach my $column (@columns) { # all columns, including the missing ones
				$values{$column} = ''; # reinitialise (will stay empty for missing columns and if have '---')
			    }
			    
			}
			else {
			    # pad those lines for which there are variables set but no match
			    my $i = @columns; # number of columns
			    $lines[$lineno++] .= "$padding$outdelim" x $i;
			}

		    }

		}
		elsif(/.+/) {

		    # The line is not empty so, by elimination, it must be the column names in the table
		    @labels = split /$indelim+/;
		    $nlabels = @labels;
		    
		    my @columns2 = @columns;
		    map {s/(.+)(=.*)/$1/} @columns2;
		    
		    if (!$table && !$tablename) { # fudge for GUI-style .tbl files
			my $i = @columns2;
			$lines[0] .= "\"$file\""; # use filename as a substitute for tablename
			$lines[0] .= "$padding$outdelim" x $i;
		    }
		    
		    # Print header line before substituting true column name for given label / format
		    if ($header eq 'matlab') { # put (short) table name
			$lines[1] .= (join "_$tablename$outdelim", @columns2)."_$tablename$outdelim";
		    }
		    elsif ($header eq 'excel') { # table name already printed
			$lines[1] .= (join $outdelim, @columns2).$outdelim;
		    }
		    elsif (!@index) { # Phase table for System or other full table (assume first)
			$lines[1] .= (join $outdelim, @columns2).$outdelim;
		    }
		    else { # Phase table for phase
			$lines[1] .= $tablename.$outdelim;
		    }

		    # Set up index array as the variable and put column label in first place
		    if (!@index) {
			unshift @columns, 'index' unless (lc $columns[0] eq 'index');
			$index[0] = 'index';
			# This space will be used for the column name in the matched table
			push @index, 'index';
		    }
		    
		    # Check that index is present (in practice it will be the first column)
		    # Case and '_0' do not matter but otherwise it must be an exact match (i.e. 'pressure' not 'P')
		    $ncolumns = @columns;
		    foreach my $column (@columns) {

			my $label = (grep {(lc $column eq lc $_) || (lc "$column\_0" eq lc $_) || (lc $column eq lc "$_\_0")
					       || (lc "$column"."1" eq lc $_) || (lc $column eq lc "$_"."1")
					       || (lc $column eq "\"$_\"") } @labels)[0];

			if ($label) {
			    $column = $label; # the column name
			}
			elsif ($column =~ /(.+)=/) {

			    $label = $column;
			    $column = $1;			    
			    $ncolumns--; # this tracks the number of 'regular' columns
			    
			    $label =~ s/$column=//;
			    $label =~ s/\{/\$values\{\'/g;
			    $label =~ s/\}/\'\}/g;
			    push @labels, $label;
			    
			    # need some kind of check in here?

			}
			else {
			    ($outfile eq '>') || 
				warn ("Warning: could not find column '$column' in table '$table', file \"$file\" at line $.!\n");
			}
			$values{$column} = ''; # initialise (will stay empty for missing columns and if have '---')
		    }

		}
		else {

		    # blank line signifies end of relevant table so jump out of loop
		    last COLUMNS;

		}
		@output = ();
	    }

	      # hit end of file, either successfully or not

	      unless ($tablename || !$table) {
		  # For solid phases, or any phase in the trace elements output, a whole table could be missing
		  ($outfile eq '>') || warn ("Warning: could not find table '$table' in file \"$file\" at line $.!\n");
		  my $i = @columns; # number of columns
		  if ($header eq 'matlab') { 
		      $lines[1] .= (join "_$table$outdelim", @columns)."_$table$outdelim";
		  }
		  else {
		      $lines[0] .= "\"$table\"";
		      $lines[1] .= (join $outdelim, @columns).$outdelim;
		  }
		  $lines[0] .= "$padding$outdelim" x $i;
	      }

	      if ($lineno > @index) {
		  ($outfile eq '>') || warn ("Warning: first table in \"$colfile\" ('$table' in file \"$file\") ".
					     "does not have entry for each recorded step (at line $.) !\n");
	      }
	      while($lineno < @index) { # length of @index array is zero if no variables are set
		  my $i = @columns; # nuber of columns
		  $lines[$lineno++] .= "$padding$outdelim" x $i; # else pad the rest of the printed table
	      }
	      ($outfile eq '>') || warn("Processed columns: @columns\n\n");

	      close FILE;

	      $lineno = 2;
	      $table = $tablename = ''; # reset for next go
	      @output = ();
	      %values = ();
	      @columns = ();
	      @labels = ();
	      
	  }

      }

	# With incorrect line endings the whole file may be one long line
	if ($. < 3) {
	    ($outfile eq '>') || warn ("Please run file-format.pl on \"$colfile\" and try again!\n");
	    last FOREACH;
	}
	close COLUMNLIST;

    }

}

# Put the header in if coming from run-alphamelts.pl
$lines[0] = ($header eq 'excel') ? $lines[0] : $header;
# Remove extra header line in matlab version
shift @lines if ($header eq 'matlab');

# Output with whatever the correct line endings for the system are
if (@ARGV) {

    print map {$_ .= "\n"} @lines;
    ($outfile eq '>') || warn "Wrote output.\n";

}
else { 

    $outfile =~ s/^(\"|\')//; # Trim any leftover quotes from start and end of file name
    $outfile =~ s/(\"|\')$//;
    # Fix any apostrophes in the name (if had single quotes or no quotes originally)
    $outfile =~ s/\'\\\'\'/\'/g || $outfile =~ s/\\\'/\'/g;
    $outfile =~ s/\\//g unless $windows; # Remove escaping from any other characters (e.g. '(' and ')')

    if (open(OUTPUT, ">$outfile")) {
	print OUTPUT map {$_ .= "\n"} @lines;
	close (OUTPUT);
	print ("Wrote output.\n");
    }
    else {
	warn "Error: could not open file \"$outfile\" : $!\n";
    }
    print ("Press return to finish.\n");
    $_ = <STDIN>;

}
