#!/usr/bin/perl -w

# run-alphamelts.pl 2.3
use strict;

my (@argv2, $in_file, $settings_file, $melts_file, $log_file, $batch_file,
    $new_file, $old_file, $out_file, $column_file, $table_file, $table_row);
my ($version, $newversion, $path, $run_path, $batch, $auto, $sets, %switch);
my ($windows, $program, $rs, $run, $debug, $delim, $line, @lines, $title);

@argv2 = ();
until (@argv2) {

    if (@ARGV) {

        @argv2 = @ARGV;

    }
    else {

        my @temp;

        print ("No command line switches specified! ");
        print ("Please enter switches now (or press return for defaults).\n");

        $_ = <STDIN>;
        chomp;

        s/\s+$//; # Trim any trailing white space
        s/\'/\"/g; # Windows can only deal with name enclosed in double quotes; *nix can take single or double

        s/\s+-/ --/g; # Assume that filenames etc. cannot begin with '-'
        @temp = split /\s+-/;
        map { push @argv2, (split /\s+/, $_, 2); } @temp;

    }

    $rs = $/; # expected record separator i.e. \r\n or \n
    $delim = '/';
    $windows = '';
    if (exists $ENV{'OS'}) {
        if ($ENV{'OS'} =~ /windows/i) {
            $windows = 1;
        }
    }
    $delim = '\\' if ($windows);

    $debug = ($windows) ? '2>NUL' : '2>/dev/null';

    $version = "2.03";
    print "Checking for updates ... $version ... ";
    $newversion = ! system("curl -m1 https://magmasource.caltech.edu/alphamelts/alphameltsdate");
    if ($newversion) {
        $newversion = `curl https://magmasource.caltech.edu/alphamelts/alphameltsdate $debug`;
        $newversion =~ s/(\d\.\d+)\.(\d+)/$1$2/;
    }
    $version =~ s/(\d\.)0*(\d+)/$1$2/;

    #          01234567890123456789012345678901234567890123456789012345678901234567890123456789
    #          This is run-alphamelts.pl X.X; use it with alphamelts X.X or updates X.X.X
    # -h prints the following help before quitting
    if (@argv2 && $argv2[0] eq '-h') {
        warn ("\n usage: run-alphamelts\.pl                                                      \n".
              "\n [-h]                 print this brief help                                     \n".
              " [-f input_file]      general file for environment variables, settings and menu \n".
              "                      input for alphamelts (can use -s, and -b or -a to split   \n".
              "                      into separate input_file, settings_file and/or batch_file)\n".
              " [-m melts_file]      .melts format template file to use for compositions etc.  \n".
              " [-n file_name]       file name base to use for new melts_file(s) (default uses \n".
              "                      original file name from the melts_file template)          \n".
              " [-t table_file]      table of compositions, P, T, and/or fO2 to be substituted \n".
              "                      into the melts_file template                              \n".
              " [-s [settings_file]] (optional) separate file with alphamelts settings to add  \n".
              "                      to melts_file(s) (if no file is specified will add lines  \n".
              "                      to suppress certain phases for backwards compatibility)   \n".
              " [-b [batch_file]]    batch processing mode using alphamelts menu commands from \n".
              "                      input_file, or batch_file if specified                    \n".
              " [-a [batch_file]]    same as '-b', except with '-t' the menu commands are      \n".
              "                      automatically repeated for each new melts_file            \n".
              "                      (ALPHAMELTS_CALC_MODE must be set to MELTS, MELTSandCO2,  \n".
              "                      MELTSandCO2_H2O, or pMELTS in input_file)                 \n".
              " [-p output_path]     path of directory for output files (if not the same as    \n".
              "                      working directory)                                        \n".
              " [-l log_file]        name of additional file to record environment, settings   \n".
              "                      and menu input (default logfile.txt)                      \n".
              " [-c column_file]     run 'column-pick.pl filename' for each filename listed in \n".
              "                      column_file before moving the output files to output_path \n".
              " [-o output_file]     filename for output of column-pick.pl (default output is  \n".
              "                      alphaMELTS_tbl.txt)                                       \n".
              " [-d]                 debug calculation issues by printing more program output  \n".
              "                      to screen                                                 \n".
              " [-x]                 process input files as usual but do not run the alphamelts\n".
              "                      executable (or use with '-p' or '-c' to move or run       \n".
              "                      column-pick on previously generated output files)         \n".
              "\n This is run-alphamelts.pl $version; use it with alphamelts $version or updates $version.X.\n\n");
        next;
    }
    elsif ($newversion) {
        $ENV{'ALPHAMELTS_VERSION'} = "$newversion";
    }

    $path = $batch = '';
    $run_path = $0; # $0 is the name (and path) of the run-alphamelts script (not the .bat)
    $run_path =~ s/run-alphamelts.pl//;
    $run = '';
    $program = 'alphamelts2'; # for backwards compatibility until alphamelts 1.9 is retired

    $in_file = $settings_file = '';
    $log_file = '';
    $table_file = $melts_file = $batch_file = $column_file = '';
    $out_file = 'alphaMELTS_tbl.txt';

    if (grep /^-/, @argv2) {

        my $temp;

        # check command line switches
        map {

            if (/^-/) {

                $temp = $_;

            }
            else {

                s/^\"//; # Trim any leftover quotes from start and end of list
                s/\"$//;
                # Fix any apostrophes in the name (if had single quotes or no quotes or double quotes originally)
                s/\"\\\"\"/\'/g || s/\\\"/\'/g || s/\"/\'/g;

                s/\\//g unless $windows; # Remove escaping from any other characters (e.g. '(' and ')')
                $switch{$temp} = $_;

            }

        } @argv2;

        $path = $switch{'-p'} if exists $switch{'-p'};
        $path .= $delim if ($path); # path and / for good measure

        $program = '' if (grep /^-x$/, @argv2);
        $debug = '' if (grep /^-[dgv]$/, @argv2);
        print "Process ID = $$\n" unless $debug;
        $run = 'gdb ' if (grep /^-g$/, @argv2);
        $run = 'valgrind --leak-check=full --show-leak-kinds=definite,possible'.
            '--track-origins=yes --expensive-definedness-checks=yes ' if (grep /^-v/, @argv2);

        $auto = grep /^-a/, @argv2;
        $batch = grep /^-b/, @argv2;
        $batch_file = $switch{'-b'} if exists $switch{'-b'};
        $batch_file = $switch{'-a'} if exists $switch{'-a'};

        $sets = grep /^-s/, @argv2;
        $settings_file = $switch{'-s'}  if exists $switch{'-s'};

        $in_file = $switch{'-f'}  if exists $switch{'-f'};
        $melts_file = $switch{'-m'}  if exists $switch{'-m'};
        $new_file = $switch{'-n'}  if exists $switch{'-n'};
        $table_file = $switch{'-t'}  if exists $switch{'-t'};
        $log_file = $switch{'-l'}  if exists $switch{'-l'};
        $column_file = $switch{'-c'}  if exists $switch{'-c'};
        $out_file = $switch{'-o'}  if exists $switch{'-o'};

        if ((grep /^-[^hfmstnaboplcdxgv]$/, @argv2) || (grep /^-..+$/, @argv2)) {
            grep {warn "RUN-ALPHAMELTS.PL ERROR: Unknown option \'$_\'\n" if $_ =~ /^-[^hfmstnaboplcdxgv]$/} @argv2;
            grep {warn "RUN-ALPHAMELTS.PL ERROR: Unknown option \'$_\'\n" if $_ =~ /^-..+$/} @argv2;
            next;
        }
        if (exists $switch{'-d'} || exists $switch{'-g'} || exists $switch{'-v'}) {
            grep {warn "RUN-ALPHAMELTS.PL ERROR: Unknown option \'$_ $switch{$_}\'\n" if exists $switch{$_}} ('-d', '-g', '-v');
            next;
        }

    }
    elsif (@argv2) {
        grep {warn "RUN-ALPHAMELTS.PL ERROR: Unknown option \'$_\'\n"} @argv2;
        next;
    }
    else {
        unshift @argv2, "\(none\)";
    }

    if ($program) {
        unless (open (LOGFILE, '>', 'logfile.txt')) {
            my $temp = ($log_file) ? ' temporary' : '';
            warn "RUN-ALPHAMELTS.PL WARNING: Cannot open (temporary) log_file logfile.txt on 1st try ($!)\n";
            unless ((chdir "$run_path") && (open (LOGFILE, '>', 'logfile.txt'))) {
                warn "RUN-ALPHAMELTS.PL WARNING: Cannot open (temporary) log_file logfile.txt on 2nd try ($!)\n";
                next;
            }
        }
    }
    if ($in_file) {
        unless (open (INPUT, '<', "$in_file")) {
            warn "RUN-ALPHAMELTS.PL ERROR: Cannot open input_file \"$in_file\" ($!)\n";
            next;
        }

        @lines = <INPUT>; # read whole file into array of lines - O.K. as it's not very big
        @lines = split /[\r\n]+/, (join '', @lines); # line endings may not be correct

        close(INPUT);

        if ($program) {

            print LOGFILE map {$_ .= "\n"} @lines;
            map {s/^\s*//; s/\=//; s/!.*//; chomp} @lines; # remove leading white space and trailing comments
            @lines = grep /.+/, @lines;

            map {

                # Actually only ALPHAMELTS_CALC_MODE, MINT etc. - others are substituted into melts_file
                if (/^ALPHAMELTS/) {
                    %ENV = (
                        %ENV, # existing environment variables
                        split # split by white space into key and value
                        )
                }

            } @lines;
        }

    }

    if ($settings_file) {
        unless (open (INPUT, '<', "$settings_file")) {
            warn "RUN-ALPHAMELTS.PL ERROR: Cannot open settings_file \"$settings_file\" ($!)\n";
            next;
        }

        @lines = (); # throw out input file lines
        @lines = <INPUT>; # read whole file into array of lines - O.K. as it's not very big
        @lines = split /[\r\n]+/, (join '', @lines); # line endings may not be correct

        close(INPUT);

        map {s/^\s*//; s/\=//; s/!.*//; chomp} @lines; # remove leading white space and trailing comments
        @lines = grep /.+/, @lines;

    }
    elsif ($sets) {
        @lines = (
            "Suppress: tridymite",
            "Suppress: cristobalite",
            "Suppress: rutile",
            "Suppress: corundum",
	    "Limit number: plagioclase 1",
	    "Limit number: alkali-feldspar 1"
	    )
    }

    if($melts_file) {

        unless (open (OLDMELTS, '<', "$melts_file")) {
            warn "RUN-ALPHAMELTS.PL ERROR: Cannot open melts_file \"$melts_file\" ($!)\n";
            next;
        }

        my @newlines = <OLDMELTS>; # read whole melts file into array of lines - O.K. as it's not very big
        @newlines = split /[\r\n]+/, (join '', @newlines); # line endings may not be correct

        map {s/^\s*//; s/!.*//; chomp} @newlines; # remove leading white space and trailing comments
        @newlines =  grep /.+/, @newlines;
        map {$_ .= "\n"} @newlines;

        close OLDMELTS;

        if ($new_file) {
            $melts_file = $new_file;
        }
        else {
            unless (rename $melts_file, "$melts_file\_bak") {
                warn "RUN-ALPHAMELTS.PL ERROR: Cannot rename melts_file \"$melts_file\" ($!)\n";
                next;
            }
        }
        unless (open (NEWMELTS, '>', "$melts_file")) {
            warn "RUN-ALPHAMELTS.PL ERROR: Cannot open melts_file \"$melts_file\" ($!)\n";
            next;
        }

        # ignore environment variables and menu input for now
        @lines = grep /:/, @lines;

        map {
            $line = "$_\n";
            print "$line"; # write line to screen

            if (/(Mode: fractionate|Suppress:)/i) {
                my $match = $1;
                if(/none/i) {
                    @newlines = grep !(/$match/i), @newlines ; # no phase names contain 'none'
                }
                else {
                    push @newlines, $line unless grep /$line/, @newlines;
                }
            }
            elsif (/(Initial Composition:|Initial Trace:|Limit number:)/i) {
                push @newlines, $line unless grep {s/$_/$line/i if $_ =~ (split /\s+/, $line)[2]} @newlines;
            }
            elsif (/Assimilant:/i) {
                push @newlines, $line unless grep {s/$_/$line/i if $_ =~ (split /\s+/, $line)[1]} @newlines;
            }
            else {
                # substitute old line in melts file for new line from input file
                # or add new line if not previously set in melts file
                push @newlines, $line unless grep {s/$_/$line/i if $_ =~ (split /\:/, $line)[0]} @newlines;
            }
        } @lines;
        @lines = ();

        print NEWMELTS @newlines; # write new melts file
        close NEWMELTS;

        if($table_file) {

            unless (open (TABLE, '<', "$table_file")) {
                warn "RUN-ALPHAMELTS.PL ERROR: Cannot open table_file \"$table_file\" ($!)\n";
                next;
            }

            @lines = <TABLE>; # read whole file into array of lines - O.K. if it's not very big
            @lines = grep /.+/, (split /[\r\n]+/, (join '', @lines)); # line endings may not be correct

            close TABLE;

            my @used_oxides = (); $table_row = 1; # anything non-zero would work
            map {

                s/^\s*//;
                if(@used_oxides) {

                    my (%values, $label, $value);

                    $new_file = $melts_file;
                    $table_row++;
                    $new_file =~ s/\./$table_row\./ || ($new_file .= $table_row);

                    open (NEWMELTS, '>', "$new_file");

                    grep {# add row number to title
                        s/\($table_file row .*/\($table_file row $table_row\)/i ||
                            s/^Title.*/$& \($table_file row $table_row\)/i
                    } @newlines;

                    $label = 0;
                    map {$values{$used_oxides[$label++]} = $_} split /[ \t,]+/;

                    foreach $label (@used_oxides) {
                        $value = $values{$label};
                        grep {# initial composition or trace (check for equality so P does not match Pb etc.)...
                            /^(Initial|Log)/i && ((lc $label eq lc ((split)[2])) ? s/\s$label.*/ $label\ $value/i
                                                  : s/\s$label:.*/ $label:\ $value/i) # ... or any other initial value
                        } @newlines;
                    }

                    print NEWMELTS @newlines; # write new melts file
                    close NEWMELTS;

                }
                else {

                    my ($label, @test);
                    @used_oxides = split /[ \t,]+/;

                    # Offset: in table file works - need to fix

                    foreach $label (@used_oxides) { #at least one match
                        @test = grep {/^(Initial|Log)/i && ((lc $label eq lc ((split)[2])) || /$label\: /i)} @newlines;
                        if(@test) {
                            $table_row = 0; #use as a counter of lines after header
                            next;
                        }
                    }
                    @used_oxides = () if ($table_row); # not yet got to header

                }
            } @lines;

        }
        @newlines = ();
    }
    elsif ($table_file || $settings_file) {
        warn "RUN-ALPHAMELTS.PL ERROR: Please provide a melts_file to use as template for the table_file\n" if $table_file;
        warn "RUN-ALPHAMELTS.PL ERROR: Please provide a melts_file to use as template for the settings_file\n" if $settings_file;
        next;
    }
    @lines = ();

    if ($program && $auto && !exists $ENV{'ALPHAMELTS_CALC_MODE'}) {
        warn "RUN-ALPHAMELTS.PL ERROR: Please put ALPHAMELTS_CALC_MODE in input_file to run in automatic mode \n";
        next;
    }
    if ($auto || $batch) {

        $old_file = ($batch_file) ? $batch_file : $in_file;
        if ($old_file) {

            my ($i);

            $new_file = "auto_batch.txt";

            unless (open (BATCH, '>', "$new_file")) {
                warn "RUN-ALPHAMELTS.PL ERROR: Cannot open batch_file \"$new_file\" ($!)\n";
                next;
            }
            unless (open (OLDFILE, '<', "$old_file")) {
                warn "RUN-ALPHAMELTS.PL ERROR: Cannot open batch_file \"$old_file\" ($!)\n";
                next;
            }

            @lines = <OLDFILE>; # read whole file into array of lines - O.K. as it's not very big
            @lines = split /[\r\n]+/, (join '', @lines); # line endings may not be correct

            close OLDFILE;

            map {s/^\s*//; s/!.*//; chomp} @lines; # remove leading white space and trailing comments
            @lines = grep /.+/, @lines;
            @lines = grep {! /(^ALPHAMELTS|:)/ } @lines;

            map {$_ .= "\n"} @lines;

            if ($table_file) {

                my $new_melts = $melts_file;
                pop @lines; # take off last line (should be 'x')

                for ($i = 1; $i <= $table_row; $i++) {

                    my $old_melts = $new_melts;
                    $new_melts = $melts_file;
                    $new_melts =~ s/\./$i\./ || ($new_melts .= $i);

                    map {s/$old_melts/$new_melts/g} @lines;
                    print BATCH @lines
                }

                print BATCH "x\n";

            }
            else {
                print BATCH @lines;
            }

            close BATCH;
            @lines = ();

        }
        else {
            warn "RUN-ALPHAMELTS.PL ERROR: Please provide a batch_file (or input_file) to run in automatic mode\n";
            next;
        }
    }

    #if (exists $ENV{'ALPHAMELTS_INTEGRATE_FILE'}) {
#       my $temp_file = $ENV{'ALPHAMELTS_INTEGRATE_FILE'};
#       $ENV{'ALPHAMELTS_INTEGRATE_FILE'} = "$path$temp_file";
    #}

    map {print "$_ $ENV{$_}\n" if /^ALPHAMELTS/;} keys %ENV;

# ABOUT TO RUN ALPHAMELTS ******************

    if ($program) {

        my ($me, @months, $year, $month, $date, $hour, $min);

        ($me = $ENV{'USERNAME'}) || ($me = $ENV{'USER'}); # system environment variable
        chomp $me;

        @months = ('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec');

        ($min, $hour, $date, $month, $year) = (localtime)[1..5];

        $year += 1900;
        $month = $months[$month];
        $hour = sprintf "%2.2d", $hour;
        $min = sprintf "%2.2d", $min;

        # record run conditions in log file
        print LOGFILE "\n\! <><> alphaMELTS $version run by $me <><> \n".
            "! Time started: $month $date $year at $hour\:$min \n\n".
            "! Command line: @argv2 \n\n";
        print LOGFILE ($auto || $batch) ? "! Batch file used: $batch_file \n" :  "! alphaMELTS menu input used: \n\n";
        close (LOGFILE);

        if ($auto || $batch) {

            # $program does not have spaces etc. but $run_path might
            ((-f "$run_path$program") && !(system "$run\"$run_path$program\" 1 $debug < auto_batch.txt")) ||
                ((-f "$run_path$program.bat") && !(system "$run\"$run_path$program\" 1 $debug < auto_batch.txt")) ||
                (!(-f "$run_path$program") && !(system "$run$program 1 $debug < auto_batch.txt")) ||
                warn "RUN-ALPHAMELTS.PL WARNING: alphamelts may have crashed!\n";

        }
        else {

            # print relevant file names before execution if using alphamelts interactively
            if ($table_file) {
                $new_file = $melts_file;
                $new_file =~ s/\./1\./ || ($new_file .= '1');
                print "\nPlease use melts files: $new_file, ";
                $new_file =~ s/1\./2\./ || s/1$/2/;
                print "$new_file etc.\n\n";
            }
            else {
                ($melts_file) ? print "\nPlease use melts file: $melts_file\n\n" : print "\n\n";
            }

            # Try $run_path first as can test for the file; then try path (e.g. for Mac double-click)
            ((-f "$run_path$program") && !(system "$run\"$run_path$program\" $debug")) ||
                (!(-f "$run_path$program") && !(system "$run\"$program\" $debug")) ||
                warn "RUN-ALPHAMELTS.PL WARNING: alphamelts may have crashed!\n";

        }
        print "\n\n";

    }

# FINISHED RUNNING ALPHAMELTS 2 **********************

    (-d $path || mkdir "$path") if ($path);

    # Test solids first so don't waste time processing old file
    unless (-f "Phase_main_tbl.txt") {
        warn "RUN-ALPHAMELTS.PL WARNING: Cannot find output file \"Phase_main_tbl.txt\" ($!)\n";
        warn "Please check that alphamelts ran properly.\n";
        next;
    }

    if (-f "liquid-model-batch.inp") { # MELTS intermediate file
        unlink "liquid-model-batch.inp";
    }

    if (open (INPUT, "<Phase_main_tbl.txt")) {
        @lines = <INPUT>;
        close (INPUT);
        unless (grep /(index\s+)(.+)(Pressure\s+)(.+)(Temperature.*)/, @lines) {
            # probably means the Solids file has already been processed
            warn "RUN-ALPHAMELTS.PL WARNING: Incorrect format in output file \"Phase_main_tbl.txt\".\n";
            warn "File has already been processed? Please check that alphamelts ran properly.\n";
            next;
        }
    }

    my $phase_file = 'Phase_main_tbl.txt';

    # now start processing output; first check file exists
    unless ((-f $phase_file) && (rename $phase_file, "$phase_file\_bak") &&
            (open (INPUT, '<', "$phase_file\_bak")) && (open (TABLE, "+>$phase_file"))) {
        warn "RUN-ALPHAMELTS.PL ERROR: Cannot open output file \"$phase_file\" ($!)\n";
        next;
    }

    my (@columns, $label, @string, $pt, $key, $val, %therm_comp, $el_list);

    $title = '';
    while (<INPUT>) {

        chomp;
        if ($title) {
            print TABLE "$title\n" if ($title); # print the title
            $title = '';
        }
        elsif (/^Title\:/) {
            $title = $_;  # print the title next time (above)
        }
        elsif (/^index/) { # get the new P and T value

            $pt = $_;
            $pt =~ s/(index\s+)(.+)(Pressure\s+)(.+)(Temperature.*)/$2$4/;
            s/(index\s+)(.+)(Pressure\s+)(.+)(Temperature\s+)//;
            @columns = split /\s+/; # open once to get element list
            $pt .= shift @columns;
            $el_list = "@columns";

        }
        else {

            s/\ oxide/_oxide/g;
            s/\ ss/_ss/g;
            @columns = split /\s+/;
            $label = shift @columns;

            if($pt) {
                # therm_comp is a hash, with the phase name as the key and an array of lines as the corresponding entry
                push @{$therm_comp{"$label"}}, "$pt @columns\n"; # add a line to the appropriate array
            }
            else {
                # probably means the Solids file has already been processed
                warn "RUN-ALPHAMELTS.PL WARNING: Incorrect format in output file \"$out_file\".\n";
                warn "File has already been processed? Please check that alphamelts ran properly.\n";
                next;
            }

        }
    }

    close INPUT;
    unless (unlink "$phase_file\_bak")  {
        warn "Could not delete backup file \"$phase_file\_bak\" ($!)\n";
        next;
    }

    @lines = ();
    push @lines, "Name: col_phase_mass.txt\n".
        "Header: Phase Masses:\n".
        "Format: prn\n".
        "Padding: '0.0'\n\n".
        "File: System_main_tbl.txt\n".
        "Table: System\n".
        "Columns: index Pressure Temperature\n\n";

    $pt = "index Pressure Temperature mass H S V Cp";

    foreach $key (sort keys %therm_comp) {
        if ($key =~ /^liquid/) { # liquid(s) first

            $val = $therm_comp{$key};
            @string = @$val;

            print TABLE "\n$key thermodynamic data and composition:\n$pt viscosity $el_list\n";
            print TABLE @string;

            push @lines, "File: Phase_main_tbl.txt\n".
                "Table: $key\n".
                "Columns: mass\n\n";

        }
    }
    foreach $key (sort keys %therm_comp) {
        unless ($key =~ /^liquid/) { # solids after

            $val = $therm_comp{$key};
            @string = @$val;

            print TABLE "\n$key thermodynamic data and composition:\n$pt formula $el_list\n";
            print TABLE @string;

            push @lines, "File: Phase_main_tbl.txt\n".
                "Table: $key\n".
                "Columns: mass\n\n";
        }
    }

    close TABLE;

    unless (open (COLFILE, '>', 'col_phase_mass.txt')) {
        warn "RUN-ALPHAMELTS.PL ERROR: Cannot open column_file \"col_phase_mass.txt\" ($!)\n";
        next;
    }

    print COLFILE @lines;

    close COLFILE;

    $program = 'column-pick.pl';
    ((-f "$run_path$program") && !(system "$run\"$run_path$program\" col_phase_mass.txt $debug > Phase_mass_tbl.txt")) ||
        (!(-f "$run_path$program") && !(system "$run$program col_phase_mass.txt $debug > Phase_mass_tbl.txt")) ||
        warn "RUN-ALPHAMELTS.PL WARNING: column-pick.pl may not have run properly!\n";

    map {s/col_phase_mass/col_phase_vol/; s/Masses/Volumes/; s/mass/V/g} @lines;

    unless (open (COLFILE, '>', 'col_phase_vol.txt')) {
        warn "RUN-ALPHAMELTS.PL ERROR: Cannot open column_file \"col_phase_vol.txt\" ($!)\n";
        next;
    }

    print COLFILE @lines;

    close COLFILE;

    ((-f "$run_path$program") && !(system "$run\"$run_path$program\" col_phase_vol.txt $debug > Phase_vol_tbl.txt")) ||
        (!(-f "$run_path$program") && !(system "$run$program col_phase_vol.txt $debug > Phase_vol_tbl.txt")) ||
        warn "RUN-ALPHAMELTS.PL WARNING: column-pick.pl may not have run properly!\n";


}
continue {

    # If forgot -c or -p switches don't have to rerun whole program - just choose menu option 'x'
    if ($column_file) {
        $program = 'column-pick.pl';
        ((-f "$run_path$program") && !(system "$run\"$run_path$program\" $column_file > $out_file")) ||
            (!(-f "$run_path$program") && !(system "$run$program $column_file")) ||
            warn "RUN-ALPHAMELTS.PL WARNING: column-pick.pl may not have run properly!\n";
    }

    $program = ($windows) ? 'copy /Y' : 'cp -f'; # perl's copy requires a module to be installed
    if ($path) {
        system "$program logfile.txt \"$path$log_file\"";

        $program = ($windows) ? 'move /Y' : 'mv -f';
        my @files = glob '*_tbl.txt *.tbl *.out';
        for my $file (@files) {
            system "$program $file \"$path\"";
        }
        system "$program logfile.txt \"$path\"";
    }
    elsif ($log_file) {
        system "$program logfile.txt $log_file";
    }

    unless (@ARGV) {
        %switch = ();
        print ("Run again (y or n)?\n");
        $_ = <STDIN>;
        chomp;
        # Mostly stays open if user accidently drags and drops before answering...
        @argv2 = () unless (/^$/ || /^n/i); # shuts if no answer or begins with n/N.
    }

}
