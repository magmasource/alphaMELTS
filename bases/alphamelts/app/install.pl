#!/usr/bin/perl -w

use strict;

my($windows, $pwd, $bin, $install, $home, $path, $priv_bin, $program, $file, @files, $in_file);

$windows = ($^O =~ /MSWin32/i) ? 1 : '';
if(exists($ENV{'OS'})) {
    if($ENV{'OS'} =~ /Windows/i) {
        $windows = 1;
    }
}

my $done = '';
until ($done) {

    $pwd = $0; # $0 is the name (and path) of the install script
    $pwd =~ s/install.pl//;
    $pwd = (($windows) ? '.\\' : './') unless ($pwd); # for 'perl ...' or Windows command line
    # for double-clicking on Mac or 'Open with' on Windows
    unless (chdir "$pwd") {
        warn "Could not change to directory $pwd: $!";
        last;
    }

    if($windows) {

        $home = "$ENV{'USERPROFILE'}";
        $home =~ s/^\"//;
        $home =~ s/\"$//;

        print "\nUSERPROFILE=$ENV{'USERPROFILE'}\n";
        print "PATH=$ENV{'PATH'}\n\n";

        $bin = "$home\\Documents\\bin";

        $file = "alphamelts2.bat run-alphamelts.pl.bat column-pick.pl.bat file-format.pl.bat";
        # where.exe seems to hang if STDERR is redirected
        @files = `for \%i in ($file) do \@find /v "echo off" \"\%~\$PATH:i\" 2> nul`;

        @files = grep {
            if (/.+/) {
                chomp;
                s/(^\-+ )(.*)/$2 \-\>/ || s/(\.exe\"*)/$1\r\n/ || s/\%\*/\r\n/;
            }
        } @files;

    }
    else {

        # $ENV{'PWD'} may be empty, if running as sudo, or be home directory, if double-clicked on Mac
        # $pwd may be './' if not double-clicked.
        $pwd = `pwd`;
        chomp $pwd;
        $home = $ENV{'HOME'};

        print "\nHOME=$ENV{'HOME'}\n";
        print "PATH=$ENV{'PATH'}\n\n";

        $bin = "$home\/bin";
        if ($^O !~ /linux/) {
            $bin = ($ENV{'PATH'} =~ /^.*:*$bin:*.*/) ? $bin : "$home\/Documents\/bin";
        }

        $file = "alphamelts2 run-alphamelts.pl column-pick.pl file-format.pl";
        @files = `which $file 2> /dev/null | xargs ls -go`;
        @files = grep s/(.* \d\d\:\d\d)(.*) \-\> (.*)/$2 \-\> $3/, @files;

    }

    if (@files) {

        print "Currently installed version:\n\n @files\n";
        print "\nContinue with installation (y or n)? ";
        $_ = <STDIN>;
        chomp;
        last if (/^$/ || /^n/i);
        print "\n";

    }

    $install = $pwd;

    print "Enter full path for installation directory, or press return to use the current directory".
        " given in brackets\n[$install]: ";
    $files[0] = <STDIN>;
    print "\nEnter full path of directory to put links in, or press return to use the default location".
        " given in brackets\n[$bin]: ";
    $files[1] = <STDIN>;

    for (my $i = 0; $i <=2; $i++) {

        $_ = $files[$i];
        unless (/^$/) {

            s/\s+$//;  # Trim any trailing white space
            s/\'/\"/g; # Windows can only deal with name enclosed in double quotes;
            s/^\"//;   # *nix can take single or double
            s/\"$//;   # Trim any leftover quotes from start and end of list

            # Fix any apostrophes in the name (if had single quotes or no quotes or double quotes originally)
            s/\"\\\"\"/\'/g || s/\\\"/\'/g || s/\"/\'/g;
            s/\\//g unless ($windows); # Remove escaping from any other characters (e.g. '(' and ')')

        }

    }

    $priv_bin = 0;
    if($windows) {
        $path = `echo ;%PATH%; | find /C /I \";$bin;\"`;
        chomp $path;
    }
    else {
        $path = ($ENV{'PATH'} =~ /^.*:*$bin:*.*/);
    }
    unless (exists $ENV{'SUDO_USER'}) {
        print "\nThe '$bin' directory appears".(($path) ? " " : " not ")."to be in your path.\n";
        unless ($path) {
            if (($^O =~ /linux/) && ($ENV{'SHELL'} =~ /bash/) && ($bin eq "$home\/bin") && !(-d "$bin")) {
                print "The system may add it to the path after installation.\n";
                $priv_bin = 1;
            }
            else {
               print "install.pl will try to add it to the path during installation.\n";
            }
        }
    }

    print "\nContinue with installation (y or n)? ";
    $_ = <STDIN>;
    chomp;
    last if (/^$/ || /^n/i);
    print "\n";

    unless ((-d "$install") || (mkdir "$install")) {
        warn "Could not make directory '$install': $!";
        last;
    }
    unless ((-d "$bin") || (mkdir "$bin")) {
        warn "Could not make directory '$bin': $!";
        last;
    }

    $program = '';
    if ($windows) {

        # back up relevant parts of registry as .reg file as go along
        # HKLM = HKEY_LOCAL_MACHINE ... settings for all users (requires admin password)
        # HKCU = HKEY_CURRENT_USER ... settings for current user (overrides local machine ones)
        # HKCR = HKEY_CLASSES_ROOT ... a combination of user and machine settings

        my $key = 'HKCU\Software\Classes\.pl';
        system "reg export $key command.reg /y" unless (system "reg query $key");
        system "reg add $key /ve /t REG_SZ /d command_auto_file /f";
        print "Command was: reg add $key /ve /t REG_SZ /d command_auto_file /f\n\n";

        $key = 'HKCU\Software\Classes\Applications\perl.exe';
        system "reg export $key perl.reg /y" unless (system "reg query $key");

        # Use whichever perl program has been assigned to open install.pl
        system "reg add $key\\shell\\open\\command /ve /t REG_SZ /d \"$^X \\\"%1\\\" %*\" /f";
        print "Command was: reg add $key\\shell\\open\\command /ve /t REG_SZ /d \"$^X \\\"%1\\\" %*\" /f\n\n";

        # Right-clicking allows user to view / edit script
        system "reg add $key\\shell\\edit\\command /ve /t REG_EXPAND_SZ /d".
            " \"\%SystemRoot\%\\system32\\Notepad.exe %1\" /f";
        print "Command was: reg add $key\\shell\\edit\\command /ve /t REG_EXPAND_SZ /d".
            " \"\%SystemRoot\%\\system32\\Notepad.exe %1\" /f\n\n";

        $key = 'HKCU\Software\Classes\command_auto_file';
        system "reg export $key command_auto_file.reg /y" unless (system "reg query $key");
        system "reg copy HKCU\\Software\\Classes\\Applications\\perl.exe $key /s /f";
        print "Command was: reg copy HKCU\\Software\\Classes\\Applications\\perl.exe $key /s /f\n\n";

        if (-f 'alphamelts_win64.exe') {
            $program = 'alphamelts_win64';
        }

        if ($program) {
            $program .= '.exe';
            unless (open(BAT, ">$bin\\alphamelts2.bat")) {
                warn "Could not open alphamelts2.bat file: $!\n";
                last;
            }
            print BAT "\@echo off\n";
            print BAT "\"$install\\$program\"\n";
            close(BAT);
            # extra in case there is a problem setting the Path
            unless (open(BAT, ">$install\\alphamelts2.bat")) {
                warn "Could not open alphamelts2.bat file: $!\n";
                last;
            }
            print BAT "\@echo off\n";
            print BAT "\"$install\\$program\"\n";
            close(BAT);
        }
        else {
            warn "Warning: installing .pl scripts only as cannot find alphamelts2 executable!";
        }

        @files = ('run-alphamelts', 'column-pick', 'file-format');
        foreach $file (@files) {
            unless (open(BAT, ">$bin\\$file\.pl\.bat")) {
                warn "Could not open $file\.pl\.bat file: $!\n";
                last;
            }
            print BAT "\@echo off\n";
            print BAT "\"$install\\$file\.pl\" %*\n";
            close(BAT);
        }

        unless ($path) {

            my $key = 'HKCU\Environment';

            if (system "reg query $key /v Path") {
                print "Adding \"$bin\" to the Path variable.\n";
                system "setx Path \"$bin\"";
                print "Command was: setx Path \"$bin\"\n\n";
            }
            else {
                system "reg export $key path.reg /y";

                print "\n<><><><> IMPORTANT: the 'Path' variable exists so not going to edit it! <><><><>\n\n".
                    "      Please ensure that the following directory is in you path...\n\n            $bin\n\n".
                    "      On Windows 11, click the Start Button and select 'Settings' or\n".
                    "      click the cog symbol to open Windows Settings. Type 'environment'\n".
                    "      in the search box and select the option 'Edit environment variables\n".
                    "      for your account' from the drop-down menu.\n\n".
                    "      Select the user's Path variable, and click 'Edit'. Select 'New',\n".
                    "      then copy and paste the following for the new value:\n\n            $bin\n\n";

            }
        }

    }
    else {

        $program = '';
        if ($^O =~ /linux/i) {
            $program = 'alphamelts_linux' if (-f 'alphamelts_linux');
        }
        elsif (!$windows) {
            if (-f 'alphamelts2_macos') {
                # Intel binary (to be phased out eventually)
                $program = 'alphamelts2_macos';
            }
            else {
                # Apple Silicon (except v2.1.0 - 2.2.1, which were fat binaries)
                $program = 'alphamelts_macos' if (-f 'alphamelts_macos');
            }
        }

        @files = ('column-pick.pl', 'file-format.pl', 'run-alphamelts.pl');
        push @files, glob '*.melts';
        push @files, glob '*.txt';
        push @files, glob 'tutorial'.(($windows) ? "\\" : "/").'*.txt';
        push @files, glob 'tutorial'.(($windows) ? "\\" : "/").'*.melts';
        $file = join ' ', @files;

        if (-x "$pwd/file-format.pl") {
            (exists $ENV{'SUDO_USER'}) ? system "sudo -u $ENV{'SUDO_USER'} \"$pwd/file-format.pl\" $file" :
            system "$^X \"$pwd/file-format.pl\" $file";
        }
        else { # previous install means local file-format not executable?
            (exists $ENV{'SUDO_USER'}) ? system "sudo -u $ENV{'SUDO_USER'} \"file-format.pl\" $file" :
            system "$^X file-format.pl $file";
        }

        $done = 1; # perl's symlink function has no force
        # On Mac force will replace a regular file with a link to itself!
        if ($program) {
            # extra in case there is a problem setting the Path
            $done = ((-f "$bin/alphamelts2") && !(-l "$bin/alphamelts2")) ? '' : $done;
            if ($done) {
                # Can deal with spaces in bin folder only if user makes it and drags in escaped name (not quoted)
                $done = (system "ln -sf \"$install/$program\" \"$bin/alphamelts2\"") ? '' : $done;
                $done = (system "ln -sf \"$install/$program\" \"$install/alphamelts2\"") ? '' : $done;
            }
        }
        else {
            warn "Warning: installing .pl scripts only as cannot find alphamelts2 executable!";
        }

        @files = ('column-pick.pl', 'file-format.pl', 'run-alphamelts.pl');
        for $file (@files) {
            $done = ((-f "$bin/$file") && !(-l "$bin/$file")) ? '' :
                    ((system "ln -sf \"$install/$file\" $bin/$file") ? '' : $done);
        }

        unless ($done) {
            warn "Warning: could not make all links! Please check that the links directory is writable.\n";
            warn "Also, please ensure that the installation and links directories are not the same.\n";
            last;
        }

        unless ($path || exists $ENV{'SUDO_USER'} || $priv_bin) {
            if ($ENV{'SHELL'} =~ /zsh/) {
                unless (open(ZSH, ">>$home/.zshrc")) {
                    warn "Could not open .zshrc file: $!";
                    last;
                }
                print ZSH "\n# added by alphaMELTS 2 installer\n";
                print ZSH "export PATH=\"$bin:\$PATH\"\n\n";
                close(ZSH);
            }
            elsif ($ENV{'SHELL'} =~ /bash/) {
                my $profile = (-f "$home/.bash_profile") ? "$home/.bash_profile" : "$home/.profile";
                unless (open(BASH, ">>$profile")) {
                    warn "Could not open $profile file: $!";
                    last;
                }
                print BASH "\n# added by alphaMELTS 2 installer\n";
                print BASH  "export PATH=\"$bin:\$PATH\"\n\n";
                close(BASH);
            }
        }

    }

    my $copy = ($windows) ? 'copy /Y' : 'cp -f'; # perl's copy requires a module to be installed

    @files = glob '*-*.pl';
    push @files, $program if ($program);

    $done = 1;
    for $file (@files) {
        if ($install ne $pwd) {
            $done = (system "$copy $file \"$install\"") ? '' : $done;
        }
        unless ($windows) {
            (chmod 0755, "$install/$file") || warn "Could not make all files executable!\n";
        }
    }
    unless ($done) {
        warn "Could not copy all executable files!\n";
        last;
    }

    print "Installation complete!";
    if (exists $ENV{'SUDO_USER'}) {
        print " Please run the install script again without 'sudo'.\n".
        "Check that the '$bin' directory is in your path." unless ($path);
    }
    else {
        print " Please log out and in and then run the install script again.\n".
        "Check that the '$bin' directory is in your path then choose 'n' to exit." unless ($path);
    }
    print "\n\n";

}

unless ($done) {
    print "\nInstallation aborted or may be incomplete!\n";
    print "Please log out and in and then and run the install script again.\n".
    "Check the '$bin' directory is in your path then choose 'n' to exit.\n" unless ($path);
}
print "Press return to finish.\n";
$_ = <STDIN>;
