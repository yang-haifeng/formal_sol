#!/usr/bin/perl

$pwd = `pwd` ;
chop($pwd) ;
print "Installing fml_sol in the directory:\n  '$pwd'\n";
$home = $ENV{"HOME"};
$bin  = $home . "/bin" ;
if(!(-e $bin)) {
    print "You must have a bin/ directory in your home directory\n" ;
    print "-----Shall I make $bin for you?\n" ;
    $input = <STDIN> ;
    print $input ;
    if($input=~/^[yY]/) {
        print "Creating $bin for you...\n" ;
        system("mkdir $bin") ;
    }   
}
$path = $ENV{"PATH"};
if(!($path =~ /$bin/)) {
    print "The $bin directory exists, but it not in the PATH environment variable\n" ;
    print "You must put the \n" ;
    print "$bin directory \n" ;
    print "in the path yourself (in the .tcshrc file if you use the tcsh shell...)\n" ;
    print "If you do it now, don't forget to type 'rehash'.\n" ;
}
print "  Creating a link 'fml_sol' in '$bin/'\n" ;
$fml_sol    = $pwd . "/fml_sol" ;
$fml_sollnk = $bin . "/fml_sol" ;
if(!(-e $fml_sollnk)) {
    print "------ Warning: file $fml_sollnk did not exist previously. You might want to type 'rehash'\n" ;
}
open(FILE,">$fml_sollnk") || die "Could not open file\n" ;
print FILE "#!/usr/bin/perl\n" ;
print FILE "system(\"$fml_sol \@ARGV\");" ;
close (FILE) ;
`chmod u+rwx $fml_sollnk` ;
