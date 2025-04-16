use English qw(-no_match_vars);

open $logfile,'<', 'CHANGELOG.md' or die "Cannot open 'CHANGELOG.md': $OS_ERROR\n";

local $RS = "## [";
<$logfile>; # skip
my $changes = <$logfile>;
chomp($changes);

my $version = '';
if ( $changes =~ /^([\d.]+)\]/ ) {
    $version = $1;
}

print STDOUT '## Release Notes for v[',$changes,"\n";

local $RS = "\n";
while (my $line = <$logfile>) {
    if ($line =~ /^\[$version\]/) {
        print STDOUT $line,"\n";
        last;
    }
}
