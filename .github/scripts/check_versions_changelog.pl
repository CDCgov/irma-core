#!/usr/bin/env perl
use English qw(-no_match_vars);
use Carp qw(croak);
use strict;
use warnings;

local $RS = undef;
open my $fh, '<', 'CHANGELOG.md' or croak "Can't open CHANGELOG.md: $OS_ERROR";
my $changelog = <$fh>;
close $fh or croak "Can't close CHANGELOG.md: $OS_ERROR";

local $RS = "\n";
my $toml_version = ( split '#', qx(cargo pkgid) )[1];
chomp($toml_version);

if ( $changelog =~ /^## \[(.*?)\] - (\S+?)$/sm ) {
    my ( $version, $date ) = ( $1, $2 );

    if ( $date ne "TBD" ) {
        if ( $date !~ /^20\d{2}-[01]\d-[0-3]\d$/smx ) {
            die "Version $version has invalid date format: $date (expected yyyy-mm-dd or TBD)\n";
        }

        if ( $version ne $toml_version ) {
            die "Cargo.toml ($toml_version) mismatches changelog ($version / $date)!\n";
        }

    } elsif ( $toml_version !~ /dev$/smx ) {
        die "Cargo.toml version ($toml_version) should have a '-dev' suffix since the changelog is: ($version / $date)!\n";
    }

    if ( $changelog !~ /<!-- Versions -->.*?^\[$version\]:/sm ) {
        die "Version $version not linked!\n";
    }
}

