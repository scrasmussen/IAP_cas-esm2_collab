package Build::ChemNamelist;

#-------------------------------------------------------------------------------------
# generates species lists for chemitry namelist settings
#-------------------------------------------------------------------------------------
# Date         Contributor      Modification
#-------------------------------------------------------------------------------------
# 26 Jan 2011  Francis Vitt     Created 
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------

use strict;
use Exporter;
use FindBin qw($Bin);
use lib "$Bin/perl5lib";
use Build::ChemPreprocess qw(get_species_list);

our @ISA = qw(Exporter);
our @EXPORT = qw(set_dep_lists);
our $VERSION = 1.00;

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
sub set_dep_lists
{
    my ( $cfgdir, $chem_pkg, $chem_src_dir, $print_lvl ) = @_;

    my ( $gas_wetdep_list, $aer_wetdep_list, $drydep_list, $mam_drydep_list ) ;

    my @species_list = get_species_list($chem_src_dir);
    if ($print_lvl>=2) {print "Chemistry species : @species_list \n" ;}

    $gas_wetdep_list = get_gas_wetdep_list( $cfgdir, $print_lvl, @species_list );

    $aer_wetdep_list = get_aer_wetdep_list( $cfgdir, $print_lvl, @species_list );

    $drydep_list = get_drydep_list( $cfgdir, $print_lvl, @species_list );

    if (($chem_pkg =~ '_mam') || ($chem_pkg == 'trop_strat_new')) {

	$mam_drydep_list = get_mamdep_list( $cfgdir, $print_lvl, @species_list );

	if (length($aer_wetdep_list)>2) {
	    $aer_wetdep_list = $aer_wetdep_list . "," . $mam_drydep_list;
	} else {
	    $aer_wetdep_list = $mam_drydep_list;
	}

    }

    return (  $gas_wetdep_list, $aer_wetdep_list, $drydep_list, $mam_drydep_list );
}

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
sub get_drydep_list
{
    my ($cfg_dir,$print_lvl,@species_list) = @_;

    my $master_file = "$cfg_dir/namelist_files/master_drydep_list.xml";

    my $list = get_dep_list($master_file,$print_lvl,@species_list);

    if ($print_lvl>=2) {print " dry dep list : $list  \n" ;}

    return ($list);
}

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
sub get_aer_wetdep_list
{
    my ($cfg_dir,$print_lvl,@species_list) = @_;

    my $master_file = "$cfg_dir/namelist_files/master_aer_wetdep_list.xml";

    my $list = get_dep_list($master_file,$print_lvl,@species_list);

    if ($print_lvl>=2) {print " aer wet dep list : $list  \n" ;}
    return ($list);
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
sub get_gas_wetdep_list
{
    my ($cfg_dir,$print_lvl,@species_list) = @_;

    my $master_file = "$cfg_dir/namelist_files/master_gas_wetdep_list.xml";

    my $list = get_dep_list($master_file,$print_lvl,@species_list);

    if ($print_lvl>=2) {print " gas wet dep list : $list  \n" ;}

    return ($list);
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
sub get_mamdep_list
{
    my ($cfg_dir,$print_lvl,@species_list) = @_;

    my $master_file = "$cfg_dir/namelist_files/master_mam_dep_list.xml";

    my $list = get_dep_list($master_file,$print_lvl,@species_list);

    if ($print_lvl>=2) {print " gas MAM dep list : $list  \n" ;}

    return ($list);
}

#-------------------------------------------------------------------------------
sub get_dep_list
{
    my ($master_file,$print_lvl,@species_list) = @_;

    if ($print_lvl>=2){ print "Using chemisty master list file $master_file \n"; }

    my @master_list = read_master_list_file($master_file);

    my $list = '';
    my $first = 1; my $pre = "";
    foreach my $name (sort @species_list) {
	foreach my $item (@master_list) {
	    if ($name eq $item) { 
		$list .= $pre .  quote_string($name) ;
                if ($first) { $pre = ","; $first = 0; }
	    }
	}
    }

    if ( length($list)<1 ) {$list = quote_string(' ') ;}

    return ($list);
}

#-------------------------------------------------------------------------------
sub read_master_list_file
{
    my ($master_file) = @_;

    require XML::Lite;

    my @master_list ;
    my $xml = XML::Lite->new($master_file);
    my $root = $xml->root_element();
    my @children = $root->get_children();
    foreach my $child (@children) {
	my $content = $child->get_content();
	my @list = split( ('\s+|\s*,+\s*') ,$content);
	foreach my $item (@list) {
	    if ( length( $item) > 0 ){
		push ( @master_list, $item );
	    }
	}
    }

    return (@master_list);
}
#-------------------------------------------------------------------------------

sub quote_string {
    my $str = shift;
    $str =~ s/^\s+//;
    $str =~ s/\s+$//;
    unless ($str =~ /^['"]/) {        #"'
        $str = "\'$str\'";
    }
    return $str;
}
1; # to appease require 
