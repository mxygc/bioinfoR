package Submitter;
use strict;
use Carp;
use SOAP::Lite;
use HTTP::Cookies;
use Exporter;

our @ISA = qw /Exporter/;
our @EXPORT_OK = qw /get_soap/;

sub get_soap {
    my ($auth_email, $input, $species, $id_type, $list_type, $list_name, $category) = @_;
    my ($check, $inputH, $input_id, $list, $all_species, @all_species, $species_idx);
    my $soap = SOAP::Lite->uri('http://service.session.sample')->proxy(
        'http://david.abcc.ncifcrf.gov/webservice/services/DAVIDWebService',
        cookie_jar => HTTP::Cookies->new( ignore_discard => 1 )
    );

    # user authentication by email address
    # For new user registration, go to http://david.abcc.ncifcrf.gov/webservice/register.htm
    $check = $soap->authenticate($auth_email)->result;
    print "User authentication: $check\n";
    die "email $auth_email not authenticated!\n" unless ( lc($check) eq "true" );
    open $inputH, "<", $input or croak "cannot open fine for input!";
    while( <$inputH> ) {
        $input_id .= $_;
    }
    close $inputH;
    #addList
    $list = $soap->addList( $input_id, $id_type, $list_name, $list_type )->result; 
    print "$list of input IDs were mapped\n";
    print "seems none mapped!\n" and exit(-1) if ( $list == 0 );
    # list all species names
    $all_species = $soap->getSpecies()->result;
    print "All possible species: $all_species\n";
    @all_species = split ",", $all_species;
    $species_idx = -1;
    for my $i (0 .. $#all_species) {
        if ( $all_species[$i] =~ m/$species/i ) {
            $species_idx = $i;
            last;
        }
    }
    # set species
    if ( $species_idx > -1 ) {
        $soap->setCurrentSpecies($species_idx);
    } else {
        $soap->setCurrentSpecies("0");
    print "no species specified, choosing the first option\n";
    }

    # confirm current species
    $species = $soap->getCurrentSpecies()->result;
    print "Current species: $species -> $all_species[$species]\n";

    # set and confirm category
    $category = $soap->setCategories($category)->result;
    print "Current categories: $category\n";
    return $soap;
}
1;
