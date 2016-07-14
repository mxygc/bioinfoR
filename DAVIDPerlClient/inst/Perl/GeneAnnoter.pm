use strict;
use Carp;

#use YAML;
use Scalar::Util qw /reftype/;
use Exporter;

our @ISA       = qw /Exporter/;
our @EXPORT_OK = qw /get_gene_annot/;

sub get_gene_annot {
    my ( $soap, $out_dir, $list_name, $category ) = @_;
    my (
        $output,            $tableReport,  @tableRecords,
        @ids,               %annotations,  @tableRecordKeys,
        @tableRecordValues, @annots,       $id,
        $term,              $terms_string, @annot_terms,
        @annot_keys,        @annot_values
    );
    $output = "$out_dir/gene_annot.txt";

    # output chart report
    open my $outputH, ">", $output or croak "cannot open file for output!";
    print $outputH "ID\tGene Name\tSpecies";
    $category =~ s/,/\t/g;
    print $outputH "$category\n";

    $tableReport  = $soap->getTableReport();
    @tableRecords = $tableReport->paramsout;
    unshift @tableRecords, $tableReport->result;
    print "total gene annotation records: ", scalar @tableRecords, "\n";
    for (@tableRecords) {
        @ids               = ();
        %annotations       = ();
        @tableRecordKeys   = keys %{$_};
        @tableRecordValues = values %{$_};
        if ( ref( $tableRecordValues[3] ) eq "ARRAY" ) {
            for ( @{ $tableRecordValues[3] } ) {
                push @ids, $_->{"array"};
            }
        }
        elsif ( ref( $tableRecordValues[3] ) eq "HASH" ) {
            push @ids, $tableRecordValues[3]->{"array"};
        }
#        print "@ids\n";
        print $outputH join ",", @ids;
        print $outputH "\t";
        print $outputH "$tableRecordValues[2]\t";
        print $outputH "$tableRecordValues[4]\t";
        $DB::single = 1;
        if ( reftype( $tableRecordValues[1] ) eq "ARRAY" ) {
            @annots = @{ $tableRecordValues[1] };
        }
        elsif ( reftype( $tableRecordValues[1] ) eq "HASH" ) {
            push @annots, $tableRecordValues[1];
        }
        # create HashTable for each annotation
        for (@annots) {
            $id           = '';
            $term         = '';
            $terms_string = '';
            @annot_terms  = ();
            @annot_keys   = keys %{$_};
            @annot_values = values %{$_};
            if ( ref( $annot_values[0] ) eq "ARRAY" ) {
                @annot_terms = @{ $annot_values[0] };
            }
            elsif ( ref( $annot_values[0] ) eq "" ) {
                push @annot_terms, $annot_values[0];
            }
            $terms_string = join ",",
              map { ( $id, $term ) = split( /\$/, $_ ); $term } @annot_terms;
            $annotations{ $annot_values[1] } = $terms_string;
        }
        print $outputH join "\t",
          map { exists $annotations{$_} ? $annotations{$_} : "" } split ",",
          $category;
        print $outputH "\n";
    }
    close $outputH;
    print "gene annot results saved to $output\n";
}
1;

