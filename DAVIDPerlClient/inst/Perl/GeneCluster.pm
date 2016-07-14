use strict;
use Carp;
use Exporter;

our @ISA       = qw /Exporter/;
our @EXPORT_OK = qw /get_gene_clust/;

sub get_gene_clust {
    my (
        $soap,                       $out_dir,
        $list_name,                  $id_type,
        $gc_similarity_term_overlap, $gc_initial_group_member,
        $gc_final_group_member,      $gc_multi_linkage,
        $gc_similarity
    ) = @_;
    my (
        $output,                        $geneClusterReport,
        @simpleGeneClusterRecords,      @simpleGeneClusterRecordKeys,
        @simpleGeneClusterRecordValues, @listRecords,
        $scoreValue,                    @listRecords_keys,
        @listRecords_values
    );
    $output = "$out_dir/gene_clust.txt";
    open my $outputH, ">", $output or croak "cannot open file for output!";

    $geneClusterReport =
      $soap->getGeneClusterReport( $gc_similarity_term_overlap,
        $gc_initial_group_member, $gc_final_group_member, $gc_multi_linkage,
        $gc_similarity );
    @simpleGeneClusterRecords = $geneClusterReport->paramsout;
    print "total gene clustering records: "
      . ( @simpleGeneClusterRecords + 1 ) . "\n";
    print "no clusters found!\n" and return (-1)
        if ( !defined( $geneClusterReport->result ) );
    @simpleGeneClusterRecordKeys   = keys %{ $geneClusterReport->result };
    @simpleGeneClusterRecordValues = values %{ $geneClusterReport->result };
    @listRecords                   = @{ $simpleGeneClusterRecordValues[0] };
    $scoreValue                    = $simpleGeneClusterRecordValues[2];

    print $outputH "Gene Group 1\tEnrichment Score:  $scoreValue\n";
    print $outputH "$id_type\tGene Name\n";

    for my $n ( 0 .. ( @listRecords - 1 ) ) {
        @listRecords_keys   = keys %{ $listRecords[$n] };
        @listRecords_values = values %{ $listRecords[$n] };
        print $outputH "$listRecords_values[2]\t$listRecords_values[1]\n";
    }

    for my $k ( 0 .. ( @simpleGeneClusterRecords - 1 ) ) {
        my $itr = $k + 2;
        @simpleGeneClusterRecordKeys = keys %{ $simpleGeneClusterRecords[$k] };
        @simpleGeneClusterRecordValues =
          values %{ $simpleGeneClusterRecords[$k] };
        $scoreValue = $simpleGeneClusterRecordValues[2];
        print $outputH "\nGene Group $itr\tEnrichment Score:  $scoreValue\n";
        print $outputH "$id_type\tGene Name\n";
        @listRecords = @{ $simpleGeneClusterRecordValues[0] };
        for my $n ( 0 .. ( @listRecords - 1 ) ) {
            @listRecords_keys   = keys %{ $listRecords[$n] };
            @listRecords_values = values %{ $listRecords[$n] };

    #print geneClusterReport "$listRecords_values[1]\t$listRecords_values[2]\n";
            print $outputH "$listRecords_values[2]\t$listRecords_values[1]\n";
        }
    }
    close $outputH;
    print "gene clustering results saved to $output\n";
}
1;
