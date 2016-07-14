use strict;
use Carp;
use Exporter;

our @ISA       = qw /Exporter/;
our @EXPORT_OK = qw /get_term_clust/;

sub get_term_clust {
    my (
        $soap,                    $out_dir,
        $list_name,               $tc_similarity_gene_overlap,
        $tc_initial_group_member, $tc_final_group_member,
        $tc_multi_linkage,        $tc_similarity
    ) = @_;
    my (
        $output,
        $termClusterReport,
        @simpleTermClusterRecords,
        @simpleTermClusterRecordKeys,
        @simpleTermClusterRecordValues,
        @chartRecords,

        %chartRecord,
        $categoryName,
        $termName,
        $listHits,
        $percent,
        $ease,
        $Genes,
        $listTotals,
        $popHits,
        $popTotals,
        $foldEnrichment,
        $bonferroni,
        $benjamini,
        $FDR
    );

    my $output = "$out_dir/term_clust.txt";

    # output chart report
    open my $outputH, ">", $output or croak "cannot open file for output!";

    $termClusterReport =
      $soap->getTermClusterReport( $tc_similarity_gene_overlap,
        $tc_initial_group_member, $tc_final_group_member, $tc_multi_linkage,
        $tc_similarity );

    @simpleTermClusterRecords = $termClusterReport->paramsout;
    print "total term clustering records: ", ( @simpleTermClusterRecords + 1 ), "\n";
    print "no clusters found! \n" and return (-1)
      if ( !defined( $termClusterReport->result ) );
    @simpleTermClusterRecordKeys   = keys %{ $termClusterReport->result };
    @simpleTermClusterRecordValues = values %{ $termClusterReport->result };
    @chartRecords                  = @{ $simpleTermClusterRecordValues[1] };
    print $outputH
"Annotation Cluster 1\tEnrichment Score:  $simpleTermClusterRecordValues[2]\n";
    print $outputH
"Category\tTerm\tCount\t%\tPvalue\tGenes\tList Total\tPop Hits\tPop Total\tFold Enrichment\tBonferroni\tBenjamini\tFDR\n";

    for my $j ( 0 .. ( @chartRecords - 1 ) ) {
        %chartRecord    = %{ $chartRecords[$j] };
        $categoryName   = $chartRecord{"categoryName"};
        $termName       = $chartRecord{"termName"};
        $listHits       = $chartRecord{"listHits"};
        $percent        = $chartRecord{"percent"};
        $ease           = $chartRecord{"ease"};
        $Genes          = $chartRecord{"geneIds"};
        $listTotals     = $chartRecord{"listTotals"};
        $popHits        = $chartRecord{"popHits"};
        $popTotals      = $chartRecord{"popTotals"};
        $foldEnrichment = $chartRecord{"foldEnrichment"};
        $bonferroni     = $chartRecord{"bonferroni"};
        $benjamini      = $chartRecord{"benjamini"};
        $FDR            = $chartRecord{"afdr"};
        print $outputH
"$categoryName\t$termName\t$listHits\t$percent\t$ease\t$Genes\t$listTotals\t$popHits\t$popTotals\t$foldEnrichment\t$bonferroni\t$benjamini\t$FDR\n";
    }
    for my $k ( 0 .. ( @simpleTermClusterRecords - 1 ) ) {
        my $itr = $k + 2;
        @simpleTermClusterRecordValues =
          values %{ $simpleTermClusterRecords[$k] };
        @chartRecords = @{ $simpleTermClusterRecordValues[1] };
        print $outputH
"\nAnnotation Cluster $itr\tEnrichment Score:  $simpleTermClusterRecordValues[2]\n";
        print $outputH
"Category\tTerm\tCount\t%\tPvalue\tGenes\tList Total\tPop Hits\tPop Total\tFold Enrichment\tBonferroni\tBenjamini\tFDR\n";
        for my $j ( 0 .. ( @chartRecords - 1 ) ) {
            %chartRecord    = %{ $chartRecords[$j] };
            $categoryName   = $chartRecord{"categoryName"};
            $termName       = $chartRecord{"termName"};
            $listHits       = $chartRecord{"listHits"};
            $percent        = $chartRecord{"percent"};
            $ease           = $chartRecord{"ease"};
            $Genes          = $chartRecord{"geneIds"};
            $listTotals     = $chartRecord{"listTotals"};
            $popHits        = $chartRecord{"popHits"};
            $popTotals      = $chartRecord{"popTotals"};
            $foldEnrichment = $chartRecord{"foldEnrichment"};
            $bonferroni     = $chartRecord{"bonferroni"};
            $benjamini      = $chartRecord{"benjamini"};
            $FDR            = $chartRecord{"afdr"};
            print $outputH
"$categoryName\t$termName\t$listHits\t$percent\t$ease\t$Genes\t$listTotals\t$popHits\t$popTotals\t$foldEnrichment\t$bonferroni\t$benjamini\t$FDR\n";
        }
    }
    close $outputH;
    print "term clustering results saved to $output\n";
}
1;
