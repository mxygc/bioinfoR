use strict;
use Carp;
use Exporter;

our @ISA         = qw /Exporter/;
our @EXPORTER_OK = qw /get_term_enrich/;

sub get_term_enrich {
    my ( $soap, $out_dir, $list_name, $te_count, $te_ease ) = @_;
    my $output = "$out_dir/term_enrich.txt";

    # output chart report
    open my $outputH, ">", $output or croak "cannot open file for output!";
    print $outputH
"Category\tTerm\tCount\t%\tPvalue\tGenes\tList Total\tPop Hits\tPop Total\tFold Enrichment\tBonferroni\tBenjamini\tFDR\n";

    # getChartReport
    my ( $chartReport, @chartRecords );
    my (
        %chartRecord, $categoryName, $termName,       $listHits,
        $percent,     $ease,         $Genes,          $listTotals,
        $popHits,     $popTotals,    $foldEnrichment, $bonferroni,
        $benjamini,   $FDR
    );
    $chartReport = $soap->getChartReport( $te_ease, $te_count );

    #    print $chartReport, "\n";
    @chartRecords = $chartReport->paramsout;

    #    print @chartRecords, "\n";
    unshift( @chartRecords, ( $chartReport->result ) );
    print "total term enrichment records: ", scalar @chartRecords, "\n";
    map {
        %chartRecord    = %{$_};
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
    } @chartRecords;
    close $outputH;
    print "term enrichment results saved to $output\n";
}
1;
