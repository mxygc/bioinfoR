use strict;
use Carp;
use Getopt::Long;
use File::Basename;
use File::Path qw /mkpath/;
use File::Spec;

use Submitter qw /get_soap/;
use TermEnricher qw /get_term_enrich/;
use TermCluster qw /get_term_clust/;
use GeneCluster qw /get_gene_clust/;
use GeneAnnoter qw /get_gene_annot/;

my $input      = undef;
my $out_dir    = undef;
my $auth_email = "vishalrp\@uci.edu";
my $task       = "term_enrich,term_clust,gene_clust,gene_annot";
my $species    = "Homo sapiens";
my $id_type    = "AFFYMETRIX_3PRIME_IVT_ID";
my $list_name  = undef;
my $list_type  = "0";

my $category =
"OMIM_DISEASE,COG_ONTOLOGY,SP_PIR_KEYWORDS,UP_SEQ_FEATURE,GOTERM_BP_FAT,GOTERM_CC_FAT,GOTERM_MF_FAT,BBID,BIOCARTA,KEGG_PATHWAY,INTERPRO,PIR_SUPERFAMILY,SMART";

my $te_count = 2;
my $te_ease  = 0.1;

my $tc_similarity_gene_overlap = 3;      #overlap
my $tc_initial_group_member    = 3;      # initialSeed
my $tc_final_group_member      = 3;      # finalSeed
my $tc_multi_linkage           = 0.5;    # linkage
my $tc_similarity              = 50;     # kappa

my $gc_similarity_term_overlap = 4;      #overlap
my $gc_initial_group_member    = 4;      # initialSeed
my $gc_final_group_member      = 4;      # finalSeed
my $gc_multi_linkage           = 0.5;    # linkage
my $gc_similarity              = 35;     # kappa

sub usage {
    print <<EOF;
Usage: perl david.pl
       --input input_id.txt
       [--out_dir output_dir]\tdefault: where the input is
       [--auth_email email_addr]\tdefault: vishalrp\@uci.edu [warning: this email is only for demostration, please register your own address at DAVID.]
       [--task term_enrich[,term_clust][,gene_clust][,gene_annot]]\tdefault: 'term_enrich,term_clust,gene_clust,gene_annot'
       [--species species]\tdefault: Homo sapiens
       [--id_type id_type]\tdefault: AFFYMETRIX_3PRIME_IVT_ID
       [--list_name list_name]\tdefault: input_filename
       [--list_type 0]\tdefault: '0' (gene list)
       [--category categories_seperated_by_comma]\tdefault: 'OMIM_DISEASE,COG_ONTOLOGY,SP_PIR_KEYWORDS,UP_SEQ_FEATURE,GOTERM_BP_FAT,GOTERM_CC_FAT,GOTERM_MF_FAT,BBID,BIOCARTA,KEGG_PATHWAY,INTERPRO,PIR_SUPERFAMILY,SMART'
       [--te_count term_enrichment_count_cutoff]\tdefault: 2
       [--te_ease term_enrichment_EASE_score_cutoff]\tdefault: 0.1
       [--tc_similarity_gene_overlap term_clustering_similarity_gene_overlap]\tdefault: 3
       [--tc_initial_group_member term_clustering_initial_group_member]\tdefault: 3
       [--tc_final_group_member term_clustering_final_group_member]\tdefault: 3
       [--tc_multi_linkage term_clustering_multiple_linkage_cutoff]\tdefault: 0.5
       [--tc_similarity term_clustering_similarity_cutoff]\tdefault: 50 (=50%)
       [--gc_similarity_gene_overlap gene_clustering_similarity_gene_overlap]\tdefault: 4
       [--gc_initial_group_member gene_clustering_initial_group_member]\tdefault: 4
       [--gc_final_group_member gene_clustering_final_group_member]\tdefault: 4
       [--gc_multi_linkage gene_clustering_multiple_linkage_cutoff]\tdefault: 0.5
       [--gc_similarity gene_clustering_similarity_cutoff]\tdefault: 35 (=35%)
EOF
}
GetOptions(
    "input=s"                      => \$input,
    "out_dir=s"                    => \$out_dir,
    "auth_email=s"                 => \$auth_email,
    "task=s"                       => \$task,
    "species=s"                    => \$species,
    "id_type=s"                    => \$id_type,
    "list_type=i"                  => \$list_type,
    "list_name=s"                  => \$list_name,
    "category=s"                   => \$category,
    "te_count=i"                   => \$te_count,
    "te_ease=f"                    => \$te_ease,
    "tc_similarity_gene_overlap=f" => \$tc_similarity_gene_overlap,
    "tc_initial_group_member=i"    => \$tc_initial_group_member,
    "tc_final_group_member=i"      => \$tc_final_group_member,
    "tc_multi_linkage=f"           => \$tc_multi_linkage,
    "tc_similarity=i"              => \$tc_similarity,
    "gc_similarity_term_overlap=f" => \$gc_similarity_term_overlap,
    "gc_initial_group_member=i"    => \$gc_initial_group_member,
    "gc_final_group_member=i"      => \$gc_final_group_member,
    "gc_multi_linkage=f"           => \$gc_multi_linkage,
    "gc_similarity=i"              => \$gc_similarity,
  )
  or do {
    print("Error in command line arguments\n");
    &usage();
    exit;
  };

croak "no input file!" if ( !defined $input );
$input = File::Spec->rel2abs($input)
  if ( !File::Spec->file_name_is_absolute($input) );
$out_dir = ( &fileparse($input) )[1] if ( !defined $out_dir );

&mkpath($out_dir);
if ( !defined $list_name ) {
    $list_name = basename($input);
    $list_name =~ s/\.[^.]*$//;
}

print "Seems no registered email address provided, using 'vishalrp\@uci.edu' temporarily. To avoid using up its quota very soon, please register and use your own.\n" if ($auth_email eq 'vishalrp@uci.edu');

my $soap = &get_soap(
    $auth_email, $input,     $species, $id_type,
    $list_type,  $list_name, $category
);

my @tasks = split ",", $task;
for my $task (@tasks) {
    if ( $task eq "term_enrich" ) {
        &get_term_enrich( $soap, $out_dir, $list_name, $te_count, $te_ease );
    }
    elsif ( $task eq "term_clust" ) {
        &get_term_clust(
            $soap,                    $out_dir,
            $list_name,               $tc_similarity_gene_overlap,
            $tc_initial_group_member, $tc_final_group_member,
            $tc_multi_linkage,        $tc_similarity
        );
    }
    elsif ( $task eq "gene_clust" ) {
        &get_gene_clust(
            $soap,                       $out_dir,
            $list_name,                  $id_type,
            $gc_similarity_term_overlap, $gc_initial_group_member,
            $gc_final_group_member,      $gc_multi_linkage,
            $gc_similarity
        );
    }
    elsif ( $task eq "gene_annot" ) {
        &get_gene_annot( $soap, $out_dir, $list_name, $category );
    }
}
