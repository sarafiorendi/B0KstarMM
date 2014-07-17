#!/usr/bin/perl

############################################################
# Simple program to generate files that run other programs #
############################################################

$q2bin  = @ARGV[0] ;
$nfiles = 100 ;

$fileOut = "myRunBin" . $q2bin . ".sh" ;
#$fileOut = "myRun.sh" ;
print "@@@ File output name: " . $fileOut . " @@@\n" ;
open(OUT, ">" . $fileOut) ;

for ($count = 0; $count < $nfiles; $count++)
{
#    $cmd = "./ExtractYield 86 CombBkgToyBin" . $q2bin . "_" . $count . ".root yesEffCorr " . $q2bin . " " . $count ;

#    $cmd = "hadd CombBkgToy_" . $count . ".root CombBkgToyBin0_" . $count . ".root CombBkgToyBin1_" . $count . ".root CombBkgToyBin2_" . $count . ".root CombBkgToyBin3_" . $count . ".root CombBkgToyBin5_" . $count . ".root CombBkgToyBin7_" . $count . ".root CombBkgToyBin8_" . $count . ".root" ;

#    $cmd = "hadd MCcocktail_1xData_All_" . $count . ".root singleCand_B0ToJPsiKst_MC_NTuple_reduced_0.root singleCand_B0ToPsi2SKst_MC_NTuple_reduced_0.root CombBkgToy_" . $count . ".root singleCand_B0ToKstMuMu_MC_NTuple_" . $count . ".root" ;

    $cmd = "Qsub -l lnxfarm -e -o cocktailMC_" . $q2bin . "_" . $count . ".log -N C" . $q2bin . $count . " .././ExtractYield 6 TestMultySamples/MCcocktail_1xData_All_" . $count . ".root yesEffCorr " . $q2bin . " " . $count ;

    print OUT "$cmd" ;
    print OUT "\n" ;
}

close(OUT) ;
