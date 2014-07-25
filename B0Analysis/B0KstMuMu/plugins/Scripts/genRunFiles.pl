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

print OUT "mkdir " . @ARGV[1] . "\n" ;
print OUT "cp ../ExtractYield " . @ARGV[1] . "\n" ;
print OUT "cd " . @ARGV[1] . "\n" ;

for ($count = 0; $count < $nfiles; $count++)
{
#    $cmd = "./ExtractYield 86 CombBkgToyBin" . $q2bin . "_" . $count . ".root yesEffCorr " . $q2bin . " " . $count ;

#    $cmd = "hadd CombBkgToy_" . $count . ".root CombBkgToyBin0_" . $count . ".root CombBkgToyBin1_" . $count . ".root CombBkgToyBin2_" . $count . ".root CombBkgToyBin3_" . $count . ".root CombBkgToyBin5_" . $count . ".root CombBkgToyBin7_" . $count . ".root CombBkgToyBin8_" . $count . ".root" ;

#    $cmd = "hadd MCcocktail_1xData_All_" . $count . ".root singleCand_B0ToJPsiKst_MC_NTuple_reduced_0.root singleCand_B0ToPsi2SKst_MC_NTuple_reduced_0.root CombBkgToy_" . $count . ".root singleCand_B0ToKstMuMu_MC_NTuple_" . $count . ".root" ;

    $cmd = "Qsub -l lnxfarm -e -o cocktailMC_" . $q2bin . "_" . $count . ".log -N COCK" . $q2bin . $count . " ./ExtractYield 6 ../TestMultySamples/MCcocktail_1xData_All_" . $count . ".root yesEffCorr " . $q2bin . " " . $count ;

    if ($count == 0)
    {
	print "Command: " . $cmd . "\n" ;
    }
    print OUT "$cmd" . "\n" ;
}

close(OUT) ;
