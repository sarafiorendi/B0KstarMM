#!/usr/bin/perl

############################################################
# Simple program to generate files that run other programs #
############################################################

###################
# Input variables #
###################
$ID      = @ARGV[0] ;
$q2bin   = @ARGV[1] ;
$dirName = @ARGV[2] ;
$nfiles  = 400 ;


$fileOut = "myRun_" . $ID . "_bin" . $q2bin . ".sh" ;
#$fileOut = "myRun.sh" ;
print "@@@ File output name: " . $fileOut . " @@@\n" ;
open(OUT, ">" . $fileOut) ;

print OUT "mkdir " . $dirName . "\n" ;
print OUT "cp ../ExtractYield " . $dirName . "\n" ;
print OUT "cp ../../python/ParameterFile.txt " . $dirName . "\n" ;
print OUT "cd " . $dirName . "\n" ;


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
