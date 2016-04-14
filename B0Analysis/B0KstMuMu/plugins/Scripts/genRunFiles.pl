#!/usr/bin/perl

############################################################
# Simple program to generate files that run other programs #
############################################################

###################
# Input variables #
###################
$q2bin   = @ARGV[0] ;
$dirName = @ARGV[1] ;
$nfiles  = 400 ;
$toyDir  = "Toys_fromMCsig" ;


$fileOut = "myRun_bin" . $q2bin . ".sh" ;
#$fileOut = "myRun.sh" ;
print "@@@ File output name: " . $fileOut . " @@@\n" ;
open(OUT, ">" . $fileOut) ;


if ($dirName ne "")
{
    print OUT "mkdir " . $dirName . "\n" ;
    print OUT "cp ../ExtractYield " . $dirName . "\n" ;
    print OUT "cp ../../python/ParameterFile.txt " . $dirName . "\n" ;
    print OUT "cd " . $dirName . "\n" ;
}


for ($count = 0; $count < $nfiles; $count++)
{
#    $cmd = "./ExtractYield 881 " . "singleCand_B0ToJPsiKst_MC_NTuple_reduced.root" . " noEffCorr " . $q2bin . " " . $count . " " . "singleCand_B0ToJPsiKst_MC_NTuple_reduced_" . $q2bin . "_" . $count . ".root" ;

#    $cmd = "hadd singleCand_B0ToJPsiKst_MC_NTuple_reduced_" . $count . ".root singleCand_B0ToJPsiKst_MC_NTuple_reduced_3_" . $count . ".root singleCand_B0ToJPsiKst_MC_NTuple_reduced_5_" . $count . ".root" ;

#    $cmd = "./ExtractYield 86 CombBkgToyBin" . $q2bin . "_" . $count . ".root yesEffCorr " . $q2bin . " " . $count ;

#    $cmd = "hadd CombBkgToy_" . $count . ".root CombBkgToyBin0_" . $count . ".root CombBkgToyBin1_" . $count . ".root CombBkgToyBin2_" . $count . ".root CombBkgToyBin3_" . $count . ".root CombBkgToyBin5_" . $count . ".root CombBkgToyBin7_" . $count . ".root CombBkgToyBin8_" . $count . ".root" ;

#    $cmd = "hadd MCcocktail_1xData_All_" . $count . ".root singleCand_B0ToKstMuMu_MC_NTuple_" . $count . ".root CombBkgToy_" . $count . ".root singleCand_B0ToJPsiKst_MC_NTuple_reduced_" . $count . ".root singleCand_B0ToPsi2SKst_MC_NTuple_reduced_" . $count . ".root" ;

    $cmd = "Qsub -l lnxfarm -e -o cocktailMC_" . $q2bin . "_" . $count . ".log -N COCK" . $q2bin . $count . " ./ExtractYield 6 ../TestMultySamples/MCcocktail_1xData_All_" . $count . ".root yesEffCorr " . $q2bin . " " . $count ;

#    $tmp = $count + "1" ;
#    $cmd = "mv " . $dirName . "_" . $tmp . "_" . $q2bin . " " . $toyDir ;


    if ($count == 0)
    {
	print "Command: " . $cmd . "\n" ;
    }
    print OUT "$cmd" . "\n" ;
}


if ($dirName ne "")
{
    print OUT "cd ..\n" ;
}


close(OUT) ;
