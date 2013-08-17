####################################################################
# Program to fit with different pre-generated efficiency functions #
####################################################################

#!/usr/bin/perl

if (@ARGV == 2)
{
    ############################
    # Configuration parameters #
    ############################
    $dirEffRndGen = "../../Efficiency/EffRndGenAnalyFiles/" ;
    $dataFile     = "Candidates/Data/singleCand_B0ToKstMuMu_DataRRPRv4v5v6v1B_NTuples_Merged.root" ;
    $parFile      = "../python/ParameterFile.txt" ;
    ############################

    print "Directory with randomly generated efficiencies: $dirEffRndGen \n" ;
    print "Data (or MC) file: $dataFile \n" ;
    print "Parameter file: $parFile \n" ;

    $listcmd = "ls " . $dirEffRndGen ;
    @list = `$listcmd` ;
    
    $q2BinIndx = @ARGV[1] ;
    
    $cmd .= "unset DISPLAY\n" ;
    $cmd .= "rm FitSystematics_q2Bin_" . $q2BinIndx . ".txt" . "\n" ;
    
    $listStart = 0 ;
    $listEnd = @list - 1 ;
    $listIndx = 0;
    foreach $file (@list[$listStart..$listEnd])
    {
	chomp $file ;
	$str = "./ExtractYield " . @ARGV[0] . " " . $dataFile . " EffCorrGenAnalyPDF " . $q2BinIndx . " " . $dirEffRndGen . $file . " " . $listIndx . " " . $parFile . "\n" ;
	$cmd .= "echo " . $str ;
	$cmd .= $str ;
	$listIndx++ ;
    }
    
    open(OUT, ">RunMultiExtractYield_q2Bin_" . $q2BinIndx . ".sh") ;
    print OUT "$cmd" ;
    print OUT "\n" ;
    close(OUT) ;
}
else
{
    print "Parameter missing:\n" ;
    print "./RunBatchSystEff.pl [FitType] [q^2 bin to fit (0 - ...)]\n" ;
}
