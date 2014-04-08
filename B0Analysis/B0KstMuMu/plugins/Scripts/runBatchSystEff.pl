####################################################################
# Program to fit with different pre-generated efficiency functions #
####################################################################

#!/usr/bin/perl

if (@ARGV == 2)
{
    system("source InitAnalysis.sh") ;

    ############################
    # Configuration parameters #
    ############################
    $dirEffRndGen = "../../efficiency/EffRndGenAnalyFiles/" ;
    $dataFile     = $ENV{'DATADIR'} . "/Data" . $ENV{'DATAYEAR'} . "/SingleCand/singleCand_B0ToKstMuMu_Data" . $ENV{'DATAYEAR'} . "ABCD_NTuples.root" ;
    $parFile      = "../../python/ParameterFile.txt" ;
    ############################

    print "Directory with data: $ENV{'DATADIR'} \n" ;
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
    $listIndx = 0 ;
    foreach $file (@list[$listStart..$listEnd])
    {
	chomp $file ;
	$str = ".././ExtractYield " . @ARGV[0] . " " . $dataFile . " EffCorrGenAnalyPDF " . $q2BinIndx . " " . $dirEffRndGen . $file . " " . $listIndx . " " . $parFile . "\n" ;
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
    print "./runBatchSystEff.pl [FitType] [q^2 bin to fit (0 - ...)]\n" ;
}
