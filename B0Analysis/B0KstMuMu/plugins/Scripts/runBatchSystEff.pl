#!/usr/bin/perl

####################################################################
# Program to fit with different pre-generated efficiency functions #
####################################################################

if (@ARGV == 3)
{
    system("source InitAnalysis.sh") ;

    ############################
    # Configuration parameters #
    ############################
    $dirEffRndGen = "../../efficiency/EffRndGenAnalyFilesSign_JPsi_Psi2S/" ;
    $dataFile     = $ENV{'DATADIR'} . "/Data" . $ENV{'DATAYEAR'} . "/SingleCand/singleCand_B0ToKstMuMu_Data" . $ENV{'DATAYEAR'} . "ABCD_NTuples.root" ;
    ############################

    print "Directory with data: $ENV{'DATADIR'} \n" ;
    print "Directory with randomly generated efficiencies: $dirEffRndGen \n" ;
    print "Data (or MC) file: $dataFile \n" ;

    $listcmd = "ls " . $dirEffRndGen ;
    @list = `$listcmd` ;
    
    $q2BinIndx = @ARGV[1] ;
    
    $cmd .= "unset DISPLAY\n" ;
    $cmd .= "rm FitSystematics_q2Bin_" . $q2BinIndx . ".txt" . "\n" ;
    
    $listStart = 0 ;
    $listEnd   = @list - 1 ;
    $listIndx  = 1 ;
    foreach $file (@list[$listStart..$listEnd])
    {
	chomp $file ;
	$execProg = ".././ExtractYield " . @ARGV[0] . " " . $dataFile . " yesEffCorrGen " . $q2BinIndx . " " . $dirEffRndGen . $file . " " . $listIndx . "\n" ;

	if (@ARGV[2] eq "false")
	{
	    $toRun = $execProg ;
	}
	else
	{
	    $toRun = "Qsub -l lnxfarm -e -o EffSys_" . $q2BinIndx . "_" . $listIndx . ".log -N ESYS" . $q2BinIndx . $listIndx . " " . $execProg ;
	}

	$cmd .= "echo " . $toRun ;
	$cmd .= $toRun ;
	$listIndx++ ;
    }

    open(OUT, ">runMultiExtractYield_q2Bin_" . $q2BinIndx . ".sh") ;
    print OUT "$cmd" ;
    print OUT "\n" ;
    close(OUT) ;
}
else
{
    print "Parameter missing:\n" ;
    print "./runBatchSystEff.pl [FitType] [q^2 bin to fit (0 - ...)] [runJobs[true / false]]\n" ;
}
