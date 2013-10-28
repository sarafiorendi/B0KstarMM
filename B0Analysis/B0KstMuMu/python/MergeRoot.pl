#!/usr/bin/perl


### GEN MC ###
#$dirName = "/nfs/data36/cms/dinardo/B0ToKstMuMu_GEN_NoFilter_MC_NTuples/" ;
#$dirName = "/nfs/data36/cms/dinardo/B0ToJPsiKst_GEN_NoFilter_MC_NTuples/" ;
#$dirName = "/nfs/data36/cms/dinardo/B0ToPsi2SKst_GEN_NoFilter_MC_NTuples/" ;
### Data ###
$dirName = "/bestman/storage/cms/store/user/dinardo/GridB0KstMuMu_Data2012A_01/" ;


$listcmd = "ls " . $dirName ;
@list = `$listcmd` ;

$nFiles = @list ;
$split = 700 ;
$nOut = int($nFiles / $split) + 1 ;

print "Total number of files: " . $nFiles . "\n" ;
print "Split into: " . $nOut . " files\n" ;

for ($count = 0; $count < $nOut; $count++)
{


### GEN MC ###
#    $fileOut = "B0ToKstMuMu_GEN_NoFilter_MC_NTuples" ;
#    $fileOut = "B0ToJPsiKst_GEN_NoFilter_MC_NTuples" ;
#    $fileOut = "B0ToPsi2SKst_GEN_NoFilter_MC_NTuples" ;
### Data ###
    $fileOut = "B0ToKstMuMu_Data2012A_NTuples_Merged" ;


    $cmd = "hadd " . $fileOut . "_" . $count . ".root " ;
    $Start = $split * $count ;
    $End = $Start + $split - 1 ;
    if ($count == ($nOut-1))
    {
	$End = $Start + ($nFiles % $split) - 1 ;
    }

    print "\nStarting file number: " . $Start . "\n" ;
    print "Ending file number: " . $End . "\n" ;

    foreach $file (@list[$Start..$End])
    {
	chomp $file ;
	$cmd .= $dirName . $file . " " ;
    }

    $fileOut .= "_" . $count . ".sh" ;

    print "\nFile output name: " . $fileOut . "\n" ;

    open(OUT, ">" . $fileOut) ;
    print OUT "$cmd" ;
    print OUT "\n" ;
    close(OUT) ;
}
