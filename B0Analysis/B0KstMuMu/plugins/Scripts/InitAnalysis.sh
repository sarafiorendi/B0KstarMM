if [ "$#" -ne 2 ] ; then
    echo @@@ Error: Parameter missing @@@
    echo "Synopsis: InitAnalysis.sh analysis_path unset_display[0=false;1=true]"
else
    export ANALYPATH=$1

    export DATAYEAR=2012
    export DATADIR=/nfs/data37/cms/dinardo/Data"$DATAYEAR"B0KstMuMuResults

    echo @@@ Analysis environment variable: $ANALYPATH @@@
    echo @@@ Directory with data: $DATADIR @@@
    echo @@@ Data year: $DATAYEAR @@@

    if [ "$2" -eq 1 ]; then
	unset DISPLAY
	echo @@@ Unset DISPLAY @@@
    fi
fi
