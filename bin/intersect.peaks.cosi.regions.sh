#!/bin/bash

#********
# USAGE *
#********

display_usage() { 
	echo -e "DESCRIPTION: provided a .bed file of cosi/cosit regions and a .bed or .bigBed file of peaks, it returns, for each region, the intersecting peaks\n"
	echo -e "\t--regions <bed file of desired regions>\n"
	echo -e "\t--peaks <bed/bigBed files of peaks> (if multiple files, must be comma-separated)\n"
	echo -e "\t--ext <bed/bigBed> (file extension of peaks; default: 'bigBed')\n"
	echo -e "\t--filter <yes/no> (i.e. whether to filter the peaks according to FC and pvalue; default: no)\n"
	echo -e "\t--FC <threshold for log2(FC)> (default: 0)\n"
	echo -e "\t--pvalue <threshold for -log10(p-value)> (default: 0)\n"
        echo -e "\t--outFile <output filename> (default: 'out.intersectBed.tsv')\n"
	echo -e "\t--outFolder <output folder> (default: cwd)\n"
	echo -e "\t--verbose <yes/no> (default: no)\n"
	echo -e "\t--keep <yes/no> (default: no)\n"
} 


if [[  $1 == "--help" ||  $1 == "-h" ]]
then
    	display_usage
        exit 0
fi


if [  $# -le 5  ]
then
	echo "ERROR: insufficient number of arguments"
    	display_usage
        exit 1
fi



#******************
# READING OPTIONS *
#******************

while [[ $# -gt 1 ]]; do

	key="$1"
	
	case $key in

	--regions)
	regions="$2"
	shift
	;;
    	
	--peaks)
	peaks="$2"
	shift
	;;

	--ext)
	ext="$2"
	shift
	;;

	--filter)
	filter="$2"
	shift
	;;
	
	--FC)
	FC="$2"
	shift
	;;

	--pvalue)
	pvalue="$2"
	shift
	;;

	--outFolder)
	outFolder="$2"
	shift
	;;

	--outFile)
	outFile="$2"
	shift
	;;
	
	--keep)
	keep="$2"
	;;

	--verbose)
	verbose="$2"
	;;
	*)
	
	;;
	esac
	shift
done


: ${ext:="bigBed"}
: ${filter:="no"}
: ${FC:=0}
: ${pvalue:=0}
: ${outFile:="out.intersectBed.tsv"}
: ${outFolder:="."}
: ${keep:="no"}
: ${verbose:="no"}




if [[ "$verbose" == "yes" ]]
then
	echo "Reading options .."
	echo "peaks bedFile =" "${peaks}"
	echo "file extension =" "${ext}"
	echo "threshold to log2(FC) =" "${FC}"
	echo "threshold to -log10(pvalue) =" "${pvalue}"
	echo "output folder: " "${outFolder}"
	echo "output file: " "${outFile}"
	echo -e "keep tmp files: " "${keep}\n"
fi




#********************************
# RETRIEVING PEAK CALLING FILES *
#********************************


if [[ "$verbose" == "yes" ]]
then
	echo -e "Retrieving peak calling files ..\n"
fi


for file in $(echo $peaks | tr "," "\n"); do
	if [[ "$ext" == "bigBed" ]]
	then
		bigBedToBed $file "$outFolder"/"$(basename $file .bigBed).bed"
		cat "$outFolder"/"$(basename $file .bigBed).bed" >> "$outFolder"/tmp
		rm "$outFolder"/"$(basename $file .bigBed).bed"
	else
		cat $file >> "$outFolder"/tmp
	fi
done



#******************
# FILTERING PEAKS *
#******************


if [[ "$filter" == "yes" ]]
then
	if [[ "$verbose" == "yes" ]]
	then
        	echo -e "Filtering peaks ..\n"
	fi
	
	cat "$outFolder"/tmp | \
	sort -k1,1 -k2,2n | \
	awk -v FC="$FC" -v pvalue="$pvalue" '
	BEGIN{FS="\t"}
	$7 >= FC && $8 >= pvalue' > "$outFolder"/filtered.peaks
	rm "$outFolder"/tmp

else
	mv "$outFolder"/tmp "$outFolder"/filtered.peaks

fi


#******************************************
# PEAK INTERSECTION WITH SELECTED REGIONS *
#******************************************


if [[ "$verbose" == "yes" ]]
then
    	echo -e "Intersecting peaks with selected regions ..\n"
fi


if [[ -e $regions ]]
then
	bedtools intersect -a $regions -b "$outFolder"/filtered.peaks -wao > "$outFolder"/"$outFile"
fi



#*********************
# REMOVING TMP FILES *
#*********************

if [[ "$keep" == "no" ]]
then
	if [[ "$verbose" == "yes" ]]
	then
		echo -e "Removing tmp files ..\n"
	fi
	
	rm "$outFolder"/filtered.peaks
fi

