#!/bin/bash

: <<'END'
  
    This file is part of NUQA.

    NUQA is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    NUQA is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with NUQA.  If not, see <http://www.gnu.org/licenses/>.
*/
END

#!/bin/bash
CWD="$(pwd)"

#Initialize parameters
files=()

#Default values:
k='21'
t='1'
dm='j'
Pre='nuqa'

while [[ $# -gt 0 ]]
do
	key="$1"

	case $key in
		-k|--kmer)
		k="$2"
		shift # past argument
		shift # past value
		;;
		-t|--threads)
		t="$2"
		shift # past argument
		shift # past value
		;;
		-d|--dist)
		dm="$2"
		shift # past argument
		shift # past value
		;;
		-p|--prefix)
		Pre="$2"
		shift # past argument
		shift # past value
		;;
		*)    # unknown option
		files+=("$1") # save it in an array for later
		shift # past argument
		;;
	esac
done




echo k     = "$k"
echo t     = "$t"
echo dm    = "$dm"
echo pre   = "$Pre"

#Determine file path and file names
#declare -a files=("$@")
filepath="${files[0]%/*}"
cd $filepath
files=("${files[@]##*/}")
files=("${files[@]%.*}")

#echo ${files[@]}

#Initialise intermediate/tmp files
declare -a counts=(0)
declare -a fastq=( "${files[@]/%/.fastq}" )
declare -a bc=( "${files[@]/%/.bc}" )
declare -a input=( "${files[@]/%/.jf}" )
declare -a inter=( "${files[@]/%/.jf.counts}" )
declare -a sorted=( "${files[@]/%/.jf.counts.sorted}" )

#Initialise parameters
fLen=${#files[@]}
Labs='Lab.txt'
suff='.txt'
n='norm'
cm='MasterCounts'
us='_'


#Assign Labels
rm $Labs
for ((l=0; l<$fLen; l++))
do
	printf ${files[l]}'\n' >> $Labs
done
wait
#: <<'END'

#jellyfish bloom counter
for ((y=0; y<$fLen; y++))
do
        jellyfish bc -m $k -s 500M -t $t -o ${bc[y]} ${fastq[y]}
done
wait

#jellyfish k-mer counting using bloom counter
for ((x=0; x<$fLen; x++))
do
        jellyfish count -m $k -s 100M -t $t --bc ${bc[x]} -o ${input[x]} ${fastq[x]}
done
wait
rm *.bc

#jellyfish dump k-mers and counts to textfile
for ((i=0; i<$fLen; i++))
do
        jellyfish dump -c ${input[i]} > ${inter[i]}
done
wait
rm *.jf

#sort files for merging
for ((j=0; j<$fLen; j++))
do
        (LC_ALL=C sort -o ${sorted[j]} ${inter[j]})
done
wait

#Find k-mer count totals for normalization
rm $n$k$suff
for ((l=0; l<$fLen; l++))
do
        cat ${sorted[$l]} | awk '{s+=$2} END {printf("%.0f\n", s)}' >> $n$k$suff
done
wait
rm *.counts

#merge files
$CWD/nuqa_merge.exe ${sorted[@]} > $cm$k$Pre$suff
#END

#distance calculations
$CWD/nuqa_distance.exe -n $n$k$suff -l $Labs -d $dm $cm$k$Pre$suff > $Pre$us$k$dm$suff

rm *.sorted



printf $Pre'_'$k$dm'.txt\nY' > $Pre'_input.txt'
#
rm outfile
phylip neighbor dnaml < $Pre'_input.txt'
cat outtree >> $Pre'_'$k$dm'.nw'
rm outtree

