#conda activate
#conda install -c bioconda itsx
itsdir=data/intermediate/ITSx
mkdir $itsdir
ITSx -i $itsdir/input_itsx.fasta -t fungi -o $itsdir
