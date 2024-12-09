#!/bin/bash

#! sbatch directives begin here ###############################
#! Name of the job:
#SBATCH -J ConvertChr
#SBATCH --account=BISC020662


#! hmem has high memory with only 8 nodes
#SBATCH --nodes=1
#SBATCH --mem=10G
#SBATCH --ntasks=1

THREADS=$1

module load lib/hdf5/1.10.6-gcc

## STEPS ########
EXTRACT4DSITES=$1
RUNPHYLOFITONCHUNKS=$2
RUNPHYLOFIT=$3
RUNEXTRACTWHOLEMAF=$4
RUNPHYLOP=$5
RUNPHASTCONS=$6
RUNPHASTCONSNONHELICONIUS=$7
RUNPHASTCONSHELICONIUS=$8
RUNPHASTCONSONALLSPECIES=$9
RUNCNEEPREP=${10}
PREPAREPHYLOACC=${11}
RUNPHYLOACC=${12}
RUNENRICHMENT=${13}
RUNPHYLOPREGDOMAINS=${14}
DATAFORPLOTS=${15}
EXTRACTUCE=${16}
#################

WORKDIR=$WORK/HeliconiniiProject/HeliconGenomeAlignmentAnnotation/
DMELPROTDB=$WORK/Databases/Flybase/Dmel.Prot.FlybaseRefseq.fasta
SWISSUNIPROT=$WORK/Databases/SwissProt/uniprot_sprot
THREADS=30
HAL=$WORKDIR/63Nymphalidae.Cactus-v1.5.hal
MAF=$WORKDIR/63Nymphalidae.Cactus-v1.5.maf
TREE=`cat $WORKDIR/IqTree2_2.ML.63speciesFinal.nwk.tree`
CATHapp=$APPS2/cath-tools-genomescan/apps
FUNFAM=$WORK/Databases/FunFam/funfam-hmm3-v4_2_0.lib

mkdir -p $WORKDIR/CNEEselectionAnalysis
cd $WORKDIR/CNEEselectionAnalysis

if [[ $EXTRACT4DSITES == 1 ]]; then

	SPECIES='Eisa Hmel'
	. $APPS2/cactus-bin-v2.5.1/cactus_env/bin/activate
	
	for SP in $SPECIES; do
	
		ASSEMBLY=`ls -1 /user/work/tk19812/HeliconiniiProject/CompleteAssemblies/AssembledGenomes/$SP.assembly.v*.fasta`
		GFF3=`ls -1 $WORKDIR/$SP.v*.annotation.CAT.gff3`
	
		if [[ ! -s $SP.bed ]]; then
			echo -e "Step 1. Get 4-fold degenerate sites based on galGal4 NCBI annotations"
			echo -e "Get the coding regions"
			gffread --no-pseudo -C -T -o $SP.gtf $GFF3
			grep -P "\tCDS\t" $SP.gtf > ${SP}_filt.gtf
		
			echo -e "Convert to GenePred format..."
			gtfToGenePred ${SP}_filt.gtf $SP.gp
			echo -e "... and to bed"
			genePredToBed $SP.gp $SP.bed
		fi
	
		if [[ ! -s ${SP}_4d.bed ]]; then
			echo -e "Extract 4D sites from a subsample of species"
			echo -e "hal4dExtract --conserved --hdf5InMemory $HAL $SP $SP.bed ${SP}_4d.bed"
			hal4dExtract --conserved --hdf5InMemory $HAL $SP $SP.bed ${SP}_4d.bed 
		else
			ls -l ${SP}_4d.bed
		fi
	done

	#echo -e "hal4dExtract is computing don't kill it"
	#wait
	echo -e "hal4dExtract is done!"
	
	for SP in $SPECIES; do
		if [[ ! -s $SP.extract_maf/neut4d_input_$SP.maf ]]; then
			mkdir -p $SP.extract_maf
			cd $SP.extract_maf
	
			cp ../${SP}_4d.bed .
			if [[ $SP == "Eisa" ]]; then
				#SPLIST="Dple,Mcin,Smor,Dpha,Djun,Herd,Hcha,Hsar,Haoe,Hdor,Hmel"
				SPLIST="Smor,Dpha,Djun,Herd,Hcha,Hsar,Haoe,Hdor,Hmel"
			elif [[ $SP == "Hmel" ]]; then
				SPLIST="Smor,Dpha,Djun,Eisa,Herd,Hcha,Hsar,Haoe,Hdor"
			fi
			hal2mafMP.py --numProc 48 --targetGenomes $SPLIST --refGenome $SP --noAncestors --noDupes --refTargets ${SP}_4d.bed $HAL neut4d_input_$SP.maf &
			#hal2maf --refGenome $SP --noAncestors --noDupes --refTargets ${SP}_4d.bed ../$HAL neut4d_input_$SP.maf
			cd ..
		fi
	done
	wait

fi

if [[ $RUNPHYLOFITONCHUNKS == 1 ]]; then

	. $APPS2/cactus-bin-v2.5.1/cactus_env/bin/activate

	SP='Eisa'

	mkdir -p $SP.extract_maf.chunks && cd $SP.extract_maf.chunks
	cp ../$SP.extract_maf/${SP}_4d.bed .

	ASSEMBLYFAI=`ls -1 $WORK/HeliconiniiProject/CompleteAssemblies/AssembledGenomes/$SP.assembly.v*.fasta.fai`
	cat $ASSEMBLYFAI | cut -f1 | shuf | while read scf; do
		if [[ -s $scf.ok ]]; then
			ls -l $scf.*.63way.mod
		elif [[ -s $scf.processing ]]; then
			rm $scf.processing
		else
			echo -e "`date` $scf is processing" > $scf.processing
			grep -w $scf ${SP}_4d.bed > $scf.ToSplit
			if [[ -s $scf.ToSplit ]]; then
				echo -e "Split $scf"
				rm -f $chunk.temp.
				split -l 50 $scf.ToSplit $scf.ToSplit.
				chunks=`ls -1 $scf.ToSplit.*`
				echo $chunks
				for chunk in $chunks; do
					out=`echo $chunk | sed s/ToSplit.//`
					if [[ ! -s $out.maf ]]; then
						echo -e "Convert $chunk"
						hal2mafMP.py --numProc 5 --refGenome $SP --noAncestors --noDupes --refTargets $chunk $HAL $out.maf &
#						hal2maf --refGenome $SP --noAncestors --noDupes --refTargets $chunk $HAL $out.maf 
					fi
				done
	
				wait

				for chunk in $chunks; do
					out=`echo $chunk | sed s/ToSplit.//`
					if [[ -s $out.maf ]]; then
						if [[ ! -s $out.63way.mod ]]; then
							NTAXA=`grep '^s' $out.maf | cut -f2 | cut -f1 -d'.' | sort -u | wc -l`
		
							if [ $NTAXA -gt 62 ]; then
								mafSorter --m $out.maf --seq $SP.$scf > $out.4d.$SP.Sorted.maf
								msa_view $out.4d.$SP.Sorted.maf --in-format MAF --out-format SS > $out.4d.$SP.ss
								time phyloFit --subst-mod SSREV --precision HIGH --tree $TREE --out-root $out.63way --sym-freqs --log $out.neut4d.log $out.4d.$SP.ss & 
							fi
						fi
					fi
				done
				
				wait
			fi
			mv $scf.processing $scf.ok
		fi
	done

	ALLMODS=`ls -1 *.63way.mod`
	for mod in $ALLMODS; do
		if [[ ! -s $mod ]]; then
			echo -e "Remove $mod"
			ls -l $mod
			rm $mod
		fi
	done

	ls *.63way.mod > 63Way.cons.txt
	phyloBoot --read-mods '*63Way.cons.txt' --output-average ave_63way.mod
fi

if [[ $RUNPHYLOFIT == 1 ]]; then

	. $APPS2/cactus-bin-v2.5.1/cactus_env/bin/activate

	SPECIES='Eisa'

	for SP in $SPECIES; do
		ASSEMBLYFAI=`ls -1 $WORK/HeliconiniiProject/CompleteAssemblies/AssembledGenomes/$SP.assembly.v*.fasta.fai`
		cd $SP.extract_maf.chunks

		echo -e "Update the neutral model"
		phyloFit --subst-mod SSREV --init-model ave_63way.mod --out-root neut_${SP}_63way.Updated --sym-freqs --log Updating_neut4d_63-way_$SP.log ../$SP.extract_maf/neut4d_63-way_input_$SP.maf &> neut_${SP}_63way-Updated_63way.out
		echo -e "Neutral model UPDATED!!!! you should be happy about that!"

		cp neut_${SP}_63way.Updated.mod neut_${SP}_final.mod

		echo -e "Step 4. Adjust background GC content to be reflective of the average GC content across the non-ancestral genomes in the alignment"
		###code added to get GC content for each genome, by sampling every 100 bp
		mkdir -p baseComp
		cd baseComp
		for TARGET in $(halStats --genomes $HAL); do
			#output is fraction_of_As fraction_of_Gs fraction_of_Cs fraction_of_Ts
			if [[ ! -s $TARGET.basecomp ]]; then
				halStats --baseComp $TARGET,100 $HAL > $TARGET.basecomp
			fi
		done
		cd ..

		echo -e "get average gc content in non-ancestral genomes and update models:"
		GC=$(cat baseComp/*.basecomp | awk '{SUM+=$2;SUM+=$3;print SUM/42}' | tail -n 1)

		echo $GC

		modFreqs neut_${SP}_final.mod $GC > neut_${SP}_corrected.mod
	done
fi

if [[ $RUNEXTRACTWHOLEMAF == 1 ]]; then

	SPECIES="Herd"

	. $APPS2/cactus-bin-v2.5.1/cactus_env/bin/activate

	for SP in $SPECIES; do
		ASSEMBLYFAI=`ls -1 $WORK/HeliconiniiProject/CompleteAssemblies/AssembledGenomes/$SP.assembly.v*.fasta.fai`
		mkdir -p $SP.extract_maf && cd $SP.extract_maf

		cat $ASSEMBLYFAI | cut -f1 | shuf | while read scf; do
			if [[ ! -s ${scf}.maf ]]; then
				if [[ ! -s ${scf}.processing ]]; then
					echo $scf
					echo "`date` $scf.processing" > $scf.processing
					echo -e "hal2mafMP.py --numProc 50 --refGenome $SP --noAncestors --noDupes --refSequence $scf $HAL $scf.maf"
					hal2mafMP.py --numProc 50 --refGenome $SP --noAncestors --noDupes --refSequence $scf $HAL $scf.maf
					mv $scf.processing $scf.ok
				elif [[ -s ${scf}.processing ]]; then
					cat ${scf}.processing
				elif [[ -s ${scf}.ok ]]; then
					echo "$scf DONE!"
				fi
			fi
		done
		cd ..
	done

fi

if [[ $RUNPHYLOP == 1 ]]; then

	#SPECIES="Eisa Hmel"
	SPECIES=Hmel

	. $APPS2/cactus-bin-v2.5.1/cactus_env/bin/activate

	echo -e "### RUN PHYLOP ##"
	#run halPhyloPMP.py with 12 processors per on each neutral version
	#use the _corrected version of each model
	for SP in $SPECIES; do
		ASSEMBLYFAI=`ls -1 $WORK/HeliconiniiProject/CompleteAssemblies/AssembledGenomes/$SP.assembly.v*.fasta.fai`
		mkdir -p $SP.PhyloP && cd $SP.PhyloP

		if [[ ! -s neut_Eisa_corrected.mod ]]; then
			cp ../Eisa.extract_maf.chunks/neut_Eisa_corrected.mod neut_Eisa_corrected.mod
		fi

		#echo -e "halPhyloPMP.py --numProc 12 $HAL $SP.neutralMod ${SP}_phyloP_ver1.wig"
		#halPhyloPMP.py --numProc 24 $HAL $SP neut_${SP}_corrected.mod ${SP}_phyloP.wig &> halPhyloPMP_$SP.log &
		cat $ASSEMBLYFAI | cut -f1 | while read scf; do
			if [[ ! -s ${SP}_${scf}.phyloP.bw ]]; then
				if [[ ! -s $scf.phyloP.processing ]]; then
					echo $scf.phyloP.processing > $scf.phyloP.processing
					#hal2mafMP.py --numProc 30 --refGenome $SP --noAncestors --noDupes --refSequence $scf $HAL $scf.maf
					if [[ -s ../$SP.extract_maf/$scf.ok ]]; then
						echo -e "phyloP --mode CONACC --wig-scores --method LRT neut_${SP}_corrected.mod ../$SP.extract_maf/$scf.maf"
						phyloP --mode CONACC --wig-scores --method LRT neut_${SP}_corrected.mod ../$SP.extract_maf/$scf.maf > ${SP}_${scf}.phyloP.wig
						wigToBigWig ${SP}_${scf}.phyloP.wig $ASSEMBLYFAI ${SP}_${scf}.phyloP.bw
						gzip ${SP}_${scf}.phyloP.wig
					fi
					rm $scf.phyloP.processing
				fi
			fi
		done

		cat *_phyloP.wig > $SP.WholeGenome_PhyloP.wig
		wigToBigWig $SP.WholeGenome_PhyloP.wig $ASSEMBLYFAI $SP.WholeGenome_PhyloP.bw

		cd ..
	done

#	for SP in $SPECIES; do	
#	#finally also run tree version
#		mkdir $SP.halTreePhyloP && cd $SP.halTreePhyloP
#		echo -e "halTreePhyloP.py --prec 5 --bigWig --numProc 24 $HAL ../$SP.neutralMod/neut_${SP}_corrected.mod"
#		halTreePhyloP.py --prec 5 --bigWig --numProc 24 $HAL neut_${SP}_corrected.mod . &> halTreePhyloP.log
#		tail halTreePhyloP.log
#		exit 
#		cd ..
#	done
#	wait
	## END RUN PHYLOP -- NOTE THESE RESULTS NEED TO BE CHECKED ###
fi


if [[ $RUNPHASTCONS == 1 ]]; then

	. $APPS2/cactus-bin-v2.5.1/cactus_env/bin/activate

	### RUNNING PHASTCONS ###
	#to run phastCons, we need to take a slightly different approach as there is no direct interface with hal
	#so the first step is to export the MAFs that we want, in this case starting with two: chicken, ostrich
	#for each MAF, we then run phastCons
	#sources: https://genome.ucsc.edu/cgi-bin/hgTrackUi?db=hg38&g=cons100way and http://compgen.cshl.edu/phast/phastCons-HOWTO.html

	SPECIES="Eisa"
	

	for SP in $SPECIES; do
		echo $SP
		ASSEMBLYFAI=`ls -1 $WORK/HeliconiniiProject/CompleteAssemblies/AssembledGenomes/$SP.assembly.v*.fasta.fai`
		BEDGENES=`ls -1 $WORKDIR/$SP.v*.annotation.CAT.bed`
		GFF3=`ls -1 $WORKDIR/$SP.v*.annotation.CAT.gff3`

		#mkdir -p $SP.PhastCons && cd $SP.PhastCons
		mkdir -p $SP.PhastCons.NonHeliconius && cd $SP.PhastCons.NonHeliconius

		if [[ ! -s neut_Eisa_corrected.mod ]]; then
			cp ../Eisa.extract_maf.chunks/neut_Eisa_corrected.mod neut_Eisa_corrected.mod
#			cp ../Eisa.PhastCons/ave*cons.mod .
		fi	

		if [[ ! -s ave.cons.mod ]] && [[ ! -s ave.noncons.mod ]]; then
			cat $ASSEMBLYFAI | cut -f1 | shuf | while read scf; do
				echo $scf
				if [[ ! -s $scf.maf.processing ]]; then
					if [[ ! -s $scf.DownSmpl.maf ]]; then
						echo -e "hal2maf on $scf"
						echo -e "hal2maf on $scf" > $scf.maf.processing
						NUMPROC=12

						#echo -e "hal2mafMP.py --numProc $NUMPROC --refGenome $SP --noAncestors --noDupes --targetGenomes Dple,Bany,Mcin,Jcoe,Smor,Pdid,Dpha,Diul,Ptel,Djun,Avcr,Avpe,Avfl,Elam,Evib,Eali,Elyb,Etal,Hcha,Herd,Haoe,Hmel --refSequence $scf $HAL $scf.DownSmpl.maf"
						#hal2mafMP.py --numProc $NUMPROC --refGenome $SP --noAncestors --noDupes --targetGenomes Dple,Bany,Mcin,Jcoe,Smor,Pdid,Dpha,Diul,Ptel,Djun,Avcr,Avpe,Avfl,Elam,Evib,Eali,Elyb,Etal,Hcha,Herd,Haoe,Hmel --refSequence $scf $HAL $scf.DownSmpl.maf

						if [[ -s $scf.DownSmpl.maf ]]; then
							mv $scf.maf.processing $scf.maf.ok
						fi 
					fi
	
					if [[ -s $scf.maf.ok ]]; then
						if [[ ! -s $scf.cons.mod ]] && [[ ! -s $scf.noncons.mod ]]; then
							echo -e "phastCons $scf.DownSmpl.maf neut_${SP}_corrected.mod --estimate-rho $scf --expected-length 12 --target-coverage 0.25"
							phastCons $scf.DownSmpl.maf neut_${SP}_corrected.mod --estimate-rho $scf --expected-length 12 --target-coverage 0.25 &> $scf.rhoEstimate_$SP.log
							echo -e "phastCons --score --most-conserved $scf.mostcons.bed $scf.DownSmpl.maf $scf.cons.mod,$scf.noncons.mod > $scf.scores.wig"
							phastCons --score --most-conserved $scf.mostcons.bed $scf.DownSmpl.maf $scf.cons.mod,$scf.noncons.mod > $scf.scores.wig
							wigToBigWig $scf.scores.wig $ASSEMBLYFAI $scf.scores.bw
							gzip $scf.scores.bw & 
						fi
					fi
				fi
			done
wait
exit
			if [[ ! -s ave.cons.mod ]]; then
				rm Eisa00*.mod
				ls *.cons.mod > cons.txt
				ls *.noncons.mod > noncons.txt
				phyloBoot --read-mods '*cons.txt' --output-average ave.cons.mod
				phyloBoot --read-mods '*noncons.txt' --output-average ave.noncons.mod
			fi
		fi
exit
		if [[ ! -s WholeGenome${SP}PhastCons.scores.bw ]]; then
		cat $ASSEMBLYFAI | cut -f1 | while read scf; do
				if [[ ! -s $scf.mostcons.bed ]]; then
					if [[ ! -s $scf.phastCons.processing ]]; then
						echo $scf
						echo $scf > $scf.phastCons.processing
						#hal2mafMP.py --numProc 12 --refGenome $SP --noAncestors --noDupes --targetGenomes Dple,Bany,Mcin,Jcoe,Smor,Pdid,Dpha,Diul,Ptel,Djun,Avcr,Avpe,Avfl,Elam,Evib,Eali,Elyb,Etal,Hcha,Herd,Haoe,Hmel --refSequence $scf $HAL $scf.DownSmpl.maf
						phastCons --score --most-conserved $scf.mostcons.bed $scf.DownSmpl.maf ave.cons.mod,ave.noncons.mod > $scf.scores.wig 
						wigToBigWig $scf.scores.wig $ASSEMBLYFAI $scf.scores.bw
						rm $scf.phastCons.processing
					fi
				fi
			done
exit
			cat $SP*.scores.wig > WholeGenome${SP}PhastCons.scores.wig
			wigToBigWig WholeGenome${SP}PhastCons.scores.wig $ASSEMBLYFAI WholeGenome${SP}PhastCons.scores.bw
		fi
	done
fi

if [[ $RUNPHASTCONSNONHELICONIUS == 1 ]]; then

    . $APPS2/cactus-bin-v2.5.1/cactus_env/bin/activate

    ### RUNNING PHASTCONS ###
    #to run phastCons, we need to take a slightly different approach as there is no direct interface with hal
    #so the first step is to export the MAFs that we want, in this case starting with two: chicken, ostrich
    #for each MAF, we then run phastCons
    #sources: https://genome.ucsc.edu/cgi-bin/hgTrackUi?db=hg38&g=cons100way and http://compgen.cshl.edu/phast/phastCons-HOWTO.html

    SP="Eisa"


	echo $SP
	ASSEMBLYFAI=`ls -1 $WORK/HeliconiniiProject/CompleteAssemblies/AssembledGenomes/$SP.assembly.v*.fasta.fai`
	BEDGENES=`ls -1 $WORKDIR/$SP.v*.annotation.CAT.bed`
	GFF3=`ls -1 $WORKDIR/$SP.v*.annotation.CAT.gff3`

	mkdir -p $SP.PhastCons.NonHeliconius && cd $SP.PhastCons.NonHeliconius

	if [[ ! -s neut_Eisa_corrected.mod ]]; then
		cp ../Eisa.extract_maf.chunks/neut_Eisa_corrected.mod neut_Eisa_corrected.mod
	fi

	if [[ ! -s ave.cons.mod ]] && [[ ! -s ave.noncons.mod ]]; then
		cat $ASSEMBLYFAI | cut -f1 | shuf | while read scf; do
			if [[ ! -s $scf.maf.processing ]]; then
				if [[ ! -s $scf.DownSmpl.maf ]]; then
					echo -e "hal2maf on $scf"
					echo -e "hal2maf on $scf" > $scf.maf.processing
					NUMPROC=12

					echo -e "hal2mafMP.py --numProc $NUMPROC --refGenome $SP --noAncestors --noDupes --targetGenomes Dple,Bany,Mcin,Jcoe,Smor,Pdid,Dpha,Diul,Ptel,Djun,Avcr,Avpe,Avfl,Elam,Evib,Eali,Elyb,Etal --refSequence $scf $HAL $scf.DownSmpl.maf"
					hal2mafMP.py --numProc $NUMPROC --refGenome $SP --noAncestors --noDupes --targetGenomes Dple,Bany,Mcin,Jcoe,Smor,Pdid,Dpha,Diul,Ptel,Djun,Avcr,Avpe,Avfl,Elam,Evib,Eali,Elyb,Etal --refSequence $scf $HAL $scf.DownSmpl.maf

					if [[ -s $scf.DownSmpl.maf ]]; then
						mv $scf.maf.processing $scf.maf.ok
					fi
				fi

				#if [[ -s $scf.maf.ok ]]; then
				#	echo -e "phastCons $scf.DownSmpl.maf neut_${SP}_corrected.mod --estimate-rho $scf --expected-length 12 --target-coverage 0.25"
				#	phastCons $scf.DownSmpl.maf neut_${SP}_corrected.mod --estimate-rho $scf --expected-length 12 --target-coverage 0.25 &> $scf.rhoEstimate_$SP.log &
				#fi
			fi
		done
exit
		rm Eisa00*.mod
		ls *.cons.mod > cons.txt
		ls *.noncons.mod > noncons.txt
		phyloBoot --read-mods '*cons.txt' --output-average ave.NonHel.cons.mod
		phyloBoot --read-mods '*noncons.txt' --output-average ave.NonHel.noncons.mod
	fi

exit
	if [[ ! -s WholeGenome${SP}PhastCons.scores.bw ]]; then
		cat $ASSEMBLYFAI | cut -f1 | while read scf; do
			if [[ ! -s $scf.mostcons.bed ]]; then
				if [[ ! -s $scf.phastCons.processing ]]; then
					echo $scf
					echo $scf > $scf.phastCons.processing
					phastCons --score --most-conserved $scf.mostcons.bed $scf.DownSmpl.maf ave.cons.mod,ave.noncons.mod > $scf.scores.wig
					wigToBigWig $scf.scores.wig $ASSEMBLYFAI $scf.scores.bw
					rm $scf.phastCons.processing
				fi
			fi
		done

		cat $SP*.scores.wig > WholeGenome${SP}PhastCons.scores.wig
		wigToBigWig WholeGenome${SP}PhastCons.scores.wig $ASSEMBLYFAI WholeGenome${SP}PhastCons.scores.bw
	fi
fi

if [[ $RUNPHASTCONSHELICONIUS == 1 ]]; then

    . $APPS2/cactus-bin-v2.5.1/cactus_env/bin/activate

    ### RUNNING PHASTCONS ###
    #to run phastCons, we need to take a slightly different approach as there is no direct interface with hal
    #so the first step is to export the MAFs that we want, in this case starting with two: chicken, ostrich
    #for each MAF, we then run phastCons
    #sources: https://genome.ucsc.edu/cgi-bin/hgTrackUi?db=hg38&g=cons100way and http://compgen.cshl.edu/phast/phastCons-HOWTO.html

    SP="Hmel"


    echo $SP
    ASSEMBLYFAI=`ls -1 $WORK/HeliconiniiProject/CompleteAssemblies/AssembledGenomes/$SP.assembly.v*.fasta.fai`
    BEDGENES=`ls -1 $WORKDIR/$SP.v*.annotation.CAT.bed`
    GFF3=`ls -1 $WORKDIR/$SP.v*.annotation.CAT.gff3`

    mkdir -p $SP.PhastCons.Heliconius && cd $SP.PhastCons.Heliconius

    if [[ ! -s neut_Eisa_corrected.mod ]]; then
        cp ../Eisa.extract_maf.chunks/neut_Eisa_corrected.mod neut_Eisa_corrected.mod
    fi  

    if [[ ! -s ave.cons.mod ]] && [[ ! -s ave.noncons.mod ]]; then
        cat $ASSEMBLYFAI | cut -f1 | shuf | while read scf; do
            if [[ ! -s $scf.maf.processing ]]; then
                if [[ ! -s $scf.DownSmpl.maf ]]; then
                    echo -e "hal2maf on $scf"
                    echo -e "hal2maf on $scf" > $scf.maf.processing
                    NUMPROC=12

                    echo -e "hal2mafMP.py --numProc $NUMPROC --refGenome $SP --noAncestors --noDupes --rootGenome Heliconius --refSequence $scf $HAL $scf.DownSmpl.maf"
                    hal2mafMP.py --numProc $NUMPROC --refGenome $SP --noAncestors --noDupes --rootGenome Heliconius --refSequence $scf $HAL $scf.DownSmpl.maf

                    if [[ -s $scf.DownSmpl.maf ]]; then
                        mv $scf.maf.processing $scf.maf.ok
                    fi  
                fi  

                #if [[ -s $scf.maf.ok ]]; then
                #   echo -e "phastCons $scf.DownSmpl.maf neut_${SP}_corrected.mod --estimate-rho $scf --expected-length 12 --target-coverage 0.25"
                #   phastCons $scf.DownSmpl.maf neut_${SP}_corrected.mod --estimate-rho $scf --expected-length 12 --target-coverage 0.25 &> $scf.rhoEstimate_$SP.log &
                #fi 
            fi  
        done
exit
        rm Eisa00*.mod
        ls *.cons.mod > cons.txt
        ls *.noncons.mod > noncons.txt
        phyloBoot --read-mods '*cons.txt' --output-average ave.NonHel.cons.mod
        phyloBoot --read-mods '*noncons.txt' --output-average ave.NonHel.noncons.mod
    fi    
        
exit
    if [[ ! -s WholeGenome${SP}PhastCons.scores.bw ]]; then
        cat $ASSEMBLYFAI | cut -f1 | while read scf; do
            if [[ ! -s $scf.mostcons.bed ]]; then
                if [[ ! -s $scf.phastCons.processing ]]; then
                    echo $scf
                    echo $scf > $scf.phastCons.processing
                    phastCons --score --most-conserved $scf.mostcons.bed $scf.DownSmpl.maf ave.cons.mod,ave.noncons.mod > $scf.scores.wig
                    wigToBigWig $scf.scores.wig $ASSEMBLYFAI $scf.scores.bw
                    rm $scf.phastCons.processing
                fi  
            fi  
        done

        cat $SP*.scores.wig > WholeGenome${SP}PhastCons.scores.wig
        wigToBigWig WholeGenome${SP}PhastCons.scores.wig $ASSEMBLYFAI WholeGenome${SP}PhastCons.scores.bw
    fi
fi

if [[ $RUNPHASTCONSONALLSPECIES == 1 ]]; then

	. $APPS2/cactus-bin-v2.5.1/cactus_env/bin/activate

	### RUNNING PHASTCONS ###
	#to run phastCons, we need to take a slightly different approach as there is no direct interface with hal
	#so the first step is to export the MAFs that we want, in this case starting with two: chicken, ostrich
	#for each MAF, we then run phastCons
	#sources: https://genome.ucsc.edu/cgi-bin/hgTrackUi?db=hg38&g=cons100way and http://compgen.cshl.edu/phast/phastCons-HOWTO.html
	
	SP="Hmel"
	
	
	echo $SP
	ASSEMBLYFAI=`ls -1 $WORK/HeliconiniiProject/CompleteAssemblies/AssembledGenomes/$SP.assembly.v*.fasta.fai`
	BEDGENES=`ls -1 $WORKDIR/$SP.v*.annotation.CAT.bed`
	GFF3=`ls -1 $WORKDIR/$SP.v*.annotation.CAT.gff3`
	
	mkdir -p $SP.PhastCons.AllSpecies && cd $SP.PhastCons.AllSpecies
	
	if [[ ! -s neut_Eisa_corrected.mod ]]; then
		cp ../neut_Eisa_corrected.mod neut_Eisa_corrected.mod
	fi

	if [[ ! -s ave.Heliconiinae.cons.mod ]] && [[ ! -s ave.Heliconiinae.noncons.mod ]]; then
#		cat $ASSEMBLYFAI | cut -f1 | shuf | while read scf; do
#			if [[ ! -s $scf.maf.processing ]]; then
#				if [[ ! -s $scf.Heliconiinae.maf ]]; then
#					echo -e "hal2maf on $scf"
#					echo -e "hal2maf on $scf" > $scf.maf.processing
#					NUMPROC=12
#
# 					echo -e "hal2mafMP.py --numProc $NUMPROC --refGenome $SP --noAncestors --noDupes --rootGenome Heliconiinae --refSequence $scf $HAL $scf.Heliconiinae.maf"
# 					hal2mafMP.py --numProc $NUMPROC --refGenome $SP --noAncestors --noDupes --rootGenome Heliconiinae --refSequence $scf $HAL $scf.Heliconiinae.maf
#
# 					if [[ -s $scf.Heliconiinae.maf ]]; then
# 						mv $scf.maf.processing $scf.maf.ok
#					fi
#				fi
#
#				if [[ -s $scf.maf.ok ]]; then
#					echo -e "phastCons $scf.Heliconiinae.maf neut_${SP}_corrected.mod --estimate-rho $scf --expected-length 12 --target-coverage 0.25"
#					phastCons $scf.Heliconiinae.maf neut_Eisa_corrected.mod --estimate-rho $scf --expected-length 12 --target-coverage 0.25 &> $scf.rhoEstimate_$SP.log
#					gzip $scf.Heliconiinae.maf
#				fi 
#			fi
#		done

		rm cons.txt noncons.txt

		awk '{ if ($2 > 100000) print $1 }' $ASSEMBLYFAI | while read scf; do
	        ls $scf.cons.mod >> cons.txt
	        ls $scf.noncons.mod >> noncons.txt
		done
	        phyloBoot --read-mods '*cons.txt' --output-average ave.Heliconiinae.cons.mod
    	    phyloBoot --read-mods '*noncons.txt' --output-average ave.Heliconiinae.noncons.mod
    fi

    if [[ ! -s WholeGenome${SP}PhastCons.scores.bw ]]; then
        cat $ASSEMBLYFAI | cut -f1 | while read scf; do
            if [[ ! -s $scf.maf.ok ]]; then
                if [[ ! -s $scf.phastCons.processing ]]; then
                    echo $scf
                    echo $scf > $scf.phastCons.processing
					rm -r ./js
					cactus-hal2maf ./js $HAL $scf.maf --refGenome $SP --noAncestors --onlyOrthologs --noDupes --refSequence $scf --chunkSize 100000
					rm -r ./js
                    phastCons --score --most-conserved $scf.mostcons.bed $scf.maf ave.Heliconiinae.cons.mod,ave.Heliconiinae.noncons.mod > $scf.scores.wig
                    wigToBigWig $scf.scores.wig $ASSEMBLYFAI $scf.scores.bw
					gzip $scf.scores.wig
                    rm $scf.phastCons.processing
                fi
            fi
        done

		cat $ASSEMBLYFAI | cut -f1 | while read scf; do
			ls -l $scf.mostcons.bed
		done

		zcat $SP*.scores.wig.gz > WholeGenome${SP}PhastCons.scores.wig
		wigToBigWig WholeGenome${SP}PhastCons.scores.wig $ASSEMBLYFAI WholeGenome${SP}PhastCons.scores.bw
    fi
fi

if [[ $RUNCNEEPREP == 1 ]]; then 

	. $APPS2/cactus-bin-v2.5.1/cactus_env/bin/activate

	#SPECIES="Eisa"
	SPECIES="Hmel"

	for SP in $SPECIES; do	
		echo $SP
		ASSEMBLYFAI=`ls -1 $WORK/HeliconiniiProject/CompleteAssemblies/AssembledGenomes/$SP.assembly.v*.fasta.fai`
		BEDGENES=`ls -1 $WORKDIR/$SP.v*.annotation.CAT.bed`
		GFF3=`ls -1 $WORKDIR/$SP.v*.annotation.CAT.gff3`

		#cd $SP.PhastCons
		cd $SP.PhastCons.AllSpecies

		if [[ ! -s ${GFF3%.gff3}.Exons.bed ]] && [[ ! -s ${GFF3%.gff3}.CDSonly.bed ]]; then
			echo -e "Extracting CDS and Exons"
			awk '{ if ( $3 == "CDS" ) print $0 }' $GFF3 | bedtools sort -i - | bedtools merge -i - > ${GFF3%.gff3}.CDSonly.bed
			awk '{ if ( $3 == "exon" ) print $0 }' $GFF3 | bedtools sort -i - | bedtools merge -i - > ${GFF3%.gff3}.Exons.bed
		fi

		if [[ ! -s All-ConsRegions_intgCNEEs.bed ]]; then
			cat $ASSEMBLYFAI | cut -f1 | while read scf; do
				if [[ ! -s $scf.Filtered.mostcons.codingCEs.bed ]]; then
	
					NCEs=`wc -l $scf.mostcons.bed | cut -f1 -d' '`
					echo $NCEs

					if [ "$NCEs" -gt 0 ]; then

						echo -e "Merge flanking regions 5bp apart for $scf"
						bedtools merge -o first,mean -c 4,5 -d 5 -i $scf.mostcons.bed | awk '{ if ( $5 > 0 ) print $0 }' > $scf.Filtered.mostcons.bed
		
						echo -e "Removing CEs shorter than 50bp"
						awk '{ if ( $3-$2+1 > 50 ) print $0}' $scf.mostcons.bed | awk '{ if ( $5 > 0 ) print $0 }' > $scf.Filtered.mostcons.bed
						wc -l $scf.mostcons.bed
						wc -l $scf.Filtered.mostcons.bed
		
						echo -e "Extracting genic CEs"
						bedtools intersect -u -a $scf.Filtered.mostcons.bed -b $GFF3 > $scf.Filtered.mostcons.gCEs.bed
						wc -l $scf.Filtered.mostcons.gCEs.bed
		
						echo -e "Extracting intergenic CNEEs within genes - (All the others)"
						ExcludeBedRecords.py $scf.Filtered.mostcons.gCEs.bed $scf.Filtered.mostcons.bed > $scf.Filtered.mostcons.intgCNEEs.bed
						wc -l $scf.Filtered.mostcons.intgCNEEs.bed
		
						echo -e "Extracting exonic CEs"
						bedtools intersect -u -a $scf.Filtered.mostcons.gCEs.bed -b ${GFF3%.gff3}.Exons.bed > $scf.Filtered.mostcons.exonicCEs.bed
						wc -l $scf.Filtered.mostcons.exonicCEs.bed
		
						echo -e "Extracting intronic CNEEs - (All the others)"
						ExcludeBedRecords.py $scf.Filtered.mostcons.exonicCEs.bed $scf.Filtered.mostcons.gCEs.bed > $scf.Filtered.mostcons.intronicCEs.bed
						wc -l $scf.Filtered.mostcons.intronicCEs.bed
					
						echo -e "Extracting UTR CEs"
						bedtools intersect -u -a $scf.Filtered.mostcons.exonicCEs.bed -b ${GFF3%.gff3}.CDSonly.bed > $scf.Filtered.mostcons.codingCEs.bed
						wc -l $scf.Filtered.mostcons.codingCEs.bed
		
						echo -e "Extracting CDS CEs - All the others"
						ExcludeBedRecords.py $scf.Filtered.mostcons.codingCEs.bed $scf.Filtered.mostcons.exonicCEs.bed > $scf.Filtered.mostcons.UTRCEs.bed
						wc -l $scf.Filtered.mostcons.UTRCEs.bed
		
						echo -e "##############################"
						sleep 10s
					fi
				fi
			done
	
			cat *.Filtered.mostcons.bed | bedtools sort -i - > All-ConsRegions_Merged.bed
			cat *.Filtered.mostcons.intgCNEEs.bed | bedtools sort -i - > All-ConsRegions_intgCNEEs.bed
			cat *.Filtered.mostcons.gCEs.bed | bedtools sort -i - > All-ConsRegions_gCEs.bed
			cat *.Filtered.mostcons.intronicCEs.bed | bedtools sort -i - > All-ConsRegions_intronicCEs.bed
			cat *.Filtered.mostcons.UTRCEs.bed | bedtools sort -i - > All-ConsRegions_UTRCEs.bed
			cat *.Filtered.mostcons.codingCEs.bed | bedtools sort -i - > All-ConsRegions_codingCEs.bed
			wc -l All-ConsRegions*.bed
		fi

		echo -e "Extract and filter out CNEEs"
		echo -e "ExtractCEs.sh $SP $GFF3 $ASSEMBLYFAI $HAL"
		#ExtractCEs.sh $SP $GFF3 $ASSEMBLYFAI $HAL 

		cat $ASSEMBLYFAI | cut -f1 | while read scf; do
			cd $scf
			cat *.CNEE.bed > ../$scf.intCNEEs.bed
#			#cp *.fasta ../AllFastas
			cd ..

			NCE=`wc -l $scf.intCNEEs.bed | cut -f1 -d' ' &`
			FASTAS=`ls -1 AllFastasIntrons/$scf*.fasta | wc -l &`
			wait

			if [ $NCE -eq $FASTAS ]; then
				echo $NCE $FASTAS
				if [[ ! -s $scf.tar.gz ]]; then
					tar cvf $scf.tar.gz $scf
					ls -l $scf.tar.gz 
				fi
				#rm -r $scf
			fi
		done	
		cd ..
	done
fi

if [[ $PREPAREPHYLOACC == 1 ]]; then

	SP='Hmel'
	#cd $SP.PhastCons
	cd $SP.PhastCons.AllSpecies

	#cat *.intCNEEs.bed | bedtools sort -i - > $SP.phyloacc.allCNEEs.bed

	#PreparePhyloAccInput.py -b $SP.phyloacc.allCNEEs.bed -t $WORKDIR/IqTree2_2.ML.63speciesFinal.nwk.tree -d AllFastasIntrons -o $SP.phyloacc.Input.allCNEE.bed -f $SP.phyloacc.Input.allCNEE.fasta > $SP.phyloacc.RenamedIntCNEEs.bed

	. ~/.bashrc

	conda activate phyloacc-env

	# Gregariuos
#	phyloacc.py -a $SP.phyloacc.Input.allCNEE.fasta -b $SP.phyloacc.Input.allCNEE.bed \
#		-m neut_Eisa_corrected.mod -o $SP.phyloacc-allCNEEs.Gregarious \
#		-t "Hwal;Hbur;Hxan;Hdor;Hheb;Haoe;Hhec;Hsap;Hhew;Hcon;Hele;Hant;Hsar;Hleu;Hric;Hdem;Hert;Evib;Djun" \
#		-c "Ptel;Diul;Dpha;Pdid;Hbes;Hnum;Hism;Hheu;Htim;Hpac;Hcyd;Hmel;Heth;Helv;Hpar;Hatt;Hhel;Hnat;Hege;Htel;Hhor;Hcly;Herd;Hpet;Hhim;Hlat;Heet;Hher;Hper;Hcha;Hhie;Etal;Elyb;Eali;Elam;Eisa;Avfl;Avcr;Avpe" \
#		-g "Dple;Bany;Mcin;Jcoe;Smor" -n 10 -batch 500 -j 5 -part . -r st

	# Heliconius stem
#	phyloacc.py -a $SP.phyloacc.Input.allCNEE.fasta -b $SP.phyloacc.Input.allCNEE.bed \
#		-m neut_Eisa_corrected.mod -o $SP.phyloacc-allCNEEs \
#		-t "Htel;Hcly;Hhor;Hhec;Hher;Hpet;Herd;Heet;Hlat;Hhim;Hcha;Hper;Hert;Hdem;Hric;Hleu;Hsar;Hant;Hele;Hcon;Hhew;Hsap;Haoe;Hheb;Hhie;Hxan;Hdor;Hege;Hbur;Hwal;Hnat;Hbes;Hnum;Hism;Heth;Hpar;Helv;Hatt;Hhel;Hmel;Hcyd;Hpac;Hheu;Htim" \
#		-c "Pdid;Dpha;Diul;Ptel;Djun;Avfl;Avpe;Avcr;Eisa;Elam;Evib;Eali;Etal;Elyb" \
#		-g "Dple;Bany;Mcin;Jcoe;Smor" -n 10 -batch 500 -j 5 -part . -r st #-l scOG.All.Astral.Rooted.tree #

	snakemake -p -s  $WORKDIR/CNEEselectionAnalysis/$SP.PhastCons.AllSpecies/$SP.phyloacc-allCNEEs.Gregarious/phyloacc-job-files/snakemake/run_phyloacc.smk \
		--configfile $WORKDIR/CNEEselectionAnalysis/$SP.PhastCons.AllSpecies/$SP.phyloacc-allCNEEs.Gregarious/phyloacc-job-files/snakemake/phyloacc-config.yaml \
		--profile    $WORKDIR/CNEEselectionAnalysis/$SP.PhastCons.AllSpecies/$SP.phyloacc-allCNEEs.Gregarious/phyloacc-job-files/snakemake/profiles/slurm_profile \
		--dryrun
fi

if [[ $RUNPHYLOACC == 1 ]]; then

	SP='Hmel'
	PAWD=$SP.PhastCons.AllSpecies
	cd $PAWD
	PHYLOACCDIR=$SP.phyloacc-allCNEEs.Gregarious

	. ~/.bashrc
	conda activate phyloacc-env

	for i in $(seq 1 297); do
		if [[ ! -s $WORKDIR/CNEEselectionAnalysis/$SP.PhastCons/phyloacc/phyloacc-job-files/phyloacc-output/$i-phyloacc-st-out/$i-phyloacc.log ]]; then
			rm $WORKDIR/CNEEselectionAnalysis/$SP.PhastCons/phyloacc/phyloacc-job-files/cfgs/$i.ok
		fi
		if [[ -s $WORKDIR/CNEEselectionAnalysis/$PAWD/$PHYLOACCDIR/phyloacc-job-files/cfgs/$i.ok ]]; then
			echo -e "`date` batch $i is DONE"
		elif [[ -s $WORKDIR/CNEEselectionAnalysis/$PAWD/$PHYLOACCDIR/phyloacc-job-files/cfgs/$i.running ]]; then
			echo -e "`date` batch $i is running"
		else
			echo -e "`date` batch $i Started!"
			echo -e "`date` batch $i" > $WORKDIR/CNEEselectionAnalysis/$PAWD/$PHYLOACCDIR/phyloacc-job-files/cfgs/$i.running
			PhyloAcc-ST $WORKDIR/CNEEselectionAnalysis/$PAWD/$PHYLOACCDIR/phyloacc-job-files/cfgs/$i-st.cfg &> $WORKDIR/CNEEselectionAnalysis/$PAWD/$PHYLOACCDIR/phyloacc-job-files/phyloacc-output/$i-phyloacc-st-out/$i-phyloacc.log
			mv $WORKDIR/CNEEselectionAnalysis/$PAWD/$PHYLOACCDIR/phyloacc-job-files/cfgs/$i.running $WORKDIR/CNEEselectionAnalysis/$PAWD/$PHYLOACCDIR/phyloacc-job-files/cfgs/$i.ok
		fi
	done

	phyloacc_post.py -i $WORKDIR/CNEEselectionAnalysis/$PAWD/$PHYLOACCDIR -o $WORKDIR/CNEEselectionAnalysis/$PAWD/$PHYLOACCDIR-Output

#	ParsePhyloACCoutput.py -i $WORKDIR/CNEEselectionAnalysis/$PAWD/$PHYLOACCDIR-Output/results/elem_lik.txt -b $SP.phyloacc.RenamedIntCNEEs.bed
fi

if [[ $RUNENRICHMENT == 1 ]]; then

	SP='Hmel'
	ASSEMBLYFAI=`ls -1 $WORK/HeliconiniiProject/CompleteAssemblies/AssembledGenomes/$SP.assembly.v*.fasta.fai`
	#BEDGENES=`ls -1 $WORK/HeliconiniiProject/HeliconGenomeAlignmentAnnotation/UPDATEannotations/$SP.v*.annotation.CAT.bed`
	GFF3=`ls -1 $WORK/HeliconiniiProject/HeliconGenomeAlignmentAnnotation/UPDATEannotations/$SP.v*.annotation.CAT.gff3`
	#CNEES=$WORKDIR/CNEEselectionAnalysis/Eisa.PhastCons/Eisa.AccCNEEs.M1.Heliconius.bed
	#ALLCNEES=$WORKDIR/CNEEselectionAnalysis/Eisa.PhastCons/Eisa.phyloacc.AllCNEEs.bed
	#aCNEES=$WORKDIR/CNEEselectionAnalysis/$SP.PhastCons.AllSpecies/$SP.phyloacc.RenamedIntCNEEs.M1.Gregarious.bed
	aCNEES=$WORKDIR/CNEEselectionAnalysis/$SP.PhastCons.AllSpecies/$SP.phyloacc.RenamedIntCNEEs.M1.StrictGregarious.bed 
	ALLCNEES=$WORKDIR/CNEEselectionAnalysis/$SP.PhastCons.AllSpecies/$SP.phyloacc.RenamedIntCNEEs.bed
	EnrichmentDIR=$SP.Enrichment.StrictCNEEs

	. $APPS2/goatools/venv/bin/activate

	OUTPUTTBL=$WORK/HeliconiniiProject/HeliconGenomeAlignmentAnnotation/CNEEselectionAnalysis/$SP.PhastCons.AllSpecies/$SP.phyloacc-allCNEEs.Gregarious-Output/results/elem_lik.txt

	if [ ! -s $aCNEES ]; then
		#ParsePhyloACCoutputGregariousProg.py -i $OUTPUTTBL -b $ALLCNEES
		ParsePhyloACCoutputGregariousProgStrict.py -i $OUTPUTTBL -b $ALLCNEES
	fi

	mkdir -p $EnrichmentDIR && cd $EnrichmentDIR

	if [[ ! -s $SP.GeneGenomicRegions.bed ]]; then
		gffread $GFF3 -g ${ASSEMBLYFAI%.fai} -y tmp && RemoveDotAsStop.py -f tmp | sed 's/transcript://' > $SP.annotation.CAT.Prot.fasta && rm tmp
		LongestIsoformPerLocus.py $SP.annotation.CAT.Prot.fasta > $SP.annotation.CAT.LongestAA.fasta

		diamond blastp --ultra-sensitive --threads $THREADS --db $SWISSUNIPROT.dmnd --max-target-seqs 10 --query $SP.annotation.CAT.LongestAA.fasta --outfmt 6 \
			--out $SP.annotation.CAT.LongestAA.Diamond.outfmt6.out

		cut -f1 $SP.annotation.CAT.LongestAA.Diamond.outfmt6.out | sort -u > QueryHits.ids

		fasta_id_extractor.v2.0.py -f $SP.annotation.CAT.LongestAA.fasta -l QueryHits.ids --method=exclude > $SP.annotation.CAT.LongestAANoHits.fasta

		$CATHapp/cath-genomescan.pl -i $SP.annotation.CAT.LongestAANoHits.fasta -l $FUNFAM -o $SP.annotation.CAT.LongestAANoHits.FFresults

		echo -e "# `date`: Retrieve FunFam alignments and GO annotations using CATH API"
		cd $SP.annotation.CAT.LongestAANoHits.FFresults && mkdir -p GOterms

		grep "^$SP" $SP.annotation.CAT.LongestAANoHits.crh | awk '{ if ( $7 < 1e-5 ) print $0 }' | while read line; do
			Trans=`echo $line | awk '{ print $1 }'`
			OG=`echo $Trans | cut -f1 -d'.'`
			if [[ ! -s $OG.consensus.CATH.GOids.tsv ]]; then
				echo $OG
				FF=`echo $line | awk '{ print $2 }'`
				$CATHapp/retrieve_FunFam_aln_GO_EC_anno_CATH-API.pl $FF v4_2_0 . > $Trans.FF.log
				GO=`cat $Trans.FF.log | grep GO`
				info=`echo $line | sed 's/ /,/g'`

				GetGOterms.py -a $GO -f $info -g /user/work/tk19812/Databases/Flybase/go.obo_20230401
			fi
		done

		cat *.CATH.GOids.tsv > $SP.annotation.CAT.LongestAANoHits_CATH.GOids.tsv
		cd ..

		GenomicReagionAssociationArule.py -f $ASSEMBLYFAI -g $GFF3 -b $SP.GeneGenomicRegions.bed
	fi

	mkdir -p Permutations
	wc -l $aCNEES
	for i in $(seq 1 10000); do
		if [[ ! -s Permutations/ReshuffledCNEEs.$i.bed ]]; then
			bedtools shuffle -chrom -chromFirst -noOverlapping -i $aCNEES -g $ASSEMBLYFAI | bedtools sort -i - > Permutations/ReshuffledCNEEs.$i.bed
			wc -l Permutations/ReshuffledCNEEs.$i.bed
		fi
	done

	#mkdir -p AllCNEEsReshulled
	#NALLCNEES=`wc -l $ALLCNEES | cut -f1 -d' '`
	#REP=0
	#for i in $(seq 1 10000); do
	#	REP=`bc <<< $REP+1`
	#	if [ ! -s AllCNEEsReshulled/ReshuffleedAllCNEEs.$REP.bed ]; then
	#		bedtools shuffle -chrom -chromFirst -noOverlapping -i $ALLCNEES -g $ASSEMBLYFAI | bedtools sort -i - > AllCNEEsReshulled/ReshuffleedAllCNEEs.$REP.bed
	#	fi
	#done
	#	REP=`bc <<< $REP+1`
	#	bedtools shuffle -chrom -chromFirst -noOverlapping -i $ALLCNEES -g $ASSEMBLYFAI | bedtools sort -i - > AllCNEEsReshulled/ReshuffleedAllCNEEs.$REP.bed &
	#	REP=`bc <<< $REP+1`
	#	bedtools shuffle -chrom -chromFirst -noOverlapping -i $ALLCNEES -g $ASSEMBLYFAI | bedtools sort -i - > AllCNEEsReshulled/ReshuffleedAllCNEEs.$REP.bed &
	#	REP=`bc <<< $REP+1`
	#	bedtools shuffle -chrom -chromFirst -noOverlapping -i $ALLCNEES -g $ASSEMBLYFAI | bedtools sort -i - > AllCNEEsReshulled/ReshuffleedAllCNEEs.$REP.bed &
	#	REP=`bc <<< $REP+1`
	#	bedtools shuffle -chrom -chromFirst -noOverlapping -i $ALLCNEES -g $ASSEMBLYFAI | bedtools sort -i - > AllCNEEsReshulled/ReshuffleedAllCNEEs.$REP.bed &
	#	REP=`bc <<< $REP+1`
	#	bedtools shuffle -chrom -chromFirst -noOverlapping -i $ALLCNEES -g $ASSEMBLYFAI | bedtools sort -i - > AllCNEEsReshulled/ReshuffleedAllCNEEs.$REP.bed &
	#	REP=`bc <<< $REP+1`
	#	bedtools shuffle -chrom -chromFirst -noOverlapping -i $ALLCNEES -g $ASSEMBLYFAI | bedtools sort -i - > AllCNEEsReshulled/ReshuffleedAllCNEEs.$REP.bed &
	#	REP=`bc <<< $REP+1`
	#	bedtools shuffle -chrom -chromFirst -noOverlapping -i $ALLCNEES -g $ASSEMBLYFAI | bedtools sort -i - > AllCNEEsReshulled/ReshuffleedAllCNEEs.$REP.bed &
	#	REP=`bc <<< $REP+1`
	#	bedtools shuffle -chrom -chromFirst -noOverlapping -i $ALLCNEES -g $ASSEMBLYFAI | bedtools sort -i - > AllCNEEsReshulled/ReshuffleedAllCNEEs.$REP.bed &
	#	REP=`bc <<< $REP+1`
	#	bedtools shuffle -chrom -chromFirst -noOverlapping -i $ALLCNEES -g $ASSEMBLYFAI | bedtools sort -i - > AllCNEEsReshulled/ReshuffleedAllCNEEs.$REP.bed &
	#	REP=`bc <<< $REP+1`
	#	bedtools shuffle -chrom -chromFirst -noOverlapping -i $ALLCNEES -g $ASSEMBLYFAI | bedtools sort -i - > AllCNEEsReshulled/ReshuffleedAllCNEEs.$REP.bed &
	#	wait
	#done

#	mkdir -p PermutationsRandomCNEEs
#	NCNEES=`wc -l $aCNEES`
#	echo -e "$ALLCNEES $NCNEES"
#	for i in $(seq 1 10000); do
#		if [[ ! -s PermutationsRandomCNEEs/RandomCNEEs.$i.bed ]]; then
#			ExtractXrandomeElementFromBED.py $ALLCNEES $NCNEES > PermutationsRandomCNEEs/RandomCNEEs.$i.bed
#			wc -l PermutationsRandomCNEEs/RandomCNEEs.$i.bed
#		fi
#	done

	########## GO Terms Enrichment #########################################################################################################################
	#mkdir -p GOtermEnrichment && cd GOtermEnrichment
	#CATH=$SP.annotation.CAT.LongestAANoHits.FFresults/$SP.annotation.CAT.LongestAANoHits_CATH.GOids.tsv
	#. $APPS2/goatools/venv/bin/activate

	#echo -e "Namespace\tGOid\tGOname\tNumbCNEEs\tGenRegFenrich\tGenReagPval\tGenomeFraction\tPermFenrich\tPermPval\tPermuProb\tREPermP\tREAvgObsDist\tNloci\tloci" > ../CNEE.GOEnrichment.tsv

	#NAMESPACE="biological_process"

	##for job in $(seq 1 1); do
	##	#sleep 10s
	##	PhyloAccGOtermEnrichment.sh $SP $aCNEES $ASSEMBLYFAI ../$CATH $NAMESPACE & 
	##done
	##wait

	#cat *.CNEE.GOEnrichment.tsv >> ../CNEE.GOEnrichment.tsv
	#cd ..
	##cat CNEE.GOEnrichment.tsv | sed 's/ /_/g' | awk '{ if ($6 < 0.05) print $0}'
	#. $APPS2/ETE3/bin/activate
	#GOtermFDRcorrection.py -t CNEE.GOEnrichment.tsv > CNEE.GOEnrichment.FDR.tsv
	########################################################################################################################################################

	########## Gene Enrichment #############################################################################################################################
	#mkdir -p GeneEnrichment
	#cd $WORKDIR/CNEEselectionAnalysis/$EnrichmentDIR
	#REP=0
	#for i in $(seq 1 1000); do
	#	for j in $(seq 1 10); do
	#		REP=`bc <<< $REP+1`
	#		echo -e "`date` Permutation $REP"
	#		NearestGeneOnPermutation.py -b PermutationsRandomCNEEs/RandomCNEEs.$REP.bed -g $GFF3 | bedtools sort -i - > GeneEnrichment/ReshuffledCNEEs.$REP.bed &
	#		REP=`bc <<< $REP+1`
	#		NearestGeneOnPermutation.py -b PermutationsRandomCNEEs/RandomCNEEs.$REP.bed -g $GFF3 | bedtools sort -i - > GeneEnrichment/ReshuffledCNEEs.$REP.bed &
	#		REP=`bc <<< $REP+1`
	#		NearestGeneOnPermutation.py -b PermutationsRandomCNEEs/RandomCNEEs.$REP.bed -g $GFF3 | bedtools sort -i - > GeneEnrichment/ReshuffledCNEEs.$REP.bed &
	#		REP=`bc <<< $REP+1`
	#		NearestGeneOnPermutation.py -b PermutationsRandomCNEEs/RandomCNEEs.$REP.bed -g $GFF3 | bedtools sort -i - > GeneEnrichment/ReshuffledCNEEs.$REP.bed &
	#		REP=`bc <<< $REP+1`
	#		NearestGeneOnPermutation.py -b PermutationsRandomCNEEs/RandomCNEEs.$REP.bed -g $GFF3 | bedtools sort -i - > GeneEnrichment/ReshuffledCNEEs.$REP.bed &
	#		REP=`bc <<< $REP+1`
	#		NearestGeneOnPermutation.py -b PermutationsRandomCNEEs/RandomCNEEs.$REP.bed -g $GFF3 | bedtools sort -i - > GeneEnrichment/ReshuffledCNEEs.$REP.bed &
	#		REP=`bc <<< $REP+1`
	#		NearestGeneOnPermutation.py -b PermutationsRandomCNEEs/RandomCNEEs.$REP.bed -g $GFF3 | bedtools sort -i - > GeneEnrichment/ReshuffledCNEEs.$REP.bed &
	#		REP=`bc <<< $REP+1`
	#		NearestGeneOnPermutation.py -b PermutationsRandomCNEEs/RandomCNEEs.$REP.bed -g $GFF3 | bedtools sort -i - > GeneEnrichment/ReshuffledCNEEs.$REP.bed &
	#		REP=`bc <<< $REP+1`
	#		NearestGeneOnPermutation.py -b PermutationsRandomCNEEs/RandomCNEEs.$REP.bed -g $GFF3 | bedtools sort -i - > GeneEnrichment/ReshuffledCNEEs.$REP.bed &
	#		REP=`bc <<< $REP+1`
	#		NearestGeneOnPermutation.py -b PermutationsRandomCNEEs/RandomCNEEs.$REP.bed -g $GFF3 | bedtools sort -i - > GeneEnrichment/ReshuffledCNEEs.$REP.bed &
	#		#NearestGeneOnPermutation.py -b Permutations/ReshuffledCNEEs.$REP.bed -g $GFF3 | bedtools sort -i - > GeneEnrichment/ReshuffledCNEEs.$REP.bed &
	#		wait
	#	done
	#done
	#
	#ComputeGeneEnrichment.py -b $aCNEES -d GeneEnrichment -n 10000 > $SP.AccCNEEs.M1.Heliconius.RandCNEEsGeneEnrichment.tsv 
	########################################################################################################################################################

	########## Spatial Enrichment ##########################################################################################################################
	mkdir -p SpatialEnrichment
	SlidingWindowsGenerator.py -f $ASSEMBLYFAI -w 100000 > $SP.100k.SlidWin.bed 
	bedtools intersect -c -a $SP.100k.SlidWin.bed -b $aCNEES > $SP.ObsSpatialEnrichment.bed
	
	REP=0
	for i in $(seq 1 1000); do
		REP=`bc <<< $REP+1`
		echo -e "`date` Permutation $REP"
		bedtools intersect -c -a $SP.100k.SlidWin.bed -b PermutationsRandomCNEEs/RandomCNEEs.$REP.bed > SpatialEnrichment/ReshuffledCNEEs.$REP.bed &
		REP=`bc <<< $REP+1`
		bedtools intersect -c -a $SP.100k.SlidWin.bed -b PermutationsRandomCNEEs/RandomCNEEs.$REP.bed > SpatialEnrichment/ReshuffledCNEEs.$REP.bed &
		REP=`bc <<< $REP+1`
		bedtools intersect -c -a $SP.100k.SlidWin.bed -b PermutationsRandomCNEEs/RandomCNEEs.$REP.bed > SpatialEnrichment/ReshuffledCNEEs.$REP.bed &
		REP=`bc <<< $REP+1`
		bedtools intersect -c -a $SP.100k.SlidWin.bed -b PermutationsRandomCNEEs/RandomCNEEs.$REP.bed > SpatialEnrichment/ReshuffledCNEEs.$REP.bed &
		REP=`bc <<< $REP+1`
		bedtools intersect -c -a $SP.100k.SlidWin.bed -b PermutationsRandomCNEEs/RandomCNEEs.$REP.bed > SpatialEnrichment/ReshuffledCNEEs.$REP.bed &
		REP=`bc <<< $REP+1`
		bedtools intersect -c -a $SP.100k.SlidWin.bed -b PermutationsRandomCNEEs/RandomCNEEs.$REP.bed > SpatialEnrichment/ReshuffledCNEEs.$REP.bed &
		REP=`bc <<< $REP+1`
		bedtools intersect -c -a $SP.100k.SlidWin.bed -b PermutationsRandomCNEEs/RandomCNEEs.$REP.bed > SpatialEnrichment/ReshuffledCNEEs.$REP.bed &
		REP=`bc <<< $REP+1`
		bedtools intersect -c -a $SP.100k.SlidWin.bed -b PermutationsRandomCNEEs/RandomCNEEs.$REP.bed > SpatialEnrichment/ReshuffledCNEEs.$REP.bed &
		REP=`bc <<< $REP+1`
		bedtools intersect -c -a $SP.100k.SlidWin.bed -b PermutationsRandomCNEEs/RandomCNEEs.$REP.bed > SpatialEnrichment/ReshuffledCNEEs.$REP.bed &
		REP=`bc <<< $REP+1`
		bedtools intersect -c -a $SP.100k.SlidWin.bed -b PermutationsRandomCNEEs/RandomCNEEs.$REP.bed > SpatialEnrichment/ReshuffledCNEEs.$REP.bed &
	done
	
	ComputeSpatialEnrichment.py -b $SP.ObsSpatialEnrichment.bed -d SpatialEnrichment -n 10000 | bedtools sort -i - > $SP.AccCNEEs.M1.Heliconius.RandSpatialEnrichment.bed
	awk '{ if ($6 < 0.05) print $0 }' $SP.AccCNEEs.M1.Heliconius.RandSpatialEnrichment.bed > $SP.AccCNEEs.M1.Heliconius.RandSpatialEnrichment.Sign.bed
	########################################################################################################################################################
exit
	########## ATAC peak Enrichment ########################################################################################################################
	mkdir -p ATACenrichment && cd ATACenrichment

	ATACs="Eisa.Head_50Overlap.all.common.sort.clean.m1.bed Eisa.Wings_50Overlap.all.common.sort.clean.m1.bed"
	ATACDIR="/user/work/tk19812/HeliconiniiProject/HeliconGenomeAlignmentAnnotation/CNEEselectionAnalysis/ATACpeaks"

	############################### ATAC enrichment on all CNEEs
#	for ATAC in $ATACs; do
#		echo $ATAC
#
#		touch $ATAC.Observed
#
#		REP=0
#		for i in $(seq 1 1000); do
#			echo -e "`date` Permutation $REP"
#			REP=`bc <<< $REP+1`
#			NALLCNEES=`wc -l ../AllCNEEsReshulled/ReshuffledAllCNEEs.$REP.bed | cut -f1 -d' '` && ATACOVL=`bedtools intersect -a ../AllCNEEsReshulled/ReshuffledAllCNEEs.$REP.bed -b $ATACDIR/$ATAC | wc -l` && bc -l <<< $ATACOVL/$NALLCNEES >> $ATAC.Observed &
#			REP=`bc <<< $REP+1`
#			NALLCNEES=`wc -l ../AllCNEEsReshulled/ReshuffledAllCNEEs.$REP.bed | cut -f1 -d' '` && ATACOVL=`bedtools intersect -a ../AllCNEEsReshulled/ReshuffledAllCNEEs.$REP.bed -b $ATACDIR/$ATAC | wc -l` && bc -l <<< $ATACOVL/$NALLCNEES >> $ATAC.Observed &
#			REP=`bc <<< $REP+1`
#			NALLCNEES=`wc -l ../AllCNEEsReshulled/ReshuffledAllCNEEs.$REP.bed | cut -f1 -d' '` && ATACOVL=`bedtools intersect -a ../AllCNEEsReshulled/ReshuffledAllCNEEs.$REP.bed -b $ATACDIR/$ATAC | wc -l` && bc -l <<< $ATACOVL/$NALLCNEES >> $ATAC.Observed &
#			REP=`bc <<< $REP+1`
#			NALLCNEES=`wc -l ../AllCNEEsReshulled/ReshuffledAllCNEEs.$REP.bed | cut -f1 -d' '` && ATACOVL=`bedtools intersect -a ../AllCNEEsReshulled/ReshuffledAllCNEEs.$REP.bed -b $ATACDIR/$ATAC | wc -l` && bc -l <<< $ATACOVL/$NALLCNEES >> $ATAC.Observed &
#			REP=`bc <<< $REP+1`
#			NALLCNEES=`wc -l ../AllCNEEsReshulled/ReshuffledAllCNEEs.$REP.bed | cut -f1 -d' '` && ATACOVL=`bedtools intersect -a ../AllCNEEsReshulled/ReshuffledAllCNEEs.$REP.bed -b $ATACDIR/$ATAC | wc -l` && bc -l <<< $ATACOVL/$NALLCNEES >> $ATAC.Observed &
#			REP=`bc <<< $REP+1`
#			NALLCNEES=`wc -l ../AllCNEEsReshulled/ReshuffledAllCNEEs.$REP.bed | cut -f1 -d' '` && ATACOVL=`bedtools intersect -a ../AllCNEEsReshulled/ReshuffledAllCNEEs.$REP.bed -b $ATACDIR/$ATAC | wc -l` && bc -l <<< $ATACOVL/$NALLCNEES >> $ATAC.Observed &
#			REP=`bc <<< $REP+1`
#			NALLCNEES=`wc -l ../AllCNEEsReshulled/ReshuffledAllCNEEs.$REP.bed | cut -f1 -d' '` && ATACOVL=`bedtools intersect -a ../AllCNEEsReshulled/ReshuffledAllCNEEs.$REP.bed -b $ATACDIR/$ATAC | wc -l` && bc -l <<< $ATACOVL/$NALLCNEES >> $ATAC.Observed &
#			REP=`bc <<< $REP+1`
#			NALLCNEES=`wc -l ../AllCNEEsReshulled/ReshuffledAllCNEEs.$REP.bed | cut -f1 -d' '` && ATACOVL=`bedtools intersect -a ../AllCNEEsReshulled/ReshuffledAllCNEEs.$REP.bed -b $ATACDIR/$ATAC | wc -l` && bc -l <<< $ATACOVL/$NALLCNEES >> $ATAC.Observed &
#			REP=`bc <<< $REP+1`
#			NALLCNEES=`wc -l ../AllCNEEsReshulled/ReshuffledAllCNEEs.$REP.bed | cut -f1 -d' '` && ATACOVL=`bedtools intersect -a ../AllCNEEsReshulled/ReshuffledAllCNEEs.$REP.bed -b $ATACDIR/$ATAC | wc -l` && bc -l <<< $ATACOVL/$NALLCNEES >> $ATAC.Observed &
#			REP=`bc <<< $REP+1`
#			NALLCNEES=`wc -l ../AllCNEEsReshulled/ReshuffledAllCNEEs.$REP.bed | cut -f1 -d' '` && ATACOVL=`bedtools intersect -a ../AllCNEEsReshulled/ReshuffledAllCNEEs.$REP.bed -b $ATACDIR/$ATAC | wc -l` && bc -l <<< $ATACOVL/$NALLCNEES >> $ATAC.Observed &
#			wait
#		done
#	done
#
	for ATAC in $ATACs; do
		NALLCNEES=`wc -l $ALLCNEES | cut -f1 -d' '`
		OBSATACOVL=`bedtools intersect -a $ALLCNEES -b $ATACDIR/$ATAC | wc -l`
		echo $NALLCNEES $OBSATACOVL
		bc -l <<< $OBSATACOVL/$NALLCNEES

		wc -l $ATACDIR/$ATAC
		awk '{ sum+=$3-$2-1 } END { print sum }' $ATACDIR/$ATAC
		wc -l $ATAC.Observed
		echo -e "ATACpeakEnrichment.py -p $ATAC.Observed -o $OBSATACOVL -t $NALLCNEES"
		ATACpeakEnrichment.py -p $ATAC.Observed -o $OBSATACOVL -t $NALLCNEES
	done
#	
	############################### ATAC enrichment on aCNEEs

#	for ATAC in $ATACs; do
#		echo $ATAC
#
#		touch $ATAC.aCNEEsObserved
#	
#		REP=0
#		for i in $(seq 1 1000); do
#			echo -e "`date` Permutation $REP"
#			REP=`bc <<< $REP+1`
#			NALLCNEES=`wc -l ../PermutationsRandomCNEEs/RandomCNEEs.$REP.bed | cut -f1 -d' '` && ATACOVL=`bedtools intersect -a ../PermutationsRandomCNEEs/RandomCNEEs.$REP.bed -b $ATACDIR/$ATAC | wc -l` && bc -l <<< $ATACOVL/$NALLCNEES >> $ATAC.aCNEEsObserved &
#			REP=`bc <<< $REP+1`
#			NALLCNEES=`wc -l ../PermutationsRandomCNEEs/RandomCNEEs.$REP.bed | cut -f1 -d' '` && ATACOVL=`bedtools intersect -a ../PermutationsRandomCNEEs/RandomCNEEs.$REP.bed -b $ATACDIR/$ATAC | wc -l` && bc -l <<< $ATACOVL/$NALLCNEES >> $ATAC.aCNEEsObserved &
#			REP=`bc <<< $REP+1`
#			NALLCNEES=`wc -l ../PermutationsRandomCNEEs/RandomCNEEs.$REP.bed | cut -f1 -d' '` && ATACOVL=`bedtools intersect -a ../PermutationsRandomCNEEs/RandomCNEEs.$REP.bed -b $ATACDIR/$ATAC | wc -l` && bc -l <<< $ATACOVL/$NALLCNEES >> $ATAC.aCNEEsObserved &
#			REP=`bc <<< $REP+1`
#			NALLCNEES=`wc -l ../PermutationsRandomCNEEs/RandomCNEEs.$REP.bed | cut -f1 -d' '` && ATACOVL=`bedtools intersect -a ../PermutationsRandomCNEEs/RandomCNEEs.$REP.bed -b $ATACDIR/$ATAC | wc -l` && bc -l <<< $ATACOVL/$NALLCNEES >> $ATAC.aCNEEsObserved &
#			REP=`bc <<< $REP+1`
#			NALLCNEES=`wc -l ../PermutationsRandomCNEEs/RandomCNEEs.$REP.bed | cut -f1 -d' '` && ATACOVL=`bedtools intersect -a ../PermutationsRandomCNEEs/RandomCNEEs.$REP.bed -b $ATACDIR/$ATAC | wc -l` && bc -l <<< $ATACOVL/$NALLCNEES >> $ATAC.aCNEEsObserved &
#			REP=`bc <<< $REP+1`
#			NALLCNEES=`wc -l ../PermutationsRandomCNEEs/RandomCNEEs.$REP.bed | cut -f1 -d' '` && ATACOVL=`bedtools intersect -a ../PermutationsRandomCNEEs/RandomCNEEs.$REP.bed -b $ATACDIR/$ATAC | wc -l` && bc -l <<< $ATACOVL/$NALLCNEES >> $ATAC.aCNEEsObserved &
#			REP=`bc <<< $REP+1`
#			NALLCNEES=`wc -l ../PermutationsRandomCNEEs/RandomCNEEs.$REP.bed | cut -f1 -d' '` && ATACOVL=`bedtools intersect -a ../PermutationsRandomCNEEs/RandomCNEEs.$REP.bed -b $ATACDIR/$ATAC | wc -l` && bc -l <<< $ATACOVL/$NALLCNEES >> $ATAC.aCNEEsObserved &
#			REP=`bc <<< $REP+1`
#			NALLCNEES=`wc -l ../PermutationsRandomCNEEs/RandomCNEEs.$REP.bed | cut -f1 -d' '` && ATACOVL=`bedtools intersect -a ../PermutationsRandomCNEEs/RandomCNEEs.$REP.bed -b $ATACDIR/$ATAC | wc -l` && bc -l <<< $ATACOVL/$NALLCNEES >> $ATAC.aCNEEsObserved &
#			REP=`bc <<< $REP+1`
#			NALLCNEES=`wc -l ../PermutationsRandomCNEEs/RandomCNEEs.$REP.bed | cut -f1 -d' '` && ATACOVL=`bedtools intersect -a ../PermutationsRandomCNEEs/RandomCNEEs.$REP.bed -b $ATACDIR/$ATAC | wc -l` && bc -l <<< $ATACOVL/$NALLCNEES >> $ATAC.aCNEEsObserved &
#			REP=`bc <<< $REP+1`
#			NALLCNEES=`wc -l ../PermutationsRandomCNEEs/RandomCNEEs.$REP.bed | cut -f1 -d' '` && ATACOVL=`bedtools intersect -a ../PermutationsRandomCNEEs/RandomCNEEs.$REP.bed -b $ATACDIR/$ATAC | wc -l` && bc -l <<< $ATACOVL/$NALLCNEES >> $ATAC.aCNEEsObserved &
#			wait
#		done
#	done

	for ATAC in $ATACs; do
		NCNEES=`wc -l $aCNEES | cut -f1 -d' '`
		OBSATACOVL=`bedtools intersect -a $aCNEES -b $ATACDIR/$ATAC | wc -l`
		echo $NCNEES $OBSATACOVL
		bc -l <<< $OBSATACOVL/$NCNEES

		wc -l $ATACDIR/$ATAC
		awk '{ sum+=$3-$2-1 } END { print sum }' $ATACDIR/$ATAC
		wc -l $ATAC.Observed
		echo -e "ATACpeakEnrichment.py -p $ATAC.Observed -o $OBSATACOVL -t $NCNEES"
		ATACpeakEnrichment.py -p $ATAC.Observed -o $OBSATACOVL -t $NCNEES
	done


	########################################################################################################################################################
fi

if [[ $RUNPHYLOPREGDOMAINS == 1 ]];then
	. $APPS2/cactus-bin-v2.5.1/cactus_env/bin/activate

	SP='Hmel'
    ASSEMBLYFAI=`ls -1 $WORK/HeliconiniiProject/CompleteAssemblies/AssembledGenomes/$SP.assembly.v*.fasta.fai`
    BEDGENES=`ls -1 $WORKDIR/$SP.v*.annotation.CAT.bed`
    GWORKDIRFF3=`ls -1 $WORKDIR/$SP.v*.annotation.CAT.gff3`
    #CNEES=$WORKDIR/CNEEselectionAnalysis/Eisa.PhastCons/$SP.AccCNEEs.M1.Heliconius.bed
    CNEES=$WORKDIR/CNEEselectionAnalysis/Hmel.PhastCons.AllCNEEs/Hmel.phyloacc.RenamedIntCNEEs.M1.StrictGregarious.bed
    #SPENRICH=$WORKDIR/CNEEselectionAnalysis/$SP.Enrichment.AllCNEEs/Eisa.AccCNEEs.M1.Heliconius.SpatialEnrichment.Sign.bed
    SPENRICH=$WORKDIR/CNEEselectionAnalysis/$SP.Enrichment.StrictCNEEs/$SP.AccCNEEs.M1.Heliconius.RandSpatialEnrichment.Sign.bed
    NEUTMOD=$WORKDIR/CNEEselectionAnalysis/$SP.PhastCons/neut_Eisa_corrected.mod
    CONSMODDIR=$WORKDIR/CNEEselectionAnalysis/$SP.PhastCons/MODs
	REGREGIONS=$WORKDIR/CNEEselectionAnalysis/$SP.Enrichment.StrictCNEEs/$SP.GeneGenomicRegions.bed

	mkdir -p $SP.GrerariousConvCopyNum && cd $SP.GrerariousConvCopyNum

	GregSpecies=Djun,Evib,Hhec,Hert,Hdem,Hric,Hleu,Hsar,Hant,Hele,Hcon,Hhew,Hsap,Haoe,Hheb,Hxan,Hdor,Hbur,Hwal
	NonGregSpecies=Avcr,Avpe,Avfl,Eisa,Elam,Eali,Elyb,Etal,Htel,Hcly,Hhor,Hher,Hpet,Herd,Heet,Hlat,Hhim,Hcha,Hper,Hhie,Hege,Hnat,Hbes,Hnum,Hism,Heth,Hpar,Helv,Hatt,Hhel,Hmel,Hcyd,Hpac,Hheu,Htim


#	if [[ ! -s Windows.10k.5ks.bed ]]; then
#		echo -e "`date`: Generate sliding windows -w 10000 -s 5000"
#		bedtools makewindows -g $ASSEMBLYFAI -w 10000 -s 5000 > Windows.10k.5ks.bed
#	fi
#
#	cut -f1 $ASSEMBLYFAI | while read chr; do
#		echo $chr
#		if [[ ! -s $chr.CopyNum.maf ]]; then
#			hal2maf --refGenome Hmel --noAncestors --refSequence $chr $HAL $chr.CopyNum.maf
#		fi
#
#		if [[ ! -s $chr.CopyNum.Ratio.10k.5ks.mean.bed ]]; then
#			MAF2copyNumberCompare.py -i $chr.CopyNum.maf -a $GregSpecies -b $NonGregSpecies > $chr.CopyNum.Ratio.BedGraph
#
#			bedGraphToBigWig $chr.CopyNum.Ratio.BedGraph $ASSEMBLYFAI $chr.CopyNum.Ratio.bw &
#		fi
#	done
#
#	cat *CopyNum.Ratio.BedGraph > AllCNVratio.BedGraph
#
#	bedtools map -a Windows.10k.5ks.bed -b AllCNVratio.BedGraph -c 4 -o mean > AllCNVratio.BedGraph.10k.5ks.mean.bed

#exit
	cd ..
	mkdir -p $SP.GrerariousConvRegDomain && cd $SP.GrerariousConvRegDomain

#	ExtractSubTreeNeutralModel.py $NEUTMOD Djun,Evib,Hhec,Hert,Hdem,Hric,Hleu,Hsar,Hant,Hele,Hcon,Hhew,Hsap,Haoe,Hheb,Hxan,Hdor,Hbur,Hwal > ${NEUTMOD%.mod}.Gregarious.mod
#	ExtractSubTreeNeutralModel.py $NEUTMOD Avcr,Avpe,Avfl,Eisa,Elam,Eali,Elyb,Etal,Htel,Hcly,Hhor,Hher,Hpet,Herd,Heet,Hlat,Hhim,Hcha,Hper,Hhie,Hege,Hnat,Hbes,Hnum,Hism,Heth,Hpar,Helv,Hatt,Hhel,Hmel,Hcyd,Hpac,Hheu,Htim > ${NEUTMOD%.mod}.NonGregarious.mod

	wc -l $REGREGIONS

	echo -e "Chr\tStart\tEnd\tGene\tNclades" > RegDomainNclades.bed
	echo -e "Locus\tGroup\tNumSites\tAccAverage\tConsAverage\tOverallAverage\tAccMedian\tConsMedian\tOverallMedian" > RegDomainPhyloPconacc.tsv

	cut -f1-4 $REGREGIONS | while read region; do
		WindFile=`echo -e $region | cut -f1-3 -d ' ' | sed 's/ /_/g'`
		CHR=`echo -e $region | cut -f1 -d ' '`
		GENE=`echo -e $region | cut -f4 -d ' '`
		mkdir -p $CHR && cd $CHR
		echo -e "$GENE $WindFile"
		echo -e "$region" | sed "s/ /\t/g" > $WindFile.bed
		#hal2maf --refGenome $SP --onlyOrthologs --noAncestors --noDupes --targetGenomes Djun,Evib,Hhec,Hert,Hdem,Hric,Hleu,Hsar,Hant,Hele,Hcon,Hhew,Hsap,Haoe,Hheb,Hxan,Hdor,Hbur,Hwal --refTargets $WindFile.bed $HAL $WindFile.Gregarious.maf &
		#hal2maf --refGenome $SP --onlyOrthologs --noAncestors --noDupes --targetGenomes Avcr,Avpe,Avfl,Eisa,Elam,Eali,Elyb,Etal,Htel,Hcly,Hhor,Hher,Hpet,Herd,Heet,Hlat,Hhim,Hcha,Hper,Hhie,Hege,Hnat,Hbes,Hnum,Hism,Heth,Hpar,Helv,Hatt,Hhel,Hmel,Hcyd,Hpac,Hheu,Htim --refTargets $WindFile.bed $HAL $WindFile.NonGregarious.maf & 
		#
		#hal2maf --refGenome Hmel --noAncestors --targetGenomes Djun,Evib,Hhec,Hert,Hdem,Hric,Hleu,Hsar,Hant,Hele,Hcon,Hhew,Hsap,Haoe,Hheb,Hxan,Hdor,Hbur,Hwal --refTargets $WindFile.bed $HAL $WindFile.CopyNumberGregarious.maf &
		#hal2maf --refGenome Hmel --noAncestors --targetGenomes Avcr,Avpe,Avfl,Eisa,Elam,Eali,Elyb,Etal,Htel,Hcly,Hhor,Hher,Hpet,Herd,Heet,Hlat,Hhim,Hcha,Hper,Hhie,Hege,Hnat,Hbes,Hnum,Hism,Heth,Hpar,Helv,Hatt,Hhel,Hmel,Hcyd,Hpac,Hheu,Htim --refTargets $WindFile.bed $HAL $WindFile.CopyNumberNonGregarious.maf
		#wait

		#phyloP --base-by-base --mode CONACC --wig-scores --method LRT $NEUTMOD --branch Djun,Evib,Hhec $WindFile.Gregarious.fasta | gzip - > $WindFile.Gregarious.phyloP.wig.gz &
		#phyloP --base-by-base --mode CONACC --wig-scores --method LRT $NEUTMOD --branch Avcr,Avpe,Avfl $WindFile.NonGregarious.fasta | gzip - > $WindFile.NonGregarious.phyloP.wig.gz
		#wait

		if [[ ! -s $WindFile.maf ]]; then
			hal2maf --refGenome $SP --onlyOrthologs --noAncestors --noDupes --refTargets $WindFile.bed $HAL $WindFile.maf
		fi

		if [[ -s $WindFile.maf ]]; then
			#GREGSPLIST=`SpeciesSelectionMAF.py $WindFile.maf Djun,Evib,Hhec,Hert,Hdem,Hric,Hleu,Hsar,Hant,Hele,Hcon,Hhew,Hsap,Haoe,Hheb,Hxan,Hdor,Hbur,Hwal`
			#NONGREGSPLIST=`SpeciesSelectionMAF.py $WindFile.maf Avcr,Avpe,Avfl,Eisa,Elam,Eali,Elyb,Etal,Htel,Hcly,Hhor,Hher,Hpet,Herd,Heet,Hlat,Hhim,Hcha,Hper,Hhie,Hege,Hnat,Hbes,Hnum,Hism,Heth,Hpar,Helv,Hatt,Hhel,Hmel,Hcyd,Hpac,Hheu,Htim`
			GREGSPLIST='Djun,Evib,Hhec,Hert,Hdem,Hric,Hleu,Hsar,Hant,Hele,Hcon,Hhew,Hsap,Haoe,Hheb,Hxan,Hdor,Hbur,Hwal'
			NONGREGSPLIST='Avfl,Eisa,Herd,Heet,Hlat,Hhim,Hcha,Hper,Hhie,Hege,Hnat,Hbes,Hnum,Hism,Heth,Hpar,Helv,Hcyd,Hpac'

#			if [[ ! -s $WindFile.NonGregarious.phyloP.wig.gz ]]; then
				phyloP --base-by-base --mode CONACC --wig-scores --method LRT $NEUTMOD --branch $GREGSPLIST $WindFile.maf | gzip - > $WindFile.Gregarious.phyloP.wig.gz &
				phyloP --base-by-base --mode CONACC --wig-scores --method LRT $NEUTMOD --branch $NONGREGSPLIST $WindFile.maf | gzip - > $WindFile.NonGregarious.phyloP.wig.gz
				wait
#			fi

			MAF2Nclades.py -g $GENE -i $WindFile.maf -s Djun:Evib:Hhec:Hert,Hdem,Hric,Hleu,Hsar,Hant,Hele,Hcon,Hhew,Hsap:Haoe:Hheb,Hxan,Hdor:Hbur,Hwal >> ../RegDomainNclades.bed
			ParsePhyloPwig.py $WindFile.Gregarious.phyloP.wig.gz Gregarious >> ../RegDomainPhyloPconacc.tsv
			ParsePhyloPwig.py $WindFile.NonGregarious.phyloP.wig.gz NonGregarious >> ../RegDomainPhyloPconacc.tsv
		else
			ls -l $WindFile.maf
		fi

		cd ..
	done


	TransposeRegDomainPhyloPconaccTbl.py -t RegDomainPhyloPconacc.tsv -b RegDomainNclades.bed -o RegDomainPhyloPconaccOut.tsv
fi



if [[ $DATAFORPLOTS == 1 ]];then

	. $APPS2/cactus-bin-v2.5.1/cactus_env/bin/activate

	SP='Hmel'
	ASSEMBLYFAI=`ls -1 $WORK/HeliconiniiProject/CompleteAssemblies/AssembledGenomes/$SP.assembly.v*.fasta.fai`
	BEDGENES=`ls -1 $WORKDIR/$SP.v*.annotation.CAT.bed`
	GWORKDIRFF3=`ls -1 $WORKDIR/$SP.v*.annotation.CAT.gff3`
	#CNEES=$WORKDIR/CNEEselectionAnalysis/Eisa.PhastCons/$SP.AccCNEEs.M1.Heliconius.bed
	CNEES=$WORKDIR/CNEEselectionAnalysis/Hmel.PhastCons.AllCNEEs/Hmel.phyloacc.RenamedIntCNEEs.M1.StrictGregarious.bed
	#SPENRICH=$WORKDIR/CNEEselectionAnalysis/$SP.Enrichment.AllCNEEs/Eisa.AccCNEEs.M1.Heliconius.SpatialEnrichment.Sign.bed
	SPENRICH=$WORKDIR/CNEEselectionAnalysis/$SP.Enrichment.StrictCNEEs/$SP.AccCNEEs.M1.Heliconius.RandSpatialEnrichment.Sign.bed
	NEUTMOD=$WORKDIR/CNEEselectionAnalysis/$SP.PhastCons/neut_Eisa_corrected.mod
	CONSMODDIR=$WORKDIR/CNEEselectionAnalysis/$SP.PhastCons/MODs

	mkdir -p $SP.DataForPlots && cd $SP.DataForPlots

	cut -f1-3 $SPENRICH | while read region; do
		WindFile=`echo -e $region | cut -f1-3 -d ' ' | sed 's/ /_/g'`
		CHR=`echo -e $region | cut -f1 -d ' '`
		echo $WindFile
		echo -e "$region" | sed "s/ /\t/g" > $WindFile.bed
		hal2maf --refGenome $SP --onlyOrthologs --noAncestors --noDupes --refTargets $WindFile.bed $HAL $WindFile.63way.maf &
		hal2maf --refGenome $SP --onlyOrthologs --noAncestors --noDupes --rootGenome Heliconius --refTargets $WindFile.bed $HAL $WindFile.Heliconius.maf &
		hal2maf --refGenome $SP --onlyOrthologs --noAncestors --noDupes --targetGenomes Dple,Bany,Mcin,Jcoe,Smor,Pdid,Dpha,Diul,Ptel,Djun,Avcr,Avpe,Avfl,Elam,Evib,Eali,Elyb,Etal --refTargets $WindFile.bed $HAL $WindFile.Non-Heliconius.maf & 
		halAlignmentDepth --refSequence $CHR --noAncestors $HAL $SP > $SP.63wayCov.wig & 
	done
	wait

	cut -f1-3 $SPENRICH | while read region; do
		WindFile=`echo -e $region | cut -f1-3 -d ' ' | sed 's/ /_/g'`
		CHR=`echo -e $region | cut -f1 -d ' '`
		echo $WindFile
		phyloP --mode CONACC --wig-scores --method LRT $NEUTMOD $WindFile.63way.maf > $WindFile.63way.CONACCphyloP.wig 
		phyloP --mode CONACC --wig-scores --method LRT $NEUTMOD $WindFile.Heliconius.maf > $WindFile.Heliconius.CONACCphyloP.wig 
		phyloP --mode CONACC --wig-scores --method LRT $NEUTMOD $WindFile.Non-Heliconius.maf > $WindFile.Non-Heliconius.CONACCphyloP.wig 

		phastCons --score $WindFile.63way.maf $CONSMODDIR/$CHR.cons.mod,$CONSMODDIR/$CHR.noncons.mod > $WindFile.63way.scores.wig 
		phastCons --score $WindFile.Heliconius.maf $CONSMODDIR/$CHR.cons.mod,$CONSMODDIR/$CHR.noncons.mod > $WindFile.Heliconius.scores.wig 
		phastCons --score $WindFile.Non-Heliconius.maf $CONSMODDIR/$CHR.cons.mod,$CONSMODDIR/$CHR.noncons.mod > $WindFile.Non-Heliconius.scores.wig 

		sed -i "s/(null)/$CHR/" $WindFile.63way.scores.wig
		sed -i "s/(null)/$CHR/" $WindFile.Heliconius.scores.wig
		sed -i "s/(null)/$CHR/" $WindFile.Non-Heliconius.scores.wig
		wait
	done

	cut -f1-3 $SPENRICH | while read region; do
		WindFile=`echo -e $region | cut -f1-3 -d ' ' | sed 's/ /_/g'`
		CHR=`echo -e $region | cut -f1 -d ' '`
		for wig in $(ls -1 *.wig); do
			wigToBigWig $wig $ASSEMBLYFAI ${wig%.wig}.bw
			gzip $wig &
		done
	done
fi

if [[ $EXTRACTUCE == 1 ]]; then

	SP='Hmel'
	ASSEMBLYFAI=`ls -1 $WORK/HeliconiniiProject/CompleteAssemblies/AssembledGenomes/$SP.assembly.v*.fasta.fai`
	BEDGENES=`ls -1 $WORKDIR/$SP.v*.annotation.CAT.bed`
	GFF3=`ls -1 $WORKDIR/$SP.v*.annotation.CAT.gff3`
	intCNEES=$WORKDIR/CNEEselectionAnalysis/$SP.PhastCons/$SP.phyloacc.RenamedIntCNEEs.bed
	intgCNEES=$WORKDIR/CNEEselectionAnalysis/$SP.PhastCons/$SP.phyloacc.RenamedIntgCNEEs.bed
	intFASTA=$WORKDIR/CNEEselectionAnalysis/$SP.PhastCons/AllFastasIntrons
	intgFASTA=$WORKDIR/CNEEselectionAnalysis/$SP.PhastCons/AllFastasIntergenic
	intNotInfCNEEs=$WORKDIR/CNEEselectionAnalysis/$SP.PhastCons/Int-No-informative-CNEE-loci.txt
	intgNotInfCNEEs=$WORKDIR/CNEEselectionAnalysis/$SP.PhastCons/Intg-No-informative-CNEE-loci.txt
	intPhyloACCDIR=$WORKDIR/CNEEselectionAnalysis/$SP.PhastCons/phyloacc-Intronic-Output/results/elem_lik.txt
	intgPhyloACCDIR=$WORKDIR/CNEEselectionAnalysis/$SP.PhastCons/phyloacc-Intergenic-Output/results/elem_lik.txt

	cd $SP.PhastCons
	
	mkdir -p NotInformativeCNEEs

#	touch UCE.txt

#	grep -v '^#' $intPhyloACCDIR | grep CNEE | cut -f2 | cut -f1 -d'.' | while read CNEE; do
#	#grep -w M0 $intPhyloACCDIR | grep -v '^#' | cut -f2 | cut -f1 -d'.' | while read CNEE; do
#		MOD=`grep -w $CNEE $intPhyloACCDIR | cut -f3`
#		region=`grep -w $CNEE $intCNEES | cut -f1-3 | sed "s/\t/_/g"`
#		FASTA=$intFASTA/$region.fasta
#		PID=`Pidentity.py $FASTA`
#
#		echo -e "int$CNEE\t$MOD\t$region\t$PID" >> UCE.txt
#		if (( $(echo "$PID > 0.80" | bc -l) )); then
#			PhyloACC.ReorderCNEEalignment.py -f $intFASTA/$region.fasta > NotInformativeCNEEs/int$CNEE.$region.fasta
#			ls -l NotInformativeCNEEs/int$CNEE.$region.fasta
#		fi
#	done

	ls -l $intgPhyloACCDIR
	grep -v '^#' $intgPhyloACCDIR | grep CNEE | cut -f2 | cut -f1 -d'.' | while read CNEE; do
	#grep -w M0 $intPhyloACCDIR | grep -v '^#' | cut -f2 | cut -f1 -d'.' | while read CNEE; do
		MOD=`grep -w $CNEE $intgPhyloACCDIR | cut -f3`
		region=`grep -w $CNEE $intgCNEES | cut -f1-3 | sed "s/\t/_/g"`
		FASTA=$intgFASTA/$region.fasta
		PID=`Pidentity.py $FASTA`

		echo -e "intg$CNEE\t$MOD\t$region\t$PID" >> UCE.txt
		if (( $(echo "$PID > 0.80" | bc -l) )); then
		PhyloACC.ReorderCNEEalignment.py -f $intFASTA/$region.fasta > NotInformativeCNEEs/intg$CNEE.$region.fasta
		ls -l NotInformativeCNEEs/intg$CNEE.$region.fasta
		fi
	done

fi








