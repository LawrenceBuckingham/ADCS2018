# Computes the information gain for each cluster in a SP500000 codebook.

function GetRunId() {
	if [ -f runId.txt ] ; then runId=$(cat runId.txt); else runId=835800000; fi
	(( runId++ ))
	echo -n $runId >runId.txt
	echo runId = ${runId}
	export runId
}

function TrecEval () {
	trec_eval_tc_compact ${homologs} ${codebook}_run_${runId}.rankings ${codebook}.trec false 101
	head --lines=1 ${codebook}.trec  >${codebook}.precision_recall_summary
	tail --lines=1 ${codebook}.trec >>${codebook}.precision_recall_summary
}

function DomainKMedoids () {
	echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
	if [ ! -f ${dbFasta} ]
	then
		echo "Training set ${dbFasta} does not exist. Creating it now."
		SplitFastaHomologs
	fi


	GetRunId
	task=DomainKMedoids
	script=${repeat}_${task}_run_${runId}.bash

	echo >${script} \
	DomainKMedoids \
		--domains ../swissprot/swissprot.domains \
		--db ${dbFasta} \
		--protos ${domainProtos} \
		--clusters ${domainClusters} \
		--kmerLength ${kmerLength} \
		--idIndex 2 \
		--classIndex 3 \
		--isCaseSensitive false \
		--threshold ${threshold} \
		--seed ${runId} \
		--numThreads ${numThreads} \
		--matrixid 62
	
	chmod +x ./${script}
	start=$(date +%s)
	echo "${task} ${repeat} run_${runId} start: ${start}"
	echo
	./${script}
	end=$(date +%s)
	echo
	echo "${task} ${repeat} run_${runId} end: ${end}"
	(( elapsed = end - start ))
	echo "${task} ${repeat} run_${runId} elapsed: ${elapsed}"
	echo
}

function GetLargestProtosByClass () {
	echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
	
	if [ ! -f ${domainProtos} ] ; then
		echo "File ${domainProtos} does not exist. Creating it now."
		DomainKMedoids
	fi

	GetRunId
	task=GetLargestProtosByClass
	script=${codebook}_run_${runId}.bash

	echo >${script} \
	GetLargestProtosByClass \
		--db ${dbFasta} \
		--protosIn ${domainProtos} \
		--clustersIn ${domainClusters} \
		--kmerLength ${kmerLength} \
		--idIndex 2 \
		--classIndex 3 \
		--protosPerClass ${protosPerClass} \
		--protosOut ${dbProtos} \
		--clustersOut ${dbClusters}
	chmod +x ./${script}
	
	chmod +x ./${script}
	start=$(date +%s)
	echo "${task} ${codebook}_run_${runId} start: ${start}"
	echo
	./${script}
	end=$(date +%s)
	echo
	echo "${task} ${codebook}_run_${runId} end: ${end}"
	(( elapsed = end - start ))
	echo "${task} ${codebook}_run_${runId} elapsed: ${elapsed}"
	echo
}
	
function EncodePartition () {
	echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

	if [ ! -f ${dbProtos} ]
	then
		echo "File ${dbProtos} does not exist. Creating it now."
		GetLargestProtosByClass
	fi

	task=EncodePartition
	script=${codebook}_${task}_run_${runId}.bash
	echo >${script} 
	
	for fileType in train test
	do
		echo AAClustSigEncode	\
			--seqFile ${partition}.${fileType}.faa \
			--protoFile ${dbProtos} \
			--wordLength ${wordLength} \
			--dist BlosumDistance \
			--matrixId 62 \
			--idindex 2 \
			--classindex 3 \
			--outFile ${codebook}.${fileType}.signatures \
			--threshold ${threshold} \
			--assignNearest ${assignNearest} \
			--numthreads ${numThreads} \
			>>${script}
	done
	chmod +x ./${script}
	
	start=$(date +%s)
	echo "EncodePartition ${codebook}_run_${runId} start: ${start}"
	echo
	./${script}
	end=$(date +%s)
	echo
	echo "EncodePartition ${codebook}_run_${runId} end: ${end}"
	(( elapsed = end - start ))
	echo "EncodePartition ${codebook}_run_${runId} elapsed: ${elapsed}"
	echo
}

function SigRankPartition () {
	echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
	
	if [ ! -f ${dbSignatures} ] ; then
		echo "File ${dbSignatures} does not exist. Creating it now."
		EncodePartition
	fi

	GetRunId

	numClusters=$(wc ${dbProtos} | awk '{print $1 / 2}')

	echo >${codebook}_run_${runId}.bash \
	AAClustSig	\
		--dbSigs ${dbSignatures} \
		--querySigs ${querySignatures} \
		--sigLength ${numClusters} \
		--outFile ${codebook}_run_${runId}.rankings \
		--numthreads ${numThreads} \
		--maxresults ${maxRecords} \
		--mode merge

	chmod +x ./${codebook}_run_${runId}.bash
		
	start=$(date +%s)
	echo "SigRankPartition ${codebook}_run_${runId} start: ${start}"
		echo
		./${codebook}_run_${runId}.bash
		end=$(date +%s)
		echo
		echo "SigRankPartition ${codebook}_run_${runId} end: ${end}"
		(( elapsed = end - start ))
		echo "SigRankPartition ${codebook}_run_${runId} elapsed: ${elapsed}"
		echo

	TrecEval
	rm ${codebook}_run_${runId}.rankings
}

function AAClusterFirst () {
	echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
	if [ ! -f ${dbClusters} ] ; then
		echo "File ${dbClusters} does not exist. Creating it now."
		GenerateClusters
	fi
	
	GetRunId

	echo >${codebook}_run_${runId}.bash \
	AAClusterFirst \
		--fastaFile ${db}.faa \
		--clusterIn ${dbClusters} \
		--protoIn ${dbProtos} \
		--clusterOut ${codebook}.clusters \
		--protoOut ${codebook}.protos \
		--idIndex 2 \
		--classIndex 3 \
		--numClusters ${numClusters} \
		--numthreads ${numThreads} \
		--wordLength ${wordLength} \
		--matrixId 62

	chmod +x ./${codebook}_run_${runId}.bash
		
	start=$(date +%s)
	echo "AAClusterFirst ${codebook}_run_${runId} start: ${start}"
		echo
		./${codebook}_run_${runId}.bash
		end=$(date +%s)
		echo
		echo "AAClusterFirst ${codebook}_run_${runId} end: ${end}"
		(( elapsed = end - start ))
		echo "AAClusterFirst ${codebook}_run_${runId} elapsed: ${elapsed}"
}

#-------------------------------------------------------------------

function GenerateClusters () {
	echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
	GetRunId

	echo >${codebook}_run_${runId}.bash \
	AAClust \
		--fastaFile ${db}.faa \
		--clusterOut ${dbClusters} \
		--protoOut ${dbProtos} \
		--idIndex 2 \
		--classIndex 3 \
		--numthreads ${numThreads} \
		--wordLength ${wordLength} \
		--matrixId 62 \
		--increment ${numThreads} \
		--seed $(date +%s) \
		--threshold ${threshold}
	chmod +x ./${codebook}_run_${runId}.bash
		
	start=$(date +%s)
	echo "GenerateClusters ${codebook}_run_${runId} start: ${start}"
	echo
	./${codebook}_run_${runId}.bash
	end=$(date +%s)
	echo
	echo "GenerateClusters ${codebook}_run_${runId} end: ${end}"
	(( elapsed = end - start ))
	echo "GenerateClusters ${codebook}_run_${runId} elapsed: ${elapsed}"
}

#-------------------------------------------------------------------

function SplitFastaHomologs () {
	echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
	GetRunId
	task=SplitFastaHomologs

	echo >${repeat}_${task}_run_${runId}.bash \
	SplitFastaHomologs \
		--fasta ${db}.faa \
		--homologs ${db}.homologs \
		--outstub ${repeat}_p \
		--idIndex 2 \
		--seed $(date +%s) \
		--parts 10
	chmod +x ./${repeat}_${task}_run_${runId}.bash
		
	start=$(date +%s)
	echo "${repeat}_${task}_run_${runId} start: ${start}"
	echo
	./${repeat}_${task}_run_${runId}.bash
	end=$(date +%s)
	echo
	echo "${repeat}_${task}_run_${runId} end: ${end}"
	(( elapsed = end - start ))
	echo "${repeat}_${task}_run_${runId} elapsed: ${elapsed}"
}

#-------------------------------------------------------------------
numThreads=96

wordLength=30
kmerLength=${wordLength}
threshold=305
defaultDistance=375

db=sp100000
maxRecords=4000

assignNearest=false

#-------------------------------------------------------------------

for (( r = 1; r <= 10; r ++ ))
do
	repeat=${db}_t${threshold}_k${wordLength}_r${r}

	if [ ! -f ${repeat}_p.01.train.faa ]
	then
		SplitFastaHomologs
	fi

for part in 01 02 03 04 05 06 07 08 09 10
do
	partition=${repeat}_p.${part}
	domainProtos=${partition}.domain.protos
	domainClusters=${partition}.domain.clusters
	homologs=${partition}.homologs
	dbFasta=${partition}.train.faa

	if [ ! -f ${domainProtos} ]
	then
		DomainKMedoids
	fi

for protosPerClass in 45 50 55 60
do
	#---------------------------------------------------
	# Set up derived quantities here.
	#---------------------------------------------------
		codebook=${partition}_P${protosPerClass}
		dbProtos=${codebook}.train.protos
		dbClusters=${codebook}.train.clusters
		dbSignatures=${codebook}.train.signatures
		querySignatures=${codebook}.test.signatures
	#--------------------------------------------------- 

	if [ ! -f ${codebook}.precision_recall_summary ]; then
		echo "Computing ${codebook}.precision_recall_summary"
		SigRankPartition
	fi
done
done
done
