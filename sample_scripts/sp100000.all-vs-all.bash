# computes MAP and precision-recall curves for sp100000 dataset at a range of codebook sizes

function GetRunId() {
	if [ -f runId.txt ] ; then runId=$(cat runId.txt); else runId=835800000; fi
	(( runId++ ))
	echo -n $runId >runId.txt
	# echo runId = ${runId}
	export runId
}

function TrecEval () {
	trec_eval_tc_compact ../${db}.homologs ~/tmp/${codebook}.run_${runId}.rankings ${codebook}.run_${runId}.trec false 101
	tail --lines=102 ${codebook}.run_${runId}.trec >${codebook}.run_${runId}.precision_recall_summary
	rm ~/tmp/${codebook}.run_${runId}.rankings
}
	
function EncodeDatabase () {
	echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
	if [ ! -f ${codebook}.protos ]; then
		echo "File ${codebook}.protos does not exist. Creating it."
		AAClusterFirst
	fi

	GetRunId

	echo >${codebook}.run_${runId}.bash \
	AAClustSigEncode	\
		--seqFile ../${db}.faa \
		--protoFile ${codebook}.protos \
		--wordLength ${wordLength} \
		--dist BlosumDistance \
		--matrixId 62 \
		--idindex 2 \
		--classindex 3 \
		--outFile ${dbSignatures} \
		--threshold ${threshold} \
		--assignNearest ${assignNearest} \
		--numthreads ${numThreads} 
	chmod +x ./${codebook}.run_${runId}.bash
		
	start=$(date +%s)
	echo "EncodeDatabase ${codebook}.run_${runId} start: ${start}"
	echo
	./${codebook}.run_${runId}.bash
	end=$(date +%s)
	echo
	echo "EncodeDatabase ${codebook}.run_${runId} end: ${end}"
	(( elapsed = end - start ))
	echo "EncodeDatabase ${codebook}.run_${runId} elapsed: ${elapsed}"
	echo ${elapsed} >${codebook}.run_${runId}.elapsed
	echo
}

function SigRankDatabase () {
	echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
	
	if [ ! -f ${dbSignatures} ] ; then
		echo "File ${dbSignatures} does not exist. Creating it."
		EncodeDatabase
	fi

	GetRunId
	
	echo >${codebook}.run_${runId}.bash \
	AAClustSig	\
		--dbSigs ${dbSignatures} \
		--querySigs ${dbSignatures} \
		--sigLength ${numClusters} \
		--outFile ~/tmp/${codebook}.run_${runId}.rankings \
		--numthreads ${numThreads} \
		--maxresults ${maxRecords} \
		--mode merge

	chmod +x ./${codebook}.run_${runId}.bash
		
	start=$(date +%s)
	echo "SigRankDatabase ${codebook}.run_${runId} start: ${start}"
		echo
		./${codebook}.run_${runId}.bash
		end=$(date +%s)
		echo
		echo "SigRankDatabase ${codebook}.run_${runId} end: ${end}"
		(( elapsed = end - start ))
		echo "SigRankDatabase ${codebook}.run_${runId} elapsed: ${elapsed}"
	echo ${elapsed} >${codebook}.run_${runId}.elapsed
		echo

	TrecEval
	# rm ~/tmp/${codebook}.run_${runId}.rankings
}

function AAClusterFirst () {
	echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
	if [ ! -f ${dbClusters} ] ; then
		echo "File ${dbClusters} does not exist. Creating it."
		GenerateClusters
	fi
	
	GetRunId

	echo >${codebook}.run_${runId}.bash \
	AAClusterFirst \
		--fastaFile ../${db}.faa \
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

	chmod +x ./${codebook}.run_${runId}.bash
		
	start=$(date +%s)
	echo "AAClusterFirst ${codebook}.run_${runId} start: ${start}"
		echo
		./${codebook}.run_${runId}.bash
		end=$(date +%s)
		echo
		echo "AAClusterFirst ${codebook}.run_${runId} end: ${end}"
		(( elapsed = end - start ))
		echo "AAClusterFirst ${codebook}.run_${runId} elapsed: ${elapsed}"
		echo ${elapsed} >${codebook}.run_${runId}.elapsed
}

function GenerateClusters () {
	echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
	GetRunId

	echo >${codebook}.run_${runId}.bash \
	AAClust \
		--fastaFile ../${db}.faa \
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
	chmod +x ./${codebook}.run_${runId}.bash
		
	start=$(date +%s)
	echo "GenerateClusters ${codebook}.run_${runId} start: ${start}"
	echo
	./${codebook}.run_${runId}.bash
	end=$(date +%s)
	echo
	echo "GenerateClusters ${codebook}.run_${runId} end: ${end}"
	(( elapsed = end - start ))
	echo "GenerateClusters ${codebook}.run_${runId} elapsed: ${elapsed}"
	echo ${elapsed} >${codebook}.run_${runId}.elapsed
}

function KmerRankDatabase () {
	GetRunId

	echo >${codebook}.run_${runId}.bash \
	KmerRank \
		--dbfile ../${db}.faa \
		--queryfile ../${db}.faa \
		--idIndex 2 \
		--classindex -1 \
		--fraglength 1 \
		--fragmode HausdorffAverageAverage \
		--kmermode BestOfBest \
		--kmerlength ${wordLength} \
		--alphabet AA \
		--matrix 62 \
		--dist BlosumDistance \
		--rankingfile ~/tmp/${codebook}.run_${runId}.rankings \
		--codebookfile ${codebook}.clusters \
		--prototypeFile ${codebook}.protos \
		--thresholdDistance ${threshold} \
		--defaultDistance ${defaultDistance} \
		--numthreads ${numThreads} \
		--maxrecords ${maxRecords} \
		--sampleSize 0 \
		--skip 1
	
	chmod +x ./${codebook}.run_${runId}.bash
	
	start=$(date +%s)
	echo "KmerRankDatabase ${codebook}.run_${runId} start: ${start}"
	echo
	./${codebook}.run_${runId}.bash
	end=$(date +%s)
	echo
	echo "KmerRankDatabase ${codebook}.run_${runId} end: ${end}"
	(( elapsed = end - start ))
	echo "KmerRankDatabase ${codebook}.run_${runId} elapsed: ${elapsed}"
	echo ${elapsed} >${codebook}.run_${runId}.elapsed
	echo

	TrecEval
	# rm ~/tmp/${codebook}.run_${runId}.rankings
}

function GetDomainCoverage () {
	echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
	GetRunId

	echo >${codebook}.run_${runId}.bash \
	GetDomainCoverage \
		--domainIn        ~/swissprot/swissprot.domains   \
		--protoIn         ${db}_t${threshold}_k${kmerLength}_c${codebookSize}.protos \
		--seqIn           ../${db}.faa             \
		--idIndex         2                                \
		--classIndex      -1                               \
		--kmerLength      30                               \
		--matrixId        62                               \
		--isCaseSensitive false                            \
		--outFile         ${db}_t${threshold}_k${kmerLength}_c${codebookSize}.domain_coverage \
		--notInDomain     ${db}_t${threshold}_k${kmerLength}_c${codebookSize}.not_in_domain
	
	chmod +x ./${codebook}.run_${runId}.bash
	start=$(date +%s)
	echo "GetDomainCoverage ${codebook}.run_${runId} start: ${start}"
	echo
	./${codebook}.run_${runId}.bash
	end=$(date +%s)
	echo
	echo "GetDomainCoverage ${codebook}.run_${runId} end: ${end}"
	(( elapsed = end - start ))
	echo "GetDomainCoverage ${codebook}.run_${runId} elapsed: ${elapsed}"
	echo ${elapsed} >${codebook}.run_${runId}.elapsed
	echo
}

#-------------------------------------------------------------------
numThreads=96

wordLength=30
threshold=305
defaultDistance=375

db=sp100000
maxRecords=4000

assignNearest=false

#-------------------------------------------------------------------

for numClusters in \
100000 \
39385 \
35148 \
30654 \
25961 \
21094 \
16075 \
10874 \
5503  \
1114  
do
	#---------------------------------------------------
	# Set up derived quantities here.
	#---------------------------------------------------
	dbProtos=../${db}_t${threshold}_k${wordLength}.protos
	dbClusters=../${db}_t${threshold}_k${wordLength}.clusters

	codebook=${db}_t${threshold}_k${wordLength}_c${numClusters}
	dbSignatures=${codebook}.signatures
	#---------------------------------------------------

	SigRankDatabase
done
