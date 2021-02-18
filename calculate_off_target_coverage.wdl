task bedtools_genomecov{
    File bam
    File bai
    String sampleid
    Int disk = 30
    Int memory = 8
    
    runtime{
        docker : "biocontainers/bedtools:v2.28.0_cv2"
        memory : "${memory} GB"
        cpu : "1"
        disks : "local-disk ${disk} HDD"
    }


    command <<<
    echo "logging resources"
        set -xeuo pipefail
        function runtimeInfo() {            
            echo [$(date)]
            echo \* CPU usage: $(top -bn 2 -d 0.01 | grep '^%Cpu' | tail -n 1 | awk '{print $2}')%
            echo \* Memory usage: $(free -m | grep Mem | awk '{ OFMT="%.0f"; print ($3/$2)*100; }')%
            echo \* Disk usage: $(df | grep cromwell_root | awk '{ print $5 }')
        }
        while true;
            do runtimeInfo;
            sleep 30;
        done &
        
        bedtools genomecov -ibam ${bam} -bg -max 4 > ${sampleid}.coverage.txt
    >>>

    output{
        File bedtools_cvg = "${sampleid}.coverage.txt"
    }
}

task exclusion_regions{
    File bedtools_file
    File exclusion_script
    String sampleid
    Int disk = 30
    Int memory = 4
    
    runtime{
        docker : "breardon/calc_mutational_burden:1.1.1"
        memory : "${memory} GB"
        cpu : "1"
        disks : "local-disk ${disk} HDD"
    }


    command <<<
    echo "logging resources"
        set -xeuo pipefail
        function runtimeInfo() {            
            echo [$(date)]
            echo \* CPU usage: $(top -bn 2 -d 0.01 | grep '^%Cpu' | tail -n 1 | awk '{print $2}')%
            echo \* Memory usage: $(free -m | grep Mem | awk '{ OFMT="%.0f"; print ($3/$2)*100; }')%
            echo \* Disk usage: $(df | grep cromwell_root | awk '{ print $5 }')
        }
        while true;
            do runtimeInfo;
            sleep 30;
        done &
        
        cp ${exclusion_script} create_exclusion_bed.py
        
        python create_exclusion_bed.py --bedtools_genomecov_file ${bedtools_file} --name ${sampleid}
    >>>

    output{
        File exclusion_bed = "${sampleid}.exclusion.bed"
    }
}

task BedToIntervalList {
	File bedFile
	String bedFileName = basename(bedFile, ".bed")
	File refFastaDict = "gs://getzlab-workflows-reference_files-oa/hg19/Homo_sapiens_assembly19.dict"
	Int memory = 18
	Int disk = 15
    Int command_mem = 16

	command <<<
		java "-Xmx${command_mem}G" -jar /usr/gitc/picard.jar BedToIntervalList \
		I=${bedFile} \
		O=${bedFileName}.interval_list \
		SD=${refFastaDict}
	>>>

	output {
		File intervalList = "${bedFileName}.interval_list"
	}

	runtime {
		docker: "broadinstitute/genomes-in-the-cloud:2.3.1-1500064817"
		memory: "${memory} GB"
		cpu: "1"
		disks: "local-disk ${disk} HDD"
	}
}

task depthOfCov {
	File bam
	File bai
	Int minBaseQuality
	Int minMappingQuality
	String sampleid
	File exclusionintervals
	File geneList = "gs://fc-f36b3dc8-85f7-4d7f-bc99-a4610229d66a/broadinstitute/GermlinePipeline/DoC/GermlinePipeline_DoC_geneTrack.canonical.sorted.cleaned.cosmic_synonyms.refSeq.txt"
	File refFasta = "gs://fc-f36b3dc8-85f7-4d7f-bc99-a4610229d66a/broadinstitute/reference/hg19/fasta/Homo_sapiens_assembly19.fasta"
	File refFastaDict = "gs://fc-f36b3dc8-85f7-4d7f-bc99-a4610229d66a/broadinstitute/reference/hg19/fasta/Homo_sapiens_assembly19.dict"
	File refFastaIndex = "gs://fc-f36b3dc8-85f7-4d7f-bc99-a4610229d66a/broadinstitute/reference/hg19/fasta/Homo_sapiens_assembly19.fasta.fai"
	Int disk = 65
	Int memory = 36
	Int preemptible = 3
    Int command_mem = 32

	command <<<
		ln -s ${bam} bamfile.bam
		ln -s ${bai} bamfile.bai

		java -Xmx${command_mem}g -jar /usr/GenomeAnalysisTK.jar \
		-R ${refFasta} \
		-T DepthOfCoverage \
		-o ${sampleid} \
		-omitBaseOutput \
		-pt sample \
		-geneList ${geneList} \
		-I bamfile.bam \
        -XL ${exclusionintervals}
		--minBaseQuality ${minBaseQuality} \
		--minMappingQuality ${minMappingQuality}
        
        cat "${sampleid}.sample_summary" | cut -f 3 | tail -n1  > sample_mean_coverage.txt
	>>>

	output {
		File sampleGeneSummary = "${sampleid}.sample_gene_summary"
		File sampleSummary = "${sampleid}.sample_summary"
        Float sampleMeanCoverage = read_float("sample_mean_coverage.txt")
		File sampleStatistics = "${sampleid}.sample_statistics"
		File sampleIntervalSummary = "${sampleid}.sample_interval_summary"
		File sampleIntervalStatistics = "${sampleid}.sample_interval_statistics"
		File sampleCumulativeCoverageProportions = "${sampleid}.sample_cumulative_coverage_proportions"
		File sampleCumulativeCoverageCounts = "${sampleid}.sample_cumulative_coverage_counts"
	}

	runtime {
		docker: "broadinstitute/gatk3:3.7-0"
		memory: "${memory} GB"
		cpu: "1"
		disks: "local-disk ${disk} HDD"
		preemptible: "${preemptible}"
	}
}

workflow calculate_off_target_coverage_workflow{
	String sampleid
    File bam
    File bai

    call bedtools_genomecov{
    	input:
        	sampleid = sampleid,
            bam = bam,
            bai = bai

    }
    
    call exclusion_regions {
    	input:
        	sampleid = sampleid,
        	bedtools_file = bedtools_genomecov.bedtools_cvg
    }
    
    call BedToIntervalList {
    	input:
        	bedFile = exclusion_regions.exclusion_bed
    }
    
    call depthOfCov {
    	input:
        	bam = bam,
            bai = bai,
            exclusionintervals = BedToIntervalList.intervalList,
            sampleid = sampleid
    }
}