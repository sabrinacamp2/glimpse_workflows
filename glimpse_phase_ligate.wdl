workflow glimpse_phase_ligate_workflow {

	Array[File] refpanel_curated
    Array[File] refpanel_curated_index
	Array[File] genetic_maps
	Array[File] genotypelikelihoods
	Array[File] genotypelikelihoods_index
	Array[File] chunks
	Array[String] chr = ['1', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '2', '20', '21', '22', '3', '4', '5', '6', '7', '8', '9', 'X']
	String sampleid
	Int memoryGb
	Int preemptible



call preparescatter {

}

scatter (idx in preparescatter.scatterIndices) {
    call glimpse_phase {
        input:
        	refpanel_curated=refpanel_curated[idx],
            refpanel_curated_index=refpanel_curated_index[idx],
            geneticmap=genetic_maps[idx],
            genotypelikelihood=genotypelikelihoods[idx],
            genotypelikelihood_index=genotypelikelihoods_index[idx],
            chr_names=chr[idx],
            chunks = chunks[idx],
            sampleid=sampleid,
            memoryGb=memoryGb,
            preemptible=preemptible
        }
    call glimpse_ligate{
    	input:
    		imputed_chunks=glimpse_phase.imputed_chunks,
    		chr_names=chr[idx],
            sampleid=sampleid,
            memoryGb=memoryGb,
            preemptible=preemptible

    }
    }

   }

task preparescatter {

	Int numsamples_minusone
    
    command {
        # Create index for split intervals
		seq 0 ${numsamples_minusone} | cat > indices.dat
        }

	output {
		Array[Int] scatterIndices=read_lines("indices.dat")
	}
        runtime {
        docker: "vanallenlab/mutect:1.1.6"
        memory: "2 GB"
        preemptible: 3
    }
}

task glimpse_phase {
	
	File chunks
	File geneticmap
	File genotypelikelihood
	File genotypelikelihood_index
	String chr_names
	String sampleid
	String dollar = "$"
	File refpanel_curated
	File refpanel_curated_index
	Int memoryGb
	Int cpu = 8
	Int preemptible
    
    Int? diskspaceGB_buffer = 20
    Int? diskspaceGB = ceil(size(geneticmap, "G") + size(genotypelikelihood, "G") + size(refpanel_curated, "G")  
                             + diskspaceGB_buffer)


	command <<<
    
    	ln -s ${refpanel_curated} refpanel.bcf
        ln -s ${refpanel_curated_index} refpanel.bcf.csi
        ln -s ${genotypelikelihood} gt.vcf.gz
        ln -s ${genotypelikelihood_index} gt.vcf.gz.csi

        while IFS="" read -r LINE || [ -n "${dollar}LINE" ];
        do
        	printf -v ID "%02d" $(echo ${dollar}LINE | cut -d" " -f1)
        	IRG=$(echo $LINE | cut -d" " -f3)
        	ORG=$(echo $LINE | cut -d" " -f4)
        	OUT=${sampleid}.${chr_names}.${dollar}{ID}.bcf
        	GLIMPSE_phase --input gt.vcf.gz --reference refpanel.bcf --map ${geneticmap} --input-region ${dollar}{IRG} --output-region ${dollar}{ORG} --output ${dollar}{OUT} --thread 8
        	bcftools index -f ${dollar}{OUT}
        done < ${chunks}
	>>>

	output {
		Array[File] imputed_chunks = glob("${sampleid}.${chr_names}*.bcf*")
	}

	runtime {
		docker: "vanallenlab/glimpse:1.0.1"
		memory: "${memoryGb} GB"
		cpu: "${cpu}"
		disks: "local-disk ${diskspaceGB} HDD"
		preemptible: "${preemptible}"
	}
}

task glimpse_ligate{

	Array[File] imputed_chunks
	String sampleid
	String chr_names
	Int memoryGb
	Int preemptible
    
    Int? diskspaceGB_buffer = 20
    Int? diskspaceGB = ceil(size(imputed_chunks, "G") + diskspaceGB_buffer)

	command {

        mv ${sep = " " imputed_chunks} .
        ls *.bcf > list.${chr_names}.txt

		
        cat list.${chr_names}.txt
		

		GLIMPSE_ligate --input list.${chr_names}.txt --output ${sampleid}.${chr_names}.merged.bcf
		bcftools index -f ${sampleid}.${chr_names}.merged.bcf

	}

	output {
		File ligated = "${sampleid}.${chr_names}.merged.bcf"
		File ligated_index = "${sampleid}.${chr_names}.merged.bcf.csi"
	}

	runtime {
		docker: "vanallenlab/glimpse:1.0.1"
		memory: "${memoryGb} GB"
		disks: "local-disk ${diskspaceGB} HDD"
		preemptible: "${preemptible}"
	}

}