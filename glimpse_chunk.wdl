workflow glimpse_chunk_workflow {

	Array[File] refpanel_sitesvcf_array
	Array[File] refpanel_sitesvcf_index_array
	Array[String] chr = ['1', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '2', '20', '21', '22', '3', '4', '5', '6', '7', '8', '9', 'X']

call preparescatter {

}

scatter (idx in preparescatter.scatterIndices) {
    call glimpse_chunk {
        input:
            sitesvcf=refpanel_sitesvcf_array[idx],
            sitesvcfindex=refpanel_sitesvcf_index_array[idx],
            chr_names=chr[idx]
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

task glimpse_chunk {
	
	File sitesvcf
	File sitesvcfindex
	String chr_names
	Int diskSpaceGb
	Int memoryGb
	Int preemptible


	command <<<
    
    	ln -s ${sitesvcf} sites.vcf.gz
        ln -s ${sitesvcfindex} sites.vcf.gz.csi

		GLIMPSE_chunk --input sites.vcf.gz --region ${chr_names} --window-size 2000000 --buffer-size 200000 --output chunks.${chr_names}.txt

	>>>

	output {
		File chunks = "chunks.${chr_names}.txt"
	}

	runtime {
		docker: "vanallenlab/glimpse:1.0.1"
		memory: "${memoryGb} GB"
		disks: "local-disk ${diskSpaceGb} HDD"
		preemptible: "${preemptible}"
	}
}