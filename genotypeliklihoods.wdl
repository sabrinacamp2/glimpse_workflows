workflow gtlikelihoods_workflow {
	Array[File] refpanel_sitestsv
    Array[File] refpanel_sitestsv_index
	Array[File] refpanel_sitesvcf
    Array[File] refpanel_sitesvcf_index
	Array[String] chr = ['1', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '2', '20', '21', '22', '3', '4', '5', '6', '7', '8', '9', 'X']


call preparescatter {

}

scatter (idx in preparescatter.scatterIndices) {
    call gtlikelihoods {
        input:
            sitesvcf=refpanel_sitesvcf[idx],
            sitesvcfindex=refpanel_sitesvcf_index[idx],
            sitestsv=refpanel_sitestsv[idx],
            sitestsvindex=refpanel_sitestsv_index[idx],
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


task gtlikelihoods {
	
	File samplebam
    File samplebam_index
	File sitesvcf
    File sitesvcfindex
	File sitestsv
    File sitestsvindex
	String chr_names
	File reffasta
	File refdict
	File refidx
	String sampleid
	Int memoryGb
	Int preemptible
    
    Int? diskspaceGB_buffer = 20
    Int? diskspaceGB = ceil(size(samplebam, "G") + size(samplebam_index, "G") + size(sitesvcf, "G")  
                             + diskspaceGB_buffer)


	command <<<
    
    	ln -s ${samplebam} sample.bam
        ln -s ${samplebam_index} sample.bai
        ln -s ${sitesvcf} sites.vcf.gz
        ln -s ${sitesvcfindex} sites.vcf.gz.csi
        ln -s ${sitestsv} sites.tsv.gz
        ln -s ${sitestsvindex} sites.tsv.gz.tbi

		## compute genotype likelihoods
		bcftools mpileup -f ${reffasta} -A -I -E -a 'FORMAT/DP' -T sites.vcf.gz -r ${chr_names} sample.bam -Ou | bcftools call -Aim -C alleles -T sites.tsv.gz -Oz -o ${sampleid}.${chr_names}.vcf.gz
		bcftools index -f ${sampleid}.${chr_names}.vcf.gz

	>>>

	output {
		File genotypelikelihood = "${sampleid}.${chr_names}.vcf.gz"
		File genotypelikelihood_index = "${sampleid}.${chr_names}.vcf.gz.csi"
	}

	runtime {
		docker: "vanallenlab/samtools_bcftools_htslib:1.9"
		memory: "${memoryGb} GB"
		disks: "local-disk ${diskspaceGB} HDD"
		preemptible: "${preemptible}"
	}
}