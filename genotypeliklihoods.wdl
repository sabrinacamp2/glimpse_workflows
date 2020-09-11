workflow gtlikelihoods_workflow {
	Array[File] refpanel_sitestsv = ${refpanel_sitestsv}
	Array[File] refpanel_sitesvcf = ${refpanel_sitesvcf}
	Array[Int] chr = ['1', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '3', '4', '5', '6', '7', '8', '9', 'X']
	
}

call preparescatter{

}

scatter (idx in preparescatter.scatterIndices) {
    call gtlikelihoods {
        input:
            sitesvcf=refpanel_sitesvcf[idx],
            sitestsv=refpanel_sitestsv[idx],
            chr_names=chr[idx]
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
	File sitesvcf
	File sitestsv
	String chr_names
	File reffasta
	File refdict
	File refidx
	String sampleid
	Int diskSpaceGb
	Int memoryGb
	Int preemptible


	command <<<

		## compute genotype likelihoods
		bcftools mpileup -f ${reffasta} -I -E -a 'FORMAT/DP' -T ${sitesvcf} -r ${chr_names} ${samplebam} -Ou | bcftools call -Aim -C alleles -T ${sitestsv} -Oz -o ${sampleid}.${chr_name}.vcf.gz
		bcftools index -f ${sampleid}.${chr_name}.vcf.gz

	>>>

	output {
		Array[File] genotypelikelihood = "${sampleid}.${chr_name}.vcf.gz"
		Array[File] genotypelikelihood_index = "thousGP.${name}.sites.vcf.gz.csi"
	}

	runtime {
		docker: "vanallenlab/samtools_bcftools_htslib:1.9"
		memory: "${memoryGb} GB"
		disks: "local-disk ${diskSpaceGb} HDD"
		preemptible: "${preemptible}"
	}
}