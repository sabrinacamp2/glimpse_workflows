workflow glimpse_concordance_workflow {

	Array[File] gnomad_af
	Array[File] gnomad_af_index
	Array[File] exome_imputed
	Array[File] exome_imputed_index
	Array[File]	genome_imputed
	Array[File] genome_imputed_index
	Array[String] chr = ['1', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '2', '20', '21', '22', '3', '4', '5', '6', '7', '8', '9', 'X']
	Int diskSpaceGb
	Int memoryGb
	Int preemptible

call preparescatter {

}

scatter (idx in preparescatter.scatterIndices) {
	call convertformat{
		input:
			gnomad_af=gnomad_af[idx],
			gnomad_af_index=gnomad_af_index[idx],
			diskSpaceGb=diskSpaceGb,
			memoryGb=memoryGb,
			preemptible=preemptible
	}
    call glimpse_concordance {
        input:
            gnomad_af=convertformat.gnomad_af_bcf,
            gnomad_af_index=convertformat.gnomad_af_bcf_index,
            chr_names=chr[idx],
            exome_imputed=exome_imputed[idx],
            exome_imputed_index=exome_imputed_index[idx],
            genome_imputed=genome_imputed[idx],
            genome_imputed_index=genome_imputed_index[idx],
            diskSpaceGb=diskSpaceGb,
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

task convertformat {
	File gnomad_af
	File gnomad_af_index
	File? exclude_samples
	String name = basename(gnomad_af, ".vcf.bgz")
	Int diskSpaceGb
	Int memoryGb
	Int preemptible


	command <<<

		ln -s ${gnomad_af} gnomad.vcf.gz
		ln -s ${gnomad_af_index} gnomad.vcf.gz.tbi

		bcftools view gnomad.vcf.gz -Ob > ${name}.bcf
		bcftools index -f ${name}.bcf

	>>>

	output {
		File gnomad_af_bcf = "${name}.bcf"
		File gnomad_af_bcf_index = "${name}.bcf.csi"
	}

	runtime {
		docker: "quay.io/biocontainers/bcftools:1.9--ha228f0b_3"
		memory: "${memoryGb} GB"
		disks: "local-disk ${diskSpaceGb} HDD"
		preemptible: "${preemptible}"
	}
}
task glimpse_concordance {
	
	File gnomad_af
	File gnomad_af_index
	File exome_imputed
	File exome_imputed_index
	File genome_imputed
	File genome_imputed_index
	String sampleid
	String chr_names
	Int diskSpaceGb
	Int memoryGb
	Int cpu
	Int preemptible


	command <<<
    
    	ln -s ${gnomad_af} gnomad.bcf
    	ln -s ${gnomad_af_index} gnomad.bcf.csi
    	ln -s ${exome_imputed} exome.bcf
    	ln -s ${exome_imputed_index} exome.bcf.csi
    	ln -s ${genome_imputed} genome.bcf
    	ln -s ${genome_imputed_index} genome.bcf.csi


		GLIMPSE_concordance --input ${chr_names} gnomad.bcf genome.bcf exome.bcf --minDP 8 --output /cromwell_root/${sampleid}.${chr_names} --minPROB 0.9999 \
		--bins 0.00000 0.00100 0.00200 0.00500 0.01000 0.05000 0.10000 0.20000 0.50000ex \
		--info_af AF_nfe --thread 4

		echo "names of output files"
		ls -lh /cromwell_root/
	>>>

	output {
		File plotting = "${sampleid}.${chr_names}.rsquare.grp.txt.gz"
	}

	runtime {
		docker: "vanallenlab/glimpse:1.0.1"
		memory: "${memoryGb} GB"
		cpu : "${cpu}"
		disks: "local-disk ${diskSpaceGb} HDD"
		preemptible: "${preemptible}"
	}
}