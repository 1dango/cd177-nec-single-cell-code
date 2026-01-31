############################################################
## Cell Ranger preprocessing for mouse scRNA-seq
## Cell Ranger v5.1.0
###########################################################
cellranger mkref \
          --genome=mm10_ensembl98 \
            --fasta=GRCm38.primary_assembly.genome.fa \
              --genes=Mus_musculus.GRCm38.98.gtf

              # NEC sample
              cellranger count \
                        --id=NEC \
                          --transcriptome=mm10_ensembl98 \
                            --fastqs=/path/to/NEC_fastq/ \
                              --sample=NEC \
                                --expect-cells=8000 \
                                  --localcores=16 \
                                    --localmem=64

                                    # CTRL sample
                                    cellranger count \
                                            --id=CTRL \
                                                --transcriptome=mm10_ensembl98 \
                                                  --fastqs=/path/to/CTRL_fastq/ \
                                                    --sample=CTRL \
                                                      --expect-cells=8000 \
                                                        --localcores=16 \
                                                          --localmem=64

