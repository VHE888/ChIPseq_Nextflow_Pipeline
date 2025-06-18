# ChIP-seq Nextflow Pipeline

This is a Nextflow-based ChIP-seq analysis pipeline. The pipeline includes quality control, alignment, peak calling, and differential binding analysis with downstream visualization.

## 📁 Structure

```
Chipseq-pipeline/
├── main.nf                    # Main Nextflow workflow
├── nextflow.config            # Resource and parameter settings
├── test.nf                    # Optional test flow
├── modules/                   # Nextflow modules
├── refs/                      # Genome reference files
├── notebook/                  # Jupyter notebooks for downstream analysis
├── samplesheets/              # Sample metadata (full, subset, replicate)
├── results/
│   ├── plots/
│   └── tables/                # Results tables 
├── qc_reports/                # MultiQC and other HTML reports
├── .gitignore
└── README.md
```


## 🔁 Workflow Steps
1. **Quality Control**  
   - FastQC + MultiQC  
   - Trimming

2. **Alignment**  
   - BWA / Bowtie2 alignment to reference genome  
   - Deduplication and filtering

3. **Peak Calling**  
   - MACS3 for peak detection  
   - Normalization for input controls

4. **Differential Analysis**  
   - DiffBind or custom R-based comparison  
   - Outputs significant peaks and annotations

5. **Visualization**  
   - Peak overlap plots  
   - Enrichment figures  

## 🚀 Run
```bash
nextflow run main.nf -profile standard
```

## 👤 Author

**Wenshou He**  
GitHub: [@VHE888](https://github.com/VHE888)  
Boston University - Bioinformatics MSc
