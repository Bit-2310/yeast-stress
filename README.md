# yeast-stress
Analysis of oxidative stress response in yeast using RNA-Seq data and automated workflows.
# Yeast Stress Analysis

## Overview
This project investigates the **oxidative stress response in yeast** using RNA-Seq data. The analysis focuses on understanding the differential gene expression and identifying key pathways and genes involved in the response to oxidative stress induced by hydrogen peroxide (**H₂O₂**).

The project employs an automated workflow using **Nextflow** and **Docker** for reproducibility and scalability.

---

## Project Workflow
### Key Steps:
1. **Data Acquisition**:
   - Dataset: PRJNA184040 from the Sequence Read Archive (SRA).
   - Tools used: SRA Toolkit for downloading raw FASTQ files.

2. **Quality Control**:
   - **Tool**: FastQC.
   - Outcome: Visualized read quality and detected potential issues.

3. **Read Trimming**:
   - **Tool**: Fastp.
   - Purpose: Removed low-quality bases and adapter sequences.

4. **Alignment**:
   - **Tool**: HISAT2.
   - Reference Genome: Saccharomyces cerevisiae genome.

5. **Gene Quantification**:
   - **Tool**: FeatureCounts.
   - Outcome: Generated raw count data for downstream analysis.

6. **Differential Gene Expression Analysis**:
   - **Tool**: edgeR.
   - Output: List of significantly differentially expressed genes (DEGs), volcano plots, and expression summaries.

7. **Visualization**:
   - **Tools**: R (ggplot2, pheatmap).
   - Visualizations:
     - Volcano plots for DEGs.
     - Heatmaps for gene expression patterns.
     - Principal Component Analysis (PCA).

---

## Results
### Key Findings:
- Identified oxidative stress-responsive genes involved in:
  - Reactive oxygen species (ROS) detoxification.
  - Cellular stress signaling pathways.
- Differential expression patterns indicate a robust transcriptional response to oxidative stress.

### Outputs:
- Gene count matrix: `results/counts/gene_counts.txt`
- Differential expression results: `results/edger/`
- Plots:
  - Volcano plots: `results/edger/volcano_plots/`
  - Heatmaps: `results/visualization/heatmaps/`

---

## Directory Structure
```plaintext
project/
│
├── data/                  # Raw and processed data
├── results/               # Outputs of analysis
│   ├── qc/                # Quality control results
│   ├── trimmed/           # Trimmed reads
│   ├── alignment/         # Aligned BAM files
│   ├── counts/            # Gene counts
│   ├── edger/             # Differential gene expression analysis
│   └── visualization/     # Plots and figures
├── src/                   # Source code (scripts and workflows)
├── docker-compose.yml     # Docker configuration
├── Dockerfile             # Docker container setup
├── nextflow.config        # Nextflow pipeline configuration
└── README.md              # Project description
```

---

## Tools and Technologies
- **Programming Languages**: Python, R
- **Workflow Management**: Nextflow
- **Containerization**: Docker
- **RNA-Seq Tools**:
  - FastQC
  - Fastp
  - HISAT2
  - FeatureCounts
  - edgeR
- **Visualization**:
  - ggplot2
  - pheatmap

---

## How to Run the Project
1. Clone the repository:
   ```bash
   git clone https://github.com/<your-username>/yeast_stress.git
   cd yeast_stress
   ```

2. Install dependencies:
   - Docker:
     ```bash
     sudo apt install docker
     ```
   - Nextflow:
     ```bash
     curl -s https://get.nextflow.io | bash
     ```

3. Run the pipeline:
   ```bash
   nextflow run main.nf
   ```

---

## References
- Dataset: [PRJNA184040](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA184040)
- Tools and documentation:
  - [Nextflow Documentation](https://www.nextflow.io/)
  - [Docker Documentation](https://docs.docker.com/)

---

## License
This project is licensed under the [MIT License](LICENSE).

---

## Acknowledgments
This project was inspired by the need to automate RNA-Seq analysis workflows for studying oxidative stress responses in yeast. Special thanks to the developers of Nextflow, Docker, and the RNA-Seq tools used.
