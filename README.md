# PatientProfiler_BRCA

## PatientProfiler: A network-based approach to personalized medicine

Deciphering the intricate mechanisms underlying reprogramming in cancer cells is a crucial challenge in oncology, as it holds the key to advance our ability to diagnose and treat cancer patients. For this reason, comprehensive and patient-specific multi-omic characterization of tumor specimens has become increasingly common in clinical practice. While these efforts have advanced our understanding of the molecular mechanisms underlying breast cancer progression, the identification of personalized therapeutic approaches remains a distant goal. The main shortcoming is the absence of a robust computational framework to integrate and interpret the available multi-dimensional data and to drive translational solutions.

To fill this gap, we developed PatientProfiler, a computational pipeline that leverages causal interaction data, annotated in our in-house manually-curated resource, SIGNOR, to address how the genetic and molecular background of single patients contributes to the establishment of a malignant phenotype. PatientProfiler is an open-source, R-based package composed of several functions that allows multi-omic data analysis and standardization, generation of patient-specific mechanistic models of signal transduction, and extraction of network-based prognostic biomarkers.

We successfully benchmarked the tool to genomic, transcriptomic, (phospho)proteomic, and clinical data derived from 122 biopsies of treatment-naïve breast cancer, available at the CPTAC portal. We identified patient-specific mechanistic models (one patient, one network) that recapitulate dysregulated signaling pathways in breast cancer. This collection of models provides valuable insights into the underlying mechanisms of tumorigenesis and disease progression. Moreover, in-depth topological exploration of these networks has allowed us to define seven subgroups of Breast cancer patients, each associated with a unique transcriptomic signature and a distinct prognostic value.

In summary, our work demonstrates that PatientProfiler is a tool for patient-specific network analysis, advancing personalized medicine towards the identification of actionable biomarkers and tailored therapeutic strategies.

Read the full paper on [biorxiv](https://www.biorxiv.org/content/10.1101/2025.01.31.635886v1.full)

## Repository

The repository is divided in Steps:

1.  [Step1: Harmonization of input data](https://html-preview.github.io/?url=https://github.com/SaccoPerfettoLab/PatientProfiler_BRCA/blob/main/Step1/Step1.html)

2.  [Step2: Protein activity inference](https://html-preview.github.io/?url=https://github.com/SaccoPerfettoLab/PatientProfiler_BRCA/blob/main/Step2/Step2.html)

3.  [Step3: Generation of mechanistic models](https://html-preview.github.io/?url=https://github.com/SaccoPerfettoLab/PatientProfiler_BRCA/blob/main/Step3/Step3.html)

4.  [Step4: Network-based startification](https://html-preview.github.io/?url=https://github.com/SaccoPerfettoLab/PatientProfiler_BRCA/blob/main/Step4/Step4.html)

5.  [Step5: Identification of biomarkers](https://html-preview.github.io/?url=https://github.com/SaccoPerfettoLab/PatientProfiler_BRCA/blob/main/Step5/Step5.html)

Next, two additional folder contain the PatientProfiler robustness assessment and benchmarking:

-   **Robustness assessment** is divided in:

    -   [Randomization analysis](https://html-preview.github.io/?url=https://github.com/SaccoPerfettoLab/PatientProfiler_BRCA/blob/main/Robustness_analysis/randomization.html)

    -   [Cross-validation](https://html-preview.github.io/?url=https://github.com/SaccoPerfettoLab/PatientProfiler_BRCA/blob/main/Robustness_analysis/crossvalidation.html)

-   [Benchmark: comparison with other methods to find signatures](https://html-preview.github.io/?url=https://github.com/SaccoPerfettoLab/PatientProfiler_BRCA/blob/main/Benchmark/Benchmark.html)

<img src="./img/Figure1-2.svg" width="720" height="800"/>
