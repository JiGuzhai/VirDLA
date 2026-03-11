# Virdla: Promoter and Enhancer Analysis Tool

## Method
Virdla(Virtual Dual-Luciferase Assay) is an AlphaGenome-powered in silico tool that predicts dual-luciferase assay activity ratios from DNA sequences input, enabling high-throughput virtual screening of regulatory elements.
The method integrates RNA-seq and ATAC-seq predictions from the AlphaGenome DNA model to quantify insulator-mediated transcriptional regulation.

### Core Algorithm

The tool analyzes DNA sequences containing two promoters (p1 and p2), two coding sequences (CDS1 and CDS2), and an insulator element positioned between them. The analysis proceeds through three main steps:

**1. Sequence Analysis:**
- Constructs full sequence: p1 + CDS1 + insulator + p2 + CDS2
- Pads sequence to 16KB with N nucleotides
- Predicts RNA-seq and ATAC-seq profiles using AlphaGenome model
- Calculates area under the curve (AUC) for each genomic region

**2. Empty Sequence Baseline:**
- Constructs sequence without p1 promoter: CDS1 + insulator + p2 + CDS2
- Predicts RNA-seq profile for CDS1 region
- Establishes baseline expression level (empty_parameter)

**3. VirReport Score Calculation:**
```
VirReport = 0.5 × (CDS1_RNA - empty_parameter) / CDS2_RNA + 0.5 × (p1_ATAC / p2_ATAC)
```

The VirReport score integrates:
- Relative CDS1 expression after accounting for baseline (RNA component)
- Chromatin accessibility ratio between promoters (ATAC component)

### Data Processing

- **Noise Filtering:** Uses Median Absolute Deviation (MAD) method to filter noise from RNA-seq and ATAC data
- **Data Validation:** Removes infinite and NaN values before analysis
- **Coordinate Calculation:** Maps genomic regions to padded sequence coordinates

### Biological Interpretation

- Higher VirReport scores indicate stronger insulation effect
- The score reflects both transcriptional output and chromatin accessibility
- Normalization by CDS2 expression accounts for overall construct activity

## Usage

### Prerequisites

1. AlphaGenome API key
2. Python environment with required packages:
   - alphagenome
   - numpy
   - argparse
   - json

### Configuration

Create a `test_box.json` file with the following parameters:

```json
{
    "ot": "EFO:0001203",
    "t": "experiment_title",
    "p1": "promoter1_sequence",
    "CDS1": "coding_sequence_1",
    "insulator": "insulator_sequence",
    "p2": "promoter2_sequence",
    "CDS2": "coding_sequence_2"
}
```

**Parameters:**
- `ot`: Ontology term for cell type (e.g., "EFO:0001203" for MCF7 cells)
- `t`: Experiment title
- `p1`: First promoter sequence (5' promoter)
- `CDS1`: First coding sequence
- `insulator`: Insulator element sequence
- `p2`: Second promoter sequence (3' promoter)
- `CDS2`: Second coding sequence

### Running the Analysis

1. **Replace API Key:**
   Open `Virdla.py` and replace `"your_api_key"` with your actual AlphaGenome API key:
   ```python
   model = dna_client.create("your_actual_api_key_here")
   ```

2. **Execute the script:**
   ```bash
   python Virdla.py
   ```

3. **Output:**
   The script will print:
   - Configuration file loaded
   - API loaded status
   - RNA-seq and ATAC AUC values for each region
   - Empty parameter (baseline CDS1 expression)
   - Final VirReport score

### Command Line Arguments

You can also override configuration parameters via command line:

```bash
python Virdla.py -ot EFO:0001203 -p1 ATG... -CDS1 ATG... -insulator GATC... -p2 GAC... -CDS2 ATG...
```

### Interpreting Results

**Example Output:**
```
Config loaded from: test_box.json
API loaded
...
Virana: 0.6789
```

**Key Metrics:**
- **RNA AUC:** Transcriptional activity in each region
- **ATAC AUC:** Chromatin accessibility in each region
- **Empty Parameter:** Baseline CDS1 expression without p1 promoter
- **VirReport Score:** Final insulation effect metric,expressed as a percentage of the second promoter(normaly CMV promoter.)

### Validation

The tool has been validated against experimental data:
- CMV promoter prediction: 100.8% CMV (accurate baseline)
- hTERT promoter predictions correlate with literature values
- Insulator effects are detectable through VirReport score changes

## Supplement

1.PubMed entry: https://pubmed.ncbi.nlm.nih.gov/16331888/
The specific 6.5% ±1% value for hTERT promoter activity in A549 cells (compared to CMV),
but in the Virdla prediction,the hTERT promoter activity is 13.2%CMV in A549 cells.
2.In the aticle"Caspase-8 Gene Therapy Using the Human Telomerase Reverse Transcriptase Promoter for Malignant Glioma Cells. "Fig2A and 2B
the Prediction list below
True,Predict
{hTERT-181}{100.0}{MCF7}{EFO_0001203}  5.48%CMV 100
{hTERT-378}{115.0}{MCF7}{EFO_0001203}  8.53%CMV 157
{hTERT-1375}{155.0}{MCF7}{EFO_0001203} 5.83%CMV 106
3.when test the CMV promoter the prediction is 100.8%CMV,so we think the Virdla prediction is not necessary to crrect the value by the influence of insulator factors.
