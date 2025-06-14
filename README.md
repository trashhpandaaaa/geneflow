# 🧬 GeneFlow: Human Cell Expression Analysis Platform

A comprehensive, production-ready toolkit for analyzing gene expression data from human osteosarcoma cells (GSE1000 dataset). This platform provides a complete pipeline from raw microarray data to publication-ready results with integrated machine learning classification.

## ✨ Features

- **🔬 Comprehensive Analysis**: Process 22K+ genes with advanced statistical methods
- **🤖 Machine Learning**: Automated classification with multiple algorithms
- **📊 Interactive Dashboard**: Streamlit-powered web interface for real-time exploration
- **📈 Publication-Ready**: High-quality plots and formatted results
- **🎓 Educational**: Step-by-step analysis options for learning

## 🚀 Quick Start

### Installation

```bash
# Clone the repository
git clone <repository-url>
cd crops-for-life

# Install dependencies
pip install -r requirements.txt
```

### Choose Your Analysis Approach

```bash
# Option 1: Full automated pipeline
python human_cell_analysis_clean.py

# Option 2: Step-by-step educational analysis
python simple_analysis_clean.py

# Option 3: Interactive web dashboard
streamlit run gene_expression_dashboard.py
```

## 📊 What This Platform Does

| Feature | Description |
|---------|-------------|
| **Analyzes** | 10 osteosarcoma cell samples across 5 surface types (amino acid conjugated vs control) |
| **Processes** | 22K+ genes with log transformation, normalization, and quality control |
| **Identifies** | Differentially expressed genes with statistical significance testing |
| **Predicts** | Treatment groups using machine learning (Random Forest, SVM, etc.) |
| **Generates** | Publication-ready plots, data exports, and comprehensive analysis reports |

## 🗂️ Dataset Information

- **Study**: GSE1000 - Osteosarcoma TE85 cell responses to amino acid surfaces
- **Platform**: Affymetrix Human Genome U133A Array
- **Conditions**: 5 surface types × 2 time points (6h, 32h)
- **Sample Size**: 10 biological samples
- **Applications**: Biomaterial design, cancer research, surface engineering

## 📁 Project Structure

```
crops-for-life/
├── 📁 dataset/
│   └── GSE1000_series_matrix.txt      # Raw microarray data
├── 📁 results/                        # Generated outputs and visualizations
├── 🐍 gene_expression_analyzer_clean.py      # Core analysis engine
├── 🐍 human_cell_analysis_clean.py           # Full automated pipeline
├── 🐍 simple_analysis_clean.py               # Educational step-by-step version
├── 🐍 gene_expression_dashboard.py           # Interactive Streamlit dashboard
├── 📄 requirements.txt                       # Python dependencies
└── 📖 README.md                              # This file
```

## 🔬 Analysis Pipeline

Our comprehensive analysis follows these steps:

1. **📥 Data Loading**
   - Parse GEO series matrix format
   - Extract 22K+ gene expression values
   - Validate data integrity

2. **🔧 Preprocessing**
   - Log2 transformation for normalization
   - Filter low-expression genes
   - Handle missing values

3. **✅ Quality Control**
   - Principal Component Analysis (PCA)
   - Sample correlation analysis
   - Outlier detection and visualization

4. **🧮 Differential Expression**
   - Statistical hypothesis testing
   - Multiple testing correction (FDR)
   - Volcano plot generation

5. **🤖 Machine Learning**
   - Train multiple classifiers
   - Cross-validation assessment
   - Feature importance analysis

6. **📊 Results Generation**
   - Export processed data to CSV
   - Generate publication-quality plots
   - Create comprehensive summary reports

## 📈 Expected Results

| Metric | Expected Range |
|--------|----------------|
| **Differential Genes** | 100-500 significant genes (FDR < 0.05) |
| **ML Classification Accuracy** | 85-95% with cross-validation |
| **Generated Outputs** | 8+ publication-ready plots and data files |
| **Processing Time** | 2-5 minutes for complete analysis |

## 💻 Code Example

```python
from gene_expression_analyzer_clean import GeneExpressionAnalyzer

# Initialize analyzer with your data
analyzer = GeneExpressionAnalyzer('dataset/GSE1000_series_matrix.txt')

# Parse the data file
analyzer.parse_geo_file()

# Get sample metadata
metadata = analyzer.get_sample_metadata()

# Perform differential expression analysis
de_results = analyzer.differential_expression_analysis(
    group1_samples=['GSM15786', 'GSM15795'],
    group2_samples=['GSM15785', 'GSM15790']
)

# Run machine learning classification
ml_results = analyzer.machine_learning_classification(
    metadata, 
    'surface_type'
)

# Generate visualizations
analyzer.generate_quality_plots()
analyzer.create_volcano_plot(de_results)
```

## 🛠️ Requirements

### System Requirements
- **Python**: 3.8 or higher
- **Memory**: Minimum 4GB RAM recommended
- **Storage**: ~100MB for data and results

### Key Dependencies
- `pandas` - Data manipulation and analysis
- `scikit-learn` - Machine learning algorithms
- `matplotlib` & `seaborn` - Static visualizations
- `plotly` - Interactive plots
- `streamlit` - Web dashboard framework
- `scipy` & `numpy` - Scientific computing

```bash
pip install -r requirements.txt  # Installs all required packages
```

## 🎯 Applications

### Research Areas
- **🔬 Biomedical Research**: Drug discovery, cancer biology, biomarker identification
- **🧪 Biomaterial Development**: Surface engineering, cell-material interactions
- **🎓 Education**: Bioinformatics training, machine learning with biological data
- **🏭 Quality Control**: Manufacturing consistency, experimental validation

### Use Cases
- Compare gene expression between treatment conditions
- Identify biomarkers for disease states
- Validate experimental hypotheses
- Train students in computational biology

## 🚨 Troubleshooting

### Common Issues and Solutions

| Problem | Solution |
|---------|----------|
| **Import Error** | Ensure you're using: `from gene_expression_analyzer_clean import GeneExpressionAnalyzer` |
| **Memory Issues** | Increase filtering thresholds to reduce the number of genes processed |
| **Dashboard Port Conflict** | Use alternative port: `streamlit run gene_expression_dashboard.py --server.port=8502` |
| **Missing Dependencies** | Run `pip install -r requirements.txt` to install all packages |
| **Data File Not Found** | Verify the dataset file path: `dataset/GSE1000_series_matrix.txt` |

### Getting Help
1. Check the code examples in the analysis scripts
2. Review this troubleshooting section
3. Examine the generated log files in the `results/` directory
4. Create an issue with detailed error information

## 📊 Output Files

The analysis generates several types of output files:

### Data Files
- `expression_data_simple.csv` - Processed expression matrix
- `sample_metadata_simple.csv` - Sample information and groupings
- `differential_expression_simple.csv` - Statistical results

### Visualizations
- `pca_analysis.png` - Principal component analysis
- `sample_correlation.png` - Sample-to-sample correlations
- `expression_distribution.png` - Data quality metrics
- `quality_control_plots.png` - Comprehensive QC dashboard

### Reports
- `analysis_summary.txt` - Complete analysis summary
- `summary_simple.txt` - Key findings and statistics

## 🤝 Contributing

We welcome contributions! Please feel free to:
- Report bugs or issues
- Suggest new features
- Submit pull requests
- Improve documentation

## 📜 License

This project is available for research and educational purposes. Please cite appropriately if used in publications.

---

**🧬 GeneFlow** - Making professional gene expression analysis accessible to everyone 📊✨
