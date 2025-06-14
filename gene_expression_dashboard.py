"""
Interactive Gene Expression Analysis Dashboard
==============================================

A Streamlit dashboard for exploring and analyzing human cell gene expression data
from the GSE1000 osteosarcoma study.

Run with: streamlit run gene_expression_dashboard.py
"""

import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
import warnings

warnings.filterwarnings('ignore')

# Page configuration
st.set_page_config(
    page_title="Gene Expression Analysis Dashboard",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS for better styling
st.markdown("""
<style>
.metric-card {
    background-color: #f0f2f6;
    padding: 1rem;
    border-radius: 10px;
    border-left: 5px solid #1f77b4;
}
.success-card {
    background-color: #d4edda;
    padding: 1rem;
    border-radius: 10px;
    border-left: 5px solid #28a745;
}
</style>
""", unsafe_allow_html=True)

@st.cache_data
def load_data():
    """Load and preprocess the gene expression data."""
    from gene_expression_analyzer_clean import GeneExpressionAnalyzer
    
    try:
        analyzer = GeneExpressionAnalyzer('dataset/GSE1000_series_matrix.txt')
        analyzer.parse_geo_file()
        
        expression_log = np.log2(analyzer.expression_data + 1)
        min_expression_threshold = 2.0
        genes_above_threshold = (expression_log.mean(axis=1) > min_expression_threshold)
        expression_filtered = expression_log.loc[genes_above_threshold]
        
        sample_metadata = analyzer.get_sample_metadata()
        
        return analyzer, expression_filtered, sample_metadata
    except Exception as e:
        st.error(f"Error loading data: {str(e)}")
        return None, None, None

def plot_pca_interactive(expression_data, sample_metadata):
    """Create interactive PCA plot."""
    X_pca = expression_data.T
    pca = PCA(n_components=3)
    pca_result = pca.fit_transform(X_pca)
    
    pca_df = pd.DataFrame({
        'PC1': pca_result[:, 0],
        'PC2': pca_result[:, 1],
        'PC3': pca_result[:, 2],
        'Sample_ID': sample_metadata['sample_id'],
        'Surface_Type': sample_metadata['surface_type'],
        'Time_Point': sample_metadata['time_point'],
        'Treatment_Group': sample_metadata['treatment_group']
    })
    
    fig = px.scatter_3d(
        pca_df, 
        x='PC1', y='PC2', z='PC3',
        color='Surface_Type',
        symbol='Time_Point',
        hover_name='Sample_ID',
        hover_data=['Treatment_Group'],
        title=f'PCA Analysis (PC1: {pca.explained_variance_ratio_[0]:.1%}, PC2: {pca.explained_variance_ratio_[1]:.1%}, PC3: {pca.explained_variance_ratio_[2]:.1%})',
        labels={
            'PC1': f'PC1 ({pca.explained_variance_ratio_[0]:.1%})',
            'PC2': f'PC2 ({pca.explained_variance_ratio_[1]:.1%})',
            'PC3': f'PC3 ({pca.explained_variance_ratio_[2]:.1%})'
        }
    )
    
    fig.update_layout(
        scene=dict(
            xaxis_title=f'PC1 ({pca.explained_variance_ratio_[0]:.1%})',
            yaxis_title=f'PC2 ({pca.explained_variance_ratio_[1]:.1%})',
            zaxis_title=f'PC3 ({pca.explained_variance_ratio_[2]:.1%})'
        )
    )
    
    return fig, pca

def plot_gene_expression_heatmap(expression_data, sample_metadata, n_genes=50):
    """Create interactive heatmap of top variable genes."""
    top_genes = expression_data.var(axis=1).nlargest(n_genes).index
    heatmap_data = expression_data.loc[top_genes]
    
    sample_labels = [f"{row['sample_id']}<br>{row['surface_type']}<br>{row['time_point']}" 
                    for _, row in sample_metadata.iterrows()]
    
    fig = go.Figure(data=go.Heatmap(
        z=heatmap_data.values,
        x=sample_labels,
        y=heatmap_data.index,
        colorscale='RdBu_r',
        zmid=0,
        hoverongaps=False,
        colorbar=dict(title="Expression Level")
    ))
    
    fig.update_layout(
        title=f'Top {n_genes} Most Variable Genes',
        xaxis_title='Samples',
        yaxis_title='Genes',
        height=max(400, n_genes * 10)
    )
    
    return fig

def plot_volcano(diff_expr_results, p_threshold=0.05, fc_threshold=1.5):
    """Create interactive volcano plot."""
    if diff_expr_results is None or len(diff_expr_results) == 0:
        return None
    
    results = diff_expr_results.copy()
    results['-log10_p'] = -np.log10(results['p_value'])
    results['significant'] = (results['p_value_adjusted'] < p_threshold) & (np.abs(results['log2_fold_change']) > np.log2(fc_threshold))
    
    fig = px.scatter(
        results,
        x='log2_fold_change',
        y='-log10_p',
        color='significant',
        hover_name='gene_id',
        hover_data=['fold_change', 'p_value', 'p_value_adjusted'],
        title='Volcano Plot: Differential Gene Expression',
        labels={
            'log2_fold_change': 'Log2 Fold Change',
            '-log10_p': '-Log10 P-value'
        },
        color_discrete_map={True: 'red', False: 'gray'}
    )
    
    fig.add_hline(y=-np.log10(p_threshold), line_dash="dash", line_color="blue")
    fig.add_vline(x=np.log2(fc_threshold), line_dash="dash", line_color="blue")
    fig.add_vline(x=-np.log2(fc_threshold), line_dash="dash", line_color="blue")
    
    return fig

def main():
    """Main dashboard application."""
    
    st.title("🧬 Human Cell Gene Expression Analysis Dashboard")
    st.markdown("### GSE1000: Osteosarcoma Cell Response to Amino Acid Surfaces")
    
    with st.spinner("Loading and preprocessing data..."):
        analyzer, expression_data, sample_metadata = load_data()
    
    if analyzer is None:
        st.error("Failed to load data. Please check if the dataset file exists.")
        return
    
    st.sidebar.header("📊 Analysis Options")
    
    st.sidebar.subheader("Dataset Overview")
    st.sidebar.info(f"""
    **Genes**: {expression_data.shape[0]:,}  
    **Samples**: {expression_data.shape[1]}  
    **Platform**: Affymetrix HG-U133A  
    **Cell Line**: Human Osteosarcoma TE85
    """)
    
    analysis_type = st.sidebar.selectbox(
        "Select Analysis Type",
        ["Data Overview", "PCA Analysis", "Gene Expression Heatmap", "Differential Expression", "Machine Learning Classification"]
    )

    if analysis_type == "Data Overview":
        st.header("📋 Dataset Overview")
        
        col1, col2, col3, col4 = st.columns(4)
        
        with col1:
            st.metric("Total Genes", f"{expression_data.shape[0]:,}")
        with col2:
            st.metric("Total Samples", expression_data.shape[1])
        with col3:
            st.metric("Surface Types", sample_metadata['surface_type'].nunique())
        with col4:
            st.metric("Time Points", sample_metadata['time_point'].nunique())
        
        st.subheader("Sample Information")
        st.dataframe(sample_metadata, use_container_width=True)
        
        st.subheader("Expression Statistics")
        stats_df = pd.DataFrame({
            'Metric': ['Mean Expression', 'Median Expression', 'Std Expression', 'Min Expression', 'Max Expression'],
            'Value': [
                f"{expression_data.values.mean():.2f}",
                f"{np.median(expression_data.values):.2f}",
                f"{expression_data.values.std():.2f}",
                f"{expression_data.values.min():.2f}",
                f"{expression_data.values.max():.2f}"
            ]
        })
        st.dataframe(stats_df, use_container_width=True)
        
        st.subheader("Expression Distribution")
        fig_hist = px.histogram(
            x=expression_data.values.flatten(),
            nbins=50,
            title="Distribution of Gene Expression Values",
            labels={'x': 'Expression Level', 'y': 'Frequency'}
        )
        st.plotly_chart(fig_hist, use_container_width=True)
    
    elif analysis_type == "PCA Analysis":
        st.header("Principal Component Analysis")
        
        col1, col2 = st.columns([3, 1])
        
        with col2:
            st.subheader("PCA Settings")
            color_by = st.selectbox("Color samples by:", ["Surface_Type", "Time_Point", "Treatment_Group"])
        
        with col1:
            fig_pca, pca = plot_pca_interactive(expression_data, sample_metadata)
            st.plotly_chart(fig_pca, use_container_width=True)
        
        st.subheader("Explained Variance")
        variance_df = pd.DataFrame({
            'Component': [f'PC{i+1}' for i in range(min(5, len(pca.explained_variance_ratio_)))],
            'Explained Variance (%)': [f"{var*100:.1f}%" for var in pca.explained_variance_ratio_[:5]],
            'Cumulative (%)': [f"{np.sum(pca.explained_variance_ratio_[:i+1])*100:.1f}%" for i in range(min(5, len(pca.explained_variance_ratio_)))]
        })
        st.dataframe(variance_df, use_container_width=True)
    
    elif analysis_type == "Gene Expression Heatmap":
        st.header("Gene Expression Heatmap")
        
        col1, col2 = st.columns([3, 1])
        
        with col2:
            st.subheader("Heatmap Settings")
            n_genes = st.slider("Number of genes to display:", 10, 100, 50)
            gene_selection = st.selectbox("Gene selection criteria:", ["Most Variable", "Highest Expression"])
        
        with col1:
            fig_heatmap = plot_gene_expression_heatmap(expression_data, sample_metadata, n_genes)
            st.plotly_chart(fig_heatmap, use_container_width=True)
        
        if gene_selection == "Most Variable":
            top_genes = expression_data.var(axis=1).nlargest(n_genes)
        else:
            top_genes = expression_data.mean(axis=1).nlargest(n_genes)
        
        st.subheader(f"📋 Top {n_genes} {gene_selection} Genes")
        gene_df = pd.DataFrame({
            'Gene ID': top_genes.index,
            'Value': top_genes.values
        })
        st.dataframe(gene_df, use_container_width=True)
    
    elif analysis_type == "Differential Expression":
        st.header("Differential Expression Analysis")
        
        col1, col2 = st.columns(2)
        
        with col1:
            st.subheader("Group 1 (Control)")
            group1_surface = st.multiselect(
                "Select surface types for Group 1:",
                sample_metadata['surface_type'].unique(),
                default=['Control_Polystyrene']
            )
        
        with col2:
            st.subheader("Group 2 (Treatment)")
            available_surfaces = [s for s in sample_metadata['surface_type'].unique() if s not in group1_surface]
            group2_surface = st.multiselect(
                "Select surface types for Group 2:",
                available_surfaces,
                default=[s for s in available_surfaces if 'Acid' in s or 'Serine' in s or 'Amino' in s]
            )
        
        if st.button("🔬 Run Differential Expression Analysis"):
            if len(group1_surface) > 0 and len(group2_surface) > 0:
                with st.spinner("Analyzing differential expression..."):
                    group1_samples = sample_metadata[sample_metadata['surface_type'].isin(group1_surface)]['sample_id'].tolist()
                    group2_samples = sample_metadata[sample_metadata['surface_type'].isin(group2_surface)]['sample_id'].tolist()
                    
                    diff_expr_results = analyzer.differential_expression_analysis(
                        group1_samples=group1_samples,
                        group2_samples=group2_samples,
                        method='ttest'
                    )
                    
                    st.success("Analysis completed!")
                    
                    col1, col2, col3, col4 = st.columns(4)
                    
                    significant_genes = diff_expr_results[diff_expr_results['p_value_adjusted'] < 0.05]
                    upregulated = significant_genes[significant_genes['log2_fold_change'] > 0]
                    downregulated = significant_genes[significant_genes['log2_fold_change'] < 0]
                    
                    with col1:
                        st.metric("Total Genes", len(diff_expr_results))
                    with col2:
                        st.metric("Significant Genes", len(significant_genes))
                    with col3:
                        st.metric("Upregulated", len(upregulated))
                    with col4:
                        st.metric("Downregulated", len(downregulated))
                    
                    fig_volcano = plot_volcano(diff_expr_results)
                    if fig_volcano:
                        st.plotly_chart(fig_volcano, use_container_width=True)
                    
                    st.subheader("📊 Top Significant Genes")
                    display_results = significant_genes.head(20)[['gene_id', 'fold_change', 'log2_fold_change', 'p_value', 'p_value_adjusted']]
                    st.dataframe(display_results, use_container_width=True)
                    
                    csv = diff_expr_results.to_csv(index=False)
                    st.download_button(
                        label="💾 Download Full Results",
                        data=csv,
                        file_name="differential_expression_results.csv",
                        mime="text/csv"
                    )
            else:
                st.warning("Please select surface types for both groups.")
    
    elif analysis_type == "Machine Learning Classification":
        st.header("Machine Learning Classification")
        
        col1, col2 = st.columns([2, 1])
        
        with col2:
            st.subheader("ML Settings")
            target_column = st.selectbox("Classify by:", ["surface_type", "time_point", "treatment_group"])
            test_size = st.slider("Test set size:", 0.1, 0.5, 0.3)
        
        if st.button("Train Models"):
            with st.spinner("Training machine learning models..."):
                try:
                    X = expression_data.T
                    y = sample_metadata[target_column]
                    
                    from sklearn.feature_selection import VarianceThreshold
                    selector = VarianceThreshold(threshold=0.1)
                    X_selected = selector.fit_transform(X)
                    
                    X_train, X_test, y_train, y_test = train_test_split(
                        X_selected, y, test_size=test_size, random_state=42, stratify=y
                    )
                    
                    scaler = StandardScaler()
                    X_train_scaled = scaler.fit_transform(X_train)
                    X_test_scaled = scaler.transform(X_test)
                    
                    from sklearn.linear_model import LogisticRegression
                    from sklearn.svm import SVC
                    from sklearn.ensemble import GradientBoostingClassifier
                    
                    models = {
                        'Random Forest': RandomForestClassifier(n_estimators=100, random_state=42),
                        'Logistic Regression': LogisticRegression(random_state=42, max_iter=1000),
                        'SVM': SVC(random_state=42, kernel='rbf'),
                        'Gradient Boosting': GradientBoostingClassifier(random_state=42)
                    }
                    
                    results = {}
                    for name, model in models.items():
                        model.fit(X_train_scaled, y_train)
                        train_score = model.score(X_train_scaled, y_train)
                        test_score = model.score(X_test_scaled, y_test)
                        
                        results[name] = {
                            'Train Accuracy': train_score,
                            'Test Accuracy': test_score
                        }
                    
                    st.success("✅ Model training completed!")
                    
                    results_df = pd.DataFrame(results).T
                    st.dataframe(results_df, use_container_width=True)
                    
                    best_model = max(results.keys(), key=lambda k: results[k]['Test Accuracy'])
                    st.info(f"🏆 Best performing model: **{best_model}** (Test Accuracy: {results[best_model]['Test Accuracy']:.3f})")
                    
                    fig_performance = px.bar(
                        results_df.reset_index(),
                        x='index',
                        y=['Train Accuracy', 'Test Accuracy'],
                        barmode='group',
                        title='Model Performance Comparison',
                        labels={'index': 'Model', 'value': 'Accuracy'}
                    )
                    st.plotly_chart(fig_performance, use_container_width=True)
                    
                except Exception as e:
                    st.error(f"Error during model training: {str(e)}")
    
    st.markdown("---")
    st.markdown("### 📄 About This Dashboard")
    st.info("""
    This dashboard provides interactive analysis of the GSE1000 gene expression dataset, 
    which studies human osteosarcoma TE85 cells grown on different amino acid conjugated surfaces. 
    The analysis includes data exploration, dimensionality reduction, differential expression analysis, 
    and machine learning classification.
    """)

if __name__ == "__main__":
    main()
