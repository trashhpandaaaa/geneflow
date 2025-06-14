#!/usr/bin/env python3

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from gene_expression_analyzer_clean import GeneExpressionAnalyzer
import warnings

warnings.filterwarnings('ignore')

def load_and_preprocess_data():
    print("🔄 Loading data...")
    
    data_file = 'dataset/GSE1000_series_matrix.txt'
    if not os.path.exists(data_file):
        raise FileNotFoundError(f"Dataset file not found: {data_file}")
    
    analyzer = GeneExpressionAnalyzer(data_file)
    analyzer.parse_geo_file()
    
    sample_metadata = analyzer.get_sample_metadata()
    
    expression_log = np.log2(analyzer.expression_data + 1)
    
    print(f"✅ Loaded {expression_log.shape[0]} genes, {expression_log.shape[1]} samples")
    
    return analyzer, expression_log, sample_metadata

def explore_data(expression_data, sample_metadata):
    print("\n📊 Data exploration...")
    
    print(f"Expression data shape: {expression_data.shape}")
    print(f"Sample metadata shape: {sample_metadata.shape}")
    
    print("\nSample information:")
    print(sample_metadata.to_string(index=False))
    
    print(f"\nSurface types: {sample_metadata['surface_type'].unique()}")
    print(f"Time points: {sample_metadata['time_point'].unique()}")
    
    print(f"\nExpression statistics:")
    print(f"Mean: {expression_data.values.mean():.2f}")
    print(f"Std: {expression_data.values.std():.2f}")
    print(f"Min: {expression_data.values.min():.2f}")
    print(f"Max: {expression_data.values.max():.2f}")

def create_simple_plots(expression_data, sample_metadata):
    print("\n📈 Creating plots...")
    
    os.makedirs('results', exist_ok=True)
    
    plt.figure(figsize=(10, 6))
    plt.hist(expression_data.values.flatten(), bins=50, alpha=0.7)
    plt.title('Expression Value Distribution')
    plt.xlabel('Log2 Expression')
    plt.ylabel('Frequency')
    plt.savefig('results/expression_distribution.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    plt.figure(figsize=(8, 6))
    corr_matrix = expression_data.corr()
    sns.heatmap(corr_matrix, cmap='coolwarm', center=0, square=True)
    plt.title('Sample Correlation Matrix')
    plt.savefig('results/sample_correlation.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    print("✅ Basic plots created")

def simple_differential_expression(analyzer, expression_data, sample_metadata):
    print("\n🧬 Differential expression analysis...")
    
    try:
        control_samples = sample_metadata[sample_metadata['surface_type'] == 'Control_Polystyrene']['sample_id'].tolist()
        treatment_samples = sample_metadata[sample_metadata['surface_type'] != 'Control_Polystyrene']['sample_id'].tolist()
        
        print(f"Control samples: {control_samples}")
        print(f"Treatment samples: {treatment_samples}")
        
        if len(control_samples) == 0 or len(treatment_samples) == 0:
            print("⚠️ Cannot perform analysis - insufficient samples in groups")
            return None
        
        analyzer.expression_data = expression_data
        
        diff_results = analyzer.differential_expression_analysis(
            group1_samples=control_samples,
            group2_samples=treatment_samples,
            method='ttest'
        )
        
        significant = diff_results[diff_results['p_value'] < 0.05]
        print(f"✅ Found {len(significant)} significant genes (p < 0.05)")
        
        if len(significant) > 0:
            print("\nTop 5 significant genes:")
            top_genes = significant.nsmallest(5, 'p_value')
            for _, row in top_genes.iterrows():
                print(f"  {row['gene_id']}: FC={row['fold_change']:.2f}, p={row['p_value']:.2e}")
        
        return diff_results
        
    except Exception as e:
        print(f"⚠️ Differential expression failed: {e}")
        return None

def simple_machine_learning(expression_data, sample_metadata):
    print("\n🤖 Machine learning classification...")
    
    try:
        from sklearn.ensemble import RandomForestClassifier
        from sklearn.model_selection import train_test_split
        from sklearn.preprocessing import StandardScaler
        from sklearn.feature_selection import VarianceThreshold
        
        X = expression_data.T
        y = sample_metadata['surface_type']
        
        print(f"Features: {X.shape[1]}, Samples: {X.shape[0]}, Classes: {y.nunique()}")
        
        selector = VarianceThreshold(threshold=0.5)
        X_selected = selector.fit_transform(X)
        print(f"After feature selection: {X_selected.shape[1]} features")
        
        X_train, X_test, y_train, y_test = train_test_split(
            X_selected, y, test_size=0.3, random_state=42, stratify=y
        )
        
        scaler = StandardScaler()
        X_train_scaled = scaler.fit_transform(X_train)
        X_test_scaled = scaler.transform(X_test)
        
        model = RandomForestClassifier(n_estimators=100, random_state=42)
        model.fit(X_train_scaled, y_train)
        
        train_score = model.score(X_train_scaled, y_train)
        test_score = model.score(X_test_scaled, y_test)
        
        print(f"✅ Random Forest Results:")
        print(f"  Training accuracy: {train_score:.3f}")
        print(f"  Test accuracy: {test_score:.3f}")
        
        return model, scaler, selector
        
    except Exception as e:
        print(f"⚠️ Machine learning failed: {e}")
        return None, None, None

def export_simple_results(expression_data, sample_metadata, diff_results=None):
    print("\n💾 Exporting results...")
    
    os.makedirs('results', exist_ok=True)
    
    try:
        expression_data.to_csv('results/expression_data_simple.csv')
        sample_metadata.to_csv('results/sample_metadata_simple.csv', index=False)
        
        if diff_results is not None:
            diff_results.to_csv('results/differential_expression_simple.csv', index=False)
        
        summary = {
            'genes': expression_data.shape[0],
            'samples': expression_data.shape[1],
            'surface_types': sample_metadata['surface_type'].nunique(),
            'significant_genes': len(diff_results[diff_results['p_value'] < 0.05]) if diff_results is not None else 0
        }
        
        with open('results/summary_simple.txt', 'w') as f:
            for key, value in summary.items():
                f.write(f"{key}: {value}\n")
        
        print("✅ Results exported to results/ directory")
        
    except Exception as e:
        print(f"⚠️ Export failed: {e}")

def main():
    print("🧬 Simple Human Cell Gene Expression Analysis")
    print("=" * 50)
    
    try:
        analyzer, expression_data, sample_metadata = load_and_preprocess_data()
        
        explore_data(expression_data, sample_metadata)
        
        create_simple_plots(expression_data, sample_metadata)
        
        diff_results = simple_differential_expression(analyzer, expression_data, sample_metadata)
        
        model, scaler, selector = simple_machine_learning(expression_data, sample_metadata)
        
        export_simple_results(expression_data, sample_metadata, diff_results)
        
        print("\n🎉 Simple analysis completed!")
        print("Check the results/ directory for output files and plots.")
        
    except Exception as e:
        print(f"❌ Analysis failed: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
