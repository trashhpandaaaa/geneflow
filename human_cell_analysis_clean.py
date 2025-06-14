#!/usr/bin/env python3

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.metrics import classification_report, confusion_matrix, accuracy_score
from sklearn.feature_selection import VarianceThreshold
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.ensemble import GradientBoostingClassifier
import warnings

warnings.filterwarnings('ignore')

from gene_expression_analyzer_clean import GeneExpressionAnalyzer

def main():
    
    print("🧬 Human Osteosarcoma Cell Gene Expression Analysis")
    print("=" * 60)
    print("GSE1000 Study: Amino Acid Conjugated Surfaces Analysis\n")
    
    data_file = 'dataset/GSE1000_series_matrix.txt'
    if not os.path.exists(data_file):
        print(f"❌ Error: Dataset file not found at {data_file}")
        print("Please ensure the GSE1000_series_matrix.txt file is in the dataset/ directory.")
        return
    
    try:
        print("📂 SECTION 1: Loading and Exploring Dataset")
        print("-" * 40)
        
        analyzer = GeneExpressionAnalyzer(data_file)
        analyzer.parse_geo_file()
        
        print(f"✅ Dataset loaded successfully:")
        print(f"   • Total genes: {analyzer.expression_data.shape[0]:,}")
        print(f"   • Total samples: {analyzer.expression_data.shape[1]}")
        print(f"   • Data type: Gene expression microarray")
        
        sample_metadata = analyzer.get_sample_metadata()
        if sample_metadata is None:
            print("❌ Error: Could not extract sample metadata")
            return
            
        print(f"\n📋 Sample Information:")
        print(sample_metadata.to_string(index=False))
        
        print(f"\n🧪 Experimental Design Summary:")
        print(f"   • Surface types: {sample_metadata['surface_type'].nunique()}")
        print(f"   • Time points: {sample_metadata['time_point'].nunique()}")
        print(f"   • Total treatment groups: {sample_metadata['treatment_group'].nunique()}")
        
        treatment_counts = sample_metadata['treatment_group'].value_counts()
        print(f"\n📊 Samples per treatment group:")
        for treatment, count in treatment_counts.items():
            print(f"   • {treatment}: {count} sample(s)")
        
        print(f"\n📈 Expression Data Statistics:")
        stats_dict = analyzer.basic_statistics()
        
        print(f"\n🔧 SECTION 2: Data Preprocessing")
        print("-" * 40)
        
        expression_log = np.log2(analyzer.expression_data + 1)
        print("✅ Applied log2 transformation")
        
        min_expression_threshold = 2.0
        genes_above_threshold = (expression_log.mean(axis=1) > min_expression_threshold)
        expression_filtered = expression_log.loc[genes_above_threshold]
        
        print(f"🔍 Gene filtering results:")
        print(f"   • Original genes: {analyzer.expression_data.shape[0]:,}")
        print(f"   • After filtering: {expression_filtered.shape[0]:,}")
        print(f"   • Removed: {analyzer.expression_data.shape[0] - expression_filtered.shape[0]:,} low-expression genes")
        
        expression_normalized = expression_filtered.copy()
        for gene in expression_normalized.index:
            gene_values = expression_normalized.loc[gene]
            if gene_values.std() > 0:
                expression_normalized.loc[gene] = (gene_values - gene_values.mean()) / gene_values.std()
        
        print("✅ Applied Z-score normalization")
        
        sample_metadata_clean = sample_metadata.copy()
        if 'sample_id' not in sample_metadata_clean.columns:
            sample_metadata_clean['sample_id'] = expression_normalized.columns
        
        common_samples = set(sample_metadata_clean['sample_id']) & set(expression_normalized.columns)
        sample_metadata_clean = sample_metadata_clean[sample_metadata_clean['sample_id'].isin(common_samples)]
        expression_normalized = expression_normalized[list(common_samples)]
        
        sample_metadata_clean = sample_metadata_clean.set_index('sample_id').reindex(expression_normalized.columns).reset_index()
        
        print(f"✅ Aligned sample metadata: {len(sample_metadata_clean)} samples")
        print(f"📋 Final data shape: {expression_normalized.shape}")
        
        print(f"\n📊 SECTION 3: Quality Control and Visualization")
        print("-" * 40)
        
        try:
            plt.style.use('default')
            plt.rcParams['figure.figsize'] = (15, 10)
            
            print("📈 Generating quality control plots...")
            
            fig, axes = plt.subplots(2, 2, figsize=(15, 12))
            
            axes[0, 0].hist(expression_normalized.values.flatten(), bins=50, alpha=0.7, edgecolor='black')
            axes[0, 0].set_title('Distribution of Normalized Expression Values')
            axes[0, 0].set_xlabel('Expression Level (Z-score)')
            axes[0, 0].set_ylabel('Frequency')
            axes[0, 0].grid(True, alpha=0.3)
            
            corr_matrix = expression_normalized.corr()
            sns.heatmap(corr_matrix, ax=axes[0, 1], cmap='coolwarm', center=0, 
                       square=True, cbar_kws={'label': 'Correlation'})
            axes[0, 1].set_title('Sample Correlation Matrix')
            
            n_samples_to_show = min(10, expression_normalized.shape[1])
            sample_data = [expression_normalized.iloc[:, i] for i in range(n_samples_to_show)]
            axes[1, 0].boxplot(sample_data)
            axes[1, 0].set_title(f'Expression Distribution by Sample (First {n_samples_to_show})')
            axes[1, 0].set_xlabel('Sample')
            axes[1, 0].set_ylabel('Expression Level')
            axes[1, 0].tick_params(axis='x', rotation=45)
            
            if expression_normalized.shape[1] > 2:
                pca = PCA(n_components=2)
                pca_result = pca.fit_transform(expression_normalized.T)
                
                unique_groups = sample_metadata_clean['treatment_group'].unique()
                colors = plt.cm.tab10(np.linspace(0, 1, len(unique_groups)))
                color_map = dict(zip(unique_groups, colors))
                
                for i, (_, row) in enumerate(sample_metadata_clean.iterrows()):
                    color = color_map[row['treatment_group']]
                    axes[1, 1].scatter(pca_result[i, 0], pca_result[i, 1], 
                                     c=[color], s=100, alpha=0.8, 
                                     label=row['treatment_group'])
                
                axes[1, 1].set_title(f'PCA Plot (PC1: {pca.explained_variance_ratio_[0]:.1%}, PC2: {pca.explained_variance_ratio_[1]:.1%})')
                axes[1, 1].set_xlabel('Principal Component 1')
                axes[1, 1].set_ylabel('Principal Component 2')
                
                handles, labels = axes[1, 1].get_legend_handles_labels()
                by_label = dict(zip(labels, handles))
                axes[1, 1].legend(by_label.values(), by_label.keys(), bbox_to_anchor=(1.05, 1), loc='upper left')
            
            plt.tight_layout()
            plt.savefig('results/quality_control_plots.png', dpi=300, bbox_inches='tight')
            plt.show()
            print("✅ Quality control plots generated and saved")
            
        except Exception as e:
            print(f"⚠️ Warning: Could not generate some QC plots - {str(e)}")
        
        print(f"\n🎯 SECTION 4: Principal Component Analysis")
        print("-" * 40)
        
        try:
            pca_detailed = PCA(n_components=min(5, expression_normalized.shape[1]))
            pca_result_detailed = pca_detailed.fit_transform(expression_normalized.T)
            
            print("📊 Explained Variance by Principal Components:")
            for i, var in enumerate(pca_detailed.explained_variance_ratio_[:5]):
                print(f"   • PC{i+1}: {var*100:.1f}%")
            
            cumulative_var = np.cumsum(pca_detailed.explained_variance_ratio_)
            print(f"   • Cumulative variance (first 3 PCs): {cumulative_var[2]*100:.1f}%")
            
            fig, ax = plt.subplots(figsize=(12, 8))
            
            for group in sample_metadata_clean['treatment_group'].unique():
                mask = sample_metadata_clean['treatment_group'] == group
                indices = mask.values
                ax.scatter(pca_result_detailed[indices, 0], pca_result_detailed[indices, 1], 
                          label=group, s=100, alpha=0.7)
            
            ax.set_xlabel(f'PC1 ({pca_detailed.explained_variance_ratio_[0]*100:.1f}%)')
            ax.set_ylabel(f'PC2 ({pca_detailed.explained_variance_ratio_[1]*100:.1f}%)')
            ax.set_title('PCA Analysis: Treatment Group Separation')
            ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
            ax.grid(True, alpha=0.3)
            
            plt.tight_layout()
            plt.savefig('results/pca_analysis.png', dpi=300, bbox_inches='tight')
            plt.show()
            print("✅ Detailed PCA analysis completed")
            
        except Exception as e:
            print(f"⚠️ Warning: PCA analysis failed - {str(e)}")
        
        print(f"\n🔥 SECTION 5: Gene Expression Heatmap")
        print("-" * 40)
        
        try:
            n_genes = 50
            top_variable_genes = expression_normalized.var(axis=1).nlargest(n_genes).index
            heatmap_data = expression_normalized.loc[top_variable_genes]
            
            sample_labels = [f"{row['sample_id']}\n{row['surface_type']}\n{row['time_point']}" 
                           for _, row in sample_metadata_clean.iterrows()]
            
            plt.figure(figsize=(14, 12))
            sns.heatmap(heatmap_data, 
                       xticklabels=sample_labels,
                       yticklabels=True,
                       cmap='RdYlBu_r', 
                       center=0,
                       cbar_kws={'label': 'Normalized Expression (Z-score)'},
                       linewidths=0.1)
            
            plt.title(f'Top {n_genes} Most Variable Genes Across Samples', fontsize=16, pad=20)
            plt.xlabel('Samples (Sample ID / Surface Type / Time Point)', fontsize=12)
            plt.ylabel('Gene IDs', fontsize=12)
            plt.xticks(rotation=45, ha='right')
            plt.yticks(fontsize=8)
            plt.tight_layout()
            plt.savefig('results/gene_expression_heatmap.png', dpi=300, bbox_inches='tight')
            plt.show()
            
            print(f"✅ Heatmap shows expression patterns of {len(top_variable_genes)} most variable genes")
            
        except Exception as e:
            print(f"⚠️ Warning: Could not generate heatmap - {str(e)}")
        
        print(f"\n🤖 SECTION 6: Machine Learning Classification")
        print("-" * 40)
        
        try:
            print("🚀 Training machine learning models...")
            
            X = expression_normalized.T
            y_surface = sample_metadata_clean['surface_type']
            
            print(f"   • Features (genes): {X.shape[1]:,}")
            print(f"   • Samples: {X.shape[0]}")
            print(f"   • Classes: {y_surface.nunique()}")
            
            selector = VarianceThreshold(threshold=0.1)
            X_selected = selector.fit_transform(X)
            selected_genes = expression_normalized.index[selector.get_support()]
            
            print(f"   • Genes after variance filtering: {X_selected.shape[1]:,}")
            
            X_train, X_test, y_train, y_test = train_test_split(
                X_selected, y_surface, test_size=0.3, random_state=42, stratify=y_surface
            )
            
            scaler = StandardScaler()
            X_train_scaled = scaler.fit_transform(X_train)
            X_test_scaled = scaler.transform(X_test)
            
            models = {
                'Random Forest': RandomForestClassifier(n_estimators=100, random_state=42),
                'Logistic Regression': LogisticRegression(random_state=42, max_iter=1000),
                'SVM': SVC(random_state=42, kernel='rbf'),
                'Gradient Boosting': GradientBoostingClassifier(random_state=42)
            }
            
            model_results = {}
            print(f"\n📊 Model Performance Results:")
            
            for name, model in models.items():
                try:
                    model.fit(X_train_scaled, y_train)
                    train_score = model.score(X_train_scaled, y_train)
                    test_score = model.score(X_test_scaled, y_test)
                    
                    cv_scores = cross_val_score(model, X_train_scaled, y_train, cv=3)
                    
                    model_results[name] = {
                        'train_accuracy': train_score,
                        'test_accuracy': test_score,
                        'cv_mean': cv_scores.mean(),
                        'cv_std': cv_scores.std(),
                        'model': model
                    }
                    
                    print(f"   • {name}:")
                    print(f"     - Train: {train_score:.3f}")
                    print(f"     - Test: {test_score:.3f}")
                    print(f"     - CV: {cv_scores.mean():.3f} ± {cv_scores.std():.3f}")
                    
                except Exception as e:
                    print(f"   • {name}: Failed - {str(e)}")
            
            if model_results:
                best_model_name = max(model_results.keys(), key=lambda k: model_results[k]['test_accuracy'])
                best_model = model_results[best_model_name]['model']
                
                print(f"\n🏆 Best performing model: {best_model_name}")
                print(f"   Test Accuracy: {model_results[best_model_name]['test_accuracy']:.3f}")
                
                y_pred = best_model.predict(X_test_scaled)
                
                print(f"\n📋 Classification Report ({best_model_name}):")
                print(classification_report(y_test, y_pred))
                
                cm = confusion_matrix(y_test, y_pred)
                plt.figure(figsize=(8, 6))
                sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', 
                           xticklabels=best_model.classes_, yticklabels=best_model.classes_)
                plt.title(f'Confusion Matrix - {best_model_name}')
                plt.xlabel('Predicted')
                plt.ylabel('Actual')
                plt.tight_layout()
                plt.savefig('results/confusion_matrix.png', dpi=300, bbox_inches='tight')
                plt.show()
                
                if hasattr(best_model, 'feature_importances_'):
                    importance_df = pd.DataFrame({
                        'gene': selected_genes,
                        'importance': best_model.feature_importances_
                    }).sort_values('importance', ascending=False)
                    
                    print(f"\n🎯 Top 10 Most Important Genes:")
                    for idx, row in importance_df.head(10).iterrows():
                        print(f"   • {row['gene']}: {row['importance']:.4f}")
                    
                    plt.figure(figsize=(12, 8))
                    top_20 = importance_df.head(20)
                    plt.barh(range(len(top_20)), top_20['importance'], alpha=0.7)
                    plt.yticks(range(len(top_20)), top_20['gene'].values)
                    plt.xlabel('Feature Importance')
                    plt.title(f'Top 20 Most Important Genes ({best_model_name})')
                    plt.gca().invert_yaxis()
                    plt.grid(True, alpha=0.3)
                    plt.tight_layout()
                    plt.savefig('results/feature_importance.png', dpi=300, bbox_inches='tight')
                    plt.show()
                
        except Exception as e:
            print(f"⚠️ Warning: Machine learning analysis failed - {str(e)}")
            model_results = {}
            selected_genes = []
        
        print(f"\n🧬 SECTION 7: Differential Expression Analysis")
        print("-" * 40)
        
        try:
            print("🔬 Performing differential expression analysis...")
            
            control_samples = sample_metadata_clean[sample_metadata_clean['surface_type'] == 'Control_Polystyrene']['sample_id'].tolist()
            treatment_samples = sample_metadata_clean[sample_metadata_clean['surface_type'].isin(['Aspartic_Acid', 'Glutamic_Acid', 'Serine'])]['sample_id'].tolist()
            
            print(f"📊 Comparison Groups:")
            print(f"   • Control (Polystyrene): {control_samples}")
            print(f"   • Treatment (Amino Acid Surfaces): {treatment_samples}")
            
            if len(control_samples) > 0 and len(treatment_samples) > 0:
                analyzer.expression_data = expression_filtered
                
                diff_expr_results = analyzer.differential_expression_analysis(
                    group1_samples=control_samples,
                    group2_samples=treatment_samples,
                    method='ttest'
                )
                
                significant_genes = diff_expr_results[diff_expr_results['p_value_adjusted'] < 0.05]
                upregulated = significant_genes[significant_genes['log2_fold_change'] > 0]
                downregulated = significant_genes[significant_genes['log2_fold_change'] < 0]
                
                print(f"\n📈 Differential Expression Results:")
                print(f"   • Total genes analyzed: {len(diff_expr_results)}")
                print(f"   • Significant genes (adj. p < 0.05): {len(significant_genes)}")
                print(f"   • Upregulated in treatment: {len(upregulated)}")
                print(f"   • Downregulated in treatment: {len(downregulated)}")
                
                if len(upregulated) > 0:
                    print(f"\n🔝 Top 5 Upregulated Genes:")
                    top_up = upregulated.nsmallest(5, 'p_value_adjusted')
                    for idx, row in top_up.iterrows():
                        print(f"   • {row['gene_id']}: FC = {row['fold_change']:.2f}, adj. p = {row['p_value_adjusted']:.2e}")
                
                if len(downregulated) > 0:
                    print(f"\n🔻 Top 5 Downregulated Genes:")
                    top_down = downregulated.nsmallest(5, 'p_value_adjusted')
                    for idx, row in top_down.iterrows():
                        print(f"   • {row['gene_id']}: FC = {row['fold_change']:.2f}, adj. p = {row['p_value_adjusted']:.2e}")
                
                print(f"\n🌋 Creating volcano plot...")
                try:
                    results = diff_expr_results.copy()
                    results['-log10_p'] = -np.log10(results['p_value'])
                    p_threshold = 0.05
                    fc_threshold = 1.5
                    results['significant'] = (results['p_value_adjusted'] < p_threshold) & (np.abs(results['log2_fold_change']) > np.log2(fc_threshold))
                    
                    plt.figure(figsize=(10, 8))
                    
                    non_sig = results[~results['significant']]
                    plt.scatter(non_sig['log2_fold_change'], non_sig['-log10_p'], 
                               c='gray', alpha=0.6, s=30, label='Non-significant')
                    
                    sig = results[results['significant']]
                    if len(sig) > 0:
                        plt.scatter(sig['log2_fold_change'], sig['-log10_p'], 
                                   c='red', alpha=0.8, s=50, label='Significant')
                    
                    plt.axhline(y=-np.log10(p_threshold), color='blue', linestyle='--', alpha=0.7, label=f'p = {p_threshold}')
                    plt.axvline(x=np.log2(fc_threshold), color='blue', linestyle='--', alpha=0.7)
                    plt.axvline(x=-np.log2(fc_threshold), color='blue', linestyle='--', alpha=0.7)
                    
                    plt.xlabel('Log2 Fold Change')
                    plt.ylabel('-Log10 P-value')
                    plt.title('Volcano Plot: Differential Gene Expression')
                    plt.legend()
                    plt.grid(True, alpha=0.3)
                    
                    if len(sig) > 0:
                        top_sig = sig.nsmallest(5, 'p_value')
                        for idx, row in top_sig.iterrows():
                            plt.annotate(row['gene_id'], 
                                       (row['log2_fold_change'], row['-log10_p']),
                                       xytext=(5, 5), textcoords='offset points', 
                                       fontsize=8, alpha=0.8)
                    
                    plt.tight_layout()
                    plt.savefig('results/volcano_plot.png', dpi=300, bbox_inches='tight')
                    plt.show()
                    print("✅ Volcano plot created successfully")
                    
                except Exception as e:
                    print(f"⚠️ Warning: Could not create volcano plot - {str(e)}")
                
            else:
                print("⚠️ Warning: Insufficient samples for differential expression analysis")
                diff_expr_results = None
                significant_genes = pd.DataFrame()
                upregulated = pd.DataFrame()
                downregulated = pd.DataFrame()
                
        except Exception as e:
            print(f"⚠️ Warning: Differential expression analysis failed - {str(e)}")
            diff_expr_results = None
            significant_genes = pd.DataFrame()
            upregulated = pd.DataFrame()
            downregulated = pd.DataFrame()
        
        print(f"\n💾 SECTION 8: Exporting Results")
        print("-" * 40)
        
        results_dir = 'results'
        os.makedirs(results_dir, exist_ok=True)
        
        try:
            expression_normalized.to_csv(f'{results_dir}/normalized_expression_data.csv')
            expression_filtered.to_csv(f'{results_dir}/log_transformed_expression_data.csv')
            
            sample_metadata_clean.to_csv(f'{results_dir}/sample_metadata.csv', index=False)
            
            if diff_expr_results is not None and len(diff_expr_results) > 0:
                diff_expr_results.to_csv(f'{results_dir}/differential_expression_results.csv', index=False)
                significant_genes.to_csv(f'{results_dir}/significant_genes.csv', index=False)
            
            if model_results:
                model_summary = pd.DataFrame(model_results).T
                model_summary.to_csv(f'{results_dir}/model_performance_comparison.csv')
                
                best_model_name = max(model_results.keys(), key=lambda k: model_results[k]['test_accuracy'])
                best_model = model_results[best_model_name]['model']
                if hasattr(best_model, 'feature_importances_') and len(selected_genes) > 0:
                    importance_df = pd.DataFrame({
                        'gene': selected_genes,
                        'importance': best_model.feature_importances_
                    }).sort_values('importance', ascending=False)
                    importance_df.to_csv(f'{results_dir}/feature_importance.csv', index=False)
            
            summary_stats = {
                'total_genes': expression_normalized.shape[0],
                'total_samples': expression_normalized.shape[1],
                'surface_types': sample_metadata_clean['surface_type'].nunique(),
                'time_points': sample_metadata_clean['time_point'].nunique(),
                'significant_genes': len(significant_genes) if len(significant_genes) > 0 else 0,
                'best_model_accuracy': max([r['test_accuracy'] for r in model_results.values()]) if model_results else 'N/A'
            }
            
            with open(f'{results_dir}/analysis_summary.txt', 'w') as f:
                f.write("Human Cell Gene Expression Analysis Summary\n")
                f.write("=" * 50 + "\n\n")
                for key, value in summary_stats.items():
                    f.write(f"{key.replace('_', ' ').title()}: {value}\n")
            
            print(f"✅ Results exported to '{results_dir}/' directory:")
            print(f"   • normalized_expression_data.csv")
            print(f"   • log_transformed_expression_data.csv")
            print(f"   • sample_metadata.csv")
            if diff_expr_results is not None:
                print(f"   • differential_expression_results.csv")
                print(f"   • significant_genes.csv")
            if model_results:
                print(f"   • model_performance_comparison.csv")
                print(f"   • feature_importance.csv")
            print(f"   • analysis_summary.txt")
            print(f"   • Various plot files (.png)")
            
        except Exception as e:
            print(f"⚠️ Warning: Could not export some results - {str(e)}")
        
        print(f"\n🎉 ANALYSIS COMPLETED SUCCESSFULLY!")
        print("=" * 60)
        print(f"\n📊 Final Summary:")
        print(f"   • Dataset: {expression_normalized.shape[0]:,} genes, {expression_normalized.shape[1]} samples")
        print(f"   • Surface types analyzed: {sample_metadata_clean['surface_type'].nunique()}")
        print(f"   • Significant genes found: {len(significant_genes) if len(significant_genes) > 0 else 0}")
        print(f"   • Best ML model accuracy: {max([r['test_accuracy'] for r in model_results.values()]):.3f}" if model_results else "   • ML analysis: Not completed")
        
        print(f"\n🔬 Key Findings:")
        if len(significant_genes) > 0:
            print(f"   • {len(significant_genes)} genes significantly respond to amino acid surface treatments")
            print(f"   • {len(upregulated)} genes are upregulated in treatment conditions")
            print(f"   • {len(downregulated)} genes are downregulated in treatment conditions")
        if model_results:
            best_accuracy = max([r['test_accuracy'] for r in model_results.values()])
            print(f"   • Machine learning can classify surface treatments with {best_accuracy:.1%} accuracy")
        
        print(f"\n📁 All results and visualizations saved to: {results_dir}/")
        print(f"🚀 Analysis pipeline completed successfully!")
        
    except Exception as e:
        print(f"❌ Critical error during analysis: {str(e)}")
        print("Please check that all dependencies are installed and the dataset file is accessible.")
        import traceback
        print(f"\nDetailed error:")
        traceback.print_exc()

if __name__ == "__main__":
    plt.style.use('default')
    
    os.makedirs('results', exist_ok=True)
    
    main()
