import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split, cross_val_score
import warnings

warnings.filterwarnings('ignore')

class GeneExpressionAnalyzer:
    
    def __init__(self, filepath):
        self.filepath = filepath
        self.metadata = {}
        self.expression_data = None
        self.sample_info = {}
        self.gene_annotations = {}
        
    def parse_geo_file(self):
        print("Parsing GEO series matrix file...")
        
        with open(self.filepath, 'r', encoding='utf-8') as file:
            lines = file.readlines()
        
        table_start = None
        for i, line in enumerate(lines):
            if line.strip() == '!series_matrix_table_begin':
                table_start = i + 1
                break
        
        if table_start is None:
            raise ValueError("Could not find expression data table in file")
        
        self._extract_metadata(lines[:table_start])
        self._extract_expression_data(lines[table_start:])
        
        print(f"Successfully parsed {self.expression_data.shape[0]} genes across {self.expression_data.shape[1]} samples")
        
    def _extract_metadata(self, metadata_lines):
        for line in metadata_lines:
            if line.startswith('!Series_'):
                key = line.split('\t')[0].replace('!', '').replace('"', '')
                value = line.split('\t')[1].strip().replace('"', '') if len(line.split('\t')) > 1 else ''
                self.metadata[key] = value
            elif line.startswith('!Sample_title'):
                titles = line.strip().split('\t')[1:]
                self.sample_info['titles'] = [t.replace('"', '') for t in titles]
            elif line.startswith('!Sample_geo_accession'):
                ids = line.strip().split('\t')[1:]
                self.sample_info['ids'] = [id.replace('"', '') for id in ids]
    
    def _extract_expression_data(self, data_lines):
        table_end = None
        for i, line in enumerate(data_lines):
            if line.strip() == '!series_matrix_table_end':
                table_end = i
                break
        
        if table_end is None:
            table_end = len(data_lines)
        
        expression_lines = data_lines[:table_end]
        data_rows = []
        gene_ids = []
        
        for line in expression_lines:
            if line.strip() and not line.startswith('!'):
                parts = line.strip().split('\t')
                if len(parts) > 1:
                    gene_id = parts[0].replace('"', '')
                    expression_values = [float(x) for x in parts[1:] if x.replace('.', '').replace('-', '').isdigit()]
                    if len(expression_values) > 0:
                        gene_ids.append(gene_id)
                        data_rows.append(expression_values)
        
        column_names = self.sample_info.get('ids', [f'Sample_{i+1}' for i in range(len(data_rows[0]))])
        self.expression_data = pd.DataFrame(data_rows, index=gene_ids, columns=column_names)
        
    def get_sample_metadata(self):
        if 'titles' not in self.sample_info:
            return None
            
        sample_df = pd.DataFrame({
            'sample_id': self.sample_info.get('ids', []),
            'condition': self.sample_info.get('titles', [])
        })
        
        sample_df['surface_type'] = sample_df['condition'].apply(self._parse_surface_type)
        sample_df['time_point'] = sample_df['condition'].apply(self._parse_time_point)
        sample_df['treatment_group'] = sample_df['condition'].apply(self._parse_treatment_group)
        
        return sample_df
    
    def _parse_surface_type(self, condition):
        condition_lower = condition.lower()
        if 'aspartic acid' in condition_lower:
            return 'Aspartic_Acid'
        elif 'glutamic acid' in condition_lower:
            return 'Glutamic_Acid'
        elif 'serine' in condition_lower:
            return 'Serine'
        elif 'amino' in condition_lower:
            return 'Amino_Terminated'
        elif 'polystyrene' in condition_lower:
            return 'Control_Polystyrene'
        else:
            return 'Unknown'
    
    def _parse_time_point(self, condition):
        if '6 hours' in condition:
            return '6h'
        elif '32 hours' in condition:
            return '32h'
        else:
            return 'Unknown'
    
    def _parse_treatment_group(self, condition):
        surface = self._parse_surface_type(condition)
        time = self._parse_time_point(condition)
        return f"{surface}_{time}"
    
    def basic_statistics(self):
        stats_dict = {
            'total_genes': self.expression_data.shape[0],
            'total_samples': self.expression_data.shape[1],
            'mean_expression': self.expression_data.values.mean(),
            'median_expression': np.median(self.expression_data.values),
            'std_expression': self.expression_data.values.std(),
            'min_expression': self.expression_data.values.min(),
            'max_expression': self.expression_data.values.max()
        }
        
        print("=== GENE EXPRESSION DATA STATISTICS ===")
        for key, value in stats_dict.items():
            print(f"{key.replace('_', ' ').title()}: {value:.2f}" if isinstance(value, float) else f"{key.replace('_', ' ').title()}: {value}")
        
        return stats_dict
    
    def quality_control_plots(self):
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        
        axes[0, 0].hist(self.expression_data.values.flatten(), bins=50, alpha=0.7, edgecolor='black')
        axes[0, 0].set_title('Distribution of Expression Values')
        axes[0, 0].set_xlabel('Expression Level')
        axes[0, 0].set_ylabel('Frequency')
        axes[0, 0].grid(True, alpha=0.3)
        
        corr_matrix = self.expression_data.corr()
        sns.heatmap(corr_matrix, ax=axes[0, 1], cmap='coolwarm', center=0, 
                   square=True, cbar_kws={'label': 'Correlation'})
        axes[0, 1].set_title('Sample Correlation Matrix')
        
        sample_means = self.expression_data.mean(axis=0)
        axes[1, 0].boxplot([self.expression_data.iloc[:, i] for i in range(min(10, self.expression_data.shape[1]))])
        axes[1, 0].set_title('Expression Distribution by Sample (First 10)')
        axes[1, 0].set_xlabel('Sample')
        axes[1, 0].set_ylabel('Expression Level')
        axes[1, 0].tick_params(axis='x', rotation=45)
        
        if self.expression_data.shape[1] > 2:
            pca = PCA(n_components=2)
            pca_result = pca.fit_transform(self.expression_data.T)
            
            axes[1, 1].scatter(pca_result[:, 0], pca_result[:, 1], alpha=0.7, s=100)
            axes[1, 1].set_title(f'PCA Plot (PC1: {pca.explained_variance_ratio_[0]:.1%}, PC2: {pca.explained_variance_ratio_[1]:.1%})')
            axes[1, 1].set_xlabel('Principal Component 1')
            axes[1, 1].set_ylabel('Principal Component 2')
            
            for i, txt in enumerate(self.expression_data.columns):
                axes[1, 1].annotate(txt, (pca_result[i, 0], pca_result[i, 1]), 
                                  xytext=(5, 5), textcoords='offset points', fontsize=8)
        
        plt.tight_layout()
        plt.show()
        
        return fig
        
    def differential_expression_analysis(self, group1_samples, group2_samples, method='ttest'):
        results = []
        
        for gene in self.expression_data.index:
            group1_values = self.expression_data.loc[gene, group1_samples]
            group2_values = self.expression_data.loc[gene, group2_samples]
            
            mean1 = group1_values.mean()
            mean2 = group2_values.mean()
            fold_change = mean2 / mean1 if mean1 != 0 else np.inf
            log2_fold_change = np.log2(fold_change) if fold_change > 0 and fold_change != np.inf else np.nan
            
            if method == 'ttest':
                stat, p_value = stats.ttest_ind(group1_values, group2_values)
            elif method == 'mannwhitney':
                stat, p_value = stats.mannwhitneyu(group1_values, group2_values, alternative='two-sided')
            else:
                raise ValueError("Method must be 'ttest' or 'mannwhitney'")
            
            results.append({
                'gene_id': gene,
                'mean_group1': mean1,
                'mean_group2': mean2,
                'fold_change': fold_change,
                'log2_fold_change': log2_fold_change,
                'p_value': p_value,
                'test_statistic': stat
            })
        
        results_df = pd.DataFrame(results)
        
        from scipy.stats import false_discovery_control
        results_df['p_value_adjusted'] = false_discovery_control(results_df['p_value'])
        
        results_df = results_df.sort_values('p_value')
        
        return results_df
    
    def volcano_plot(self, diff_expr_results, p_threshold=0.05, fc_threshold=1.5):
        fig, ax = plt.subplots(figsize=(10, 8))
        
        results = diff_expr_results.copy()
        results['-log10_p'] = -np.log10(results['p_value'])
        results['significant'] = (results['p_value_adjusted'] < p_threshold) & (np.abs(results['log2_fold_change']) > np.log2(fc_threshold))
        
        colors = ['gray' if not sig else 'red' for sig in results['significant']]
        scatter = ax.scatter(results['log2_fold_change'], results['-log10_p'], 
                           c=colors, alpha=0.6, s=30)
        
        ax.axhline(y=-np.log10(p_threshold), color='blue', linestyle='--', alpha=0.7, label=f'p = {p_threshold}')
        ax.axvline(x=np.log2(fc_threshold), color='blue', linestyle='--', alpha=0.7, label=f'FC = {fc_threshold}')
        ax.axvline(x=-np.log2(fc_threshold), color='blue', linestyle='--', alpha=0.7)
        
        ax.set_xlabel('Log2 Fold Change')
        ax.set_ylabel('-Log10 P-value')
        ax.set_title('Volcano Plot: Differential Gene Expression')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        top_genes = results[results['significant']].nsmallest(10, 'p_value')
        for idx, row in top_genes.iterrows():
            ax.annotate(row['gene_id'], 
                       (row['log2_fold_change'], row['-log10_p']),
                       xytext=(5, 5), textcoords='offset points', 
                       fontsize=8, alpha=0.8)
        
        plt.tight_layout()
        plt.show()
        
        return fig
    
    def machine_learning_classification(self, sample_metadata, target_column='treatment_group'):
        X = self.expression_data.T
        y = sample_metadata[target_column]
        
        common_samples = set(X.index) & set(sample_metadata['sample_id'])
        X = X.loc[list(common_samples)]
        y = sample_metadata[sample_metadata['sample_id'].isin(common_samples)][target_column]
        
        y = y.reindex(X.index)
        
        from sklearn.feature_selection import VarianceThreshold
        selector = VarianceThreshold(threshold=0.1)
        X_selected = selector.fit_transform(X)
        selected_genes = X.columns[selector.get_support()]
        
        print(f"Selected {len(selected_genes)} genes out of {X.shape[1]} for modeling")
        
        X_train, X_test, y_train, y_test = train_test_split(
            X_selected, y, test_size=0.3, random_state=42, stratify=y
        )
        
        scaler = StandardScaler()
        X_train_scaled = scaler.fit_transform(X_train)
        X_test_scaled = scaler.transform(X_test)
        
        rf = RandomForestClassifier(n_estimators=100, random_state=42)
        rf.fit(X_train_scaled, y_train)
        
        train_score = rf.score(X_train_scaled, y_train)
        test_score = rf.score(X_test_scaled, y_test)
        cv_scores = cross_val_score(rf, X_train_scaled, y_train, cv=5)
        
        feature_importance = pd.DataFrame({
            'gene': selected_genes,
            'importance': rf.feature_importances_
        }).sort_values('importance', ascending=False)
        
        results = {
            'train_accuracy': train_score,
            'test_accuracy': test_score,
            'cv_mean': cv_scores.mean(),
            'cv_std': cv_scores.std(),
            'feature_importance': feature_importance,
            'model': rf,
            'scaler': scaler,
            'selected_genes': selected_genes
        }
        
        print(f"=== MACHINE LEARNING RESULTS ===")
        print(f"Training Accuracy: {train_score:.3f}")
        print(f"Test Accuracy: {test_score:.3f}")
        print(f"Cross-validation: {cv_scores.mean():.3f} ± {cv_scores.std():.3f}")
        
        return results
    
    def plot_top_genes(self, n_genes=20, by='variance'):
        if by == 'variance':
            top_genes = self.expression_data.var(axis=1).nlargest(n_genes).index
            title = f'Top {n_genes} Most Variable Genes'
        elif by == 'mean':
            top_genes = self.expression_data.mean(axis=1).nlargest(n_genes).index
            title = f'Top {n_genes} Highest Expressed Genes'
        else:
            raise ValueError("'by' parameter must be 'variance' or 'mean'")
        
        plt.figure(figsize=(12, 8))
        sns.heatmap(self.expression_data.loc[top_genes], 
                   cmap='RdYlBu_r', center=0, 
                   xticklabels=True, yticklabels=True,
                   cbar_kws={'label': 'Expression Level'})
        plt.title(title)
        plt.xlabel('Samples')
        plt.ylabel('Genes')
        plt.tight_layout()
        plt.show()
    
    def export_results(self, output_dir='results'):
        import os
        os.makedirs(output_dir, exist_ok=True)
        
        self.expression_data.to_csv(os.path.join(output_dir, 'expression_data.csv'))
        
        sample_metadata = self.get_sample_metadata()
        if sample_metadata is not None:
            sample_metadata.to_csv(os.path.join(output_dir, 'sample_metadata.csv'), index=False)
        
        stats = self.basic_statistics()
        with open(os.path.join(output_dir, 'basic_statistics.txt'), 'w') as f:
            for key, value in stats.items():
                f.write(f"{key}: {value}\n")
        
        print(f"Results exported to {output_dir}/")

if __name__ == "__main__":
    print("Gene Expression Analysis System")
    print("===============================")
    print("This module provides tools for analyzing gene expression data.")
    print("To use: analyzer = GeneExpressionAnalyzer('path/to/your/geo_file.txt')")
    print("Then run: analyzer.parse_geo_file() to start analysis.")
