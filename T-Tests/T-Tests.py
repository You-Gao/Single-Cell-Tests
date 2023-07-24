import pandas as pd
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns
import time 
import os 
import sys

class Analysis:
    def __init__(self, lower_quartile, upper_quartile, p_value, data_path = r"C:\Users\commo\OneDrive - University of Virginia\School\STEM\BME\Fallahi Sichani Lab Work\Gene Co-Expression KO Analysis" , filter = 1):
        
        # Parameters
        self.lower_quartile = float(lower_quartile)
        self.upper_quartile = float(upper_quartile)
        self.p_value = float(p_value)
        self.filter = filter
        
        # Data
        self.data_path = data_path + '\\'
        self.save_path = "Results" + "_" + str(lower_quartile) + "_" + str(upper_quartile) + "_" + str(p_value) + "\\"
        self.mrna_data = pd.read_csv(self.data_path + 'Depmap Melanoma mRNA.csv')
        self.pancancer_data = pd.read_csv(self.data_path + 'Depmap Pancancer mRNA.csv')
        self.ko_data = pd.read_csv(self.data_path + 'Depmap Melanoma Gene Dependency.csv')
        self.cell_data = pd.read_csv(self.data_path + 'Depmap Melanoma Cells.csv')
        self.chromatin_genes = pd.read_csv(self.data_path + 'Chromatin Genes.csv')
        
        # Calculated Data
        self.pancancer_quartiles = None
        self.mrna_t_test_results = None
        self.chromatin_gene_expression_encoding = None
        self.p_value_matrix = None
        self.effect_matrix = None
        self.significant_matrix = None
        
    
    def pancancer_quartiles_set(self, pancancer_data):
        self.pancancer_quartiles = pancancer_data.quantile([self.lower_quartile, self.upper_quartile])
        
    def mrna_data_filter(self):
        mRNA_genes = self.mrna_data.columns[1:].tolist()
        chromatin_genes = self.chromatin_genes['Gene'].tolist()
        intersection = list(set(chromatin_genes).intersection(mRNA_genes))
        
        mrna_data = self.mrna_data[intersection]
        pancancer_data = self.pancancer_data[intersection]
        
        self.pancancer_quartiles_set(pancancer_data)
    
    def mrna_ttests(self):
        results = []
        chromatin_gene_groups = pd.DataFrame()
        
        for gene in self.pancancer_quartiles.columns:
            # Index by gene column and then quantile row index
            low_expression = self.pancancer_quartiles[gene][self.lower_quartile]
            
            # Testing for > 25th percentile or < 75th percentile
            high_expression = self.pancancer_quartiles[gene][self.upper_quartile]
            
            # Create a filter and then index by the gene column
            low_expression_cells = self.mrna_data[self.mrna_data[gene] <= low_expression][['Cell Line',gene]]
            high_expression_cells = self.mrna_data[self.mrna_data[gene] > high_expression][['Cell Line',gene]]
            
            # T-Test of Gene mRNA Levels in Low and High Expression Cells
            ttest = stats.ttest_ind(low_expression_cells[gene], high_expression_cells[gene],  alternative = 'two-sided')
            
            # Append the results to the results list
            # Format: (Gene, Low Samples, Low Expression Mean, High Samples High Expression Mean, T-Stat, P-Value)
            results.append((gene, len(low_expression_cells), low_expression_cells[gene].mean(), len(high_expression_cells), high_expression_cells[gene].mean(), ttest[0], ttest[1]))
            
            # Add a column to the low expression cells called "Group" and assign it the value "low"
            low_expression_cells['Group'] = 'low'
            high_expression_cells['Group'] = 'high'
            
            # Sort the cells by cell line so that they are in the same order as the gene_groups dataframe
            expression = pd.concat([low_expression_cells, high_expression_cells])
            expression.sort_values(by=['Cell Line'], inplace=True)
            
            # Append the low and high expression cells to the gene_groups list
            chromatin_gene_groups = pd.concat([chromatin_gene_groups, expression['Group']], axis=1)
            chromatin_gene_groups.rename(columns={'Group': gene}, inplace=True)
        chromatin_ttest_df = pd.DataFrame(results, columns=['Gene','Low Samples','Low Expression Mean','High Samples','High Expression Mean','T-Stat','P-Value'])
        chromatin_ttest_df.sort_values(by=['Low Samples','High Samples'], ascending=[False,True], inplace=True)
        chromatin_ttest_df.reset_index(drop=True, inplace=True)
        
        self.mrna_t_test_results = chromatin_ttest_df
        self.chromatin_gene_expression_encoding  = chromatin_gene_groups
        
    def sample_filter(self):
        # Filter the groups by ones with both low and high sample sizes > 10
        self.mrna_t_test_results = self.mrna_t_test_results[(self.mrna_t_test_results['Low Samples'] > 2) & (self.mrna_t_test_results['High Samples'] > 2)]

        # Retrieve the Genes
        chromatin_ttest_df_genes_filter = self.mrna_t_test_results['Gene']

        # Filter the gene groups by the filtered genes
        self.chromatin_gene_expression_encoding  = self.chromatin_gene_expression_encoding[chromatin_ttest_df_genes_filter]
    
    def two_sided_ttest(self, filter):
        if filter:
            self.sample_filter()
        
        self.chromatin_gene_expression_encoding = self.chromatin_gene_expression_encoding.reindex(sorted(self.chromatin_gene_expression_encoding.columns), axis=1)
        
        vip_genes = pd.Series(list(self.chromatin_gene_expression_encoding.columns))
        vip_genes.sort_values(inplace=True)
        vip_genes.reset_index(drop=True, inplace=True)
        
        self.p_value_matrix = np.zeros((self.chromatin_gene_expression_encoding.shape[1], self.chromatin_gene_expression_encoding.shape[1]))
        self.effect_matrix = np.zeros((self.chromatin_gene_expression_encoding.shape[1], self.chromatin_gene_expression_encoding.shape[1]))
        
        for gene in vip_genes:
            low_expression = self.chromatin_gene_expression_encoding[self.chromatin_gene_expression_encoding[gene] == 'low'].index
            high_expression = self.chromatin_gene_expression_encoding[self.chromatin_gene_expression_encoding[gene] == 'high'].index
            
            for gene2 in vip_genes:
                low_expression_dependency = self.ko_data.loc[low_expression, gene2]
                high_expression_dependency = self.ko_data.loc[high_expression, gene2]
                
                # T-Test of Gene Dependency in Low and High Expression Cells
                ttest = stats.ranksums(low_expression_dependency, high_expression_dependency,  alternative = 'two-sided')
    
                # Difference in Gene Dependency in Low and High Expression Cells
                self.effect_matrix[vip_genes[vip_genes == gene].index[0], vip_genes[vip_genes == gene2].index[0]] = np.mean(low_expression_dependency) - np.mean(high_expression_dependency)
                
                # Add the t-statistic to the gene effect matrix
                self.p_value_matrix[vip_genes[vip_genes == gene].index[0], vip_genes[vip_genes == gene2].index[0]] = ttest[1]
                
        self.p_value_matrix = pd.DataFrame(self.p_value_matrix, columns=vip_genes, index=vip_genes)
        self.p_value_matrix = -np.log10(self.p_value_matrix)
        
        self.effect_matrix = pd.DataFrame(self.effect_matrix, columns=vip_genes, index=vip_genes)
                
    def significant_fliter(self, p_value):
        LOG_P_VALUE = -np.log10(p_value)
        self.significant_matrix = self.p_value_matrix[self.p_value_matrix > LOG_P_VALUE]
        
    def figures_factory(self, p_value):
        LOG_P_VALUE = -np.log10(p_value)
    
        # ---------------------- #
        pos_significant_gene_tuples = []
        
        

        for gene in self.significant_matrix.columns:
            for gene2 in self.significant_matrix.index:
                if (self.p_value_matrix[gene][gene2] > LOG_P_VALUE) and (self.effect_matrix[gene][gene2] >= .2):
                    pos_significant_gene_tuples.append((gene, gene2))
                else:
                    continue
                
        pos_significant_gene_tuples = list(pos_significant_gene_tuples)
        pos_significant_gene_tuples_df = pd.DataFrame(pos_significant_gene_tuples, columns=['Gene 1', 'Gene 2'])
        
        # ---------------------- #
        neg_significant_gene_tuples = []

        for gene in self.significant_matrix.columns:
            for gene2 in self.significant_matrix.index:
                if (self.p_value_matrix[gene][gene2] > LOG_P_VALUE) and (self.effect_matrix[gene][gene2] <= -.2):
                    neg_significant_gene_tuples.append((gene, gene2))
                else:
                    continue

        neg_significant_gene_tuples = list(neg_significant_gene_tuples)
        print(len(neg_significant_gene_tuples))

        neg_significant_gene_tuples_df = pd.DataFrame(neg_significant_gene_tuples, columns=['Gene 1', 'Gene 2'])
        neg_significant_gene_tuples_df.head(5)
        
        # ---------------------- #
        sns.set_style('white')
        
        plt.figure(figsize=(10,10))

        for gene in self.p_value_matrix.columns:
            plt.scatter(self.effect_matrix[gene], self.p_value_matrix[gene], alpha=0.5, color='whitesmoke')

        # Axis Labels
        plt.xlabel('Difference in Gene Dependency (Y (XLow) - Y (XHigh)', size = 15)
        plt.ylabel('Log of p-value of Gene Dependency in Low vs. High Expression Cells', size = 15)
        plt.xticks(size=10)
        plt.yticks(size=10)

        plt.hlines(y = LOG_P_VALUE, xmin = -.7, xmax = .7, color='grey', linestyles='dashed')
        plt.vlines(x = .2, ymin = 0, ymax = 8, color='grey', linestyles='dashed')
        plt.vlines(x = -.2, ymin = 0, ymax = 8, color='grey', linestyles='dashed')

        # Color genes in pos_significant_gene_tuples_df
        for i, genes in pos_significant_gene_tuples_df.iterrows():
                plt.scatter(self.effect_matrix[genes[0]][genes[1]], self.p_value_matrix[genes[0]][genes[1]], alpha=0.5, color='red')
        # Color genes in neg_significant_gene_tuples_df
        for i, genes in neg_significant_gene_tuples_df.iterrows():
                plt.scatter(self.effect_matrix[genes[0]][genes[1]], self.p_value_matrix[genes[0]][genes[1]], alpha=0.5, color='blue')
                
        name = 'Volcano Plot of Significant Negative and Positive Gene Pairs'
        plt.savefig(self.save_path + name + '.png', dpi=300, bbox_inches='tight')
        
        # ---------------------- #
        
        plt.figure(figsize=(10,10))

        for gene in self.p_value_matrix.columns:
            plt.scatter(self.effect_matrix[gene], self.p_value_matrix[gene], alpha=0.5, color='whitesmoke')

        # Axis Labels
        plt.title('Significant Gene Combinations (Effect Size > 0.2 and p-value < 0.001)')
        plt.xlabel('Difference in Gene Dependency (Y (XLow) - Y (XHigh)', size = 15)
        plt.ylabel('Log of p-value of Gene Dependency in Low vs. High Expression Cells', size = 15)
        plt.xticks(size=10)
        plt.yticks(size=10)

        plt.hlines(y = LOG_P_VALUE, xmin = -.7, xmax = .7, color='grey', linestyles='dashed')
        plt.vlines(x = .2, ymin = 0, ymax = 8, color='grey', linestyles='dashed')
        plt.vlines(x = -.2, ymin = 0, ymax = 8, color='grey', linestyles='dashed')
        
                # Color genes in pos_significant_gene_tuples_df
        for i, genes in pos_significant_gene_tuples_df.iterrows():
                plt.scatter(self.effect_matrix[genes[0]][genes[1]], self.p_value_matrix[genes[0]][genes[1]], alpha=0.7, color='red')
        # Color genes in neg_significant_gene_tuples_df
        for i, genes in pos_significant_gene_tuples_df.iterrows():
                plt.scatter(self.effect_matrix[genes[1]][genes[0]], self.p_value_matrix[genes[1]][genes[0]], alpha=0.7, color='blue')
                
        
        # Red and Blue in legend
        plt.scatter(0,0, alpha=0.7, color='red', label='P(X,Y)')
        plt.scatter(0,0, alpha=0.7, color='blue', label='P(Y,X)')
        plt.legend(loc='upper left', fontsize=15)

        name = 'Positive Significant Gene Combinations'
        plt.savefig(self.save_path + name + '.png', dpi=300, bbox_inches='tight')
        
        # ---------------------- #
        
        plt.figure(figsize=(10,10))

        for gene in self.p_value_matrix.columns:
            plt.scatter(self.effect_matrix[gene], self.p_value_matrix[gene], alpha=0.5, color='whitesmoke')

        # Axis Labels
        plt.title('Significant Gene Combinations (Effect Size < -0.2 and p-value < 0.001)')
        plt.xlabel('Difference in Gene Dependency (Y (XLow) - Y (XHigh)', size = 15)
        plt.ylabel('Log of p-value of Gene Dependency in Low vs. High Expression Cells', size = 15)
        plt.xticks(size=10)
        plt.yticks(size=10)

        plt.hlines(y = LOG_P_VALUE, xmin = -.7, xmax = .7, color='grey', linestyles='dashed')
        plt.vlines(x = .2, ymin = 0, ymax = 8, color='grey', linestyles='dashed')
        plt.vlines(x = -.2, ymin = 0, ymax = 8, color='grey', linestyles='dashed')
        
        # Color genes in pos_significant_gene_tuples_df
        for i, genes in neg_significant_gene_tuples_df.iterrows():
                plt.scatter(self.effect_matrix[genes[1]][genes[0]], self.p_value_matrix[genes[1]][genes[0]], alpha=0.7, color='green')
        # Color genes in neg_significant_gene_tuples_df
        for i, genes in neg_significant_gene_tuples_df.iterrows():
                plt.scatter(self.effect_matrix[genes[0]][genes[1]], self.p_value_matrix[genes[0]][genes[1]], alpha=0.7, color='purple')
                
        name = 'Negative Significant Gene Combinations'
        plt.savefig(self.save_path + name + '.png', dpi=300, bbox_inches='tight')
        
        # ---------------------- #
        # Boolean array if it contains any value
        pos_y = pos_significant_gene_tuples_df['Gene 1'].value_counts().index
        pos_x = pos_significant_gene_tuples_df['Gene 2'].value_counts().index

        pos_heatmap_df = pd.DataFrame(index=pos_x, columns=pos_y)

        pos_heatmap_df.sort_index(inplace=True)
        pos_heatmap_df.sort_index(axis=1, inplace=True)

        for i, genes in pos_significant_gene_tuples_df.iterrows():
            pos_heatmap_df.loc[genes[1]][genes[0]] = 1
        pos_heatmap_df.fillna(0, inplace=True)
        
        neg_y = neg_significant_gene_tuples_df['Gene 1'].value_counts().index
        neg_x = neg_significant_gene_tuples_df['Gene 2'].value_counts().index

        neg_heatmap_df = pd.DataFrame(index=neg_x, columns=neg_y)
        
        neg_heatmap_df.sort_index(inplace=True)
        neg_heatmap_df.sort_index(axis=1, inplace=True)

        for i, genes in neg_significant_gene_tuples_df.iterrows():
            neg_heatmap_df.loc[genes[1]][genes[0]] = 1
        
        
        # ---------------------- #
        
        pos_heatmap_df_bool = pos_heatmap_df.astype(bool)
        # drop if sum is less than 1
        pos_heatmap_expression_sum = pos_heatmap_df_bool.sum(axis=1).sort_values(ascending=False)
        
        plt.figure(figsize=(7,5))

        try:
            pos_heatmap_expression_sum[pos_heatmap_expression_sum > 1].plot(kind='bar')
        except:
            print('No significant genes for mRNA expression')
            


        plt.title('Frequency of Significant Genes for mRNA Expression (>1)', size = 15)
        plt.ylabel('Frequency', size = 15)
        plt.xlabel('Genes for mRNA Expression', size = 15)
        plt.xticks(size=10)
        plt.yticks(size=10)

        plt.ylim(0,8)

        name = 'Positive Frequency of Significant Genes for mRNA Expression'
        plt.savefig(self.save_path + name + '.png', dpi=300, bbox_inches='tight')

        # ---------------------- #
        pos_heatmap_dependency_sum = pos_heatmap_df_bool.sum(axis=0).sort_values(ascending=False)
        
        plt.figure(figsize=(7,5))

        try:
            pos_heatmap_dependency_sum[pos_heatmap_dependency_sum > 1].plot(kind='bar')
        except:
            print('No significant genes for gene dependency')
            

        plt.title('Frequency of Significant Genes for Gene Dependency (>1)', size = 15)
        plt.ylabel('Frequency', size = 15)
        plt.xlabel('Genes for Dependency', size = 15)
        plt.xticks(size=10)
        plt.yticks(size=10)

        name = 'Positive Frequency of Significant Genes for Gene Dependency'
        plt.savefig(self.save_path + name + '.png', dpi=300, bbox_inches='tight')
        
        # ---------------------- #
        
        # Boolean array if it contains any value
        neg_heatmap_df_bool = neg_heatmap_df.astype(bool)
        
        plt.figure(figsize=(7,5))

        # drop if sum is less than 1
        neg_heatmap_expression_sum = neg_heatmap_df_bool.sum(axis=1).sort_values(ascending=False)
        try:
            neg_heatmap_expression_sum[neg_heatmap_expression_sum > 1].plot(kind='bar')
        except:
            print('No significant genes for mRNA expression')
            
        plt.title('Frequency of Significant Genes for mRNA Expression (>1)', size = 15)
        plt.ylabel('Frequency', size = 15)
        plt.xlabel('Genes for mRNA Expression', size = 15)
        plt.xticks(size=10)
        plt.yticks(size=10)

        plt.ylim(0,10)

        name = 'Negative Frequency of Significant Genes for mRNA Expression'
        plt.savefig(self.save_path + name + '.png', dpi=300, bbox_inches='tight')

        plt.figure(figsize=(7,5))
        neg_heatmap_dependency_sum = neg_heatmap_df_bool.sum(axis=0).sort_values(ascending=False)

        try:
            neg_heatmap_dependency_sum[neg_heatmap_dependency_sum > 1].plot(kind='bar')
        except:
            print('No significant genes for gene dependency')

        plt.title('Frequency of Significant Genes for Gene Dependency (>1)', size = 15)
        plt.ylabel('Frequency', size = 15)
        plt.xlabel('Genes for Dependency', size = 15)
        plt.xticks(size=10)
        plt.yticks(size=10)

        plt.ylim(0,10)

        name = 'Negative Frequency of Significant Genes for Gene Dependency'
        plt.savefig(self.save_path + name + '.png', dpi=300, bbox_inches='tight')
        
                
        
    def run(self):
        start = time.time()
        print('Running...')
        self.mrna_data_filter()
        self.mrna_ttests()
        self.two_sided_ttest(self.filter)
        self.significant_fliter(self.p_value)
        self.figures_factory(self.p_value)
        end = time.time()
        print(f'Took {end - start} seconds to run')
        

 
            
            
            
            
            
        


if __name__ == '__main__':
    with open('param_list.txt', 'r') as f:
        params = f.readlines()
        for param_line in params:
            args = param_line.split()
            print("Running with parameters: {}".format(args))
            if not os.path.exists("Results" + "_" + str(args[0]) + "_" + str(args[1]) + "_" + str(args[2])):
                os.mkdir("Results" + "_" + str(args[0]) + "_" + str(args[1]) + "_" + str(args[2]))
            analysis = Analysis(args[0], args[1],  args[2])
            analysis.run()
            
            
        
    
    