#! /usr/bin/python3

import os
import itertools
import click
import pandas as pd
import numpy as np

@click.command()
@click.option('-i', '--input', required=True, type=click.Path(exists=True, dir_okay=False), help='input file which contains the information for calculate output')
def cli(input):
    """ HJ Kim's Tool - get genotype combination and get frequency by population """

    Get_Combined_Freq(input)

        


class Get_Combined_Freq:
    """trait별 사용하는 SNP의 ref, alt 서열과 alt의 freq정보를 이용하여 
       1. 나올수 있는 모든 유전형의 조합을 만든다.
       2. 해당하는 유전형 조합의 빈도가 인종별로 얼마인지 구한다.
       현재 인종은 'KRGDB', 'African', 'Latino', 'East_Asian', 'Finnish', 'European'을 한정한다.
       dataframe columns= ['trait','rsid', 'ref', 'alt', 'KRGDB', 'African', 'Latino', 'East_Asian', 'Finnish', 'European']
    """
	
    
    numeric_cols = ['KRGDB', 'African', 'Latino', 'East_Asian', 'Finnish', 'European']

    def __init__(self, input):
        
        self.df = pd.read_csv(input, sep = '\t')
        self.outdir = os.path.dirname(input)
        self.fileneme = os.path.basename(input).split('.')[0]

        try:
            self.df['ref_homo'] = self.df['ref'] + self.df['ref']
            self.df['hetero'] = self.df['ref'] + self.df['alt']
            self.df['alt_homo'] = self.df['alt'] + self.df['alt']
        except:
            self.df['ref_homo'] = self.df['REF'] + self.df['REF']
            self.df['hetero'] = self.df['REF'] + self.df['ALT']
            self.df['alt_homo'] = self.df['ALT'] + self.df['ALT']   
        
        self.df_drop = self.df.drop_duplicates(subset = ['rsid'])
        
        if 'trait' in self.df.columns:
            for trait in self.df['trait'].unique():
                group = self.df.groupby("trait").get_group(trait).reset_index(drop = True)
                gt_comb = self.make_genotype_combi(group)
                result = self.calculate_freq(gt_comb)
                output = os.path.join(self.outdir, f'{self.fileneme}_{trait}.txt')
                result.to_csv(output, sep = '\t', index = False)
                
        else:
            gt_comb = self.make_genotype_combi(self.df)
            result = self.calculate_freq(gt_comb)
            output = os.path.join(self.outdir, f'{self.fileneme}_output.txt')
            result.to_csv(output, sep = '\t', index = False)


    @staticmethod
    def make_genotype_combi(df: pd.DataFrame) -> pd.DataFrame:
        """rsid, ref, alt를 이용하여 가능한 모든 유전형의 조합을 계산한다."""
        #rsid 의 조합에 의한 genotype combination생성
        genotype = list()
        snp = list()
        for i in df.index:
            rsid = df.loc[i, 'rsid']
            gt = [df.loc[i, 'ref_homo'], df.loc[i, 'hetero'], df.loc[i, 'alt_homo']]
            snp.append(rsid)
            genotype.append(gt)
        
        combi = [list(row) for row in  list(itertools.product(*genotype))]
        gt_combi = pd.DataFrame(combi, columns = snp)
    
        return gt_combi
        
    
    def calculate_freq(self, combined: pd.DataFrame) -> pd.DataFrame:
        ## 유전형의 freq값 곱해주기 
    
        df_freq = self.df.set_index('rsid')
        df_final =  combined.copy()
        
        for population in self.numeric_cols:
            tmp = combined.copy()
            for idx, rsid in enumerate(tmp.columns):
                ref_homo = df_freq.loc[rsid, 'ref_homo']
                hetero = df_freq.loc[rsid, 'hetero']
                alt_homo = df_freq.loc[rsid, 'alt_homo']
                
                ref_homo_freq = (1 - df_freq.loc[rsid, population]) * (1 - df_freq.loc[rsid, population])
                hetero_freq = (df_freq.loc[rsid, population]) * (1 - df_freq.loc[rsid, population])
                alt_homo_freq = (df_freq.loc[rsid, population])
                
                col_name = f'{rsid}_{population}'
                
                tmp[col_name] = tmp[rsid].replace(to_replace=ref_homo, value = ref_homo_freq, regex = True)
                tmp[col_name] = tmp[col_name].replace(to_replace=hetero, value = hetero_freq, regex = True)
                tmp[col_name] = tmp[col_name].replace(to_replace=alt_homo, value = alt_homo_freq, regex = True)
            
            df2 = tmp.filter(like=population)
            df2[population] = 1
            for col in tmp.filter(like=population).columns:
                df2[population] = df2[population] * df2[col]
            df_final= pd.concat([df_final, df2[[population]]], axis = 1)
        
        pd.options.display.float_format = '{:.10f}'.format
        df_final.loc['total'] = df_final.select_dtypes(pd.np.number).sum()
        
        return df_final
       

if __name__ == '__main__':
    cli()