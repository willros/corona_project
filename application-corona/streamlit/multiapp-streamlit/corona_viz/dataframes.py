import pandas as pd
from Bio import SeqIO

class CoronaDataframe:
    '''Class that handles manipulation and generation of dataframes
    needed for plotting'''
    
    def __init__(self, query: str, ref: str, mutations: str) -> None:
        
        # Only include region between 256:29674 (same as pangolin)
        self.query = SeqIO.read(query, 'fasta').seq[256:29674]
        self.ref = SeqIO.read(ref, 'fasta').seq[256:29674]
        self.mutations = mutations
    
    
    def make_alignment_df(self) -> pd.DataFrame:
        '''Turns the sequence of the aligned fasta file into a df used to plot 
        the missing (N), deletions and mutations as a genome map'''
    
        reference = [x for x in self.ref]
        query = [x for x in self.query]
        converted_values = self._convert_alignment_value()

        alignment_df = pd.DataFrame({'nt': query, 'id': 'query', 
                                 'position': range(256, len(query) + 256), 
                                 'value': converted_values, 
                                 'ref': reference})
        return alignment_df
    
    
    def make_mutations_df(self) -> pd.DataFrame:
        '''Transforms the df of mutations to a suitable format'''
        
        mutations_df = pd.read_csv(self.mutations)
        mutations_df['unique'] = [1 if x == 1 else 0 for x in mutations_df.number_mutations]
        mutations_df['kind'] = [2 if '-' in str(x) else y for x,y in zip(mutations_df.mutation, 
                                                                         mutations_df.unique)]
        mutations_df['position'] = mutations_df['mutation'].str.extract('(\d+)').astype(int)

        return mutations_df
    
    
    def make_query_df(self) -> pd.DataFrame:
        '''Turns a list of mutations for the query into a df in the right format'''
        
        mutations_list = self._extract_snp()
        query_df = pd.DataFrame({'pango': 'QUERY', 'mutation': mutations_list})
        query_df['position'] = query_df['mutation'].str.extract('(\d+)').astype(int)
    
        return query_df
    
    
    def make_mutations_inspection_df(self) -> pd.DataFrame:
        '''Returns a dataframe with all mutations and lineages 
        for every given position in the query mutation set'''
        mutations = self.make_mutations_df()
        query_df = self.make_query_df()
        mutations = mutations[mutations['position'].isin(query_df.position)]
        query_df = query_df[query_df['position'].isin(mutations.position)]


        mutations_list = mutations.groupby('position')[['mutation', 'pango']].aggregate(list).reset_index()
        mutations_list['mutation_list'] = [list(zip(x, y)) for x,y in zip(mutations_list.pango, mutations_list.mutation)]
        mutations_list.drop(columns=['pango', 'mutation'], inplace=True)

        merged = query_df.merge(mutations_list)
        merged['pango_with_mutation'] = merged['mutation_list'].apply(len)
        merged['color'] = [1 if x.endswith('N') else 2 if x.endswith('-') else 3 for x in merged.mutation]

        return merged
    
    def _convert_alignment_value(self) -> list:
        '''Converts alignments to numbers depending on it is N, 
        deletion, substitution or the same'''
        values = []
        for q,r in zip(self.query, self.ref):
            if q != r:
                if q == 'N':
                    values.append(-1)
                elif q == '-':
                    values.append(-2)
                else:
                    values.append(1)
            else:
                values.append(0)  
        return values
    
    #### MÅSTE ADDERA +1 till index senare!!!!!!!!!!! #######
    def _extract_snp(self) -> list:
        '''Helper function to extract information about mutation'''

        return [f'{a}{index}{b}' for index,(a,b) in 
                enumerate(zip(self.ref, self.query), 256) if a != b]