import streamlit as st
from Bio import SeqIO
import shutil
import os
import time
from corona_viz.plotting import CoronaPlot
from corona_viz.nextclade import Nextclade

st.markdown('*RUN NEXTCLADE*')
fasta = st.sidebar.file_uploader("FASTA file to run Nextclade on")
if fasta is None:
    st.write('**Upload fasta file with sequence**')
    
    
else:
    with open('input.fasta', 'w+') as f:
        print(fasta.read().decode(), file=f)

    nextclade = Nextclade()
    nextclade.run_nextclade()
    nextclade.make_meta_df()
    st.write(nextclade.meta_df)
    #time.sleep(30)
    #shutil.rmtree(nextclade.OUTDIR)
    #os.remove('input.fasta')
    corona = CoronaPlot(nextclade.QUERY, nextclade.REF, 'mutations.csv')


    st.markdown('**COMPARE**')
    compare_to = corona.dataframes.make_mutations_df()
    pango = st.selectbox('Select pango lineage to compare with', compare_to['pango'].unique())
    st.markdown("***")
    compare_plot = corona.plot_compare(pango)
    st.altair_chart(compare_plot)
    st.markdown("***")

    st.markdown('*INSPECTION*')
    inspect = corona.plot_mutations_inspection()
    st.altair_chart(inspect)

    st.markdown('*HEATMAP*')
    heatmap = corona.plot_mutation_heatmap()
    st.altair_chart(heatmap)
    
    st.markdown('**ALIGNMENT**')
    a, b = corona.plot_alignment_map()
    st.altair_chart((a & b))
