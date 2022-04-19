
import os
import streamlit as st
import numpy as np
from PIL import  Image

from multipage import MultiPage
from pages import upload, compare, inspect, alignment, info, heatmap

app = MultiPage()

# Main page
img = 'corona.png'
display = Image.open(img)
col1, col2 = st.columns(2)
col1.image(display, width = 200)
col2.title('Corona viz')

app.add_page('Upload Data', upload.app)
app.add_page('Compare against corona lineages', compare.app)
app.add_page('Inspect mutations', inspect.app)
app.add_page('Alignment information', alignment.app)
app.add_page('Sequence information', info.app)
app.add_page('Heatmap', heatmap.app)


app.run()


