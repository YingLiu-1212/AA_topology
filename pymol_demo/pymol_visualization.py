# -*- coding: utf-8 -*-
from pymol import cmd, stored
import os

def hex_to_rgb(hex_color):
    hex_color = hex_color.lstrip('#')
    if len(hex_color) == 3:
        hex_color = ''.join([c*2 for c in hex_color])
    r = int(hex_color[0:2], 16) / 255.0
    g = int(hex_color[2:4], 16) / 255.0
    b = int(hex_color[4:6], 16) / 255.0
    return [r, g, b]
    
cmd.set_color("RdBu_1", hex_to_rgb("#67001F")) 
cmd.set_color("RdBu_2", hex_to_rgb("#B2182B")) 
cmd.set_color("RdBu_3", hex_to_rgb("#D6604D")) 
cmd.set_color("RdBu_4", hex_to_rgb("#F4A582")) 
cmd.set_color("RdBu_5", hex_to_rgb("#FDDBC7")) 
cmd.set_color("RdBu_6", hex_to_rgb("#F7F7F7")) 
cmd.set_color("RdBu_7", hex_to_rgb("#D1E5F0")) 
cmd.set_color("RdBu_8", hex_to_rgb("#92C5DE")) 
cmd.set_color("RdBu_9", hex_to_rgb("#4393C3")) 
cmd.set_color("RdBu_10", hex_to_rgb("#2166AC")) 
cmd.set_color("RdBu_11", hex_to_rgb("#053061")) 
 
RdBu = "RdBu_1 RdBu_2 RdBu_3 RdBu_4 RdBu_5 RdBu_6 RdBu_7 RdBu_8 RdBu_9 RdBu_10 RdBu_11"

column_index = {
    'pdb_chain': 0,
    'pdb_id': 1,
    'chain': 2,
    'chain_len': 3,
    'id': 4,
    'resi': 5,
    'AA': 6,
    'CF10': 7,
    'CF10QS': 8,
    'LD15': 9,
    'LD15QS': 10,
    'CF10QS_FPI': 11
}

# ===================== User-configurable parameters =====================

gene_name = 'TP53'
pdb_id = '1ycs'
chain_id = 'A'
plot_var = 'CF10QS'

project_dir = 'D:/Git/AA_topology/pymol_demo/'

MIN_VALUE = -1
MAX_VALUE = 1  

PDB_FILE = project_dir + pdb_id + ".pdb"            
TOPO_FILE = project_dir +"AA_TopoAttr_" + pdb_id + ".txt"        
TARGET_CHAIN = chain_id               
COLOR_SCHEME = RdBu      
OUTPUT_IMAGE = project_dir + pdb_id + "_" + plot_var + "_" + chain_id + ".png"   
# ========================================================

cmd.load(PDB_FILE)
cmd.hide("everything")

cmd.alter("all", "use_value=None")

with open(TOPO_FILE) as f:
    for line_num, line in enumerate(f):

        if line_num == 0:
            continue 
            
        line = line.strip()  
        
        line = line.strip()     
        parts = line.split('\t')
        
        chain = parts[column_index['chain']]
        resi = parts[column_index['resi']]
        use_val = float(parts[column_index[plot_var]])
        selection = f"chain {chain} and resi {resi}"
        
        cmd.alter(selection, f"use_value = {use_val}; b = {use_val}")
        #print(f"setting {chain} {resi}: CF={use_val}")


target_sel = f"chain {TARGET_CHAIN}"

cmd.spectrum("b", COLOR_SCHEME, target_sel, minimum=MIN_VALUE, maximum=MAX_VALUE)

cmd.hide("everything")
cmd.show("cartoon", target_sel)
cmd.set("cartoon_loop_radius", 0.5)  

cmd.set("cartoon_smooth_loops", 1)  
cmd.set("cartoon_flat_sheets", 0)  
cmd.set("cartoon_tube_radius", 0.5)  
cmd.set("cartoon_oval_length", 1.2)  
cmd.set("cartoon_oval_width", 0.4)  
cmd.set("cartoon_sampling", 2)  
cmd.set("light_count",6)
cmd.set("specular",0)
cmd.set("ambient",0.5)


cmd.bg_color("white")
cmd.set("light_count", 6)
cmd.zoom(target_sel, buffer=10)

cmd.refresh()
cmd.rebuild()
cmd.ray(2400, 2400)
cmd.png(OUTPUT_IMAGE)

