# Read a GasLib .net-file and construct corresponding tikz code
# based on code by Rico Raber

from __future__ import division
import sys
import xml.etree.ElementTree as ET
import networkx as nx

def printFunc(tikz_string,edge,lowerbounds,upperbounds,width,color,zerocolor,valve,outputname):
    fromname=edge[0].replace("#","_");
    toname=edge[1].replace("#","_");
    id = edge[2]
    printname=id.replace("_","")
    lb = float(lowerbounds[id])
    ub = float(upperbounds[id])
    if lb >= 0.0:
        if lb == ub:
            if lb != 0.0:
                if outputname:
                    tikz_string.append("\\draw[{}, {}, >->] ({}) -- ({}) node[pos=0.5,font=\\fontsize{{5}}{{5}}\\selectfont] {{ {} }};".format(width, color, fromname, toname, printname))
                else:
                    tikz_string.append("\\draw[{}, {}, >->] ({}) -- ({});".format(width, color, fromname, toname))
            elif not valve:
                if outputname:
                    tikz_string.append("\\draw[{}, {}, >->] ({}) -- ({}) node[pos=0.5,font=\\fontsize{{5}}{{5}}\\selectfont] {{ {} }};".format(width, zerocolor, fromname, toname, printname))
                else:
                    tikz_string.append("\\draw[{}, {}, >->] ({}) -- ({});".format(width, zerocolor, fromname, toname))
        else:
            if outputname:
                tikz_string.append("\\draw[{}, {}, ->] ({}) -- ({}) node[pos=0.5,font=\\fontsize{{5}}{{5}}\\selectfont] {{ {} }};".format(width, color, fromname, toname, printname))       
            else:
                tikz_string.append("\\draw[{}, {}, ->] ({}) -- ({});".format(width, color, fromname, toname))
    elif ub <= 0.0:
        if lb == ub:
            if lb != 0.0:
                if outputname:
                    tikz_string.append("\\draw[{}, {}, >->] ({}) -- ({}) node[pos=0.5,font=\\fontsize{{5}}{{5}}\\selectfont] {{ {} }};".format(width, color, toname, fromname, printname))
                else:
                    tikz_string.append("\\draw[{}, {}, >->] ({}) -- ({});".format(width, color, toname, fromname))
            elif not valve:
                if outputname:
                    tikz_string.append("\\draw[{}, {}, >->] ({}) -- ({}) node[pos=0.5,font=\\fontsize{{5}}{{5}}\\selectfont] {{ {} }};".format(width, zerocolor, toname, fromname, printname))
                else:
                    tikz_string.append("\\draw[{}, {}, >->] ({}) -- ({});".format(width, zerocolor, toname, fromname))
        else:
            if outputname:
                tikz_string.append("\\draw[{}, {}, ->] ({}) -- ({}) node[pos=0.5,font=\\fontsize{{5}}{{5}}\\selectfont] {{ {} }};".format(width, color, toname, fromname, printname))
            else:
                tikz_string.append("\\draw[{}, {}, ->] ({}) -- ({});".format(width, color, toname, fromname))
    else:
        if outputname:
            tikz_string.append("\\draw[{}, {}, dotted] ({}) -- ({}) node[pos=0.5,font=\\fontsize{{5}}{{5}}\\selectfont] {{ {} }};".format(width, color, fromname, toname, printname))
        else:
            tikz_string.append("\\draw[{}, {}, dotted] ({}) -- ({});".format(width, color, fromname, toname))

    return tikz_string


# ------------------------------------------------------------------------------
# decide whether we scale to a small picture
scalesmall = False

# read network and outfile
if len(sys.argv) == 3:
    path = sys.argv[1]
    outfile = sys.argv[2]
else:
    print("usage <net file> <out-file>")
    exit()

filename = path.split('/')[-1][:-4].strip()
print("Parsing file <%s>" % path)
tree = ET.parse(path)
root = tree.getroot()

# Namespaces: first map moves all tags without prefix to
# "http://gaslib.zib.de/Gas" namespace. The second introduces a map to "gas".
namespaces = {'': 'http://gaslib.zib.de/Gas', 'gas': 'http://gaslib.zib.de/Gas'}

nodes = {'gas:source': [], 'gas:sink': [], 'gas:innode': []}
edges = {'gas:pipe':[], 'gas:compressorStation':[], 'gas:controlValve':[], 'gas:resistor':[], 'gas:valve':[], 'gas:shortPipe':[]}

# read all nodes from xml file
for nodetype in nodes:
    for node in root[1].findall(nodetype, namespaces):
        nodes[nodetype].append({'id': node.get('id'), 'x': float(node.get('x')), 'y': float(node.get('y'))})

# read all edges from xml file
for edgetype in edges:
    for edge in root[2].findall(edgetype, namespaces):
        e = (edge.get('from'), edge.get('to'), edge.get('id'))
        edges[edgetype].append(e)

# ---------------------------------------------------------------------
# read out file
f = open(outfile, "r")

lowerbounds = {}
upperbounds = {}
line = f.readline()
while line:
    if line.startswith('flow bounds on arc'):
        data = line.split()
        pipename = data[4]
        lb = data[6]
        ub = data[8]
        lowerbounds[pipename] = lb
        upperbounds[pipename] = ub
    line = f.readline()

f.close()


# ---------------------------------------------------------------------
# determine bounding box
x_max = max([node['x'] for nodetype in nodes for node in nodes[nodetype]])
x_min = min([node['x'] for nodetype in nodes for node in nodes[nodetype]])
x_span = x_max - x_min

y_max = max([node['y'] for nodetype in nodes for node in nodes[nodetype]])
y_min = min([node['y'] for nodetype in nodes for node in nodes[nodetype]])
y_span = y_max - y_min

# default settings
if scalesmall:
    x_scale = 10
    y_scale = 7
else:
    x_scale = 163
    y_scale = 114

if scalesmall:
    inner_sep_source = "0.4pt"
    inner_sep_sink= "0.2pt"
    inner_sep_innode= "0.2pt"
else:
    inner_sep_source = "1.5pt"
    inner_sep_sink= "1pt"
    inner_sep_innode= "0.5pt"

if scalesmall:
    pipe_width = "thin"
    compressor_width = "thin"
    control_valve_width = "thin"
    resistor_width = "thin"
    valve_width = "thin"
    short_pipe_width = "thin"
else:
    pipe_width = "thick"
    compressor_width = "thick"
    control_valve_width = "thick"
    resistor_width = "thick"
    valve_width = "thick"
    short_pipe_width = "thick"

# define colors for nodes and edges
source_fill_color = "blue"
sink_fill_color = "white!80!blue"
innode_fill_color = "black"

pipe_color = "black"
zero_pipe_color = "black!10!white"
compressor_color = "red"
zero_compressor_color = "red!30!white"
control_valve_color = "orange"
zero_control_valve_color = "orange!30!white"
resistor_color = "cyan"
zero_resistor_color = "cyan!30!white"
valve_color = "black!30!green"
zero_valve_color = "black!30!green!30!white"
short_pipe_color = "brown"
zero_short_pipe_color = "brown!30!white"

# scale coordinates
for nodetype in nodes:
    for node in nodes[nodetype]:
        x = node['x']
        x *= x_scale/x_span
        x = int(x * 10000) / 10000
        node['x'] = x
        y = node['y']
        y *= y_scale/y_span
        y = int(y * 10000) / 10000
        node['y'] = y

# start generating tex files
# lualatex might be needed for larger files; fix bug with the first package
tikz_string = ["\\RequirePackage{luatex85}"]
tikz_string.append("\\documentclass[border=3pt]{standalone}")
tikz_string.append("\\usepackage{lmodern}")
tikz_string.append("\\usepackage[T1]{fontenc}")
tikz_string.append("\\usepackage{tikz}")
tikz_string.append("\\usepgflibrary{fpu}")
tikz_string.append("\\usetikzlibrary{intersections,patterns}")
tikz_string.append("\\usetikzlibrary{circuits}")
tikz_string.append("\\usetikzlibrary{shapes}")
tikz_string.append("\\usetikzlibrary{arrows.meta}")
tikz_string.append("\\pgfdeclarelayer{bg}")
tikz_string.append("\\pgfsetlayers{bg,main}\n")

tikz_string.append("\\begin{document}")

tikz_string.append("\\begin{tikzpicture}[circuit,innode/.style={{shape=circle,fill=black,scale=0.6}},sink/.style={{draw,shape=circle,line width=2pt}},source/.style={{draw,shape=circle,font={{S}},line width=2pt}}]")

# draw nodes
for node in nodes['gas:innode']:
    name=node['id'].replace("#","\\_")
    tikz_string.append("\\node[fill = {}, circle, inner sep = {}, minimum size = 0pt] ({}) at ({},{}) {{}};".format(innode_fill_color, inner_sep_innode, name, node['x'], node['y']))
for node in nodes['gas:sink']:
    name=node['id'].replace("#","_")
    tikz_string.append("\\node[fill = {}, circle, inner sep = {}, minimum size = 0pt] ({}) at ({},{}) {{}};".format(sink_fill_color, inner_sep_sink, name, node['x'], node['y']))
for node in nodes['gas:source']:
    name=node['id'].replace("#","_")
    tikz_string.append("\\node[fill = {},  circle, inner sep = {}, minimum size = 0pt] ({}) at ({},{}) {{}};".format(source_fill_color, inner_sep_source, name, node['x'], node['y']))

# draw edges
tikz_string.append("\\begin{pgfonlayer}{bg}")
for edge in edges['gas:pipe']:
    tikz_string = printFunc(tikz_string, edge, lowerbounds, upperbounds, pipe_width, pipe_color, zero_pipe_color, False, False)
for edge in edges['gas:shortPipe']:
    tikz_string = printFunc(tikz_string, edge, lowerbounds, upperbounds, short_pipe_width, short_pipe_color, zero_short_pipe_color, False, False)
for edge in edges['gas:compressorStation']:
    tikz_string = printFunc(tikz_string, edge, lowerbounds, upperbounds, compressor_width, compressor_color, zero_compressor_color, False, False)
for edge in edges['gas:controlValve']:
    tikz_string = printFunc(tikz_string, edge, lowerbounds, upperbounds, control_valve_width, control_valve_color, zero_control_valve_color, False, True)
for edge in edges['gas:resistor']:
    tikz_string = printFunc(tikz_string, edge, lowerbounds, upperbounds, resistor_width, resistor_color, zero_resistor_color, False, False)
for edge in edges['gas:valve']:
    tikz_string = printFunc(tikz_string, edge, lowerbounds, upperbounds, valve_width, valve_color, zero_valve_color, True, False)

tikz_string.append("\\end{pgfonlayer}")
tikz_string.append("\\end{tikzpicture}")

tikz_string.append("\end{document}\n")

# create string and write to file
tikz_string = "\n".join(tikz_string)
with open(filename + ".tex", "w") as text_file:
    text_file.write(tikz_string)

print("Writing tex file to <%s>" % (filename + ".tex"))
