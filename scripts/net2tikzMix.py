# Read a GasLib .net-file and construct corresponding tikz code
# based on code by Rico Raber
# editded by Linus Strauch and Pascal Börner

from __future__ import division
import sys
import xml.etree.ElementTree as ET
import networkx as nx
import math

if len(sys.argv) == 3:
    path = sys.argv[1]
    path_lsg = sys.argv[2]
    mol_or_mass = 'mass'
    untergrenze_kantenfluss_in_prozent = 0.01
    print("Default: mass percentage used")
    print("Default: arcs with more than 0.01 percent of max flow are drwan")
elif len(sys.argv) == 4:
    path = sys.argv[1]
    path_lsg = sys.argv[2]
    mol_or_mass = sys.argv[3]
elif len(sys.argv) == 5:
    path = sys.argv[1]
    path_lsg = sys.argv[2]
    mol_or_mass = sys.argv[3]
    untergrenze_kantenfluss_in_prozent = float(sys.argv[4])
else:
    print("Error: Need arguments: net file | solution file | mol or mass % | lb flow % ")
    sys.exit(1)
if mol_or_mass != 'Mol' and mol_or_mass != 'mol' and mol_or_mass != 'mass' and mol_or_mass != 'mass':
    print('Please use correct classification for flow percentage')
    exit()



filename = path.split('/')[-1][:-4].strip()
print("Parsing file <%s>" % path)
tree = ET.parse(path)
root = tree.getroot()

filename_lsg = path_lsg.split('/')[-1][:-4].strip()
print("Parsing file <%s>" % path_lsg)
tree_lsg = ET.parse(path_lsg)
root_lsg = tree_lsg.getroot()

#Introducing functions to calculate the color of nodes and edges


def get_node_mixingratio(node_id):
    for entry in root_lsg[0][0]:
        if node_id == entry.get('id'):
            return float(entry[-1].get('value'))
    for entry in root_lsg[0][1]:
        if node_id == entry.get('id'):
            return float(entry[-1].get('value'))
    for entry in root_lsg[0][2]:
        if node_id == entry.get('id'):
            return float(entry[-1].get('value'))


# in both functions below the position of the mole / mass percentage might need to be adapted -2 is correct if only mass percentage is written in the solu file
# and the entry in the solu file is lambda
def get_edge_mixingratio(node_from, node_to):
    for pipe in root_lsg[1][0]:
        if node_from==pipe.get('from') and node_to==pipe.get('to'):
            if mol_or_mass == 'mass':
                return float(pipe[-2].get('value'))     #Aktivieren für Mol
            else:
                return float( pipe[-1].get('value'))           #Aktivieren für Mass
        
def get_cs_mixingratio(node_from, node_to):
    for pipe in root_lsg[1][2]:
        if node_from==pipe.get('from') and node_to==pipe.get('to'):
            if mol_or_mass == 'mass':
                return float(pipe[-2].get('value'))     #Aktivieren für Mol
            else:
                return float( pipe[-1].get('value'))     #Aktivieren für Mass

def value_to_color(value):
    gelb = round(value * 255)
    if gelb >= 128:
        rot = gelb - 2*(255-gelb)
        green = 255
        blau = 0
    else:
        green = 2 * gelb
        blau = 255-green
        rot = 0

    return '{rgb,255:red,'+ str(int(rot))+';green,' + str(int(green)) +';blue,'+str(int(blau))+'}'

def get_flow(node_from, node_to):
    for pipe in root_lsg[1][0]:
        if node_from==pipe.get('from') and node_to==pipe.get('to'):
            return float(pipe.find('flow').get('value'))

def get_edge_direction(node_from, node_to):
    for pipe in root_lsg[1][0]:
        if node_from==pipe.get('from') and node_to==pipe.get('to'):
            if int(pipe.find('negFlowBinvar').get('value')):
                return('<-')
            else:
                return('->')

# In the following we will introduce functions to calculate which edges and nodes are having more then 0 flow. More than zero means more then untergrenze_kantenfluss_in_prozent of the total flow
sum_flow = 0
for source in root_lsg[0][0]:
    sum_flow = sum_flow + float(source.find('flow').get('value'))

def get_edge_percentage_of_flow(node_from, node_to):
    for pipe in root_lsg[1][0]:
        if node_from == pipe.get('from') and node_to == pipe.get('to'):
            
            percentage_of_flow = float(pipe.find('flow').get('value')) / sum_flow
            return percentage_of_flow

def get_node_percentage_of_flow(node_id):
    for entry in root_lsg[0][0]:
        if node_id == entry.get('id'):
            percentage_of_flow = float(entry.find('flow').get('value')) / sum_flow
            return percentage_of_flow
    for entry in root_lsg[0][1]:
        if node_id == entry.get('id'):
            percentage_of_flow = float(entry.find('flow').get('value')) / sum_flow
            return percentage_of_flow
    for entry in root_lsg[0][2]:
        if node_id == entry.get('id'):
            percentage_of_flow = float(entry.find('flow').get('value')) / sum_flow
            return percentage_of_flow







# Namespaces: first map moves all tags without prefix to
# "http://gaslib.zib.de/Gas" namespace. The second introduces a map to "gas".
namespaces = {'': 'http://gaslib.zib.de/Gas', 'gas': 'http://gaslib.zib.de/Gas'}

nodes = {'gas:source': [], 'gas:sink': [], 'gas:innode': []}
edges = {'gas:pipe':[], 'gas:compressorStation':[], 'gas:controlValve':[], 'gas:resistor':[], 'gas:valve':[], 'gas:shortPipe':[]}

node_count = 0
# read all nodes from xml file
for nodetype in nodes:
    for node in root[1].findall(nodetype, namespaces):
        nodes[nodetype].append({'id': node.get('id'), 'x': float(node.get('x')), 'y': float(node.get('y'))})
        node_count = node_count + 1

# read all edges from xml file
for edgetype in edges:
    for edge in root[2].findall(edgetype, namespaces):
        e = (edge.get('from'), edge.get('to'))
        edges[edgetype].append(e)

#Try calculating max flow on arc
        q=0 #initialize flow as zero
for edge in edges['gas:pipe']:
    fromname=edge[0].replace("#","_");
    toname=edge[1].replace("#","_");    
    if abs (get_flow(fromname, toname)) > q:
        q=abs (get_flow(fromname, toname))

def lineWidth(node_from, node_to):
    flow = abs(get_flow(node_from, node_to))
    width= 1.5*flow /q
    if width < 0.001:
        width = 0.0
    return 'line width='+str(width)+'mm'




# determine bounding box
x_max = max([node['x'] for nodetype in nodes for node in nodes[nodetype]])
x_min = min([node['x'] for nodetype in nodes for node in nodes[nodetype]])
x_span = x_max - x_min

y_max = max([node['y'] for nodetype in nodes for node in nodes[nodetype]])
y_min = min([node['y'] for nodetype in nodes for node in nodes[nodetype]])
y_span = y_max - y_min

# default settings
x_scale = max(2 * math.sqrt(node_count), 10)
y_scale = max(7/5 * math.sqrt(node_count), 7)
inner_sep_source = "1.5pt"
inner_sep_sink= "1pt"
inner_sep_innode= "1.5pt"

pipe_width = "ultra thick"
compressor_width = "thick"
control_valve_width = "thin"
resistor_width = "thin"
valve_width = "thin"
short_pipe_width = "thin"

# define colors for nodes and edges
source_fill_color = "blue"
sink_fill_color = "red"
innode_fill_color = "black"

pipe_color = "black"
compressor_color = "cyan"
control_valve_color = "orange"
resistor_color = "cyan"
valve_color = "black!30!green"
short_pipe_color = "black"

# scale coordinates
for nodetype in nodes:
    for node in nodes[nodetype]:
        node['x'] *= x_scale/x_span
        node['y'] *= y_scale/y_span

# Finding all the nodes that need to be drawn
wichtige_knoten = set()
for edge in edges['gas:pipe']:
    fromname = edge[0].replace("#", "_");
    toname = edge[1].replace("#", "_");
    if abs(get_edge_percentage_of_flow(fromname, toname)) > untergrenze_kantenfluss_in_prozent:
        wichtige_knoten.add(fromname)
        wichtige_knoten.add(toname)
for edge in edges['gas:compressorStation'] + edges['gas:controlValve'] + edges['gas:resistor'] + edges['gas:valve'] + edges['gas:shortPipe']:
    fromname = edge[0].replace("#", "_");
    toname = edge[1].replace("#", "_");
    wichtige_knoten.add(fromname)
    wichtige_knoten.add(toname)

# start generating tex files
tikz_string = ["\\documentclass[border=3pt]{standalone}"]
tikz_string.append("\\usepackage[T1]{fontenc}")
tikz_string.append("\\usepackage{tikz}")
tikz_string.append("\\usetikzlibrary{intersections,patterns}")
tikz_string.append("\\usetikzlibrary{circuits}")
tikz_string.append("\\usetikzlibrary{shapes}")
tikz_string.append("\\pgfdeclarelayer{bg}")
tikz_string.append("\\pgfsetlayers{bg,main}\n")

tikz_string.append("\\begin{document}")

tikz_string.append("\\begin{tikzpicture}[circuit,innode/.style={{shape=circle,fill=black,scale=0.6}},sink/.style={{draw,shape=circle,line width=2pt}},source/.style={{draw,shape=circle,font={{S}},line width=2pt}}]")

# draw nodes
for node in nodes['gas:innode']:
    name=node['id'].replace("#","_")
    if name in wichtige_knoten:
        tikz_string.append("\\node[fill = {}, circle, inner sep = {}, minimum size = 0pt] ({}) at ({},{}) {{}};".format(value_to_color(get_node_mixingratio(name)), inner_sep_innode, name, node['x'],node['y']))
for node in nodes['gas:sink']:
    name=node['id'].replace("#","_")
    if name in wichtige_knoten:
        tikz_string.append("\\node[fill = {}, star, star points=5, inner sep = {}, minimum size = 0pt] ({}) at ({},{}) {{}};".format('red', inner_sep_sink, name, node['x'],node['y']))
for node in nodes['gas:source']:
    name=node['id'].replace("#","_")
    if name in wichtige_knoten:
        tikz_string.append("\\node[fill = {},  diamond, inner sep = {}, minimum size = 0pt] ({}) at ({},{}) {{}};".format('green', inner_sep_source, name, node['x'],node['y']))

# draw edges
tikz_string.append("\\begin{pgfonlayer}{bg}")
for edge in edges['gas:pipe']:
    fromname=edge[0].replace("#","_");
    toname=edge[1].replace("#","_");
    if abs(get_edge_percentage_of_flow(fromname, toname)) < untergrenze_kantenfluss_in_prozent:
        pass
    elif get_edge_mixingratio(fromname, toname) > 0.99 or get_edge_mixingratio(fromname, toname) < 0.01:
        tikz_string.append("\\draw[{}, {}, color ={}] ({}) -- ({});".format(get_edge_direction(fromname, toname), lineWidth(fromname,toname), (value_to_color(get_edge_mixingratio(fromname, toname))), fromname, toname))
    else:
        tikz_string.append("\\draw[{}, {}, color ={}] ({}) -- node[sloped,font=\\small,above] {{{}}}  ({});".format(get_edge_direction(fromname, toname),
                                                                            lineWidth(fromname, toname), (
                                                                                value_to_color(
                                                                                    get_edge_mixingratio(fromname,
                                                                                                         toname))),
                                                                            fromname, round(get_edge_mixingratio(fromname, toname), 2), toname))

for edge in edges['gas:compressorStation']:
    fromname=edge[0].replace("#","_");
    toname=edge[1].replace("#","_");
    tikz_string.append("\\draw[{}, color ={}] ({}) -- ({});".format(compressor_width, compressor_color, fromname,toname)) #hier geht auch getcolor aber muss noch in der c datei angepasst werden
for edge in edges['gas:controlValve']:
    fromname=edge[0].replace("#","_");
    toname=edge[1].replace("#","_");
    tikz_string.append("\\draw[{}, {}] ({}) -- ({});".format(control_valve_width, control_valve_color, fromname,toname))
for edge in edges['gas:resistor']:
    fromname=edge[0].replace("#","_");
    toname=edge[1].replace("#","_");
    tikz_string.append("\\draw[{}, {}] ({}) -- ({});".format(resistor_width, resistor_color, fromname,toname))
for edge in edges['gas:valve']:
    fromname=edge[0].replace("#","_");
    toname=edge[1].replace("#","_");
    tikz_string.append("\\draw[{}, {}] ({}) -- ({});".format(valve_width, valve_color, fromname,toname))
for edge in edges['gas:shortPipe']:
    fromname=edge[0].replace("#","_");
    toname=edge[1].replace("#","_");
    tikz_string.append("\\draw[{}, {}] ({}) -- ({});".format(short_pipe_width, short_pipe_color, fromname,toname))

tikz_string.append("\\end{pgfonlayer}")
tikz_string.append("\\end{tikzpicture}")

tikz_string.append("\\end{document}\n")

# create string and write to file
tikz_string = "\n".join(tikz_string)
with open(filename+filename_lsg + ".tex", "w") as text_file:
    text_file.write(tikz_string)

print("Writing tex file to <%s>" % (filename+filename_lsg + ".tex"))
