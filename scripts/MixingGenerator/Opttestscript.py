from __future__ import division
import sys
import xml.etree.ElementTree as ET
import networkx as nx
import math

# opttest datei einlesen und die gewünschten dateien in ein array schreiben
opt_test_files=[]
with open('optimal.test', 'r') as file:
    # Lies eine Zeile
    line = file.readline()
    while line:
        # Verarbeite die Zeile
        # Lies die nächste Zeile
        abschnitte = line.rsplit('/',1)
        temp=abschnitte[1][:-1]
        opt_test_files.append(temp)
        line = file.readline()


#neuen nodetag bauen
node_tag = ET.Element("node")
node_tag.set('type', 'gasProperties')

molar1_tag = ET.Element("molarMass1")
molar1_tag.set('value', '18.5674')
molar1_tag.set('unit', 'kg_per_kmol')

molar2_tag = ET.Element("molarMass2")
molar2_tag.set('value', '2.016')
molar2_tag.set('unit', 'kg_per_kmol')
#- Append Element to Parent Element.
node_tag.append(molar1_tag)
node_tag.append(molar2_tag)


# Für jeden path in opt_test werden die folgenden Schritte durchgeführt
for path in opt_test_files:
    print(path)
    tree = ET.parse(path)
    root = tree.getroot()
    #Der neue Knoten wird hinzugefügt
    root[0].insert(0, node_tag)
    #Die flow values werden ausgelesen und die beiden größten in einem array fixiert
    help_array = []
    for index in range(2,len(root[0])):
        if root[0][index].get('type') == 'entry':
            help_array.append(float(root[0][index][-1].get('value')))
    largest_values = sorted(help_array, reverse=True)[:2]
    #Die molaren Massen werden hinzugefügt
    for index in range(2,len(root[0])):
        if root[0][index].get('type') == 'entry':
            if float(root[0][index][-1].get('value')) == largest_values[0] or float(root[0][index][-1].get('value')) == largest_values[1]:
                molar_tag = ET.Element('molarMass')
                molar_tag.set('value', '16.9')
                molar_tag.set('unit', 'kg_per_kmol')
                root[0][index].append(molar_tag)
            else:
                molar_tag = ET.Element('molarMass')
                molar_tag.set('value', '18.5674')
                molar_tag.set('unit', 'kg_per_kmol')
                root[0][index].append(molar_tag)
    #Der modifizierte Tree wird in eine xml Datei geschrieben
    subfolder_path = 'modifiziert/'
    new_tree = ET.ElementTree(root)
    new_root = new_tree.getroot()
    for elem in new_root.iter():
        elem.tag = elem.tag.split('}')[-1]
    new_tree.write(subfolder_path + path)






