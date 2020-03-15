from matplotlib import pyplot as plt
from math import pi
from itertools import accumulate
import operator
import re
from collections import namedtuple
class dummynode:
  def __init__(self,**kwargs):
   self.__dict__=kwargs
  def __lt__(self,other):
    return self.x<other.x
  @property
  def isleaf(self):
    return len(self.children)==0
  @property
  def y(self):
    node=self
    while not node.isleaf: #find the bottommost and topmost leaf belonging to the node
      node=node.children[0]
    minleafid=node.leafid #the leafid is used as the y-position
    node=self
    while not node.isleaf:
      node=node.children[-1]
    maxleafid=node.leafid
    return (maxleafid+minleafid)/2.0 #return the mid-point between the top and bottom leaf
  @property
  def nodes(self):
    yield self
    for child in self.children:
      for node in child.nodes:
        yield node
  @property
  def leaves(self):
    for node in self.nodes:
      if node.isleaf:
        yield node
  @property
  def nonleaves(self):
    for node in self.nodes:
      if not node.isleaf:
        yield node

def lineplot(node):
  yield (node.x,node.y)
  for child in node.children:
    yield(node.x,child.y)
    for xy in lineplot(child):
      yield xy
    yield(node.x,child.y)
    yield (node.x,node.y)
    
def lineplot_polar(node,max_dtheta=pi/100,theta_scale=None):
  if not theta_scale:
    theta_scale=2*pi/sum([1 for a in node.leaves])
  yield (node.y*theta_scale,node.x)
  for child in node.children:
    n=max(1,round(abs(node.y-child.y)*theta_scale/max_dtheta))
    for i in range(n):
      theta_r=((node.y*(n-i-1)+child.y*(i+1))/n*theta_scale,node.x)
      yield theta_r
    for theta_r in lineplot_polar(child,max_dtheta=max_dtheta,theta_scale=theta_scale):
      yield theta_r
    for i in range(n):
      theta_r=((node.y*i+child.y*(n-i))/n*theta_scale,node.x)
      yield theta_r
    yield (node.y*theta_scale,node.x)
                  
def parse(newick): 
#based on stackoverflow answer:
#https://stackoverflow.com/questions/51373300/how-to-convert-newick-tree-format-to-a-tree-like-hierarchical-object/51375562#51375562
#by https://stackoverflow.com/users/5459839/trincot 

#    tokens = re.findall(r"([^:;,()\s]*)(?:\s*:\s*([\d.]+)\s*)?([,);])|(\S)", newick+";")
    tokens = re.findall(r"([^:;,()'\s]*|'[^']*')(?:\s*:\s*([\d.]+)\s*)?([,);])|(\S)", newick+";")
# re syntax "(?: ...)" creates a non-capturing group.
    def recurse(parent=None): # one node
        thisnode=dummynode()
        children = []
        name, length, delim, ch = tokens.pop(0)
        if ch == "(":
            while ch in "(,":
                node, ch= recurse(thisnode)
                children.append(node)
            name, length, delim, ch = tokens.pop(0)
        length=float(length) if length else 0.0
        thisnode.__init__(name=name.strip("'").replace('_',' '), length=length, 
                parent=parent, children=children)
        if parent:
          return thisnode, delim
        else:
          return thisnode
    tree=recurse()
    tree.x=0.0
    maxx=0.0
    for node in tree.nodes:
      if node.parent:
        node.x=node.parent.x+node.length
        maxx=max(maxx,node.x)
    for node in tree.nodes:
      node.x-=maxx
    for i,node in enumerate(tree.leaves):
       node.leafid=i+1
    return tree    
dictEntry=namedtuple('dictEntry','lt en de')    
with open('species_dict.txt','r',encoding='utf-8')as f:
  speciesDict=[dictEntry(latin,english,deutsch) for latin,english,deutsch in re.findall('([^:]*):\s*([^:]*):\s*([^\n]*)\n',f.read())]
  translation={name:e for e in speciesDict for name in e}
  
def translate(name,language):
  t=translation.get(name)
  tl=t._asdict().get(language) if t else None
  return tl if tl else name
  
# Example use:
for filename in ['Aves_species.nwk.txt','Aves_genus.nwk.txt','Aves_family.nwk.txt','Aves_order.nwk.txt']:
  with open(filename,'r') as f:
    tree=parse(f.read())
  sortednodes=sorted(tree.nodes)
  time=[node.x for node in sortednodes]
  delta=[len(node.children)-1 for node in sortednodes]
  species=list(accumulate(delta))
  print('Filename: {0:20s} Number of leaves: {1:4d} (vs. {3:2d} {2:4.1f} million years ago)'.format(filename,sum([1 for node in tree.leaves]),66.0,sum([(node.parent.x<a) and (node.x>a) for node in tree.nodes for a in (-66.043,) if node.parent])))
  plt.plot(time,species,label=filename)

plt.legend(loc='upper left')
plt.xlabel('million years before present')
plt.ylabel('number of species, ...')
plt.ylim(ymin=0.0)
plt.show()
plt.ylim(ymin=1.0)
plt.semilogy()
plt.show()
plt.close()
#print(plt.rcParams['figure.figsize'])
for filename in ['Aves_family.nwk.txt', 'Euteleostomi_family.nwk.txt','myTree.nwk.txt']:
  with open(filename,'r') as f:
    tree=parse(f.read())
  leafcount=len(list(tree.leaves))
  root=dummynode(id=-1,x=tree.x*1.05,name='root',children=[tree])
  plt.plot([a[0] for a in lineplot(root)],[a[1] for a in lineplot(root)],marker=None,label=filename)
  if leafcount<30:
    for node in tree.leaves:
      plt.text(node.x,node.y,' '+translate(node.name,'de'),va='center',ha='left')
  plt.legend(loc='upper left')
  plt.xlim(xmax=0)
  plt.ylim(ymin=0,ymax=leafcount+1)
  plt.xlabel('million years before present')
  plt.ylabel('species')
  plt.show()
  plt.close()
#  plt.figure(figsize=(8,8.7))
  xy=[(a[0],(a[1]-tree.x)) for a in lineplot_polar(tree)]
  plt.polar([a[0] for a in xy],[a[1] for a in xy],marker=None,label=filename)
  if leafcount<30*pi:
    for node in tree.leaves:
      theta=node.y/leafcount*pi*2
      flip=1 if ((theta/(2*pi)+0.25)%1)>0.5 else 0
      plt.text(theta,-tree.x,' '+translate(node.name,'de')+(' .' if flip else''),rotation=(theta/pi+flip)*180,rotation_mode='anchor',va='center',ha=('left','right')[flip])

  plt.legend(loc='upper left')
#  plt.xlim(xmax=0)
  plt.ylim(ymin=0,ymax=-tree.x)
#  plt.ylabel('million years before present')
#  plt.xlabel('species')
  timeticks=plt.yticks()[0]
  plt.yticks([-tree.x-t for t in timeticks[:-1]],['{:.0f}'.format(-t) for t in timeticks[:-1]])
  plt.xticks(plt.xticks()[0],())
  plt.show()
  plt.close()
