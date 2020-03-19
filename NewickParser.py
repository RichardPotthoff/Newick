from matplotlib import pyplot as plt
from math import pi,sin,cos,tan,atan2
import math
from itertools import accumulate
import operator
import re
from collections import namedtuple
class Tree:
  def __init__(self,newick=None,parent=None,name=None,length=math.nan,x=math.nan,children=None,id=-1,parent_id=-1):
    self.children=children if children else []
    self.parent=parent
    self.id=id
    self.parent_id=parent_id
    self.leafid=-1
    self.name=name
    self.length=length
    self.x=x
#   self.__dict__=kwargs
    if newick:
      self.parse(newick)
  def __lt__(self,other):
    return self.x<other.x
  def __repr__(self):
#    args=[]
#    if self.name: args.append('name={:s}'.format(self.name.__repr__()))
#    if not math.isnan(self.length): args.append('length={:s}'.format(self.length.__repr__()))
#    if len(self.children)!=0: args.append('children={:s}'.format(self.children.__repr__()))
#    return "Tree({:s})".format(', '.join(args))
    return 'Tree(newick="{:s}")'.format(self.newick)
  @property
  def newick(self):
    s=''
    if len(self.children)!=0: s+='({:s})'.format(','.join([c.newick for c in self.children]))
    if self.name: s+=self.name.__repr__() if re.search("[:;,()'\"]",self.name) else self.name.replace(' ','_')
    if not math.isnan(self.length): s+=':{:s}'.format(self.length.__repr__())
    if not self.parent: s+=';'
    return s
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

  def lineplot(self):
    yield (self.x,self.y)
    for child in self.children:
      yield(self.x,child.y)
      for xy in child.lineplot():
        yield xy
      yield(self.x,child.y)
      yield (self.x,self.y)
    
  def lineplot_polar(self,max_dtheta=pi/100,theta_scale=None):
    if not theta_scale:
      theta_scale=2*pi/sum([1 for a in self.leaves])
    yield (self.y*theta_scale,self.x)
    for child in self.children:
      n=max(1,round(abs(self.y-child.y)*theta_scale/max_dtheta))
      for i in range(n):
        theta_r=((self.y*(n-i-1)+child.y*(i+1))/n*theta_scale,self.x)
        yield theta_r
      for theta_r in child.lineplot_polar(max_dtheta=max_dtheta,theta_scale=theta_scale):
        yield theta_r
      for i in range(n):
        theta_r=((self.y*i+child.y*(n-i))/n*theta_scale,self.x)
        yield theta_r
      yield (self.y*theta_scale,self.x)
                  
  def parse(self,newick): 
  #based on stackoverflow answer:
  #https://stackoverflow.com/questions/51373300/how-to-convert-newick-tree-format-to-a-tree-like-hierarchical-object/51375562#51375562
  #by https://stackoverflow.com/users/5459839/trincot 
  
  #    tokens = re.findall(r"""([^:;,()\s]*)(?:\s*:\s*([\d.]+)\s*)?([,);])|(\S)", newick+";")
    tokens = re.findall(r"""
    ([^:;,()'"\s]*|'[^']*'|"[^"]*")                              #group1: name or 'name'(possibly empty)
    (?:\s*:\s*                                                   #group2: length (sepatator ":")
    ([+-]?(?:\d+(?:[.]\d*)?(?:e[+-]?\d+)?|[.]\d+(?:e[+-]?\d+)?)) #group2: floating point number (captured)
    \s*)?                                                        #group2 is optional
    ([,);])|                                                     #group3: delimiter or
    (\S)                                                         #group4: other character, e.g. "("
    """, newick+";",re.VERBOSE)
# re syntax "(?: ...)" creates a non-capturing group.
    def recurse(thisnode): # one node
      name, length, delim, ch = tokens.pop(0)
      if ch == "(":
          while ch in "(,":
              node, ch= recurse(Tree(parent=thisnode))
              thisnode.children.append(node)
          name, length, delim, ch = tokens.pop(0)
      length=float(length) if length else math.nan
      thisnode.name=name.strip("'").replace('_',' ')
      thisnode.length=length
      return thisnode, delim
    recurse(self)
    self.x=0.0
    maxx=0.0
    for id,node in enumerate(self.nodes):
      node.id=id
      if  node.parent:
        node.x=node.parent.x+node.length
        maxx=max(maxx,node.x)
    for node in self.nodes:
      node.x-=maxx
    for i,node in enumerate(self.leaves):
       node.leafid=i+1    

def limitslope(alpha,limit=pi/4):
  if abs(limit)>(pi/2-1e-6):
    return alpha
  s=sin(alpha)
  c=cos(alpha)
  if abs(limit)<1e-6:
    return 0 if c>0 else pi
  t=tan(abs(limit))
  return atan2(s,c+(1/t if c>0 else -1/t))
  
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
    tree=Tree(f.read())
  sortednodes=sorted(tree.nodes)
  time=[node.x for node in sortednodes]
  delta=[len(node.children)-1 for node in sortednodes]
  species=list(accumulate(delta))
#  print('Filename: {0:20s} Number of leaves: {1:4d} (vs. {3:2d} {2:4.1f} million years ago)'.format(filename,sum([1 for node in tree.leaves]),66.0,sum([(node.parent.x<a) and (node.x>a) for node in tree.nodes for a in (-66.043,) if node.parent])))
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
  plt.rcdefaults()
  plt.figure()
  with open(filename,'r') as f:
    tree=Tree(f.read())
  leafcount=len(list(tree.leaves))
  root=Tree(x=1.05*tree.x,children=[tree])
  plt.plot([a[0] for a in root.lineplot()],[a[1] for a in root.lineplot()],marker=None,label=filename)
  if leafcount<30:
    for node in tree.leaves:
      plt.text(node.x,node.y,' '+translate(node.name,'de'),va='center',ha='left',size='small')
  plt.legend(loc='upper left')
  plt.xlim(xmax=0)
  plt.ylim(ymin=0,ymax=leafcount+1)
  plt.xlabel('million years before present')
  plt.ylabel('species')
  plt.show()
  plt.close()
  plt.rcdefaults()
  plt.rc('ytick', labelsize='small')
  fig=plt.figure()#figsize=(8,8.7))
  dtheta=6.7*2*pi/leafcount
  xy=[(a[0]+dtheta,(a[1]-tree.x)) for a in tree.lineplot_polar()]
  plt.polar([a[0] for a in xy],[a[1] for a in xy],marker=None,label=filename)
  if leafcount<30*pi:
    for node in tree.leaves:
      theta=node.y/leafcount*pi*2+dtheta
      flip=1 if ((theta/(2*pi)+0.25)%1)>0.5 else 0
      plt.text(theta,-tree.x,' '+translate(node.name,'de')+(' .' if flip else''),rotation=(limitslope(theta,45*pi/180)/pi+flip)*180,rotation_mode='anchor',va='center',ha=('left','right')[flip],fontsize='small') 
  else:
    plt.legend(loc='upper left')
#  plt.xlim(xmax=0)
  plt.ylim(ymin=0,ymax=-tree.x)
#  plt.ylabel('million years before present')
#  plt.xlabel('species')
  timeticks=plt.yticks()[0]
  #set_rlabel_position(85)
  fig.axes[0].set_rlabel_position(90)
  fig.axes[0].tick_params(labelleft=False, labelright=True,
               labeltop=False, labelbottom=False)
  plt.yticks([-tree.x-t for t in timeticks[:-1]],['{:.0f}'.format(-t) for t in timeticks[:-1]],va='center',ha='center',bbox=dict(color='white',alpha=0.5))
  plt.xticks(plt.xticks()[0],())
  plt.show()
  plt.close()
#print(tree)
