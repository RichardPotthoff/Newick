import sqlite3
import re
from collections import namedtuple
import json
conn=sqlite3.connect('file:../ITIS.sqlite?mode=ro')
c=conn.cursor()
tablenames=list(c.execute("SELECT name FROM sqlite_master WHERE type='table';"))
def print_tableinfo(conn=conn,tablenames=tablenames):
  c=conn.cursor()
  for name in tablenames:
    print('Table: {:20s}'.format(name[0]))
    fields=list(c.execute('PRAGMA table_info({:s});'.format(name[0])))
    for i in (1,2):
      print('|'.join(['{:18s}'.format(x[i]) for x in fields]))
    print()
c.execute("SELECT distinct language FROM vernaculars;")
languages=c.fetchall()
def languagesTest():
  for language in languages:
    for x in c.execute("SELECT language, count() FROM vernaculars where language=?;",language):
      print('language:{:s} count:{:d}'.format(*x))
    
dictEntry=namedtuple('dictEntry','lt en de')    
with open('species_dict.txt','r',encoding='utf-8')as f:
  speciesDict=[dictEntry(latin,english,deutsch) for latin,english,deutsch in re.findall('([^:]*):\s*([^:]*):\s*([^\n]*)\n',f.read())]
  translation={name:e for e in speciesDict for name in e}
  
def translate(name,language):
  t=translation.get(name)
  tl=t._asdict().get(language) if t else None
  return tl if tl else name

def translateTest(speciesDict=speciesDict):
  for x in speciesDict:
    result=c.execute("SELECT tsn, completename FROM longnames where completename=?;",(x.lt,))
    print(x.lt,end=' ')
    for tsn,longname in result.fetchall():
      print(tsn,longname,end=' ')
      synonyms=c.execute("SELECT vernacular_name,language,approved_ind FROM vernaculars where tsn=? AND approved_ind='Y';",(tsn,))
      for x in synonyms.fetchall():
        print(', {:s}({:s}/{:s})'.format(*x),end='')
    print()
#languages=c.fetchall()
#print('vernacular_name count: {:d}'.format(len(languages)))
languages=['Scientific']+[language for language, in languages]
ld={language:i for i,language in enumerate(languages)}
vd={tsn:[] for tsn, in conn.cursor().execute("SELECT distinct tsn FROM vernaculars;")}
for tsn,vernacular_name,language,approved_ind in conn.cursor().execute("SELECT tsn,vernacular_name,language,approved_ind FROM vernaculars;"):
  vd[tsn].append((vernacular_name,ld[language],approved_ind))
dbl=[((complete_name,0,'Y'),*(vd.get(tsn) if vd.get(tsn) else [])) for tsn,complete_name in conn.cursor().execute("SELECT tsn, completename FROM longnames;")]

dbd={cn[0][0]:(cn[0],) for cn in dbl}
for x in dbl:
  dbd[x[0][0]]=(*dbd[x[0][0]],*x[1:])#merge data for duplicate keys
dbl1=list(dbd.values())

with open('ITIS.json','w') as f:
  json.dump({'Description':'Languages is an array which field#2 in "db" refers to','Languages':languages,'db':dbl1},f)
