import json
class SpeciesDB(object):
  def __init__(self,filename='ITIS.json'):
    self.filename=filename
    with open(self.filename,'r',encoding='utf-8') as f:
      data=json.load(f)
    self.languages=data['Languages']
    self.ld={name.casefold():i for i,name in enumerate(self.languages)}
    self.ld['Deutsch'.casefold()]=self.ld['German'.casefold()]
    self.ld['française'.casefold()]=self.ld['French'.casefold()]
    self.ld['español'.casefold()]=self.ld['Spanish'.casefold()]
    self.dbl=data['db']
    id_unspecified=self.ld['unspecified'.casefold()]
    id_english=self.ld['English'.casefold()]
    for species in self.dbl:
      for name in species:
        if name[1]==id_unspecified:
          name[1]=id_english    
    self.dbd={name[0].casefold():dataset for dataset in self.dbl for name in dataset[1:]}
    self.dbd.update({name[0].casefold():dataset for dataset in self.dbl for name in dataset[:1]})


  def translate(self,name,language='Scientific'):  
    species=self.dbd.get(name.casefold())
    if not species:
      return '"{:s}" is not in the database.'.format(name)
    language_id=self.ld.get(language.casefold())
    if language_id==None:
      return 'Language "{:s}" is not in the database.'.format(language)
    for x in reversed(species):#compare last entry first: appended entries override previous ones
      if x[1]==language_id:return x[0]
    return self.dbd[name][0][0] #return the scientific name if no translation is found for the target language
  def addName(self,speciesname,newlanguage,newname,approved_ind='N'):
    species=self.dbd.get(speciesname.casefold())
    if not species: 
      return False #the species does not exist. Use addSpecies to add a new species
    existing_entry=self.dbd.get(newname.casefold())
    if existing_entry:
      if existing_entry[0][0].casefold()==newname.casefold():
        return False #you cannot use an existing scientific name
    language_id=self.ld.get(newlanguage.casefold())
    new_entry=[newname,language_id,approved_ind]
    for x in species:
      if x==new_entry:return True
    species.append(new_entry)
    self.dbd[newname.casefold()]=species
  def addSpecies(self,scientific_name,approved_ind='N'):
    if self.dbd.get(scientific_name.casefold()):
      return
    species=[[scientific_name,0,approved_ind]]
    if not species: return False
    self.dbl.append(species)
    self.dbd[scientific_name.casefold()]=species

if __name__=='__main__':
  from collections import namedtuple
  import re
  dictEntry=namedtuple('dictEntry','lt en de')    
  with open('species_dict.txt','r',encoding='utf-8')as f:
    speciesDict=[dictEntry(latin,english,deutsch) for latin,english,deutsch in re.findall('([^:]*):\s*([^:]*):\s*([^\n]*)\n',f.read())]
    translation={name:e for e in speciesDict for name in e}
  db=SpeciesDB()
  for s in speciesDict:
    db.addName(s.lt,'English',s.en)
    db.addName(s.lt,'German',s.de)
  
  db.addSpecies('Unicornis')
  db.addName('Unicornis','Deutsch','Einhorn')
  db.addName('Unicornis','English','unicorn') #this does not work, because "Unicorn" is already a scientific name in the ITIS database (Platnick & Brescovit, 1995)
  
  print(db.translate('Wolf'))
  print(db.translate('Wolf','english'))
  print(db.translate('Homo sapiens','English'))
  print(db.translate('Homo sapiens','German'))
  print(db.translate('Canada Goose','German'))
