import json
with open('ITIS.json','r',encoding='utf-8') as f:
   data=json.load(f)
 
languages=data['Languages']
ld={name:i for i,name in enumerate(languages)}
dbl=data['db']
dbd={name[0]:dataset for dataset in dbl for name in dataset}

def translate(name,language='Scientific'):
  language_id=ld[language]
  species=dbd.get(name)
  if not species:
    return '"{:s}" is not in the database.'.format(name)
  for x in species:
    if x[1]==language_id:return x[0]
  return dbd[name][0][0] #return the scientific name if no translation is found for the target language
    
print(translate('Wolf'))
print(translate('Homo sapiens','English'))
