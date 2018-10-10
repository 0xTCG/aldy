# 786

import re, os, sys, collections

from bs4 import BeautifulSoup  
from common import log

class autodict(dict):  
   def __missing__(self, key):
      value = self[key] = type(self)()
      return value

def get_karolinska_database(gene, parser, force=False, url=None):   
   log.info('Karolinska: Calculating Karolinska database')
   with open('pages/{}.htm'.format(gene, 'rb')) as f:
      page = f.read()
   result = parse_karolinska_html(page, parser)
   return result

def parse_karolinska_html(page, func):
   soup = BeautifulSoup(page, 'lxml')
   data = autodict()
   dbsnp = collections.defaultdict(lambda: '*')
   
   tables = soup.select('table table')
   if len(tables) == 0 or True:
      tables = soup.select('table')
   for table in tables:
      for row in table.find_all('tr'): # skip header
         items = func(row, dbsnp)
         if len(items) < 3 or len(filter(None, items[2])) == 0 or items[0] == 'Allele' or items[0] == '': 
            continue
         items[0] = items[0].split()[0]
         if 'X' in items[0] or 'x' in items[0]:
            log.warn('Karolinska: Ignoring {}', items[0])
            continue 
         if items[0] in data:
            log.warn('Karolinska: {} already exists, overwriting it', items[0])
            log.debug('Karolinska: overwriting {} with {}', items[0], ','.join(map(str, items[2])))
         data[items[0]].update({
            'mutations': items[2],
            'phenotype': {
               'invivo': items[-3],
               'invitro': items[-2],
            }
         })
   return data
