# 786

import json, uuid

class NoIndent(object):
   def __init__(self, value):
      self.value = value

class NoIndentEncoder(json.JSONEncoder): 
   """ http://stackoverflow.com/a/25935321 """
   
   def __init__(self, *args, **kwargs):
      super(NoIndentEncoder, self).__init__(*args, **kwargs)
      self.kwargs = dict(kwargs)
      del self.kwargs['indent']
      self._replacement_map = {}

   def default(self, o):
      if isinstance(o, NoIndent):
         key = uuid.uuid4().hex
         self._replacement_map[key] = json.dumps(o.value, **self.kwargs)
         return "@@%s@@" % (key,)
      else:
         return super(NoIndentEncoder, self).default(o)

   def encode(self, o):
      result = super(NoIndentEncoder, self).encode(o)
      for k, v in self._replacement_map.iteritems():
         result = result.replace('"@@%s@@"' % (k,), v)
      return result
