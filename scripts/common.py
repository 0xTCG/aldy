#!/usr/bin/env python
# 786

import os, re, time
import pprint
import logbook, logbook.more

def colorize(text, color='green'):
   return logbook._termcolors.colorize(color, text)

log = logbook.Logger('Cypiripi')
LOG_FORMAT = '{record.message}'

sh = logbook.more.ColorizedStderrHandler(format_string=LOG_FORMAT, level='TRACE')
sh.push_application()
