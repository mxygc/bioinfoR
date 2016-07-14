#!/usr/local/bin/python2.7
# encoding: utf-8


import re
#import urllib
import urllib2
from bs4 import BeautifulSoup
import xlsxwriter
from _ast import Pass
from __future__ import division, print_function, unicode_literals
import os
import sys

DEBUG = 1
TESTRUN = 0
PROFILE = 0

__version__ = 0.1
__date__ = '2014-05-23'
__updated__ = '2014-05-23'
__all__ = [DEBUG, TESTRUN, PROFILE, __version__, __date__, __updated__]

class SABtractor:
    
    def __init__(self, gset, xlsx):
        self.gset = gset
        self.xlsx = xlsx
        self.url_base = "http://www.sabiosciences.com/rt_pcr_product/HTML"    
        self.header = {
                       'Host': 'www.sabiosciences.com',
                       'User-Agent': 'Mozilla/5.0 (Windows NT 6.2; rv:29.0) Gecko/20100101 Firefox/29.0',
                       'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8',
                       'Accept-Encoding': 'gzip, deflate',
                       'Connection': 'keep-alive'
        }
        self.patterns = [
                         re.compile(r'Gene Table(?P<table>.+)Gene Table', re.DOTALL),                                          #patt0
                         re.compile(r'(?:<)([^<>]+?):([^<>]+)(?:>?)'),                                                         #patt1
                         re.compile(r'(^.+?):(.*?$)', re.MULTILINE), #it is not multi-line matching, so no need to use ^ or $  #patt2
                         re.compile(r'<[^<>]*?>'),                                                                             #patt3
                         re.compile(r'(^\s+)|(\s+$)', re.DOTALL),                                                              #patt4
                         re.compile(r',')                                                                                      #patt5
                         ]
        
    def __del__(self):
        self.workbook.close()
        self.fin.close()

    def init(self):
        if not os.path.exists(self.gset):
            sys.stderr.write('input file not exists')
        dirname = os.path.dirname(self.xlsx)
        basename = os.path.basename(self.xlsx)
        if dirname == '':
            dirname = '.'
        self.dirname = os.path.join(dirname, re.sub("\.xlsx$", '', basename))
        if not os.path.exists(dirname):
            os.makedirs(dirname)
        try:
            self.workbook = xlsxwriter.Workbook(self.xlsx)
            self.worksheet = self.workbook.add_worksheet('README')
            self.worksheet.write_row(0, 0, ['set_name', 'set_id'])
        except Exception, e:
            if DEBUG or TESTRUN:
                raise(e)
        try:
            self.fin = open(self.gset, "r")
        except Exception, e:
            if DEBUG or TESTRUN:
                raise(e)
            
    def extract(self):           
        l = self.fin.readline().rstrip()
        i = 0
        while l:
            i = i + 1
            self.worksheet.write_row(i, 0, l.split("\t"))
            l = self.fin.readline().rstrip()
        del i
        self.fin.seek(0)
        l = self.fin.readline().rstrip()
        while l:
            set_name, set_id = l.split("\t")
            print("pricessing set %s => %s..." % (set_name, set_id))
            webpage = self.get_webpage(self, set_id)
            geneset = self.get_geneset(self, webpage)
            self.align_geneset(self, geneset, '')
            text = self.fmt_geneset(self, geneset, '\t')
            self.save_geneset(self, text, set_id, geneset)
            l = self.fin.readline().rstrip()

    def get_webpage(self, set_id):
        url = self.url_base + '/' + set_id + '.html'
        if DEBUG or TESTRUN:
            sys.stderr.write(url)
        req = urllib2.Request(url, headers = self.header)
        con = urllib2.urlopen(req)
        doc = con.read()
        con.close()
        return doc
        
    def get_geneset(self, webpage): 
        soup = BeautifulSoup(webpage)
        [s.extract() for s in soup('script')]
        [s.extract() for s in soup('noscript')]
        text = soup.prettify()
        table = self.patterns[0].search(text).group("table")
        table = re.sub("\r", "", table)
        table = re.sub("\n", "", table)
        table = re.sub(r',\s*<br\/>\s*', ',', table)
        table = re.sub(r'<br\/>', "\n", table)
        table = re.sub(r'<\/p>', "\n", table)
        table = re.sub("\.\s+", "\n", table)
        table = re.sub(r'&amp;', '&', table)
        
        (table, repn) = self.patterns[1].subn(r'<\1 - \2>', table)
        while (repn):
            (table, repn) = self.patterns[1].subn(r'<\1 - \2>', table)
        geneset = []
        t1 = ''
        for m in self.patterns[2].finditer(table):
            if m is not None:
                t = m.group(1)
                t2 = ''
                g = m.group(2)
                g = self.patterns[3].sub('', g)
                g = self.patterns[4].sub('', g)
                if g == '':
                    t = self.patterns[3].sub('', t)
                    t = self.patterns[4].sub('', t)
                    t1 = t
                else:
                    if (re.search("bold|<strong>", t) is not None):
                        t = self.patterns[3].sub('', t)
                        t = self.patterns[4].sub('', t)
                        t1 = t
                    elif (re.search("underline|<u>", t) is not None):
                        t = self.patterns[3].sub('', t)
                        t = self.patterns[4].sub('', t)        
                        t2 = t
                    else:
                        t = self.patterns[3].sub('', t)
                        t = self.patterns[4].sub('', t)        
                        t2 = t
                    genes = self.patterns[5].split(g)
                    if len(genes) < 1:
                        sys.stderr.write("in T1:%s\tT2:%s, there is not gene found!" % (t1, t2))
                        continue
                    genes = map(lambda g: self.patterns[4].sub('', g), genes)
                    geneset.append([t1, t2, genes])
        return geneset
    
    def align_geneset(self, geneset, fill = ''):
        nc = len(geneset)
        nrows = map(lambda gs: len(gs[2]), geneset)
        nr = max(nrows)
        na = map(lambda n: nr - n, nrows)
        map(lambda i: geneset[i][2].extend([fill] * na[i]), range(nc))
    
    def fmt_geneset(self, geneset, sep = ","):
        text = ""
        nc = len(geneset)
        nr = len(geneset[0][2])
        text = text + sep.join([gs[0] for gs in geneset]) + '\n'
        text = text + sep.join([gs[1] for gs in geneset]) + '\n'
        text = text + '\n'.join([sep.join([geneset[j][2][i] for j in range(nc)]) for i in range(nr)])
        return text
    
    def save_geneset_xlsx(self, set_id, geneset):
        nc = len(geneset)
        nr = len(geneset[0][2])
        self.worksheet = self.workbook.add_worksheet(set_id)
        # strange in the following, geneset does not need ".encode('utf-8')"
        [self.worksheet.write(0, i, geneset[i][0]) for i in range(nc)]
        [self.worksheet.write(1, i, geneset[i][1]) for i in range(nc)]
        [[self.worksheet.write(i + 2, j, geneset[j][2][i]) for j in range(nc)] for i in range(nr)]
    
    def save_geneset_txt(self, txt, text):
        fout = open(txt, "w")
        fout.write(text.encode('utf-8'))
        fout.close()
    
    def save_geneset(self, text, set_id, geneset):
        txt = self.basename + '/' + set_id + '.txt'
        if geneset is not None:
            self.save_geneset_xlsx(self.workbook, set_id, geneset)
        else:
            sys.stderr.write('null geneset')
        self.save_geneset_txt(txt, text)