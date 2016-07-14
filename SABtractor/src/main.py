#!/usr/local/bin/python2.7
# encoding: utf-8
'''
This is a simple gene sets extractor from SABioscience, the source website 

@author:     foehn <foehn@foxmail.com>
@license:    GPL-v2
@contact:    foehn@foxmail.com
'''

import os, sys
from SABtractor import __version__, __updated__, __date__, \
                       DEBUG, TESTRUN, PROFILE, SABtractor
from argparse import ArgumentParser, RawDescriptionHelpFormatter

class CLIError(Exception):
    '''Generic exception to raise and log different fatal errors.'''
    def __init__(self, msg):
        super(CLIError).__init__(type(self))
        self.msg = "E: %s" % msg
    def __str__(self):
        return self.msg
    def __unicode__(self):
        return self.msg

def main(argv = None): # IGNORE:C0111
    '''Command line options.'''

    if argv is None:
        argv = sys.argv
    else:
        sys.argv.extend(argv)

    program_name = os.path.basename(sys.argv[0])
    program_version = "v%s" % __version__
    program_build_date = str(__updated__)
    program_version_message = '%%(prog)s %s (%s)' % (program_version, 
                                                     program_build_date)
    program_shortdesc = __import__('__main__').__doc__.split("\n")[1]
    program_license = '''%s

  Created by foehn on %s.
  Copyright foehn 2014. All rights reserved.

  Licensed under the GNU General Public License, version 2 (GPL-2.0)
  https://opensource.org/licenses/GPL-2.0

  Distributed on an "AS IS" basis without warranties
  or conditions of any kind, either express or implied.
  
  
''' % (program_shortdesc, str(__date__))

    try:
        # Setup argument parser
        parser = ArgumentParser(description = program_license, 
                                formatter_class = RawDescriptionHelpFormatter)
        parser.add_argument('-i', '--input', type = str, dest = 'gset', 
                            nargs = 1, help = 'name of the input file containing \
                                               genesets to be extracted')
        parser.add_argument('-o', '--output', type = str, dest = 'xlsx',
                            nargs = 1, default = 'sab_genes.xlsx', help = '\
                            output file name [default: %(default)s]')
        parser.add_argument("-v", "--verbose", dest = "verbose", action = "count", 
                            help = "set verbosity level [default: %(default)s]")
        parser.add_argument('-V', '--version', action = 'version', 
                            version = program_version_message)
        parser.add_argument(dest = 'paths', help = 'paths to folder(s) with \
                            source file(s) [default: %(default)s]', 
                            metavar = 'path', nargs = '+')

        # Process arguments
        args = parser.parse_args()
        paths = args.paths
        verbose = args.verbose
        gset = args.gset
        xlsx = args.xlsx
        
        if verbose > 0:
            print("Verbose mode on")

        for inpath in paths:
            ### do something with inpath ###
            print(inpath)
        return 0
    except KeyboardInterrupt:
        ### handle keyboard interrupt ###
        return 0
    except Exception, e:
        if DEBUG or TESTRUN:
            raise(e)
        indent = len(program_name) * " "
        sys.stderr.write(program_name + ": " + repr(e) + "\n")
        sys.stderr.write(indent + "  for help use --help")
        return 2
    sab= SABtractor(gset, xlsx)
    sab.init()
    sab.extract()

    
    
if __name__ == "__main__":
    if DEBUG:
        sys.argv.append("-h")
        sys.argv.append("-v")
    if TESTRUN:
        import doctest
        doctest.testmod()
    if PROFILE:
        import cProfile
        import pstats
        profile_filename = 'profile.txt'
        cProfile.run('main()', profile_filename)
        statsfile = open("profile_stats.txt", "wb")
        p = pstats.Stats(profile_filename, stream=statsfile)
        stats = p.strip_dirs().sort_stats('cumulative')
        stats.print_stats()
        statsfile.close()
        sys.exit(0)
    sys.exit(main())