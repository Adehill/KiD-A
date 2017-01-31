#!/usr/bin/env python
#
# Make dependencies from fortran source
#

import sys, getopt
import os
import glob
import re

re_module=re.compile(r'^\s*module\s+([^\s!]+)',re.I|re.MULTILINE)
re_include=re.compile(r'^\s*include\s+["\']([^"\']+)["\']',re.I|re.MULTILINE)
re_use=re.compile(r'^\s*use\s+([^\s,!]+)',re.I|re.MULTILINE)
re_dependsOn=re.compile(r'^\s*! DEPENDS ON:\s+([^\s,!]+)',re.I|re.MULTILINE)

def makeDependencies(vpath=None, outfile='.includes', objaffix=''):
    fdo=open(outfile, 'w')
    files=glob.glob('[a-zA-Z]*.f[9]*[05]*')
    if vpath is not None:
        for path in vpath:
            files=files + glob.glob('%s/[a-zA-Z]*.f[9]*[05]*' % path)
    filenames={}
    objects={}
    incs={}
    uses={}
    dependsOn={}
    for file in files:
        fd=open(file, 'r')
        code=fd.read()
        fd.close()
        #Search for modules
        m=re_module.search(code)
        if m:
            mod=m.groups()[0].lower()
            filenames[mod]=file
            objects[mod]=re.sub('\.f.*$','.o',file)
    for file in files:
        fd=open(file, 'r')
        code=fd.read()
        fd.close()
        obj=re.sub('\.f.*$','.o',file)
        incs[obj]=[]
        uses[obj]=[]
        dependsOn[obj]=[]
        #search for includes
        m=re_include.findall(code)
        if m:
            for i in m:
                incs[obj].append(i)
        #search for uses
        m=re_use.findall(code)
        if m:
            for i in m:
                uses[obj].append(i)
        #write dependencies to file
        fdo.write('%s%s: %s ' % (objaffix, obj, file))
        for i in uses[obj]+incs[obj]+dependsOn[obj]:
            if i.lower() in filenames.keys():
                fdo.write('%s%s ' % (objaffix,objects[i.lower()]))
            for path in vpath:
                if path+'/'+i.lower() in filenames.keys():
                    fdo.write('%s%s ' % (objaffix,objects[path+i.lower()]))
#            else:
#                print 'Error: dependency not found for "%s" in %s.' % (i, file)
        fdo.write('\n\n')
    fdo.close()
        

if __name__=="__main__":
    
    opts, args = getopt.getopt(sys.argv[1:], "v:a:")

    vpath=[]
    for o, a in opts:
        if o == "-v":
            vpath.append(a)
        if o == "-a":
            objaffix=a

    makeDependencies(vpath,objaffix=a)
