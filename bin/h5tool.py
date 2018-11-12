# Copyright (c) 2011-2013, Pacific Biosciences of California, Inc.
#
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
# * Redistributions of source code must retain the above copyright
#   notice, this list of conditions and the following disclaimer.
#
# * Redistributions in binary form must reproduce the above copyright
#   notice, this list of conditions and the following disclaimer in the
#   documentation and/or other materials provided with the distribution.
#
# * Neither the name of Pacific Biosciences nor the names of its
#   contributors may be used to endorse or promote products derived from
#   this software without specific prior written permission.
#
# NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
# GRANTED BY THIS LICENSE.  THIS SOFTWARE IS PROVIDED BY PACIFIC
# BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
# WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS CONTRIBUTORS
# BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
# BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
# WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
# OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN
# IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import os
import sys
import argparse
import shutil
import datetime
import logging
import tempfile

import h5py as H5
import numpy as NP

from pbcore.io import CmpH5Reader

from pbh5tools.PBH5ToolsException import PBH5ToolsException
from pbh5tools.CmpH5Format import CmpH5Format
from pbh5tools.CmpH5Utils import *
from pbh5tools.Metrics import *

def main ():
    ### parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('inCmp', metavar='input.cmp.h5')
    parser.add_argument('whitelist', metavar='whitelist.txt')
    parser.add_argument('--outFile',
                        default = "out.cmp.h5",
                        dest='outCmp', metavar='out.cmp.h5',
                        help = "Either a pattern string or a filename")
    args = parser.parse_args()

    ### read whitelist
    wl = set(line.strip() for line in open(args.whitelist))

    ### load aligned_reads.cmp.h5
    cmpH5 = CmpH5Reader(args.inCmp)

    ### alignment indices
    idxs = []

    for (i,aln) in enumerate(cmpH5):
        qname = "{:s}/{:d}".format(aln.movieInfo.Name,aln.HoleNumber)
        if qname in wl:
            idxs.append(aln.AlnID)

    ### save subset of alignments
    doSelect(args.inCmp,args.outCmp,idxs)

def doSelect(inCmpFile, outCmpFile, idxs):
    """Take an input cmp.h5 file and a vector of indices into the
    AlnIndex and create a new cmp.h5 file from those alignments."""
    def trimDataset(groupName, alnIdxID, inCmp, outCmp, fmt, idName = 'ID'):
        ids = outCmp[fmt.ALN_INDEX][:,alnIdxID]
        nds = '/'.join([groupName, idName])
        msk = NP.array([x in ids for x in inCmp[nds].value]) # got to be an NP.array
        for dsName in inCmp[groupName].keys():
            copyDataset('/'.join([groupName, dsName]), inCmp, outCmp,
                        msk, fmt)

    def copyGroup(groupName, inCmp, outCmp):
        if groupName in inCmp:
            outCmp.copy(inCmp[groupName], groupName)

    try:
        inCmp  = H5.File(inCmpFile, 'r')
        outCmp = H5.File(outCmpFile, 'w') # fail if it exists.
        idxs   = NP.array(idxs)
        fmt    = CmpH5Format(inCmp)

        if not (NP.max(idxs) < inCmp[fmt.ALN_INDEX].shape[0] and
                NP.min(idxs) >= 0):
            raise PBH5ToolsException("Invalid idxs specified, must be within [0, %d)" %
                                     inCmp[fmt.ALN_INDEX].shape[0])

        # copy over the AlnIndex and other AlnInfo elements
        # correpsonding to idxs to new file.
        for dsName in inCmp[fmt.ALN_INFO].keys():
            copyDataset('/'.join([fmt.ALN_INFO, dsName]), inCmp, outCmp, idxs, fmt)

        # reset the ALN_ID
        outCmp[fmt.ALN_INDEX][:,fmt.ID] = \
            NP.array(range(1, outCmp[fmt.ALN_INDEX].shape[0] + 1))

        # trim the other datasets
        trimDataset(fmt.ALN_GROUP, fmt.ALN_ID, inCmp, outCmp, fmt)
        # trimDataset(fmt.REF_GROUP, fmt.REF_ID, inCmp, outCmp, fmt)
        # trimDataset(fmt.MOVIE_INFO, fmt.MOVIE_ID, inCmp, outCmp, fmt)
        # copy Ref,Movie dataset whole
        for groupName in [fmt.REF_GROUP,fmt.MOVIE_INFO]:
            for dsName in inCmp[groupName].keys():
                copyDataset('/'.join([groupName,dsName]), inCmp, outCmp, None, fmt)

        # other groups will go over whole hog
        copyGroup(fmt.FILE_LOG, inCmp, outCmp)
        copyGroup(fmt.REF_INFO, inCmp, outCmp)
        copyGroup(fmt.BARCODE_INFO, inCmp, outCmp)

        # now we copy over the actual data
        for i in xrange(0, outCmp[fmt.ALN_GROUP_ID].shape[0]):
            # figure out what reads are in this group.
            agID = outCmp[fmt.ALN_GROUP_ID][i]
            agPT = outCmp[fmt.ALN_GROUP_PATH][i]
            alnIdx = outCmp[fmt.ALN_INDEX].value
            whReads = NP.where(agID == alnIdx[:,fmt.ALN_ID])[0]
            offBegin = alnIdx[whReads, fmt.OFFSET_BEGIN]
            offEnd = alnIdx[whReads, fmt.OFFSET_END]
            totalSize = NP.sum((offEnd - offBegin) + 1) # 0 in between

            for dsName in inCmp[agPT].keys():
                fullPath = '/'.join([agPT, dsName])
                newDs = outCmp.create_dataset(fullPath, shape = (totalSize,),
                                              dtype = inCmp[fullPath].dtype)
                origDs = inCmp[fullPath]
                cs = 0
                for j in xrange(0, len(whReads)):
                    newEnd = cs + offEnd[j] - offBegin[j]
                    newDs[cs:newEnd] = origDs[offBegin[j]:offEnd[j]]
                    outCmp[fmt.ALN_INDEX][whReads[j],fmt.OFFSET_BEGIN] = cs
                    outCmp[fmt.ALN_INDEX][whReads[j],fmt.OFFSET_END] = newEnd
                    cs = newEnd

        # copy over the top-level attributes
        copyAttributes(inCmp, outCmp)

        # remove the offset table
        deleteIfExists(outCmp, fmt.REF_OFFSET_TABLE)
        deleteAttrIfExists(outCmp, fmt.INDEX_ATTR)

        # close the sucker
        logging.debug("Closing output cmp.h5 file.")
        outCmp.close()

    except Exception, e:
        logging.exception(e)
        try:
            os.remove(outCmpFile)
        except:
            pass
        raise e

if __name__ == "__main__":
    main()
