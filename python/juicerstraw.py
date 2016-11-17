#this is the module file
#Reads the genome name from the hic header and the program will output  x, y, and count
#Can only take in a .hic file
from __future__ import absolute_import, division, print_function, unicode_literals
import sys
import struct
# import urllib
import wget
import zlib
blockMap = dict()
# global version
version=0
#read function
def readcstr(f):
    # buf = bytearray()
    buf = b""
    while True:
        b = f.read(1)
        if b is None or b == b"\0":
            # return str(buf,encoding="utf-8", errors="strict")
            return buf.decode("utf-8")
        else:
            buf += b
            # buf.append(b)

#FUN(req, chr1, chr2, [c1pos1, c1pos2, c2pos1, c2pos2]) Return [master, chr1ind, chr2ind]
def readHeader(req, chr1, chr2, posilist):
    magic_string = struct.unpack('<3s', req.read(3))[0]
    req.read(1)
    if (magic_string != b"HIC"):
        print('This does not appear to be a HiC file magic string is incorrect')
        return -1
    global version
    version = struct.unpack('<i',req.read(4))[0]
    if (version < 6):
        print("Version {0} no longer supported".format(str(version)))
        return -1
    # print('HiC version:' + '  {0}'.format(str(version)))
    master = struct.unpack('<q',req.read(8))[0]
    genome = b""
    c=req.read(1)
    while (c != b'\0'):
        genome += c
        c=req.read(1)
    # print('Genome ID:','  {0}'.format(str(genome, encoding="utf-8", errors="strict")))
    # read and throw away attribute dictionary (stats+graphs)
    nattributes = struct.unpack('<i',req.read(4))[0]
    for x in range(nattributes):
      key = readcstr(req)
      value = readcstr(req)
      # print('   Key:{0}'.format(key))
      # print('   Value:{0}'.format(value))
    nChrs = struct.unpack('<i',req.read(4))[0]
    # print("Chromosomes: ")
    found1 = False
    found2 = False
    for i in range(0, nChrs):
      name = readcstr(req)
    #   print(str(name))
      length = struct.unpack('<i',req.read(4))[0]
      if (name==chr1):
          found1=True
          chr1ind=i
          if (posilist[0]==-100):
              posilist[0]=0
              posilist[1]=length
      if (name==chr2):
          found2=True
          chr2ind=i
          if (posilist[2]==-100):
              posilist[2]=0
              posilist[3]=length
    if ((not found1) or (not found2)):
      print("One of the chromosomes wasn't found in the file. Check that the chromosome name matches the genome.\n")
      return -1
    return [master, chr1ind, chr2ind, posilist[0], posilist[1], posilist[2], posilist[3]]

#FUN(fin, master, c1, c2, norm, unit, resolution) Retrun [myFilePos, c1NormEntry, c2NormEntry]
def readFooter(req, master, c1, c2, norm, unit, resolution):
    c1NormEntry=dict()
    c2NormEntry=dict()
    req.seek(master)
    nBytes = struct.unpack('<i', req.read(4))[0]
    key = str(c1) + "_" + str(c2)
    nEntries = struct.unpack('<i', req.read(4))[0]
    found = False
    for i in range(nEntries):
        stri = readcstr(req)
        fpos = struct.unpack('<q', req.read(8))[0]
        sizeinbytes = struct.unpack('<i', req.read(4))[0]
        if (stri==key):
            myFilePos = fpos
            found=True
    if (not found):
        print("File doesn't have the given chr_chr map\n")
    if (norm=="NONE"):
        return [myFilePos, 0, 0]
    nExpectedValues = struct.unpack('<i',req.read(4))[0]
    for i in range(nExpectedValues):
        str_ = readcstr(req)
        binSize = struct.unpack('<i',req.read(4))[0]
        nValues = struct.unpack('<i',req.read(4))[0]
        for j in range(nValues):
            v = struct.unpack('<d',req.read(8))[0]
        nNormalizationFactors = struct.unpack('<i',req.read(4))[0]
        for j in range(nNormalizationFactors):
            chrIdx = struct.unpack('<i',req.read(4))[0]
            v = struct.unpack('<d',req.read(8))[0]
    nExpectedValues = struct.unpack('<i',req.read(4))[0]
    for i in range(nExpectedValues):
        str_ = readcstr(req)
        str_ = readcstr(req)
        binSize = struct.unpack('<i',req.read(4))[0]
        nValues = struct.unpack('<i',req.read(4))[0]
        for j in range(nValues):
            v = struct.unpack('<d',req.read(8))[0]
        nNormalizationFactors = struct.unpack('<i',req.read(4))[0]
        for j in range(nNormalizationFactors):
            chrIdx = struct.unpack('<i',req.read(4))[0]
            v = struct.unpack('<d',req.read(8))[0]
    nEntries = struct.unpack('<i',req.read(4))[0]
    found1=False
    found2=False
    for i in range(nEntries):
        normtype = readcstr(req)
        chrIdx = struct.unpack('<i',req.read(4))[0]
        unit1 = readcstr(req)
        resolution1 = struct.unpack('<i',req.read(4))[0]
        filePosition = struct.unpack('<q',req.read(8))[0]
        sizeInBytes = struct.unpack('<i',req.read(4))[0]
        if (chrIdx==c1 and normtype==norm and unit1==unit and resolution1==resolution):
            # print("****"+str(filePosition))
            c1NormEntry['position']=filePosition
            c1NormEntry['size']=sizeInBytes
            found1=True
        if (chrIdx==c2 and normtype==norm and unit1==unit and resolution1==resolution):
            c2NormEntry['position']=filePosition
            c2NormEntry['size']=sizeInBytes
            found2=True
    if ((not found1) or (not found2)):
        print("File did not contain {0} normalization vectors for one or both chromosomes at {1} {2}\n".format(norm, resolution, unit))
        return -1
    return [myFilePos, c1NormEntry, c2NormEntry]

#FUN(fin, unit, resolution) Return [storeBlockData, myBlockBinCount, myBlockColumnCount]
def readMatrixZoomData(req, myunit, mybinsize):
    unit = readcstr(req)
    temp = struct.unpack('<i',req.read(4))[0]
    temp2 = struct.unpack('<f',req.read(4))[0]
    temp2 = struct.unpack('<f',req.read(4))[0]
    temp2 = struct.unpack('<f',req.read(4))[0]
    temp2 = struct.unpack('<f',req.read(4))[0]
    binSize = struct.unpack('<i',req.read(4))[0]
    blockBinCount = struct.unpack('<i',req.read(4))[0]
    blockColumnCount = struct.unpack('<i',req.read(4))[0]
    storeBlockData = False
    #for the initial
    myBlockBinCount = -1
    myBlockColumnCount = -1
    if (myunit==unit and mybinsize==binSize):
        myBlockBinCount=blockBinCount
        myBlockColumnCount=blockColumnCount
        storeBlockData=True
    nBlocks = struct.unpack('<i',req.read(4))[0]
    for b in range(nBlocks):
        blockNumber = struct.unpack('<i',req.read(4))[0]
        filePosition = struct.unpack('<q',req.read(8))[0]
        blockSizeInBytes = struct.unpack('<i',req.read(4))[0]
        entry=dict()
        entry['size'] = blockSizeInBytes
        entry['position'] = filePosition
        if (storeBlockData):
            blockMap[blockNumber] = entry
    return [storeBlockData, myBlockBinCount, myBlockColumnCount]

#FUN(fin, myFilePos, unit, binsize) Return [blockBinCount, blockColumnCount]
def readMatrix(req, myFilePos, unit, binsize):
    # print(str(req), "\t", str(myFilePos),"\t", str(unit), "\t", str(binsize))
    req.seek(myFilePos)
    c1 = struct.unpack('<i',req.read(4))[0]
    c2 = struct.unpack('<i',req.read(4))[0]
    nRes = struct.unpack('<i',req.read(4))[0]
    # print(str(c1),"\t",str(c2),"\t",str(nRes),"\n")
    i = 0
    found = False
    blockBinCount = -1
    blockColumnCount = -1
    while (i<nRes and (not found)):
        list1 = readMatrixZoomData(req, unit, binsize)
        found = list1[0]
        if(list1[1]!=-1 and list1[2]!=-1):
            blockBinCount = list1[1]
            blockColumnCount = list1[2]
        i=i+1
    if (not found):
        print("Error finding block data\n")
        return -1
    return [blockBinCount, blockColumnCount]

#FUN(regionIndices, blockBinCount, blockColumnCount, intra) Return blocksSet
def getBlockNumbersForRegionFromBinPosition(regionIndices, blockBinCount, blockColumnCount, intra):
    col1=int(regionIndices[0]/blockBinCount)
    col2=int((regionIndices[1]+1)/blockBinCount)
    row1=int(regionIndices[2]/blockBinCount)
    row2=int((regionIndices[3]+1)/blockBinCount)
    blocksSet=set()
    # print(str(col1)+"\t"+str(col2)+"\t"+str(row1)+"\t"+str(row2))
    for r in range(row1, row2+1):
        for c in range(col1, col2+1):
            blockNumber=r*blockColumnCount+c
            blocksSet.add(blockNumber)
    if (intra):
        for r in range(col1, col2+1):
            for c in range(row1, row2+1):
                blockNumber=r*blockColumnCount+c
                blocksSet.add(blockNumber)
    # print(str(blocksSet))
    return blocksSet

#FUN(fin, blockNumber)
def readBlock(req, blockNumber):
    idx=dict()
    if(blockNumber in blockMap):
        idx=blockMap[blockNumber]
    else:
        idx['size']=0
        idx['position']=0
    if (idx['size']==0):
        return []
    req.seek(idx['position'])
    compressedBytes = req.read(idx['size'])
    uncompressedBytes = zlib.decompress(compressedBytes)
    nRecords = struct.unpack('<i',uncompressedBytes[0:4])[0]
    v=[]
    global version
    if (version < 7):
        for i in range(nRecords):
            binX = struct.unpack('<i',uncompressedBytes[(12*i):(12*i+4)])[0]
            binY = struct.unpack('<i',uncompressedBytes[(12*i+4):(12*i+8)])[0]
            counts = struct.unpack('<f',uncompressedBytes[(12*i+8):(12*i+12)])[0]
            record = dict()
            record['binX'] = binX
            record['binY'] = binY
            record['counts'] = counts
            v.append(record)
    else:
        binXOffset = struct.unpack('<i',uncompressedBytes[4:8])[0]
        binYOffset = struct.unpack('<i',uncompressedBytes[8:12])[0]
        useShort = struct.unpack('<b',uncompressedBytes[12:13])[0]
        type_ = struct.unpack('<b',uncompressedBytes[13:14])[0]
        index=0
        if (type_==1):
            rowCount = struct.unpack('<h',uncompressedBytes[14:16])[0]
            temp = 16
            for i in range(rowCount):
                y = struct.unpack('<h',uncompressedBytes[temp:(temp+2)])[0]
                temp=temp+2
                binY = y + binYOffset
                colCount = struct.unpack('<h',uncompressedBytes[temp:(temp+2)])[0]
                temp=temp+2
                for j in range(colCount):
                    x = struct.unpack('<h',uncompressedBytes[temp:(temp+2)])[0]
                    temp=temp+2
                    binX = binXOffset + x
                    if (useShort==0):
                        c = struct.unpack('<h',uncompressedBytes[temp:(temp+2)])[0]
                        temp=temp+2
                        counts = c
                    else:
                        counts = struct.unpack('<f',uncompressedBytes[temp:(temp+4)])[0]
                        temp=temp+4
                    record = dict()
                    record['binX'] = binX
                    record['binY'] = binY
                    record['counts'] = counts
                    v.append(record)
                    index = index + 1
        elif (type_== 2):
            temp=14
            nPts = struct.unpack('<i',uncompressedBytes[temp:(temp+4)])[0]
            temp=temp+4
            w = struct.unpack('<h',uncompressedBytes[temp:(temp+2)])[0]
            temp=temp+2
            for i in range(nPts):
                row=i/w
                col=i-row*w
                bin1=int(binXOffset+col)
                bin2=int(binYOffset+row)
                if (useShort==0):
                    c = struct.unpack('<h',uncompressedBytes[temp:(temp+2)])[0]
                    temp=temp+2
                    if (c != -32768):
                        record = dict()
                        record['binX'] = bin1
                        record['binY'] = bin2
                        record['counts'] = c
                        v.append(record)
                        index = index + 1
                else:
                    counts = struct.unpack('<f',uncompressedBytes[temp:(temp+4)])[0]
                    temp=temp+4
                    if (countsnot != 0x7fc00000):
                        record = dict()
                        record['binX'] = bin1
                        record['binY'] = bin2
                        record['counts'] = counts
                        v.append(record)
                        index = index + 1
    return v

#FUN(fin, entry) Return Norm
def readNormalizationVector(req, entry):
    req.seek(entry['position'])
    nValues = struct.unpack('<i',req.read(4))[0]
    value = []
    for i in range(nValues):
        d = struct.unpack('<d',req.read(8))[0]
        value.append(d)
    return value

#FUN(norm, fname, binsize, chr1loc, chr2loc, unit) Return [x, y, counts]
def straw(norm, req, binsize, chr1loc, chr2loc, unit):
    if (not (norm=="NONE" or norm=="VC" or norm=="VC_SQRT" or norm=="KR")):
        print("Norm specified incorrectly, must be one of <NONE/VC/VC_SQRT/KR>\nUsage: straw <NONE/VC/VC_SQRT/KR> <hicFile(s)> <chr1>[:x1:x2] <chr2>[:y1:y2] <BP/FRAG> <binsize>\n")
        return -1
    if (not (unit=="BP" or unit=="FRAG")):
        print("Norm specified incorrectly, must be one of <BP/FRAG>\nUsage: straw <NONE/VC/VC_SQRT/KR> <hicFile(s)> <chr1>[:x1:x2] <chr2>[:y1:y2] <BP/FRAG> <binsize>\n")
        return -1
    c1pos1=-100
    c1pos2=-100
    c2pos1=-100
    c2pos2=-100
    chr1_arra = chr1loc.split(":")
    chr2_arra = chr2loc.split(":")
    chr1=chr1_arra[0]
    chr2=chr2_arra[0]
    if(len(chr1_arra)==3):
        c1pos1=chr1_arra[1]
        c1pos2=chr1_arra[2]
    if(len(chr2_arra)==3):
        c2pos1=chr2_arra[1]
        c2pos2=chr2_arra[2]
    list1 = readHeader(req, chr1, chr2, [c1pos1, c1pos2, c2pos1, c2pos2])
    if(list1==-1):
        return -1;
    master=list1[0]
    chr1ind=list1[1]
    chr2ind=list1[2]
    c1pos1=int(list1[3])
    c1pos2=int(list1[4])
    c2pos1=int(list1[5])
    c2pos2=int(list1[6])
    c1=min(chr1ind,chr2ind)
    c2=max(chr1ind,chr2ind)
    origRegionIndices=[]
    regionIndices=[]
    if (chr1ind > chr2ind):
        origRegionIndices.append(c2pos1)
        origRegionIndices.append(c2pos2)
        origRegionIndices.append(c1pos1)
        origRegionIndices.append(c1pos2)
        regionIndices.append(int(c2pos1/binsize))
        regionIndices.append(int(c2pos2/binsize))
        regionIndices.append(int(c1pos1/binsize))
        regionIndices.append(int(c1pos2/binsize))
    else:
        origRegionIndices.append(c1pos1)
        origRegionIndices.append(c1pos2)
        origRegionIndices.append(c2pos1)
        origRegionIndices.append(c2pos2)
        regionIndices.append(int(c1pos1/binsize))
        regionIndices.append(int(c1pos2/binsize))
        regionIndices.append(int(c2pos1/binsize))
        regionIndices.append(int(c2pos2/binsize))
    list1 = readFooter(req, master, c1, c2, norm, unit, binsize)
    if(list1==-1):
        return -1
    myFilePos=list1[0]
    c1NormEntry=list1[1]
    c2NormEntry=list1[2]
    if (norm != "NONE"):
        c1Norm = readNormalizationVector(req, c1NormEntry)
        c2Norm = readNormalizationVector(req, c2NormEntry)
    list1 = readMatrix(req, myFilePos, unit, binsize)
    # if(list1==-1):
    #     return -1
    blockBinCount=list1[0]
    blockColumnCount=list1[1]
    blockNumbers = getBlockNumbersForRegionFromBinPosition(regionIndices, blockBinCount, blockColumnCount, c1==c2)
    yActual=[]
    xActual=[]
    counts=[]
    for i_set in (blockNumbers):
        records=readBlock(req, i_set)
        # print(str(records))
        for j in range(len(records)):
            rec=records[j]
            x=rec['binX']*binsize
            y=rec['binY']*binsize
            c=rec['counts']
            if (norm != "NONE"):
                a=c1Norm[rec['binX']]*c2Norm[rec['binY']]
                if (a!=0.0):
                    c=c/(c1Norm[rec['binX']]*c2Norm[rec['binY']])
                else:
                    c="inf"
            if ((x>=origRegionIndices[0] and x<=origRegionIndices[1] and y>=origRegionIndices[2] and y<=origRegionIndices[3]) or ((c1==c2) and y>=origRegionIndices[0] and y<=origRegionIndices[1] and x>= origRegionIndices[2] and x<=origRegionIndices[3])):
	            xActual.append(x)
	            yActual.append(y)
	            counts.append(c)
    return [xActual, yActual, counts]

def juicerstraw(norm,infile,chr1loc,chr2loc,unit,binsize):
    binsize = int(binsize)
    magic_string = ""
    try: req=open(infile, 'rb')
    except:
        filename=wget.download(infile)
        req=open(filename,'rb')
    result = straw(norm, req, binsize, chr1loc, chr2loc, unit)
    if (result == -1):
        return;
    return result
