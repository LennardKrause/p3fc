def kappa_to_euler(k_omg, kappa, alpha, k_phi):
    '''
     converts kappa to eulerian geometry
     needed for .sfrm header
     r_: in radiant 
     k_: kappa geometry
     e_: euler geometry
     alternative delta:
     r_delta_ = np.arcsin(np.cos(np.deg2rad(alpha)) * np.sin(np.deg2rad(kappa) / 2.0) / np.cos(r_e_chi / 2.0))
    '''
    import numpy as np
    r_k_omg = np.deg2rad(k_omg)
    r_e_chi = 2.0 * np.arcsin(np.sin(np.deg2rad(kappa) / 2.0) * np.sin(np.deg2rad(alpha)))
    r_k_phi = np.deg2rad(k_phi)
    r_delta = np.arccos(np.cos(np.deg2rad(kappa) / 2.0) / np.cos(r_e_chi / 2.0))
    e_chi = np.round(np.rad2deg(r_e_chi), 5)
    e_omg = np.round(np.rad2deg(r_k_omg + r_delta), 5)
    e_phi = np.round(np.rad2deg(r_k_phi + r_delta), 5)
    return e_omg, e_chi, e_phi

def read_photon2_raw(fname, dim1, dim2, bytecode):
    '''
     Read a PHOTON-II raw image file
      - endianness unchecked
      - no header, pure data
    '''
    import numpy as np
    # translate the bytecode to the bytes per pixel
    bpp = len(np.array(0, bytecode).tostring())
    # determine the image size
    size = dim1 * dim2 * bpp
    # open the file
    with open(fname, 'rb') as f:
        # read the image (bytestream)
        rawData = f.read(size)
    # reshape the image into 2d array (dim1, dim2)
    # dtype = bytecode
    data = np.frombuffer(rawData, bytecode).reshape((dim1, dim2))
    return data

def read_sfrm(fname):
    import re
    import numpy as np
    '''
     Read Bruker .sfrm frame
     - header is returned as continuous stream
     - information read from header 
       - detector dimensions (NROWS, NCOLS)
       - bytes per pixel of image (NPIXELB)
       - number of pixels in 16 and 32 bit overflowtables (NOVERFL)
     - data is returned as uint32 2D-Array
    '''
    #def chunkstring(string, length):
    #    '''
    #     return header as list of tuples
    #      - splits once at ':'
    #      - keys and values are stripped strings
    #      - values with more than 1 entry are un-splitted
    #    '''
    #    return list(tuple(map(lambda i: i.strip(), string[0+i:length+i].split(':', 1))) for i in range(0, len(string), length)) 
    #header_list = chunkstring(header, 80)
    with open(fname, 'rb') as f:
        # read the first 512 bytes
        # find keyword 'HDRBLKS' 
        header_0 = f.read(512).decode()
        # header consists of HDRBLKS x 512 byte blocks
        header_blocks = int(re.findall(r'\s*HDRBLKS\s*:\s*(\d+)', header_0)[0])
        # read the remaining header
        header = header_0 + f.read(header_blocks * 512 - 512).decode()
        # extract frame info:
        # - rows, cols (NROWS, NCOLS)
        # - bytes-per-pixel of image (NPIXELB)
        # - length of 16 and 32 bit overflow tables (NOVERFL)
        nrows = int(re.findall(r'\s*NROWS\s*:\s*(\d+)', header)[0])
        ncols = int(re.findall(r'\s*NCOLS\s*:\s*(\d+)', header)[0])
        npixb = int(re.findall(r'\s*NPIXELB\s*:\s*(\d+)', header)[0])
        nov16, nov32 = list(map(int, re.findall(r'\s*NOVERFL\s*:\s*-*\d+\s+(\d+)\s+(\d+)', header)[0]))
        # calculate the size of the image
        im_size = nrows * ncols * npixb
        # bytes-per-pixel to datatype
        bpp2dt = [None, np.uint8, np.uint16, None, np.uint32]
        # reshape data, set datatype to np.uint32
        data = np.frombuffer(f.read(im_size), bpp2dt[npixb]).reshape((nrows, ncols)).astype(np.uint32)
        # read the 16 bit overflow table
        # table is padded to a multiple of 16 bytes
        read_16 = int(np.ceil(nov16 * 2 / 16)) * 16
        # read the table, trim the trailing zeros
        table_16 = np.trim_zeros(np.frombuffer(f.read(read_16), np.uint16))
        # read the 32 bit overflow table
        # table is padded to a multiple of 16 bytes
        read_32 = int(np.ceil(nov32 * 4 / 16)) * 16
        # read the table, trim the trailing zeros
        table_32 = np.trim_zeros(np.frombuffer(f.read(read_32), np.uint32))
        # assign values from 16 bit overflow table
        data[data == 255] = table_16
        # assign values from 32 bit overflow table
        data[data == 65535] = table_32
        return header, data

def decByteOffset_np(stream, dtype="int64"):
    '''
    The following code is taken from the FabIO package:
    Version: fabio-0.9.0
    Home-page: http://github.com/silx-kit/fabio
    Author: Henning Sorensen, Erik Knudsen, Jon Wright, Regis Perdreau,
            Jérôme Kieffer, Gael Goret, Brian Pauw, Valentin Valls
    '''
    """
    Analyze a stream of char with any length of exception:
                2, 4, or 8 bytes integers

    @param stream: string representing the compressed data
    @param size: the size of the output array (of longInts)
    @return: 1D-ndarray
    """
    import numpy as np
    listnpa = []
    key16 = b"\x80"
    key32 = b"\x00\x80"
    key64 = b"\x00\x00\x00\x80"
    shift = 1
    while True:
        idx = stream.find(key16)
        if idx == -1:
            listnpa.append(np.frombuffer(stream, dtype="int8"))
            break
        listnpa.append(np.frombuffer(stream[:idx], dtype="int8"))

        if stream[idx + 1:idx + 3] == key32:
            if stream[idx + 3:idx + 7] == key64:
                # 64 bits int
                res = np.frombuffer(stream[idx + 7:idx + 15], dtype="int64")
                listnpa.append(res)
                shift = 15
            else:
                # 32 bits int
                res = np.frombuffer(stream[idx + 3:idx + 7], dtype="int32")
                listnpa.append(res)
                shift = 7
        else:  # int16
            res = np.frombuffer(stream[idx + 1:idx + 3], dtype="int16")
            listnpa.append(res)
            shift = 3
        stream = stream[idx + shift:]
    return np.ascontiguousarray(np.hstack(listnpa), dtype).cumsum()

def read_pilatus_cbf(fname, *args):
    '''
     
    '''
    import numpy as np
    import re
    with open(fname, 'rb') as f:
        stream = f.read()
    start = stream.find(b'\x0c\x1a\x04\xd5') +4
    head = str(stream[:start])
    size = int(re.search(r'X-Binary-Size:\s+(\d+)', head).group(1))
    dim1 = int(re.search(r'X-Binary-Size-Fastest-Dimension:\s+(\d+)', head).group(1))
    dim2 = int(re.search(r'X-Binary-Size-Second-Dimension:\s+(\d+)', head).group(1))
    data = decByteOffset_np(stream[start:start+size]).reshape((dim2, dim1))
    return head, data

def read_pilatus_tif(fname, rows, cols, offset, bytecode):
    '''
     
    '''
    import numpy as np
    # translate the bytecode to the bytes per pixel
    bpp = len(np.array(0, bytecode).tostring())
    # determine the image size
    size = rows * cols * bpp
    # open the file
    with open(fname, 'rb') as f:
        # read the header
        h = f.read(offset)    
        # read the image (bytestream)
        rawData = f.read(size)
    header = str(h)
    # reshape the image into 2d array (rows, cols)
    # dtype = bytecode
    data = np.frombuffer(rawData, bytecode).reshape((rows, cols))
    return header, data

def read_pilatus_tif_gz(fname, rows, cols, offset, bytecode):
    '''
     
    '''
    import gzip
    import numpy as np
    # translate the bytecode to the bytes per pixel
    bpp = len(np.array(0, bytecode).tostring())
    # determine the image size
    size = rows * cols * bpp
    # open the file
    with gzip.open(fname, 'rb') as f:
        # read the header
        h = f.read(offset)    
        # read the image (bytestream)
        rawData = f.read(size)
    header = str(h)
    # reshape the image into 2d array (rows, cols)
    # dtype = bytecode
    data = np.frombuffer(rawData, bytecode).reshape((rows, cols))
    return header, data

def pilatus_pad(data, fill=-2, pad=8):
    import numpy as np
    # get the frame saint ready 
    # - multiple of 'pad' pixels
    # - pad with 'fill'
    # - return the offsets to calculate a new beamcenter
    (rows, cols) = data.shape
    pad_rows = int(np.ceil(rows / pad) * pad)
    pad_cols = int(np.ceil(cols / pad) * pad)
    offset_rows = (pad_rows - rows) // 2
    offset_cols = (pad_cols - cols) // 2
    padded = np.zeros((pad_rows, pad_cols), dtype=np.int32)
    padded.fill(fill)
    padded[offset_rows:offset_rows + rows, offset_cols:offset_cols + cols] = data
    return padded, offset_rows, offset_cols

def bruker_header():
    '''
     default Bruker header
    '''
    import collections
    import numpy as np
    
    header = collections.OrderedDict()
    header['FORMAT']  = np.array([100], dtype=np.int64)                       # Frame Format -- 86=SAXI, 100=Bruker
    header['VERSION'] = np.array([18], dtype=np.int64)                        # Header version number
    header['HDRBLKS'] = np.array([15], dtype=np.int64)                        # Header size in 512-byte blocks
    header['TYPE']    = ['Some Frame']                                        # String indicating kind of data in the frame
    header['SITE']    = ['Some Site']                                         # Site name
    header['MODEL']   = ['?']                                                 # Diffractometer model
    header['USER']    = ['USER']                                              # Username
    header['SAMPLE']  = ['']                                                  # Sample ID
    header['SETNAME'] = ['']                                                  # Basic data set name
    header['RUN']     = np.array([1], dtype=np.int64)                         # Run number within the data set
    header['SAMPNUM'] = np.array([1], dtype=np.int64)                         # Specimen number within the data set
    header['TITLE']   = ['', '', '', '', '', '', '', '', '']                  # User comments (8 lines)
    header['NCOUNTS'] = np.array([-9999, 0], dtype=np.int64)                  # Total frame counts, Reference detector counts
    header['NOVERFL'] = np.array([-1, 0, 0], dtype=np.int64)                  # SAXI Format: Number of overflows
                                                                              # Bruker Format: #Underflows; #16-bit overfl; #32-bit overfl
    header['MINIMUM'] = np.array([-9999], dtype=np.int64)                     # Minimum pixel value
    header['MAXIMUM'] = np.array([-9999], dtype=np.int64)                     # Maximum pixel value
    header['NONTIME'] = np.array([-2], dtype=np.int64)                        # Number of on-time events
    header['NLATE']   = np.array([0], dtype=np.int64)                         # Number of late events for multiwire data
    header['FILENAM'] = ['unknown.sfrm']                                      # (Original) frame filename
    header['CREATED'] = ['01-Jan-2000 01:01:01']                              # Date and time of creation
    header['CUMULAT'] = np.array([20.0], dtype=np.float64)                    # Accumulated exposure time in real hours
    header['ELAPSDR'] = np.array([10.0, 10.0], dtype=np.float64)              # Requested time for this frame in seconds
    header['ELAPSDA'] = np.array([10.0, 10.0], dtype=np.float64)              # Actual time for this frame in seconds
    header['OSCILLA'] = np.array([0], dtype=np.int64)                         # Nonzero if acquired by oscillation
    header['NSTEPS']  = np.array([1], dtype=np.int64)                         # steps or oscillations in this frame
    header['RANGE']   =  np.array([1.0], dtype=np.float64)                    # Magnitude of scan range in decimal degrees
    header['START']   = np.array([0.0], dtype=np.float64)                     # Starting scan angle value, decimal deg
    header['INCREME'] = np.array([1.0], dtype=np.float64)                     # Signed scan angle increment between frames
    header['NUMBER']  = np.array([1], dtype=np.int64)                         # Number of this frame in series (zero-based)
    header['NFRAMES'] = np.array([1], dtype=np.int64)                         # Number of frames in the series
    header['ANGLES']  = np.array([0.0, 0.0, 0.0, 0.0], dtype=np.float64)      # Diffractometer setting angles, deg. (2Th, omg, phi, chi)
    header['NOVER64'] = np.array([0, 0, 0], dtype=np.int64)                   # Number of pixels > 64K
    header['NPIXELB'] = np.array([1, 2], dtype=np.int64)                      # Number of bytes/pixel; Number of bytes per underflow entry
    header['NROWS']   = np.array([512, 1], dtype=np.int64)                    # Number of rows in frame; number of mosaic tiles in Y; dZ/dY value
                                                                              # for each mosaic tile, X varying fastest
    header['NCOLS']   = np.array([512, 1], dtype=np.int64)                    # Number of pixels per row; number of mosaic tiles in X; dZ/dX
                                                                              # value for each mosaic tile, X varying fastest
    header['WORDORD'] = np.array([0], dtype=np.int64)                         # Order of bytes in word; always zero (0=LSB first)
    header['LONGORD'] = np.array([0], dtype=np.int64)                         # Order of words in a longword; always zero (0=LSW first
    header['TARGET']  = ['Mo']                                                # X-ray target material)
    header['SOURCEK'] = np.array([0.0], dtype=np.float64)                     # X-ray source kV
    header['SOURCEM'] = np.array([0.0], dtype=np.float64)                     # Source milliamps
    header['FILTER']  = ['?']                                                 # Text describing filter/monochromator setting
    header['CELL']    = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0], dtype=np.float64) # Cell constants, 2 lines  (A,B,C,Alpha,Beta,Gamma)
    header['MATRIX']  = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], dtype=np.float64) # Orientation matrix, 3 lines
    header['LOWTEMP'] = np.array([1, -17300, -6000], dtype=np.int64)          # Low temp flag; experiment temperature*100; detector temp*100
    header['ZOOM']    = np.array([0.0, 0.0, 1.0], dtype=np.float64)           # Image zoom Xc, Yc, Mag
    header['CENTER']  = np.array([256.0, 256.0, 256.0, 256.0], dtype=np.float64) # X, Y of direct beam at 2-theta = 0
    header['DISTANC'] = np.array([5.0], dtype=np.float64)                     # Sample-detector distance, cm
    header['TRAILER'] = np.array([0], dtype=np.int64)                         # Byte pointer to trailer info (unused; obsolete)
    header['COMPRES'] = ['none']                                              # Text describing compression method if any
    header['LINEAR']  = np.array([1.0, 0.0], dtype=np.float64)                # Linear scale, offset for pixel values
    header['PHD']     = np.array([0.0, 0.0], dtype=np.float64)                # Discriminator settings
    header['PREAMP']  = np.array([0], dtype=np.int64)                         # Preamp gain setting
    header['CORRECT'] = ['UNKNOWN']                                           # Flood correction filename
    header['WARPFIL'] = ['UNKNOWN']                                           # Spatial correction filename
    header['WAVELEN'] = np.array([0.1, 0.1, 0.1], dtype=np.float64)           # Wavelengths (average, a1, a2)
    header['MAXXY']   = np.array([1, 1], dtype=np.int64)                      # X,Y pixel # of maximum counts
    header['AXIS']    = np.array([2], dtype=np.int64)                         # Scan axis (1=2-theta, 2=omega, 3=phi, 4=chi)
    header['ENDING']  = np.array([0.0, 0.5, 0.0, 0.0], dtype=np.float64)      # Setting angles read at end of scan
    header['DETPAR']  = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0], dtype=np.float64) # Detector position corrections (Xc,Yc,Dist,Pitch,Roll,Yaw)
    header['LUT']     = ['lut']                                               # Recommended display lookup table
    header['DISPLIM'] = np.array([0.0, 0.0], dtype=np.float64)                # Recommended display contrast window settings
    header['PROGRAM'] = ['Python Image Conversion']                           # Name and version of program writing frame
    header['ROTATE']  = np.array([0], dtype=np.int64)                         # Nonzero if acquired by rotation (GADDS)
    header['BITMASK'] = ['$NULL']                                             # File name of active pixel mask (GADDS)
    header['OCTMASK'] = np.array([0, 0, 0, 0, 0, 0, 0, 0], dtype=np.int64)    # Octagon mask parameters (GADDS) #min x, min x+y, min y, max x-y, max x, max x+y, max y, max y-x
    header['ESDCELL'] = np.array([0.001, 0.001, 0.001, 0.02, 0.02, 0.02], dtype=np.float64) # Cell ESD's, 2 lines (A,B,C,Alpha,Beta,Gamma)
    header['DETTYPE'] = ['Unknown']                                           # Detector type
    header['NEXP']    = np.array([1, 0, 0, 0, 0], dtype=np.int64)             # Number exposures in this frame; CCD bias level*100,;
                                                                              # Baseline offset (usually 32); CCD orientation; Overscan Flag
    header['CCDPARM'] = np.array([0.0, 0.0, 0.0, 0.0, 0.0], dtype=np.float64) # CCD parameters for computing pixel ESDs; readnoise, e/ADU, e/photon, bias, full scale
    header['CHEM']    = ['?']                                                 # Chemical formula
    header['MORPH']   = ['?']                                                 # CIFTAB string for crystal morphology
    header['CCOLOR']  = ['?']                                                 # CIFTAB string for crystal color
    header['CSIZE']   = ['?']                                                 # String w/ 3 CIFTAB sizes, density, temp
    header['DNSMET']  = ['?']                                                 # CIFTAB string for density method
    header['DARK']    = ['NONE']                                              # Dark current frame name
    header['AUTORNG'] = np.array([0.0, 0.0, 0.0, 0.0, 1.0], dtype=np.float64) # Autorange gain, time, scale, offset, full scale
    header['ZEROADJ'] = np.array([0.0, 0.0, 0.0, 0.0], dtype=np.float64)      # Adjustments to goniometer angle zeros (tth, omg, phi, chi)
    header['XTRANS']  = np.array([0.0, 0.0, 0.0], dtype=np.float64)           # Crystal XYZ translations
    header['HKL&XY']  = np.array([0.0, 0.0, 0.0, 0.0, 0.0], dtype=np.float64) # HKL and pixel XY for reciprocal space (GADDS)
    header['AXES2']   = np.array([0.0, 0.0, 0.0, 0.0], dtype=np.float64)      # Diffractometer setting linear axes (4 ea) (GADDS)
    header['ENDING2'] = np.array([0.0, 0.0, 0.0, 0.0], dtype=np.float64)      # Actual goniometer axes @ end of frame (GADDS)
    header['FILTER2'] = np.array([0.0, 0.0, 0.0, 1.0], dtype=np.float64)      # Monochromator 2-theta, roll (both deg)
    header['LEPTOS']  = ['']
    header['CFR']     = ['']
    return header
    
def write_bruker_frame(fname, fheader, fdata):
    '''
     write a bruker image
    '''
    import numpy as np
    
    ########################
    ## write_bruker_frame ##
    ##     FUNCTIONS      ##
    ########################
    def pad_table(table, bpp):
        '''
         pads a table with zeros to a multiple of 16 bytes
        '''
        padded = np.zeros(int(np.ceil(table.size * abs(bpp) / 16)) * 16 // abs(bpp)).astype(_BPP_TO_DT[bpp])
        padded[:table.size] = table
        return padded
        
    def format_bruker_header(fheader):
        '''
         
        '''
        format_dict = {(1,   'int64'): '{:<71d} ',
                       (2,   'int64'): '{:<35d} {:<35d} ',
                       (3,   'int64'): '{:<23d} {:<23d} {:<23d} ',
                       (4,   'int64'): '{:<17d} {:<17d} {:<17d} {:<17d} ',
                       (5,   'int64'): '{:<13d} {:<13d} {:<13d} {:<13d} {:<13d}   ',
                       (6,   'int64'): '{:<11d} {:<11d} {:<11d} {:<11d} {:<11d} {:<11d} ',
                       (1,   'int32'): '{:<71d} ',
                       (2,   'int32'): '{:<35d} {:<35d} ',
                       (3,   'int32'): '{:<23d} {:<23d} {:<23d} ',
                       (4,   'int32'): '{:<17d} {:<17d} {:<17d} {:<17d} ',
                       (5,   'int32'): '{:<13d} {:<13d} {:<13d} {:<13d} {:<13d}   ',
                       (6,   'int32'): '{:<11d} {:<11d} {:<11d} {:<11d} {:<11d} {:<11d} ',
                       (1, 'float64'): '{:<71f} ',
                       (2, 'float64'): '{:<35f} {:<35f} ',
                       (3, 'float64'): '{:<23f} {:<23f} {:<23f} ',
                       (4, 'float64'): '{:<17f} {:<17f} {:<17f} {:<17f} ',
                       (5, 'float64'): '{:<13f} {:<13f} {:<13f} {:<13f} {:<15f} '}
    
        headers = []
        for name, entry in fheader.items():
            # TITLE has multiple lines
            if name == 'TITLE':
                name = '{:<7}:'.format(name)
                number = len(entry)
                for line in range(8):
                    if number < line:
                        headers.append(''.join((name, '{:<72}'.format(entry[line]))))
                    else:
                        headers.append(''.join((name, '{:<72}'.format(' '))))
                continue
    
            # DETTYPE Mixes Entry Types
            if name == 'DETTYPE':
                name = '{:<7}:'.format(name)
                string = '{:<20s} {:<11f} {:<11f} {:<1d} {:<11f} {:<10f} {:<1d} '.format(*entry)
                headers.append(''.join((name, string)))
                continue
            
            # format the name
            name = '{:<7}:'.format(name)
            
            # pad entries
            if type(entry) == list or type(entry) == str:
                headers.append(''.join(name + '{:<72}'.format(entry[0])))
                continue
            
            # fill empty fields
            if entry.shape[0] == 0:
                headers.append(name + '{:72}'.format(' '))
                continue
            
            # if line has too many entries e.g.
            # OCTMASK(8): np.int64
            # CELL(6), MATRIX(9), DETPAR(6), ESDCELL(6): np.float64
            # write the first 6 (np.int64) / 5 (np.float64) entries
            # and the remainder in a new line/entry
            if entry.shape[0] > 6 and entry.dtype == np.int64:
                while entry.shape[0] > 6:
                    format_string = format_dict[(6, str(entry.dtype))]
                    headers.append(''.join(name + format_string.format(*entry[:6])))
                    entry = entry[6:]
            elif entry.shape[0] > 5 and entry.dtype == np.float64:
                while entry.shape[0] > 5:
                    format_string = format_dict[(5, str(entry.dtype))]
                    headers.append(''.join(name + format_string.format(*entry[:5])))
                    entry = entry[5:]
            
            # format line
            format_string = format_dict[(entry.shape[0], str(entry.dtype))]
            headers.append(''.join(name + format_string.format(*entry)))
    
        # add header ending
        if headers[-1][:3] == 'CFR':
            headers = headers[:-1]
        padding = 512 - (len(headers) * 80 % 512)
        end = '\x1a\x04'
        if padding <= 80:
            start = 'CFR: HDR: IMG: '
            padding -= len(start) + 2
            dots = ''.join(['.'] * padding)
            headers.append(start + dots + end)
        else:
            while padding > 80:
                headers.append(end + ''.join(['.'] * 78))
                padding -= 80
            if padding != 0:
                headers.append(end + ''.join(['.'] * (padding - 2)))
        return ''.join(headers)
    ########################
    ## write_bruker_frame ##
    ##   FUNCTIONS END    ##
    ########################
    
    # assign bytes per pixel to numpy integers
    # int8   Byte (-128 to 127)
    # int16  Integer (-32768 to 32767)
    # int32  Integer (-2147483648 to 2147483647)
    # uint8  Unsigned integer (0 to 255)
    # uint16 Unsigned integer (0 to 65535)
    # uint32 Unsigned integer (0 to 4294967295)
    _BPP_TO_DT = {1: np.uint8,
                  2: np.uint16,
                  4: np.uint32,
                 -1: np.int8,
                 -2: np.int16,
                 -4: np.int32}
    
    # read the bytes per pixel
    # frame data (bpp), underflow table (bpp_u)
    bpp, bpp_u = fheader['NPIXELB']
    
    # generate underflow table
    # does not work as APEXII reads the data as uint8/16/32!
    if fheader['NOVERFL'][0] >= 0:
        data_underflow = fdata[fdata <= 0]
        fheader['NOVERFL'][0] = data_underflow.shape[0]
        table_underflow = pad_table(data_underflow, -1 * bpp_u)
        fdata[fdata < 0] = 0

    # generate 32 bit overflow table
    if bpp < 4:
        data_over_uint16 = fdata[fdata >= 65535]
        table_data_uint32 = pad_table(data_over_uint16, 4)
        fheader['NOVERFL'][2] = data_over_uint16.shape[0]
        fdata[fdata >= 65535] = 65535

    # generate 16 bit overflow table
    if bpp < 2:
        data_over_uint8 = fdata[fdata >= 255]
        table_data_uint16 = pad_table(data_over_uint8, 2)
        fheader['NOVERFL'][1] = data_over_uint8.shape[0]
        fdata[fdata >= 255] = 255

    # shrink data to desired bpp
    fdata = fdata.astype(_BPP_TO_DT[bpp])
    
    # write frame
    with open(fname, 'wb') as brukerFrame:
        brukerFrame.write(format_bruker_header(fheader).encode('ASCII'))
        brukerFrame.write(fdata.tobytes())
        if fheader['NOVERFL'][0] >= 0:
            brukerFrame.write(table_underflow.tobytes())
        if bpp < 2 and fheader['NOVERFL'][1] > 0:
            brukerFrame.write(table_data_uint16.tobytes())
        if bpp < 4 and fheader['NOVERFL'][2] > 0:
            brukerFrame.write(table_data_uint32.tobytes())

def fix_bad_pixel(data, flag, bad_int=-2, sat_val=2**20):
    '''
     a bunch of different (unpolished!) ideas on how to deal with bad pixels,
       just for internal testing, no real application as the bad pixels are
       masked using saints active pixel mask
    '''
    import numpy as np
    # set bad pixels to sat_val
    if flag == 's':
        data[data == bad_int] = sat_val
    # set bad pixels to the 3x3 average
    # ignoring adjacent bad pixels
    # a1: iteratively until no more than
    # 4 bad pixels are within a 3x3 matrix
    # -> average over at least 5 pixels!
    elif flag == 'a1':
        A = list(np.argwhere(data == bad_int))
        while len(A) > 0:
            for n,i in enumerate(A):
                x, y = i
                # slice 3x3 array around bad pix
                # x and y must not get negative!
                temp_ar = data[max(0, x-1):x+2, max(0, y-1):y+2].copy()
                # if more than 4 bad pix in 3x3,
                # skip this pix, move on and 
                # iterate until only 4 bad pix
                if (temp_ar == bad_int).sum() >= 5:
                    continue
                # set bad pix to average
                # neglect bad pixels
                else:
                    temp_av = int(np.average(temp_ar[temp_ar >= 0]))
                    data[x, y] = temp_av
                    A.pop(n)
                    break
    # a2/a3: only use the averageing for bad pixels
    # where there is no adjacent bad pixel, set the
    # remainder to either zero (a2) or sat_val (a3)
    elif flag == 'a2' or flag == 'a3':
        A = np.argwhere(data == bad_int)
        for x, y in A:
            # slice 3x3 array around bad pix
            # x and y must not get negative!
            temp_ar = data[max(0, x-1):x+2, max(0, y-1):y+2]
            # if more bad pix in 3x3
            # skip it and set it later to zero
            if (temp_ar == bad_int).sum() > 1:
                continue
            # set bad pix to mean
            else:
                temp_av = int(np.average(temp_ar[temp_ar >= 0]))
                data[x, y] = temp_av
        # set adjacent bad pix to val
        if flag == 'a2':
            data[data == bad_int] = 0
        elif flag == 'a3':
            data[data == bad_int] = sat_val
    # set bad pixels to zero
    elif flag == 'z':
        data[data == bad_int] = 0
    return data

def get_run_info(basename):
    # try to get the run and frame number from the filename
    # any_name_runNum_frmNum is assumed.
    # maybe except for index error as well!
    try:
        _split = basename.split('_')
        frmTmp = _split.pop()
        frmLen = len(frmTmp)
        frmNum = int(frmTmp)
        runNum = int(_split.pop())
        stem = '_'.join(_split)
        return stem, runNum, frmNum, frmLen
    except ValueError:
        pass
    # SP-8 naming convention
    frmNum = int(basename[-3:])
    runNum = int(basename[-5:-3])
    stem = basename[:-6]
    return stem, runNum, frmNum, 3
    
def convert_frame_APS_Bruker(fname, path_sfrm, rows=1043, cols=981, offset=4096, overwrite=True, beamflux=None):
    '''
    
    '''
    import os, re
    import numpy as np
    from datetime import datetime as dt
    
    # split path, name and extension
    path_to, frame_name = os.path.split(fname)
    basename, ext = os.path.splitext(frame_name)
    frame_stem, frame_run, frame_num, _ = get_run_info(basename)
    
    # output file format: some_name_rr_ffff.sfrm
    outName = os.path.join(path_to, path_sfrm, '{}_{:>02}_{:>04}.sfrm'.format(frame_stem, frame_run, frame_num))

    # check if file exists and overwrite flag
    if os.path.exists(outName) and overwrite == False:
        return False
    
    # read in the frame
    header, data = read_pilatus_tif(fname, rows, cols, offset, np.int32)
    
    # get the frame saint ready 
    # - pad with zeros
    data, offset_rows, offset_cols = pilatus_pad(data)
    
    # the frame has to be rotated by 90 degrees
    data = np.rot90(data, k=1, axes=(1, 0))
    
    # the dead areas are flagged -1
    data[data == -1] = 0
    # bad pixels are flagged -2
    data[data == -2] = 0
    
    # scale the data to avoid underflow tables
    # should yield zero for Pilatus3 images!
    baseline_offset = -1 * data.min()
    data += baseline_offset
    
    # extract scan info from tif header
    scan_flx = float(re.search(r'Flux\s+(\d+\.\d+)', header).groups()[0])
    scan_ext = float(re.search(r'Exposure_time\s+(\d+\.\d+)\s+s', header).groups()[0])
    scan_exp = float(re.search(r'Exposure_period\s+(\d+\.\d+)\s+s', header).groups()[0])
    goni_dxt = float(re.search(r'Detector_distance\s+(\d+\.\d+)\s+m', header).groups()[0]) * 1000.0
    source_w = float(re.search(r'Wavelength\s+(\d+\.\d+)\s+A', header).groups()[0])
    goni_omg = float(re.search(r'Omega\s+(-*\d+\.\d+)\s+deg.', header).groups()[0])
    goni_kap = float(re.search(r'Kappa\s+(-*\d+\.\d+)\s+deg.', header).groups()[0])
    goni_phi = float(re.search(r'Phi\s+(-*\d+\.\d+)\s+deg.', header).groups()[0])
    goni_alp = float(re.search(r'Alpha\s+(-*\d+\.\d+)\s+deg.', header).groups()[0])
    scan_inc = float(re.search(r'Phi_increment\s+(-*\d+\.\d+)\s+deg.', header).groups()[0])
    temp     = re.search(r'Beam_xy\s+\((\d+\.\d+),\s+(\d+\.\d+)\)\s+pixels', header).groups()
    p_x, p_y = float(temp[0]), float(temp[1])
    
    # convert Kappa to Euler geometry
    goni_omg, goni_chi, goni_phi = kappa_to_euler(goni_omg, goni_kap, goni_alp, goni_phi)
    
    # adjust the beam center to the rotation
    beam_x = p_y + offset_rows
    beam_y = cols - p_x + offset_cols
    
    # APS to Bruker conversion:
    scan_inc = -scan_inc
    goni_omg =  90.0 + goni_omg
    goni_phi = 360.0 - goni_phi
    goni_tth = 0.0
    
    # Phi is the scan axis!
    scan_sta = goni_phi
    scan_end = goni_phi + scan_inc
    
    # calculate detector pixel per cm
    # this is normalized to a 512x512 detector format
    # PILATUS3-1M pixel size is 0.172 mm 
    pix_per_512 = round((10.0 / 0.172) * (512.0 / cols), 6)
    
    if beamflux is not None:
        try:
            scan_flx = beamflux[frame_run][frame_num -1]
        except IndexError:
            print('WARNING: Beamflux not found for {}!'.format(basename))
    
    # default bruker header
    header = bruker_header()
    
    # fill known header items
    header['NCOLS']      = [data.shape[1]]                           # Number of pixels per row; number of mosaic tiles in X; dZ/dX
    header['NROWS']      = [data.shape[0]]                           # Number of rows in frame; number of mosaic tiles in Y; dZ/dY value
    header['CENTER'][:]  = [beam_x, beam_y, beam_x, beam_y]          # 
    header['CCDPARM'][:] = [0.00, 1.00, 1.00, 1.00, 1169523]
    header['DETPAR'][:]  = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    header['DETTYPE'][:] = ['PILATUS3-1M', pix_per_512, 0.00, 0, 0.001, 0.0, 0]
    header['SITE']       = ['ANL/APS/15ID-D']                        # Site name
    header['MODEL']      = ['Synchrotron']                           # Diffractometer model
    header['TARGET']     = ['Undulator']                             # X-ray target material)
    header['USER']       = ['USER']                                  # Username
    header['SOURCEK']    = ['?']                                     # X-ray source kV
    header['SOURCEM']    = ['?']                                     # Source milliamps
    header['WAVELEN'][:] = [source_w, source_w, source_w]            # Wavelengths (average, a1, a2)
    header['FILENAM']    = [basename]
    header['CUMULAT']    = [scan_exp]                                # Accumulated exposure time in real hours
    header['ELAPSDR']    = [scan_ext]                                # Requested time for this frame in seconds
    header['ELAPSDA']    = [scan_exp]                                # Actual time for this frame in seconds
    header['START'][:]   = scan_sta                                  # Starting scan angle value, decimal deg
    header['ANGLES'][:]  = [goni_tth, goni_omg, scan_sta, goni_chi]  # Diffractometer setting angles, deg. (2Th, omg, phi, chi)
    header['ENDING'][:]  = [goni_tth, goni_omg, scan_end, goni_chi]  # Setting angles read at end of scan
    header['TYPE']       = ['Generic Phi Scan']                      # String indicating kind of data in the frame
    header['DISTANC']    = [goni_dxt / 10.0]                         # Sample-detector distance, cm
    header['RANGE']      = [abs(scan_inc)]                           # Magnitude of scan range in decimal degrees
    header['INCREME']    = [scan_inc]                                # Signed scan angle increment between frames
    header['NUMBER']     = [frame_num]                               # Number of this frame in series (zero-based)
    header['NFRAMES']    = ['?']                                     # Number of frames in the series
    header['AXIS'][:]    = [3]                                       # Scan axis (1=2-theta, 2=omega, 3=phi, 4=chi)
    header['LOWTEMP'][:] = [1, int((-273.15 + 20.0) * 100.0), -6000] # Low temp flag; experiment temperature*100; detector temp*100
    header['NEXP'][2]    = baseline_offset
    header['MAXXY']      = np.array(np.where(data == data.max()), float)[:, 0]
    header['MAXIMUM']    = [np.max(data)]
    header['MINIMUM']    = [np.min(data)]
    header['NCOUNTS'][:] = [data.sum(), scan_flx]
    header['NPIXELB'][:] = [1, 1]                                    # bytes/pixel in main image, bytes/pixel in underflow table
    header['NOVER64'][:] = [data[data > 64000].shape[0], 0, 0]
    header['NSTEPS']     = [1]                                       # steps or oscillations in this frame
    header['COMPRES']    = ['NONE']                                  # compression scheme if any
    header['TRAILER']    = [0]                                       # byte pointer to trailer info
    header['LINEAR'][:]  = [1.00, 0.00]     
    header['PHD'][:]     = [1.00, 0.00]
    header['OCTMASK'][:] = [0, 0, 0, 1023, 1023, 2046, 1023, 1023]
    header['DISPLIM'][:] = [0.0, 63.0]                               # Recommended display contrast window settings
    header['FILTER2'][:] = [90.0, 0.0, 0.0, 1.0]                     # Monochromator 2-theta, roll (both deg)
    header['CREATED']    = [dt.fromtimestamp(os.path.getmtime(fname)).strftime('%Y-%m-%d %H:%M:%S')]# use creation time of raw data!
    
    # write the frame
    write_bruker_frame(outName, header, data)
    return True

def convert_frame_SP8_Bruker(fname, path_sfrm, tth_corr=0.0, rows=1043, cols=981, offset=4096, overwrite=True, source_w=None):
    '''
     
    '''
    import os, re
    import numpy as np
    from datetime import datetime as dt
    
    # split path, name and extension
    path_to, frame_name = os.path.split(fname)
    basename, ext = os.path.splitext(frame_name)
    frame_stem, frame_run, frame_num, _ = get_run_info(basename)
    
    # output file format: some_name_rr_ffff.sfrm
    outName = os.path.join(path_to, path_sfrm, '{}_{:>02}_{:>04}.sfrm'.format(frame_stem, frame_run, frame_num))

    # check if file exists and overwrite flag
    if os.path.exists(outName) and overwrite == False:
        return False
    
    # read in the frame
    _, data = read_pilatus_tif(fname, rows, cols, offset, np.int32)
    
    # get the frame saint ready
    # - pad with zeros
    data, offset_rows, offset_cols = pilatus_pad(data)
    
    # the frame has to be rotated by 90 degrees
    data = np.rot90(data, k=1, axes=(1, 0))

    # sloppy simulation of a detector drift
    #data = np.roll(data, (frame_run//2, frame_run), axis=(0,1))
    
    # the dead areas are flagged -1
    data[data == -1] = 0
    
    # bad pixels are flagged -2
    data[data == -2] = 0
    
    # scale the data to avoid underflow tables
    # should yield zero for Pilatus3 images!
    baseline_offset = -1 * data.min()
    data += baseline_offset
    
    # info file name
    infFile = os.path.join(path_to, basename + '.inf')
    
    # check if info file exists
    if not os.path.isfile(infFile):
        print('ERROR: Info file is missing for: {}'.format(frame_name))
        return False
    
    # extract header information
    with open(infFile) as rFile:
        infoFile = rFile.read()
        #det_dim_x = int(re.search(r'SIZE1\s*=\s*(\d+)\s*;', infoFile).groups()[0])
        #det_dim_y = int(re.search(r'SIZE2\s*=\s*(\d+)\s*;', infoFile).groups()[0])
        #det_size_x, det_size_y = [float(i) for i in re.search(r'CCD_DETECTOR_SIZE\s*=\s*(\d+\.\d+)\s*(\d+\.\d+)\s*;', infoFile).groups()]
        det_beam_x, det_beam_y = [float(i) for i in re.search(r'CCD_SPATIAL_BEAM_POSITION\s*=\s*(-*\d+\.\d+)\s*(-*\d+\.\d+)\s*;', infoFile).groups()]
        det_maxv = int(re.search(r'SATURATED_VALUE\s*=\s*(\d+)\s*;', infoFile).groups()[0])
        if source_w is None:
            source_w = float(re.search(r'SCAN_WAVELENGTH\s*=\s*(\d+\.\d+)\s*;', infoFile).groups()[0])
        source_a = float(re.search(r'SOURCE_AMPERAGE\s*=\s*(\d+\.\d+)\s*mA\s*;', infoFile).groups()[0])
        source_v = float(re.search(r'SOURCE_VOLTAGE\s*=\s*(\d+\.\d+)\s*GeV\s*;', infoFile).groups()[0])
        goni_omg, goni_chi, goni_phi = [float(i) for i in re.search(r'CRYSTAL_GONIO_VALUES\s*=\s*(-*\d+\.\d+)\s*(-*\d+\.\d+)\s*(-*\d+\.\d+)\s*;', infoFile).groups()]
        goni_tth, goni_dxt = [float(i) for i in re.search(r'SCAN_DET_RELZERO\s*=\s*-*\d+\.\d+\s*(-*\d+\.\d+)\s*(-*\d+\.\d+)\s*;', infoFile).groups()]
        scan_rax = str(re.search(r'ROTATION_AXIS_NAME\s*=\s*(\w+)\s*;', infoFile).groups()[0])
        scan_num = float(re.search(r'SCAN_SEQ_INFO\s*=\s*\d+\s*\d+\s*(\d+)\s*;', infoFile).groups()[0])
        scan_sta, scan_end, scan_inc, scan_exp = [float(i) for i in re.search(r'SCAN_ROTATION\s*=\s*(-*\d+\.\d+)\s*(-*\d+\.\d+)\s*(-*\d+\.\d+)\s*(-*\d+\.\d+)\s*-*\d+\.\d+\s*-*\d+\.\d+\s*-*\d+\.\d+\s*-*\d+\.\d+\s*-*\d+\.\d+\s*-*\d+\.\d+\s*;', infoFile).groups()]
    
    # For some reason the distance is missing for some runs.
    # At SPring-8 the detector distance 'cannot' be changed.
    # So, we hard-code 130.0 on missing entry here!
    if goni_dxt == 0.0:
        goni_dxt = 130.0
    
    # initial frame dimensions are needed to calculate
    # the beamcenter of the reshaped frame and
    # adjust the beam center to the rotation
    beam_x = det_beam_y + offset_rows
    beam_y = cols - det_beam_x + offset_cols
    
    # 2-th were misaligned (pre 2019 data)
    goni_tth = goni_tth + (goni_tth * tth_corr)
    
    # SP8 to Bruker conversion:
    goni_chi = -goni_chi
    goni_phi = -goni_phi

    # calculate detector pixel per cm
    # this is normalized to a 512x512 detector format
    # PILATUS3-1M pixel size is 0.172 mm 
    pix_per_512 = round((10.0 / 0.172) * (512.0 / cols), 6)
    
    # convert SP8 to Bruker angles
    ax_name_to_num = {'Omega':2, 'Phi':3}
    axis_start = [goni_tth, goni_omg, goni_phi, goni_chi]
    axis_start[ax_name_to_num[scan_rax]-1] = scan_sta
    axis_end = [goni_tth, goni_omg, goni_phi, goni_chi]
    axis_end[ax_name_to_num[scan_rax]-1] = scan_end
    
    # default bruker header
    header = bruker_header()
    
    # fill known header items
    header['NROWS'][:]   = [data.shape[0], 2]                                   # Number of rows in frame; number of mosaic tiles in Y; dZ/dY value
    header['NCOLS'][:]   = [data.shape[1], 5]                                   # Number of pixels per row; number of mosaic tiles in X; dZ/dX
    header['CENTER'][:]  = [beam_x, beam_y, beam_x, beam_y]                     # adjust the beam center for the filling/cutting of the frame
    header['CCDPARM'][:] = [1.00, 1.00, 1.00, 0.00, det_maxv]                   # readnoise, electronsperadu, electronsperphoton, bruker_bias, bruker_fullscale
    header['DETPAR'][:]  = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]                       # Detector position corrections (Xc, Yc, Dist, Pitch, Roll, Yaw)
    header['DETTYPE'][:] = ['PILATUS3-1M', pix_per_512, 0.001, 0, 0.001, 0.001, 1]
    header['SITE']       = ['SPring-8/BL02B1']                                  # Site name
    header['MODEL']      = ['Synchrotron']                                      # Diffractometer model
    header['TARGET']     = ['Bending Magnet']                                   # X-ray target material)
    header['USER']       = ['USER']                                             # Username
    header['SOURCEK']    = [source_v]                                           # X-ray source kV
    header['SOURCEM']    = [source_a]                                           # Source milliamps
    header['WAVELEN'][:] = [source_w, source_w, source_w]                       # Wavelengths (average, a1, a2)
    header['FILENAM']    = [basename]
    header['CUMULAT']    = [scan_exp]                                           # Accumulated exposure time in real hours
    header['ELAPSDR']    = [scan_exp]                                           # Requested time for this frame in seconds
    header['ELAPSDA']    = [scan_exp]                                           # Actual time for this frame in seconds
    header['START'][:]   = scan_sta                                             # Starting scan angle value, decimal deg
    header['ANGLES'][:]  = axis_start                                           # Diffractometer setting angles, deg. (2Th, omg, phi, chi)
    header['ENDING'][:]  = axis_end                                             # Setting angles read at end of scan
    header['TYPE']       = ['Generic {} Scan'.format(scan_rax)]                 # String indicating kind of data in the frame
    header['DISTANC']    = [float(goni_dxt) / 10.0]                             # Sample-detector distance, cm
    header['RANGE']      = [abs(scan_inc)]                                      # Magnitude of scan range in decimal degrees
    header['INCREME']    = [scan_inc]                                           # Signed scan angle increment between frames
    header['NUMBER']     = [frame_num]                                          # Number of this frame in series (zero-based)
    header['NFRAMES']    = [int(scan_num)]                                      # Number of frames in the series
    header['AXIS'][:]    = [ax_name_to_num[scan_rax]]                           # Scan axis (1=2-theta, 2=omega, 3=phi, 4=chi)
    header['LOWTEMP'][:] = [1, int((-273.15 + 20.0) * 100.0), -6000]            # Low temp flag; experiment temperature*100; detector temp*100
    header['NEXP'][2]    = baseline_offset
    header['MAXXY']      = np.array(np.where(data == data.max()), float)[:, 0]
    header['MAXIMUM']    = [np.max(data)]
    header['MINIMUM']    = [np.min(data)]
    header['NCOUNTS'][:] = [data.sum(), 0]
    header['NOVER64'][:] = [data[data > 64000].shape[0], 0, 0]
    header['NSTEPS']     = [1]                                                  # steps or oscillations in this frame
    header['NPIXELB'][:] = [1, 1]                                               # bytes/pixel in main image, bytes/pixel in underflow table
    header['COMPRES']    = ['NONE']                                             # compression scheme if any
    header['TRAILER']    = [-1]                                                 # byte pointer to trailer info
    header['LINEAR'][:]  = [1.00, 0.00]                                         # bruker_linearscale, bruker_linearoffset
    header['PHD'][:]     = [1.00, 0.10]                                         # Phosphor efficiency, phosphor thickness
    header['PREAMP']     = [1]                                                  # Preamp gain setting
    header['CORRECT']    = ['INTERNAL']                                         # Flood correction filename
    header['DARK']       = ['INTERNAL']                                         # Dark current frame name
    header['WARPFIL']    = ['LINEAR']                                           # Spatial correction filename
    ox = data.shape[1]
    oy = data.shape[0]
    header['OCTMASK'][:] = [0, 0, 0, ox-1, ox-1, ox+oy-1, oy-1, oy-1]           # Octagon mask parameters (GADDS) #min x, min x+y, min y, max x-y, max x, max x+y, max y, max y-x
    header['DISPLIM'][:] = [0.0, 100.0]                                         # Recommended display contrast window settings
    header['FILTER2'][:] = [90.0, 0.0, 0.0, 1.0]                                # Monochromator 2-theta, roll (both deg)
    header['CREATED']    = [dt.fromtimestamp(os.path.getmtime(fname)).strftime('%Y-%m-%d %H:%M:%S')]# use creation time of raw data!
    
    # write the frame
    write_bruker_frame(outName, header, data)
    return True

def convert_frame_DLS_Bruker(fname, path_sfrm, rows=1679, cols=1475, offset=0, overwrite=True):
    '''
    
    '''
    import os, re
    import numpy as np
    from datetime import datetime as dt
    
    # split path, name and extension
    path_to, frame_name = os.path.split(fname)
    basename, ext = os.path.splitext(frame_name)
    frame_stem, frame_run, frame_num, _ = get_run_info(basename)
    
    # output file format: some_name_rr_ffff.sfrm
    outName = os.path.join(path_to, path_sfrm, '{}_{:>02}_{:>04}.sfrm'.format(frame_stem, frame_run, frame_num))

    # check if file exists and overwrite flag
    if os.path.exists(outName) and overwrite == False:
        return False
    
    # read in the frame
    header, data = read_pilatus_cbf(fname)
    
    # get the frame saint ready 
    # - multiple of 128x128 pixels
    # - pad with zeros
    data, offset_rows, offset_cols = pilatus_pad(data)
    
    # the dead areas are flagged -1
    data[data == -1] = 0
    # bad pixels are flagged -2
    data[data == -2] = 0
    
    # scale the data to avoid underflow tables
    baseline_offset = -1 * data.min()
    data += baseline_offset
    
    # extract scan info from tif header
    sca_ext = float(re.search(r'Exposure_time\s+(\d+\.\d+)\s+s', header).groups()[0])
    sca_exp = float(re.search(r'Exposure_period\s+(\d+\.\d+)\s+s', header).groups()[0])
    gon_dxt = float(re.search(r'Detector_distance\s+(\d+\.\d+)\s+m', header).groups()[0]) * 1000.0
    src_wav = float(re.search(r'Wavelength\s+(\d+\.\d+)\s+A', header).groups()[0])
    sta_phi = float(re.search(r'Phi\s+(-*\d+\.\d+)\s+deg.', header).groups()[0])
    inc_phi = float(re.search(r'Phi_increment\s+(-*\d+\.\d+)\s+deg.', header).groups()[0])
    sta_chi = float(re.search(r'Chi\s+(-*\d+\.\d+)\s+deg.', header).groups()[0])
    inc_chi = float(re.search(r'Chi_increment\s+(-*\d+\.\d+)\s+deg.', header).groups()[0])
    sta_omg = float(re.search(r'Omega\s+(-*\d+\.\d+)\s+deg.', header).groups()[0])
    inc_omg = float(re.search(r'Omega_increment\s+(-*\d+\.\d+)\s+deg.', header).groups()[0])
    sta_tth = float(re.search(r'Detector_2theta\s+(-*\d+\.\d+)\s+deg.', header).groups()[0])
    pil_x   = float(re.search(r'Beam_xy\s+\((\d+\.\d+),\s+(\d+\.\d+)\)\s+pixels', header).groups()[0])
    pil_y   = float(re.search(r'Beam_xy\s+\((\d+\.\d+),\s+(\d+\.\d+)\)\s+pixels', header).groups()[1])
    
    # initial frame dimensions are needed to calculate
    # the beamcenter of the reshaped frame and
    # adjust the beam center to the Bruker flip
    # numpy array starts in the upper left corner
    # Bruker starts lower left
    beam_y = rows - pil_y + offset_rows
    beam_x = pil_x + offset_cols
    
    # DLS to Bruker conversion:
    sta_tth = round(-sta_tth, 4)
    sta_omg = round(180.0 - sta_omg, 4)
    inc_omg = round(-inc_omg, 4)
    sta_chi = round(-sta_chi, 4)
    inc_chi = round(-inc_chi, 4)
    
    # ending positions
    end_phi = round(sta_phi + inc_phi, 4)
    end_chi = round(sta_chi + inc_chi, 4)
    end_omg = round(sta_omg + inc_omg, 4)
    end_tth = round(sta_tth, 4)
    
    # Scan axis (1=2-theta, 2=omega, 3=phi, 4=chi)
    ang_nam = ['Omega',  'Phi' ,  'Chi' ]
    ang_inc = [inc_omg, inc_phi, inc_chi]
    ang_sta = [sta_omg, sta_phi, sta_chi]
    sca_nam, sca_axs, sca_sta, sca_inc = [(ang_nam[i], int(i+2), ang_sta[i], round(v,4)) for i,v in enumerate(ang_inc) if v != 0.0][0]
    
    # calculate detector pixel per cm
    # this is normalized to a 512x512 detector format
    # PILATUS3 pixel size is 0.172 mm 
    pix_per_512 = round((10.0 / 0.172) * (512.0 / cols), 6)
    
    # default bruker header
    header = bruker_header()
    
    # fill known header items
    header['NCOLS']      = [data.shape[1]]                           # Number of pixels per row; number of mosaic tiles in X; dZ/dX
    header['NROWS']      = [data.shape[0]]                           # Number of rows in frame; number of mosaic tiles in Y; dZ/dY value
    header['CENTER'][:]  = [beam_x, beam_y, beam_x, beam_y]          # 
    header['CCDPARM'][:] = [0.00, 1.00, 1.00, 1.00, 1169523]
    header['DETPAR'][:]  = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    header['DETTYPE'][:] = ['PILATUS3-2M', pix_per_512, 0.00, 0, 0.001, 0.0, 0]
    header['SITE']       = ['DLS/I19-1']                             # Site name
    header['MODEL']      = ['Synchrotron']                           # Diffractometer model
    header['TARGET']     = ['Undulator']                             # X-ray target material)
    header['USER']       = ['?']                                     # Username
    header['SOURCEK']    = ['?']                                     # X-ray source kV
    header['SOURCEM']    = ['?']                                     # Source milliamps
    header['WAVELEN'][:] = [src_wav, src_wav, src_wav]               # Wavelengths (average, a1, a2)
    header['FILENAM']    = [basename]
    header['CUMULAT']    = [sca_exp]                                 # Accumulated exposure time in real hours
    header['ELAPSDR']    = [sca_ext]                                 # Requested time for this frame in seconds
    header['ELAPSDA']    = [sca_exp]                                 # Actual time for this frame in seconds
    header['START'][:]   = sca_sta                                   # Starting scan angle value, decimal deg
    header['ANGLES'][:]  = [sta_tth, sta_omg, sta_phi, sta_chi]      # Diffractometer setting angles, deg. (2Th, omg, phi, chi)
    header['ENDING'][:]  = [end_tth, end_omg, end_phi, end_chi]      # Setting angles read at end of scan
    header['TYPE']       = ['Generic {} Scan'.format(sca_nam)]       # String indicating kind of data in the frame
    header['DISTANC']    = [float(gon_dxt) / 10.0]                   # Sample-detector distance, cm
    header['RANGE']      = [abs(sca_inc)]                            # Magnitude of scan range in decimal degrees
    header['INCREME']    = [sca_inc]                                 # Signed scan angle increment between frames
    header['NUMBER']     = [frame_num]                               # Number of this frame in series (zero-based)
    header['NFRAMES']    = ['?']                                     # Number of frames in the series
    header['AXIS'][:]    = [sca_axs]                                 # Scan axis (1=2-theta, 2=omega, 3=phi, 4=chi)
    header['LOWTEMP'][:] = [1, 0, 0]                                 # Low temp flag; experiment temperature*100; detector temp*100
    header['NEXP'][2]    = baseline_offset
    header['MAXXY']      = np.array(np.where(data == data.max()), float)[:, 0]
    header['MAXIMUM']    = [np.max(data)]
    header['MINIMUM']    = [np.min(data)]
    header['NCOUNTS'][:] = [data.sum(), 0]
    header['NPIXELB'][:] = [1, 1]                                    # bytes/pixel in main image, bytes/pixel in underflow table
    header['NOVER64'][:] = [data[data > 64000].shape[0], 0, 0]
    header['NSTEPS']     = [1]                                       # steps or oscillations in this frame
    header['COMPRES']    = ['NONE']                                  # compression scheme if any
    header['TRAILER']    = [0]                                       # byte pointer to trailer info
    header['LINEAR'][:]  = [1.00, 0.00]     
    header['PHD'][:]     = [1.00, 0.00]
    header['OCTMASK'][:] = [0, 0, 0, 1023, 1023, 2046, 1023, 1023]
    header['DISPLIM'][:] = [0.0, 63.0]                               # Recommended display contrast window settings
    header['FILTER2'][:] = [90.0, 0.0, 0.0, 1.0]                     # Monochromator 2-theta, roll (both deg)
    header['CREATED']    = [dt.fromtimestamp(os.path.getmtime(fname)).strftime('%Y-%m-%d %H:%M:%S')]# use creation time of raw data!
    
    # write the frame
    write_bruker_frame(outName, header, data)
    return True
